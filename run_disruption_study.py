#!/usr/bin/env python3
"""
Concentration of Disruptive Research Among Scientists
=====================================================

Hypothesis: Only a very small percentage of scientists produce very disruptive research.

Pipeline:
  1. Download Park et al. 2023 Zenodo dataset (CD5 disruption scores)
  2. Sample STEM authors from OpenAlex
  3. Fetch author works (1980-2010)
  4. Match works to Zenodo CD5 via DOI
  5. Compute author-level disruption metrics
  6. Distribution analysis (Lorenz, Gini, power-law)
  7. Demographic profiling (logistic + OLS regressions)
  8. Field-level heterogeneity
  9. Temporal trends
 10. Validation (A/B disruption cross-check on subsample)
 11. Report and figures
"""

import csv
import gzip
import hashlib
import json
import math
import os
import random
import re
import subprocess
import sys
import tarfile
import time
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple
from urllib.parse import quote

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.api as sm

# ── Configuration ──────────────────────────────────────────────────────────────

ROOT = Path(__file__).resolve().parent
OUTPUT_DIR = Path(os.environ.get("DISRUPT_OUTPUT_DIR", str(ROOT / "output")))
CACHE_DIR = Path(os.environ.get("DISRUPT_CACHE_DIR", str(ROOT / "cache")))
DATA_DIR = Path(os.environ.get("DISRUPT_DATA_DIR", str(ROOT / "data")))

RANDOM_SEED = 42
PUB_START_YEAR = 1980
PUB_END_YEAR = 2010
TARGET_AUTHORS = int(os.environ.get("DISRUPT_TARGET_AUTHORS", "10000"))
VALIDATION_SAMPLE_SIZE = int(os.environ.get("DISRUPT_VALIDATION_N", "2000"))
FORWARD_WINDOW_YEARS = 10
MAX_CITERS_PER_WORK = 400
MIN_CITERS_FOR_D = 5
TOP_PERCENTILE = 0.05  # top 5% = "disruptive"

ZENODO_URL = "https://zenodo.org/records/7258379/files/nature_disruption_open_access.tar.gz"
OPENALEX_AUTHOR_URL = "https://api.openalex.org/authors"
OPENALEX_WORK_URL = "https://api.openalex.org/works"

# STEM fields in OpenAlex domain taxonomy
STEM_DOMAINS = {
    "Physical Sciences",
    "Life Sciences",
    "Engineering",
    "Mathematics",
    "Computer Science",
    "Health Sciences",
    "Environmental Science",
}

# Broader STEM subfields for field-level analysis
STEM_SUBFIELDS = {
    "Physics", "Chemistry", "Biology", "Engineering", "Computer Science",
    "Mathematics", "Materials Science", "Earth Sciences", "Environmental Science",
    "Agricultural Sciences", "Medicine", "Biochemistry", "Neuroscience",
}

# ── Curl Cache Client (reused from run_design_a.py) ───────────────────────────

@dataclass
class CurlCacheClient:
    cache_dir: Path
    sleep_seconds: float = 0.04
    connect_timeout_seconds: int = 15
    max_time_seconds: int = 45

    def _cache_path(self, url: str, params: Dict[str, Any]) -> Path:
        norm = url + "?" + "&".join(f"{k}={params[k]}" for k in sorted(params))
        key = hashlib.sha1(norm.encode("utf-8")).hexdigest()
        return self.cache_dir / f"{key}.json"

    def get_json(self, url: str, params: Dict[str, Any], retries: int = 8) -> Dict[str, Any]:
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        cache_path = self._cache_path(url, params)
        if cache_path.exists():
            try:
                return json.loads(cache_path.read_text())
            except json.JSONDecodeError:
                cache_path.unlink(missing_ok=True)

        cmd = [
            "curl", "-sG",
            "--connect-timeout", str(self.connect_timeout_seconds),
            "--max-time", str(self.max_time_seconds),
            "--user-agent", "disruption-concentration-study/1.0",
            url,
        ]
        for key, value in params.items():
            cmd.extend(["--data-urlencode", f"{key}={value}"])

        error: Optional[str] = None
        for attempt in range(retries):
            proc = subprocess.run(cmd, capture_output=True, text=True)
            raw = proc.stdout.strip()
            if proc.returncode == 0 and raw:
                try:
                    parsed = json.loads(raw)
                    tmp = cache_path.with_suffix(".tmp")
                    tmp.write_text(json.dumps(parsed))
                    tmp.replace(cache_path)
                    time.sleep(self.sleep_seconds)
                    return parsed
                except json.JSONDecodeError:
                    error = raw[:400]
            else:
                err_txt = proc.stderr.strip() or raw[:400]
                error = f"rc={proc.returncode} {err_txt}".strip()
            time.sleep(0.8 * (attempt + 1))

        raise RuntimeError(f"curl failed for {url} params={params}: {error}")


def short_openalex_id(url_or_id: str) -> str:
    if not url_or_id:
        return ""
    match = re.search(r"([AW]\d+)$", str(url_or_id))
    if match:
        return match.group(1)
    return str(url_or_id).replace("https://openalex.org/", "")


def fetch_openalex_pages(
    client: CurlCacheClient,
    url: str,
    params: Dict[str, Any],
    max_items: Optional[int] = None,
) -> List[Dict[str, Any]]:
    all_rows: List[Dict[str, Any]] = []
    cursor = "*"
    while True:
        q = dict(params)
        q["per-page"] = 200
        q["cursor"] = cursor
        try:
            data = client.get_json(url, q)
        except Exception:
            break
        rows = data.get("results", []) or []
        all_rows.extend(rows)
        if max_items is not None and len(all_rows) >= max_items:
            return all_rows[:max_items]
        cursor = data.get("meta", {}).get("next_cursor")
        if not cursor:
            break
    return all_rows


def earliest_publication_year(candidate: Dict[str, Any]) -> Optional[int]:
    counts = candidate.get("counts_by_year", []) or []
    years = [item.get("year") for item in counts if item.get("works_count", 0) > 0 and item.get("year")]
    if not years:
        return None
    return int(min(years))


# ── Step 1: Download and Load Zenodo Population Benchmark ─────────────────────

def download_zenodo(data_dir: Path) -> Path:
    """Download the Park et al. 2023 disruption dataset from Zenodo."""
    archive_path = data_dir / "nature_disruption_open_access.tar.gz"
    if archive_path.exists():
        print(f"  Archive already downloaded: {archive_path}", flush=True)
        return archive_path

    print(f"  Downloading from Zenodo (~1.5 GB)...", flush=True)
    cmd = [
        "curl", "-L", "-o", str(archive_path),
        "--connect-timeout", "30",
        "--max-time", "1800",
        "--retry", "3",
        "--progress-bar",
        ZENODO_URL,
    ]
    proc = subprocess.run(cmd, capture_output=False)
    if proc.returncode != 0:
        archive_path.unlink(missing_ok=True)
        raise RuntimeError(f"Failed to download Zenodo archive (rc={proc.returncode})")
    print(f"  Downloaded: {archive_path.stat().st_size / 1e9:.2f} GB", flush=True)
    return archive_path


def extract_zenodo(archive_path: Path, data_dir: Path) -> Path:
    """Extract the Zenodo archive."""
    extracted_marker = data_dir / ".extracted"
    if extracted_marker.exists():
        print("  Archive already extracted.", flush=True)
    else:
        print("  Extracting archive...", flush=True)
        with tarfile.open(archive_path, "r:gz") as tar:
            tar.extractall(path=data_dir)
        extracted_marker.write_text("done")
    return data_dir


def load_zenodo_population_benchmark(data_dir: Path) -> Optional[pd.DataFrame]:
    """Load CD5 scores from Zenodo as a population-level benchmark.

    The Park et al. 2023 replication data contains year + cd_5 + source
    but NO per-paper identifiers (no DOIs, WoS IDs, etc.). So we use it
    only as a population distribution benchmark, not for per-paper matching.
    """
    # The unified CD5 data is in unified_cdindex_df.csv.gz
    cd_file = data_dir / "nature_disruption_open_access" / "data" / "analytical" / "unified_cdindex_df.csv.gz"
    if not cd_file.exists():
        # Search for it
        candidates = list(data_dir.rglob("unified_cdindex_df.csv.gz"))
        if not candidates:
            candidates = list(data_dir.rglob("*cdindex*.csv.gz"))
        if not candidates:
            print("  WARNING: Could not find unified_cdindex_df.csv.gz", flush=True)
            return None
        cd_file = candidates[0]

    print(f"  Loading population benchmark from {cd_file.name}...", flush=True)
    try:
        df = pd.read_csv(cd_file, compression='gzip')
    except Exception as e:
        print(f"  Error loading: {e}", flush=True)
        return None

    print(f"  Loaded {len(df):,} rows, columns: {list(df.columns)}", flush=True)

    # Expected columns: year, cd_5, source
    if 'cd_5' not in df.columns:
        print("  WARNING: No cd_5 column found.", flush=True)
        return None

    # Filter to papers only (exclude patents)
    papers = df[df['source'].str.contains('WoS|Web of Science|Paper', case=False, na=False)].copy()
    if papers.empty:
        # If no clear paper source, use everything
        print("  No 'Paper' source found, using all sources.", flush=True)
        papers = df.copy()

    papers = papers.rename(columns={'cd_5': 'cd5'})
    papers['cd5'] = pd.to_numeric(papers['cd5'], errors='coerce')
    papers = papers.dropna(subset=['cd5'])

    print(f"  Population benchmark: {len(papers):,} records with CD5 scores", flush=True)
    print(f"  Sources: {papers['source'].value_counts().to_dict()}", flush=True)
    print(f"  CD5 distribution: mean={papers['cd5'].mean():.4f}, "
          f"median={papers['cd5'].median():.4f}, "
          f"std={papers['cd5'].std():.4f}", flush=True)
    print(f"  NOTE: No per-paper IDs (DOIs) available — using as population benchmark only.", flush=True)

    return papers


# ── Step 2: Sample STEM Authors from OpenAlex ─────────────────────────────────

def _is_stem_author(concepts: List[Dict[str, Any]]) -> bool:
    """Check if an author's top concepts indicate STEM affiliation."""
    stem_keywords = {
        "physics", "chemistry", "biology", "mathematics", "computer science",
        "engineering", "materials science", "medicine", "geology", "ecology",
        "astronomy", "biochemistry", "genetics", "neuroscience", "pharmacology",
        "immunology", "microbiology", "statistics", "environmental science",
        "oceanography", "meteorology", "anatomy", "physiology", "pathology",
        "radiology", "surgery", "dentistry", "nursing", "veterinary",
        "agricultural", "botany", "zoology", "organic chemistry",
        "inorganic chemistry", "nuclear", "quantum", "optics", "thermodynamics",
        "mechanical engineering", "electrical engineering", "civil engineering",
        "chemical engineering", "biomedical", "nanotechnology", "robotics",
        "artificial intelligence", "machine learning", "algorithm",
    }
    for concept in concepts[:5]:
        name = concept.get("display_name", "").lower()
        if any(kw in name for kw in stem_keywords):
            return True
    return False


def sample_stem_authors(client: CurlCacheClient, target_n: int) -> pd.DataFrame:
    """Sample STEM authors from OpenAlex.

    Strategy: broad sample with works_count/cited_by_count filters,
    then filter STEM client-side using x_concepts (the API-side concept
    filter is unreliable for the authors endpoint).
    """
    print(f"  Target: {target_n:,} unique STEM authors", flush=True)

    all_authors: Dict[str, Dict[str, Any]] = {}
    # STEM pass rate is ~12%, so we need ~8x oversampling.
    # Each seed gives ~200 authors, ~24 pass STEM filter.
    max_seeds = max(200, (target_n // 20) + 50)

    for seed_idx in range(max_seeds):
        if len(all_authors) >= target_n:
            break

        seed_val = RANDOM_SEED + seed_idx * 7

        params = {
            "filter": "works_count:>5,cited_by_count:>10",
            "select": "id,display_name,works_count,cited_by_count,summary_stats,last_known_institutions,x_concepts,counts_by_year",
            "sample": 200,
            "seed": seed_val,
        }

        try:
            data = client.get_json(OPENALEX_AUTHOR_URL, params)
        except Exception as e:
            print(f"    Seed {seed_idx} failed: {e}", flush=True)
            continue

        results = data.get("results", []) or []
        for author in results:
            aid = short_openalex_id(author.get("id", ""))
            if aid in all_authors:
                continue

            concepts = author.get("x_concepts", []) or []

            # Client-side STEM filter
            if not _is_stem_author(concepts):
                continue

            summary = author.get("summary_stats", {}) or {}
            institutions = author.get("last_known_institutions", []) or []
            inst = institutions[0] if institutions else {}
            top_concepts = [c.get("display_name", "") for c in concepts[:5]]

            all_authors[aid] = {
                "author_id": aid,
                "display_name": author.get("display_name", ""),
                "works_count": int(author.get("works_count", 0) or 0),
                "cited_by_count": int(author.get("cited_by_count", 0) or 0),
                "h_index": summary.get("h_index", np.nan),
                "i10_index": summary.get("i10_index", np.nan),
                "institution_name": inst.get("display_name", ""),
                "institution_type": inst.get("type", ""),
                "institution_country": inst.get("country_code", ""),
                "top_concept": top_concepts[0] if top_concepts else "",
                "first_pub_year": earliest_publication_year(author),
                "counts_by_year": json.dumps(author.get("counts_by_year", [])),
            }

        if (seed_idx + 1) % 10 == 0 or len(all_authors) >= target_n:
            print(f"    After seed {seed_idx + 1}: {len(all_authors):,} unique STEM authors "
                  f"(from {(seed_idx + 1) * 200} sampled)", flush=True)

    df = pd.DataFrame(list(all_authors.values())[:target_n])
    print(f"  Sampled {len(df):,} STEM authors", flush=True)
    return df


# ── Step 3: Fetch Author Works ────────────────────────────────────────────────

def fetch_author_works_for_study(
    client: CurlCacheClient, author_id: str
) -> List[Dict[str, Any]]:
    """Fetch all articles/reviews for an author published 1980-2010."""
    rows = fetch_openalex_pages(
        client,
        OPENALEX_WORK_URL,
        {
            "filter": (
                f"author.id:{author_id},"
                f"from_publication_date:{PUB_START_YEAR}-01-01,"
                f"to_publication_date:{PUB_END_YEAR}-12-31,"
                f"type:article|review"
            ),
            "select": "id,doi,publication_year,type,referenced_works,cited_by_count,authorships,primary_topic",
        },
        max_items=500,
    )

    eligible = []
    for w in rows:
        doi = w.get("doi") or ""
        if doi:
            doi = doi.replace("https://doi.org/", "").strip().lower()

        authorships = w.get("authorships", []) or []
        team_size = len(authorships)

        # Determine author position
        author_position = "unknown"
        for aship in authorships:
            a_id = short_openalex_id((aship.get("author") or {}).get("id", ""))
            if a_id == author_id:
                author_position = aship.get("author_position", "unknown")
                break

        primary_topic = w.get("primary_topic") or {}
        field = ((primary_topic.get("field") or {}).get("display_name")) or "Unknown"
        subfield = ((primary_topic.get("subfield") or {}).get("display_name")) or "Unknown"
        domain = ((primary_topic.get("domain") or {}).get("display_name")) or "Unknown"

        ref_count = len(w.get("referenced_works", []) or [])

        eligible.append({
            "work_id": short_openalex_id(w.get("id", "")),
            "doi": doi if doi else np.nan,
            "publication_year": int(w.get("publication_year", 0)),
            "cited_by_count": int(w.get("cited_by_count", 0) or 0),
            "team_size": team_size if team_size > 0 else np.nan,
            "author_position": author_position,
            "field": field,
            "subfield": subfield,
            "domain": domain,
            "ref_count": ref_count,
            "referenced_works": [short_openalex_id(x) for x in (w.get("referenced_works", []) or [])],
        })

    return eligible


# ── Step 4: Match Works to Zenodo CD5 ─────────────────────────────────────────

def match_works_to_zenodo(
    works_df: pd.DataFrame, zenodo_df: pd.DataFrame
) -> pd.DataFrame:
    """Join OpenAlex works with Zenodo CD5 scores on DOI."""
    works_with_doi = works_df.dropna(subset=["doi"]).copy()
    zenodo_with_doi = zenodo_df.dropna(subset=["doi"]).copy()

    # Normalize DOIs
    works_with_doi["doi_norm"] = works_with_doi["doi"].str.strip().str.lower()
    zenodo_with_doi["doi_norm"] = zenodo_with_doi["doi"].str.strip().str.lower()

    merged = works_with_doi.merge(
        zenodo_with_doi[["doi_norm", "cd5"]],
        on="doi_norm",
        how="left",
    )

    matched_count = merged["cd5"].notna().sum()
    total_with_doi = len(works_with_doi)
    print(f"  Works with DOIs: {total_with_doi:,}", flush=True)
    print(f"  Matched to Zenodo CD5: {matched_count:,} ({100*matched_count/max(1,total_with_doi):.1f}%)", flush=True)

    # Merge back to get all works (those without DOIs won't have cd5)
    result = works_df.copy()
    if "doi_norm" not in result.columns:
        result["doi_norm"] = result["doi"].astype(str).str.strip().str.lower()
        result.loc[result["doi_norm"].isin(["nan", "", "none"]), "doi_norm"] = np.nan
    cd5_map = merged.dropna(subset=["cd5"]).set_index("doi_norm")["cd5"]
    result["cd5"] = result["doi_norm"].map(cd5_map)

    return result


# ── Step 5: Author-Level Disruption Metrics ───────────────────────────────────

def compute_author_disruption_metrics(
    works_df: pd.DataFrame, top_threshold: float
) -> pd.DataFrame:
    """Compute author-level disruption metrics from paper-level CD5 scores."""
    scored = works_df.dropna(subset=["cd5"]).copy()
    if scored.empty:
        return pd.DataFrame()

    # Compute the top-percentile threshold
    cd5_threshold = scored["cd5"].quantile(1 - top_threshold)
    scored["is_top5pct"] = (scored["cd5"] >= cd5_threshold).astype(int)

    print(f"  CD5 threshold for top {100*top_threshold:.0f}%: {cd5_threshold:.4f}", flush=True)
    print(f"  Papers above threshold: {scored['is_top5pct'].sum():,}", flush=True)

    agg = scored.groupby("author_id").agg(
        n_scored_papers=("work_id", "count"),
        max_cd5=("cd5", "max"),
        mean_cd5=("cd5", "mean"),
        median_cd5=("cd5", "median"),
        std_cd5=("cd5", "std"),
        n_disruptive_papers=("is_top5pct", "sum"),
        any_top5=("is_top5pct", "max"),
        mean_team_size=("team_size", "mean"),
        mean_cited_by=("cited_by_count", "mean"),
        max_cited_by=("cited_by_count", "max"),
        mean_ref_count=("ref_count", "mean"),
        pub_year_min=("publication_year", "min"),
        pub_year_max=("publication_year", "max"),
        modal_field=("field", lambda x: x.mode().iloc[0] if len(x.mode()) > 0 else "Unknown"),
    ).reset_index()

    agg["frac_disruptive"] = agg["n_disruptive_papers"] / agg["n_scored_papers"]
    agg["career_length"] = agg["pub_year_max"] - agg["pub_year_min"] + 1

    return agg


# ── Step 6: Distribution Analysis ─────────────────────────────────────────────

def lorenz_curve(values: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Compute Lorenz curve coordinates."""
    sorted_vals = np.sort(values)
    n = len(sorted_vals)
    cumulative = np.cumsum(sorted_vals)
    total = cumulative[-1]
    x = np.arange(1, n + 1) / n
    y = cumulative / total if total > 0 else cumulative
    # Prepend origin
    x = np.insert(x, 0, 0)
    y = np.insert(y, 0, 0)
    return x, y


def gini_coefficient(values: np.ndarray) -> float:
    """Compute Gini coefficient."""
    if len(values) == 0:
        return np.nan
    sorted_vals = np.sort(values)
    n = len(sorted_vals)
    index = np.arange(1, n + 1)
    return (2 * np.sum(index * sorted_vals) - (n + 1) * np.sum(sorted_vals)) / (n * np.sum(sorted_vals))


def bootstrap_gini_ci(values: np.ndarray, n_boot: int = 1000, alpha: float = 0.05) -> Tuple[float, float, float]:
    """Bootstrap confidence interval for Gini coefficient."""
    rng = np.random.RandomState(RANDOM_SEED)
    ginis = []
    for _ in range(n_boot):
        sample = rng.choice(values, size=len(values), replace=True)
        g = gini_coefficient(sample)
        if not np.isnan(g):
            ginis.append(g)
    ginis = np.array(ginis)
    point = gini_coefficient(values)
    lo = np.percentile(ginis, 100 * alpha / 2)
    hi = np.percentile(ginis, 100 * (1 - alpha / 2))
    return point, lo, hi


def concentration_ratios(author_metrics: pd.DataFrame, col: str = "n_disruptive_papers") -> Dict[str, float]:
    """What % of disruptive papers come from top X% of authors."""
    sorted_df = author_metrics.sort_values(col, ascending=False).reset_index(drop=True)
    total_disruptive = sorted_df[col].sum()
    if total_disruptive == 0:
        return {"top1pct": np.nan, "top5pct": np.nan, "top10pct": np.nan}
    n = len(sorted_df)
    return {
        "top1pct": sorted_df[col].iloc[:max(1, int(0.01 * n))].sum() / total_disruptive,
        "top5pct": sorted_df[col].iloc[:max(1, int(0.05 * n))].sum() / total_disruptive,
        "top10pct": sorted_df[col].iloc[:max(1, int(0.10 * n))].sum() / total_disruptive,
    }


def fit_power_law_tail(values: np.ndarray, xmin_percentile: float = 0.90) -> Dict[str, Any]:
    """
    Fit power-law tail using MLE (Clauset-Shalizi-Newman method, simplified).
    Returns estimated alpha and KS statistic.
    """
    vals = values[values > 0]
    if len(vals) < 50:
        return {"alpha": np.nan, "xmin": np.nan, "ks_stat": np.nan, "n_tail": 0}

    xmin = np.percentile(vals, 100 * xmin_percentile)
    tail = vals[vals >= xmin]
    n_tail = len(tail)
    if n_tail < 10:
        return {"alpha": np.nan, "xmin": xmin, "ks_stat": np.nan, "n_tail": n_tail}

    # MLE for continuous power law: alpha = 1 + n / sum(ln(x/xmin))
    alpha = 1 + n_tail / np.sum(np.log(tail / xmin))

    # KS statistic
    sorted_tail = np.sort(tail)
    empirical_cdf = np.arange(1, n_tail + 1) / n_tail
    theoretical_cdf = 1 - (sorted_tail / xmin) ** (1 - alpha)
    ks_stat = np.max(np.abs(empirical_cdf - theoretical_cdf))

    return {"alpha": alpha, "xmin": xmin, "ks_stat": ks_stat, "n_tail": n_tail}


def run_distribution_analysis(
    author_metrics: pd.DataFrame, output_dir: Path
) -> Dict[str, Any]:
    """Run full distribution analysis: Lorenz, Gini, concentration, power-law."""
    results = {}

    # ── Disruption concentration ──
    # Use n_disruptive_papers for concentration analysis
    disrupt_vals = author_metrics["n_disruptive_papers"].values.astype(float)

    # Gini on disruptive output
    gini_point, gini_lo, gini_hi = bootstrap_gini_ci(disrupt_vals)
    results["gini_disruption"] = {"point": gini_point, "ci_lo": gini_lo, "ci_hi": gini_hi}
    print(f"  Gini (disruptive papers): {gini_point:.4f} [{gini_lo:.4f}, {gini_hi:.4f}]", flush=True)

    # Concentration ratios
    conc = concentration_ratios(author_metrics, "n_disruptive_papers")
    results["concentration_disruption"] = conc
    print(f"  Top 1% produce {100*conc['top1pct']:.1f}% of disruptive papers", flush=True)
    print(f"  Top 5% produce {100*conc['top5pct']:.1f}% of disruptive papers", flush=True)
    print(f"  Top 10% produce {100*conc['top10pct']:.1f}% of disruptive papers", flush=True)

    # ── Citation concentration (benchmark) ──
    cite_vals = author_metrics["max_cited_by"].values.astype(float)
    gini_cite, gini_cite_lo, gini_cite_hi = bootstrap_gini_ci(cite_vals)
    results["gini_citations"] = {"point": gini_cite, "ci_lo": gini_cite_lo, "ci_hi": gini_cite_hi}
    print(f"  Gini (max citations): {gini_cite:.4f} [{gini_cite_lo:.4f}, {gini_cite_hi:.4f}]", flush=True)

    # ── Power-law tail fit ──
    max_cd5_vals = author_metrics["max_cd5"].dropna().values
    # Shift to positive for power-law fitting
    shifted = max_cd5_vals - max_cd5_vals.min() + 0.001
    pl_fit = fit_power_law_tail(shifted)
    results["power_law_fit"] = pl_fit

    # ── Lorenz curve plot ──
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # Disruption Lorenz
    lx, ly = lorenz_curve(disrupt_vals)
    axes[0].plot(lx, ly, 'b-', linewidth=2, label=f'Disruption (Gini={gini_point:.3f})')
    axes[0].plot([0, 1], [0, 1], 'k--', alpha=0.5, label='Perfect equality')
    axes[0].set_xlabel("Cumulative % of scientists")
    axes[0].set_ylabel("Cumulative % of disruptive papers")
    axes[0].set_title("Lorenz Curve: Disruptive Research Output")
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)

    # Citation Lorenz for comparison
    cx, cy = lorenz_curve(cite_vals)
    axes[1].plot(cx, cy, 'r-', linewidth=2, label=f'Citations (Gini={gini_cite:.3f})')
    axes[1].plot(lx, ly, 'b--', linewidth=1.5, alpha=0.7, label=f'Disruption (Gini={gini_point:.3f})')
    axes[1].plot([0, 1], [0, 1], 'k--', alpha=0.5, label='Perfect equality')
    axes[1].set_xlabel("Cumulative % of scientists")
    axes[1].set_ylabel("Cumulative share")
    axes[1].set_title("Lorenz Comparison: Disruption vs. Citations")
    axes[1].legend()
    axes[1].grid(True, alpha=0.3)

    plt.tight_layout()
    fig.savefig(output_dir / "lorenz_curves.png", dpi=150)
    plt.close(fig)

    # ── Disruption distribution histogram ──
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    axes[0].hist(author_metrics["max_cd5"].dropna(), bins=80, color="steelblue", edgecolor="white", alpha=0.8)
    axes[0].set_xlabel("Max CD5 Score")
    axes[0].set_ylabel("Number of Authors")
    axes[0].set_title("Distribution of Peak Disruption (max CD5)")
    axes[0].axvline(x=0, color='red', linestyle='--', alpha=0.5)

    axes[1].hist(author_metrics["n_disruptive_papers"], bins=range(0, int(author_metrics["n_disruptive_papers"].max()) + 2),
                 color="coral", edgecolor="white", alpha=0.8)
    axes[1].set_xlabel("Number of Top-5% Disruptive Papers")
    axes[1].set_ylabel("Number of Authors")
    axes[1].set_title("Distribution of Disruptive Paper Count")
    axes[1].set_yscale("log")

    plt.tight_layout()
    fig.savefig(output_dir / "disruption_distribution.png", dpi=150)
    plt.close(fig)

    return results


# ── Step 7: Demographic Profiling ─────────────────────────────────────────────

def run_demographic_analysis(
    author_metrics: pd.DataFrame, author_demo: pd.DataFrame, output_dir: Path
) -> Dict[str, Any]:
    """Descriptive comparison and regression models for disruptor profiling."""
    # Merge metrics with demographics
    df = author_metrics.merge(author_demo, on="author_id", how="inner")

    if df.empty:
        return {"error": "No data after merge"}

    results = {}

    # ── Descriptive comparison ──
    df["is_disruptor"] = df["any_top5"].astype(int)

    desc_cols = ["mean_team_size", "n_scored_papers", "career_length", "h_index",
                 "mean_cd5", "max_cd5", "mean_cited_by"]
    available_cols = [c for c in desc_cols if c in df.columns]

    desc_table = df.groupby("is_disruptor")[available_cols].agg(["mean", "median", "std"])
    # Flatten multi-level columns for JSON serialization
    results["descriptive_table"] = {
        f"{col}_{stat}": desc_table[(col, stat)].to_dict()
        for col, stat in desc_table.columns
    }

    print("  Disruptor vs Non-disruptor comparison:", flush=True)
    for col in available_cols:
        d_mean = df.loc[df["is_disruptor"] == 1, col].mean()
        nd_mean = df.loc[df["is_disruptor"] == 0, col].mean()
        print(f"    {col}: disruptors={d_mean:.2f}, non-disruptors={nd_mean:.2f}", flush=True)

    # ── Logistic regression: Pr(disruptor) ──
    model_df = df.dropna(subset=["mean_team_size", "n_scored_papers", "career_length"]).copy()
    if len(model_df) < 50:
        results["logistic_model"] = {"error": f"Too few observations ({len(model_df)})"}
        results["ols_model"] = {"error": f"Too few observations ({len(model_df)})"}
        return results

    model_df["papers_per_year"] = model_df["n_scored_papers"] / model_df["career_length"].clip(lower=1)
    model_df["log_h_index"] = np.log1p(model_df.get("h_index", 0))

    # Build X matrix
    x_cols = ["mean_team_size", "papers_per_year", "career_length"]
    if "h_index" in model_df.columns and model_df["h_index"].notna().sum() > 50:
        x_cols.append("log_h_index")

    # Add field fixed effects (top fields only)
    if "modal_field" in model_df.columns:
        top_fields = model_df["modal_field"].value_counts().head(8).index.tolist()
        for f in top_fields[1:]:  # skip most common (reference)
            model_df[f"field_{f}"] = (model_df["modal_field"] == f).astype(int)
            x_cols.append(f"field_{f}")

    X = model_df[x_cols].copy()
    X = X.fillna(X.median())
    X = sm.add_constant(X, has_constant="add")

    try:
        y_logit = model_df["is_disruptor"].astype(float)
        logit_res = sm.Logit(y_logit, X).fit(disp=False, maxiter=100)
        results["logistic_model"] = {
            "n": int(len(model_df)),
            "pseudo_r2": float(logit_res.prsquared),
            "coefficients": {k: {"coef": float(v), "p": float(logit_res.pvalues[k]),
                                  "odds_ratio": float(np.exp(v))}
                             for k, v in logit_res.params.items()},
        }
        print(f"  Logistic model: n={len(model_df)}, pseudo-R2={logit_res.prsquared:.4f}", flush=True)
    except Exception as e:
        results["logistic_model"] = {"error": str(e)}
        print(f"  Logistic model failed: {e}", flush=True)

    # ── OLS: max_cd5 ~ covariates ──
    try:
        y_ols = model_df["max_cd5"].astype(float)
        ols_res = sm.OLS(y_ols, X).fit(cov_type="HC3")
        results["ols_model"] = {
            "n": int(len(model_df)),
            "r2": float(ols_res.rsquared),
            "adj_r2": float(ols_res.rsquared_adj),
            "coefficients": {k: {"coef": float(v), "p": float(ols_res.pvalues[k])}
                             for k, v in ols_res.params.items()},
        }
        print(f"  OLS model: n={len(model_df)}, R2={ols_res.rsquared:.4f}", flush=True)
    except Exception as e:
        results["ols_model"] = {"error": str(e)}
        print(f"  OLS model failed: {e}", flush=True)

    return results


# ── Step 8: Field-Level Heterogeneity ─────────────────────────────────────────

def field_heterogeneity(author_metrics: pd.DataFrame, output_dir: Path) -> Dict[str, Any]:
    """Compute Gini coefficients by STEM subfield."""
    results = {}
    field_ginis = []

    for field, group in author_metrics.groupby("modal_field"):
        if len(group) < 30:
            continue
        vals = group["n_disruptive_papers"].values.astype(float)
        g = gini_coefficient(vals)
        field_ginis.append({"field": field, "gini": g, "n_authors": len(group),
                           "mean_max_cd5": group["max_cd5"].mean(),
                           "pct_disruptors": group["any_top5"].mean()})

    field_df = pd.DataFrame(field_ginis).sort_values("gini", ascending=False)
    results["field_ginis"] = field_df.to_dict("records")

    if not field_df.empty:
        print("  Field-level Gini coefficients:", flush=True)
        for _, row in field_df.head(15).iterrows():
            print(f"    {row['field']}: Gini={row['gini']:.4f} (n={row['n_authors']}, "
                  f"{100*row['pct_disruptors']:.1f}% disruptors)", flush=True)

        # Plot
        plot_df = field_df.head(15).sort_values("gini")
        fig, ax = plt.subplots(figsize=(10, 6))
        bars = ax.barh(plot_df["field"], plot_df["gini"], color="steelblue", edgecolor="white")
        ax.set_xlabel("Gini Coefficient (Disruptive Papers)")
        ax.set_title("Concentration of Disruptive Research by STEM Field")
        ax.axvline(x=field_df["gini"].median(), color="red", linestyle="--", alpha=0.5, label="Median")
        ax.legend()
        plt.tight_layout()
        fig.savefig(output_dir / "field_gini_comparison.png", dpi=150)
        plt.close(fig)

    return results


# ── Step 9: Temporal Trends ───────────────────────────────────────────────────

def temporal_trends(works_df: pd.DataFrame, output_dir: Path) -> Dict[str, Any]:
    """Compute Gini in 5-year windows."""
    scored = works_df.dropna(subset=["cd5"]).copy()
    if scored.empty:
        return {"error": "No scored works"}

    cd5_threshold = scored["cd5"].quantile(1 - TOP_PERCENTILE)

    windows = []
    for start in range(1980, 2011, 5):
        end = start + 4
        window = scored[(scored["publication_year"] >= start) & (scored["publication_year"] <= end)]
        if len(window) < 50:
            continue

        window_top = (window["cd5"] >= cd5_threshold).astype(int)

        # Author-level aggregation within window
        author_disrupt = window.groupby("author_id").agg(
            n_papers=("work_id", "count"),
            n_disruptive=("cd5", lambda x: (x >= cd5_threshold).sum()),
            mean_cd5=("cd5", "mean"),
        ).reset_index()

        g = gini_coefficient(author_disrupt["n_disruptive"].values.astype(float))
        frac_disruptors = (author_disrupt["n_disruptive"] > 0).mean()

        windows.append({
            "window": f"{start}-{end}",
            "start_year": start,
            "n_papers": len(window),
            "n_authors": len(author_disrupt),
            "gini": g,
            "frac_disruptors": frac_disruptors,
            "mean_cd5": window["cd5"].mean(),
        })

    results = {"windows": windows}

    if windows:
        wdf = pd.DataFrame(windows)
        print("  Temporal trends:", flush=True)
        for _, row in wdf.iterrows():
            print(f"    {row['window']}: Gini={row['gini']:.4f}, "
                  f"disruptor%={100*row['frac_disruptors']:.1f}%, n={row['n_authors']}", flush=True)

        # Plot
        fig, axes = plt.subplots(1, 2, figsize=(14, 5))

        axes[0].plot(wdf["start_year"], wdf["gini"], "bo-", linewidth=2, markersize=8)
        axes[0].set_xlabel("Window Start Year")
        axes[0].set_ylabel("Gini Coefficient")
        axes[0].set_title("Concentration of Disruption Over Time")
        axes[0].grid(True, alpha=0.3)

        axes[1].plot(wdf["start_year"], wdf["frac_disruptors"] * 100, "rs-", linewidth=2, markersize=8)
        axes[1].set_xlabel("Window Start Year")
        axes[1].set_ylabel("% of Authors with Disruptive Papers")
        axes[1].set_title("Disruptor Fraction Over Time")
        axes[1].grid(True, alpha=0.3)

        plt.tight_layout()
        fig.savefig(output_dir / "temporal_trends.png", dpi=150)
        plt.close(fig)

    return results


# ── Step 10: Validation ───────────────────────────────────────────────────────

def compute_ab_disruption(client: CurlCacheClient, work: Dict[str, Any]) -> Dict[str, Any]:
    """Compute A/B disruption index using OpenAlex citation data."""
    focal_id = work["work_id"]
    pub_year = int(work["publication_year"])
    refs = set(work.get("referenced_works", []))
    if not refs:
        return {"a_count": 0, "b_count": 0, "d_ab10": np.nan, "window_citers": 0}

    to_year = pub_year + FORWARD_WINDOW_YEARS
    citers = fetch_openalex_pages(
        client,
        OPENALEX_WORK_URL,
        {
            "filter": f"cites:{focal_id},from_publication_date:{pub_year}-01-01,to_publication_date:{to_year}-12-31",
            "select": "id,referenced_works",
        },
        max_items=MAX_CITERS_PER_WORK,
    )

    a_count = 0
    b_count = 0
    for citer in citers:
        citer_refs = {short_openalex_id(x) for x in (citer.get("referenced_works", []) or [])}
        if refs & citer_refs:
            b_count += 1
        else:
            a_count += 1
    total = a_count + b_count
    d_ab10 = np.nan
    if total >= MIN_CITERS_FOR_D:
        d_ab10 = (a_count - b_count) / total
    return {"a_count": a_count, "b_count": b_count, "d_ab10": d_ab10, "window_citers": total}


def run_validation_benchmark(
    works_df: pd.DataFrame, zenodo_benchmark: Optional[pd.DataFrame], output_dir: Path
) -> Dict[str, Any]:
    """Validate our disruption distribution against Park et al. population benchmark.

    Since Zenodo data has no per-paper IDs, we compare distributional properties:
    - Summary stats (mean, median, std, skew)
    - Percentile comparison
    - KS test for distributional similarity
    """
    scored = works_df.dropna(subset=["cd5"]).copy()
    if scored.empty:
        return {"error": "No scored works for validation"}

    results = {
        "n_our_papers": len(scored),
        "our_mean": float(scored["cd5"].mean()),
        "our_median": float(scored["cd5"].median()),
        "our_std": float(scored["cd5"].std()),
        "our_skew": float(scored["cd5"].skew()),
        "our_pct_positive": float((scored["cd5"] > 0).mean()),
        "our_pct_strongly_disruptive": float((scored["cd5"] > 0.5).mean()),
    }

    print(f"  Our sample: n={len(scored)}, mean={results['our_mean']:.4f}, "
          f"median={results['our_median']:.4f}, std={results['our_std']:.4f}", flush=True)
    print(f"  % positive disruption: {100*results['our_pct_positive']:.1f}%", flush=True)
    print(f"  % strongly disruptive (>0.5): {100*results['our_pct_strongly_disruptive']:.1f}%", flush=True)

    if zenodo_benchmark is not None and not zenodo_benchmark.empty:
        from scipy import stats

        zen = zenodo_benchmark["cd5"].dropna()
        results["zenodo_n"] = len(zen)
        results["zenodo_mean"] = float(zen.mean())
        results["zenodo_median"] = float(zen.median())
        results["zenodo_std"] = float(zen.std())
        results["zenodo_pct_positive"] = float((zen > 0).mean())

        # KS test comparing distributions
        # Subsample zenodo for tractability
        zen_sample = zen.sample(n=min(50000, len(zen)), random_state=RANDOM_SEED)
        ks_stat, ks_p = stats.ks_2samp(scored["cd5"].values, zen_sample.values)
        results["ks_stat"] = float(ks_stat)
        results["ks_p"] = float(ks_p)

        print(f"  Zenodo benchmark: n={len(zen):,}, mean={zen.mean():.4f}, "
              f"median={zen.median():.4f}", flush=True)
        print(f"  KS test: stat={ks_stat:.4f}, p={ks_p:.4g}", flush=True)

        # Distribution comparison plot
        fig, axes = plt.subplots(1, 2, figsize=(14, 5))

        # Overlapping histograms
        bins = np.linspace(-1, 1, 80)
        axes[0].hist(scored["cd5"], bins=bins, density=True, alpha=0.6,
                     color="steelblue", label=f"Our sample (n={len(scored):,})")
        axes[0].hist(zen_sample, bins=bins, density=True, alpha=0.4,
                     color="coral", label=f"Zenodo benchmark (n={len(zen_sample):,})")
        axes[0].set_xlabel("Disruption Score")
        axes[0].set_ylabel("Density")
        axes[0].set_title("Distribution Comparison: Our Sample vs. Zenodo")
        axes[0].legend()
        axes[0].grid(True, alpha=0.3)

        # QQ plot
        our_quantiles = np.percentile(scored["cd5"], np.arange(1, 100))
        zen_quantiles = np.percentile(zen_sample, np.arange(1, 100))
        axes[1].scatter(zen_quantiles, our_quantiles, s=15, alpha=0.7)
        lims = [min(zen_quantiles.min(), our_quantiles.min()),
                max(zen_quantiles.max(), our_quantiles.max())]
        axes[1].plot(lims, lims, 'r--', alpha=0.5)
        axes[1].set_xlabel("Zenodo Benchmark Quantiles")
        axes[1].set_ylabel("Our Sample Quantiles")
        axes[1].set_title("Q-Q Plot: Disruption Scores")
        axes[1].grid(True, alpha=0.3)

        plt.tight_layout()
        fig.savefig(output_dir / "validation_benchmark.png", dpi=150)
        plt.close(fig)
    else:
        print("  No Zenodo benchmark available for comparison.", flush=True)
        results["zenodo_available"] = False

    return results


# ── Step 11: Report ───────────────────────────────────────────────────────────

def write_report(
    author_demo: pd.DataFrame,
    author_metrics: pd.DataFrame,
    works_df: pd.DataFrame,
    dist_results: Dict[str, Any],
    demo_results: Dict[str, Any],
    field_results: Dict[str, Any],
    temporal_results: Dict[str, Any],
    validation_results: Dict[str, Any],
    output_dir: Path,
) -> None:
    """Generate markdown report."""
    lines: List[str] = []
    lines.append("# Concentration of Disruptive Research Among Scientists")
    lines.append("")
    lines.append("## Hypothesis")
    lines.append("Only a very small percentage of scientists produce very disruptive research.")
    lines.append("")

    # ── Data Summary ──
    lines.append("## Data Summary")
    n_authors = len(author_metrics)
    n_works = len(works_df)
    scored_works = works_df["cd5"].notna().sum() if "cd5" in works_df.columns else 0
    lines.append(f"- **Authors sampled**: {n_authors:,}")
    lines.append(f"- **Total works fetched**: {n_works:,}")
    lines.append(f"- **Works with CD5 scores**: {scored_works:,}")
    lines.append(f"- **Publication window**: {PUB_START_YEAR}–{PUB_END_YEAR}")
    lines.append(f"- **Disruption threshold**: Top {100*TOP_PERCENTILE:.0f}% of CD5 distribution")
    lines.append("")

    # ── Distribution Results ──
    lines.append("## Distribution Analysis")
    lines.append("")

    gini_d = dist_results.get("gini_disruption", {})
    if gini_d:
        lines.append(f"### Gini Coefficient (Disruptive Papers)")
        lines.append(f"- **Point estimate**: {gini_d.get('point', 'N/A'):.4f}")
        lines.append(f"- **95% CI**: [{gini_d.get('ci_lo', 'N/A'):.4f}, {gini_d.get('ci_hi', 'N/A'):.4f}]")
        lines.append("")

    gini_c = dist_results.get("gini_citations", {})
    if gini_c:
        lines.append(f"### Gini Coefficient (Citations, benchmark)")
        lines.append(f"- **Point estimate**: {gini_c.get('point', 'N/A'):.4f}")
        lines.append(f"- **95% CI**: [{gini_c.get('ci_lo', 'N/A'):.4f}, {gini_c.get('ci_hi', 'N/A'):.4f}]")
        lines.append("")

    conc = dist_results.get("concentration_disruption", {})
    if conc:
        lines.append("### Concentration Ratios")
        lines.append(f"- Top 1% of scientists produce **{100*conc.get('top1pct', 0):.1f}%** of disruptive papers")
        lines.append(f"- Top 5% of scientists produce **{100*conc.get('top5pct', 0):.1f}%** of disruptive papers")
        lines.append(f"- Top 10% of scientists produce **{100*conc.get('top10pct', 0):.1f}%** of disruptive papers")
        lines.append("")

    pl = dist_results.get("power_law_fit", {})
    if pl and not np.isnan(pl.get("alpha", np.nan)):
        lines.append("### Power-Law Tail Fit")
        lines.append(f"- Alpha: {pl['alpha']:.3f}")
        lines.append(f"- x_min: {pl['xmin']:.4f}")
        lines.append(f"- KS statistic: {pl['ks_stat']:.4f}")
        lines.append(f"- Tail observations: {pl['n_tail']}")
        lines.append("")

    lines.append("### Figures")
    lines.append("- ![Lorenz Curves](lorenz_curves.png)")
    lines.append("- ![Disruption Distribution](disruption_distribution.png)")
    lines.append("")

    # ── Demographic Profiling ──
    lines.append("## Demographic Profiling")
    lines.append("")

    logit = demo_results.get("logistic_model", {})
    if "error" not in logit:
        lines.append("### Logistic Regression: Pr(Disruptor)")
        lines.append(f"- N = {logit.get('n', 'N/A')}")
        lines.append(f"- Pseudo-R² = {logit.get('pseudo_r2', 0):.4f}")
        lines.append("")
        lines.append("| Variable | Coefficient | Odds Ratio | p-value |")
        lines.append("|----------|------------|------------|---------|")
        for var, vals in logit.get("coefficients", {}).items():
            if var == "const":
                continue
            lines.append(f"| {var} | {vals['coef']:.4f} | {vals.get('odds_ratio', 'N/A'):.3f} | {vals['p']:.4g} |")
        lines.append("")
    else:
        lines.append(f"Logistic model: {logit.get('error', 'failed')}")
        lines.append("")

    ols = demo_results.get("ols_model", {})
    if "error" not in ols:
        lines.append("### OLS Regression: max_cd5")
        lines.append(f"- N = {ols.get('n', 'N/A')}")
        lines.append(f"- R² = {ols.get('r2', 0):.4f}")
        lines.append(f"- Adjusted R² = {ols.get('adj_r2', 0):.4f}")
        lines.append("")
        lines.append("| Variable | Coefficient | p-value |")
        lines.append("|----------|------------|---------|")
        for var, vals in ols.get("coefficients", {}).items():
            if var == "const":
                continue
            lines.append(f"| {var} | {vals['coef']:.4f} | {vals['p']:.4g} |")
        lines.append("")

    # ── Field Heterogeneity ──
    lines.append("## Field-Level Heterogeneity")
    lines.append("")
    field_ginis = field_results.get("field_ginis", [])
    if field_ginis:
        lines.append("| Field | Gini | N Authors | % Disruptors |")
        lines.append("|-------|------|-----------|-------------|")
        for fg in field_ginis[:15]:
            lines.append(f"| {fg['field']} | {fg['gini']:.4f} | {fg['n_authors']} | {100*fg['pct_disruptors']:.1f}% |")
        lines.append("")
        lines.append("![Field Gini Comparison](field_gini_comparison.png)")
        lines.append("")

    # ── Temporal Trends ──
    lines.append("## Temporal Trends")
    lines.append("")
    windows = temporal_results.get("windows", [])
    if windows:
        lines.append("| Window | Gini | % Disruptors | N Authors | N Papers |")
        lines.append("|--------|------|-------------|-----------|----------|")
        for w in windows:
            lines.append(f"| {w['window']} | {w['gini']:.4f} | {100*w['frac_disruptors']:.1f}% | {w['n_authors']} | {w['n_papers']} |")
        lines.append("")
        lines.append("![Temporal Trends](temporal_trends.png)")
        lines.append("")

    # ── Validation ──
    lines.append("## Validation (Distribution Comparison)")
    lines.append("")
    if "error" not in validation_results:
        lines.append(f"### Our Sample")
        lines.append(f"- N papers scored: {validation_results.get('n_our_papers', 'N/A')}")
        lines.append(f"- Mean disruption: {validation_results.get('our_mean', 0):.4f}")
        lines.append(f"- Median disruption: {validation_results.get('our_median', 0):.4f}")
        lines.append(f"- % positive disruption: {100*validation_results.get('our_pct_positive', 0):.1f}%")
        lines.append("")
        if "zenodo_n" in validation_results:
            lines.append(f"### Zenodo Benchmark (Park et al. 2023)")
            lines.append(f"- N records: {validation_results.get('zenodo_n', 'N/A'):,}")
            lines.append(f"- Mean disruption: {validation_results.get('zenodo_mean', 0):.4f}")
            lines.append(f"- Median disruption: {validation_results.get('zenodo_median', 0):.4f}")
            lines.append(f"- KS statistic: {validation_results.get('ks_stat', 'N/A'):.4f} "
                         f"(p={validation_results.get('ks_p', 'N/A'):.4g})")
            lines.append("")
            lines.append("![Validation Benchmark](validation_benchmark.png)")
    else:
        lines.append(f"Validation: {validation_results.get('error', 'not run')}")
    lines.append("")

    # ── Conclusion ──
    lines.append("## Conclusion")
    lines.append("")
    if gini_d and gini_d.get("point", 0) > 0.7:
        lines.append("The evidence strongly supports the hypothesis that disruptive research is highly "
                     "concentrated among a small fraction of scientists. The Gini coefficient for disruptive "
                     f"output ({gini_d['point']:.3f}) indicates extreme inequality in the production of "
                     "paradigm-shifting work.")
    elif gini_d and gini_d.get("point", 0) > 0.5:
        lines.append("The evidence supports moderate concentration of disruptive research output. "
                     f"The Gini coefficient ({gini_d['point']:.3f}) indicates substantial inequality, "
                     "though less extreme than pure citation-based measures might suggest.")
    else:
        lines.append("The distribution of disruptive research shows less concentration than hypothesized. "
                     "Further analysis with larger samples or alternative disruption measures may be warranted.")
    lines.append("")

    lines.append("## Notes")
    lines.append("- Disruption measured via A/B disruption index computed from OpenAlex citation data (10-year forward window).")
    lines.append("- Park et al. (2023) Zenodo dataset used as population-level benchmark (no per-paper IDs available for matching).")
    lines.append("- Author demographics from OpenAlex API; STEM filtered client-side via top concepts.")
    lines.append("- STEM fields only (Physical Sciences, Life Sciences, Engineering, Mathematics, CS).")
    lines.append(f"- Publication window: {PUB_START_YEAR}–{PUB_END_YEAR} to ensure adequate forward citation window.")

    report_path = output_dir / "disruption_concentration_report.md"
    report_path.write_text("\n".join(lines))
    print(f"  Report written to {report_path}", flush=True)


# ── Fallback: Pure OpenAlex Approach ──────────────────────────────────────────

def fallback_compute_disruption_openalex(
    client: CurlCacheClient, works_df: pd.DataFrame, max_papers: int = 5000
) -> pd.DataFrame:
    """Compute A/B disruption for papers when Zenodo data is unavailable."""
    print(f"  FALLBACK: Computing A/B disruption via OpenAlex for up to {max_papers} papers...", flush=True)

    eligible = works_df[works_df["ref_count"] > 0].copy()
    if len(eligible) > max_papers:
        eligible = eligible.sample(n=max_papers, random_state=RANDOM_SEED)

    results = []
    for i, (_, row) in enumerate(eligible.iterrows(), start=1):
        work = {
            "work_id": row["work_id"],
            "publication_year": row["publication_year"],
            "referenced_works": row.get("referenced_works", []),
        }
        if isinstance(work["referenced_works"], str):
            try:
                work["referenced_works"] = json.loads(work["referenced_works"])
            except (json.JSONDecodeError, TypeError):
                work["referenced_works"] = []

        metric = compute_ab_disruption(client, work)
        results.append({"work_id": row["work_id"], "cd5": metric["d_ab10"]})

        if i % 100 == 0:
            valid = sum(1 for r in results if not np.isnan(r.get("cd5", np.nan)))
            print(f"    Computed {i}/{len(eligible)} ({valid} valid)...", flush=True)

    result_df = pd.DataFrame(results)
    merged = works_df.merge(result_df, on="work_id", how="left", suffixes=("_old", ""))
    if "cd5_old" in merged.columns:
        merged.drop(columns=["cd5_old"], inplace=True)
    return merged


# ── Main Pipeline ─────────────────────────────────────────────────────────────

def main() -> None:
    random.seed(RANDOM_SEED)
    np.random.seed(RANDOM_SEED)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    DATA_DIR.mkdir(parents=True, exist_ok=True)

    client = CurlCacheClient(CACHE_DIR)

    # ── Step 1: Download and load Zenodo population benchmark ──
    print("=" * 70, flush=True)
    print("Step 1/11: Download and load Zenodo population benchmark", flush=True)
    print("=" * 70, flush=True)

    zenodo_benchmark = None
    try:
        archive_path = download_zenodo(DATA_DIR)
        extract_zenodo(archive_path, DATA_DIR)
        zenodo_benchmark = load_zenodo_population_benchmark(DATA_DIR)
    except Exception as e:
        print(f"  Zenodo download/load failed: {e}", flush=True)

    if zenodo_benchmark is not None:
        zenodo_benchmark.head(10000).to_csv(OUTPUT_DIR / "zenodo_cd5_benchmark_sample.csv", index=False)
        print(f"  Zenodo benchmark loaded: {len(zenodo_benchmark):,} records.", flush=True)
        print("  NOTE: No per-paper IDs in Zenodo data. Using OpenAlex A/B disruption for paper-level scores.", flush=True)
    else:
        print("  Zenodo benchmark not available. Will proceed without population benchmark.", flush=True)

    # ── Step 2: Sample STEM authors ──
    print("=" * 70, flush=True)
    print("Step 2/11: Sample STEM authors from OpenAlex", flush=True)
    print("=" * 70, flush=True)

    # Check for cached author data
    author_cache_path = OUTPUT_DIR / "sampled_authors.csv"
    if author_cache_path.exists():
        print("  Loading cached author sample...", flush=True)
        author_demo = pd.read_csv(author_cache_path)
        print(f"  Loaded {len(author_demo):,} authors from cache.", flush=True)
    else:
        author_demo = sample_stem_authors(client, TARGET_AUTHORS)
        author_demo.to_csv(author_cache_path, index=False)

    # ── Step 3: Fetch author works ──
    print("=" * 70, flush=True)
    print("Step 3/11: Fetch author works (1980-2010)", flush=True)
    print("=" * 70, flush=True)

    works_cache_path = OUTPUT_DIR / "all_works.csv"
    if works_cache_path.exists():
        print("  Loading cached works...", flush=True)
        works_df = pd.read_csv(works_cache_path)
        print(f"  Loaded {len(works_df):,} works from cache.", flush=True)
    else:
        all_works: List[Dict[str, Any]] = []
        total_authors = len(author_demo)
        for i, (_, author) in enumerate(author_demo.iterrows(), start=1):
            aid = author["author_id"]
            works = fetch_author_works_for_study(client, aid)
            for w in works:
                w["author_id"] = aid
                # Store referenced_works as JSON string for CSV serialization
                if "referenced_works" in w:
                    w["referenced_works_json"] = json.dumps(w.pop("referenced_works"))
            all_works.extend(works)
            if i % 500 == 0 or i == total_authors:
                print(f"  Fetched works for {i}/{total_authors} authors ({len(all_works):,} works total)", flush=True)

        works_df = pd.DataFrame(all_works)
        if not works_df.empty:
            works_df.to_csv(works_cache_path, index=False)
        print(f"  Total works: {len(works_df):,}", flush=True)

    if works_df.empty:
        print("FATAL: No works fetched. Exiting.", flush=True)
        sys.exit(1)

    # ── Step 4: Compute A/B disruption scores ──
    print("=" * 70, flush=True)
    print("Step 4/11: Compute A/B disruption scores via OpenAlex", flush=True)
    print("=" * 70, flush=True)

    scored_cache_path = OUTPUT_DIR / "works_with_cd5.csv"
    if scored_cache_path.exists():
        print("  Loading cached disruption scores...", flush=True)
        works_df = pd.read_csv(scored_cache_path)
        print(f"  Loaded {works_df['cd5'].notna().sum():,} scored works from cache.", flush=True)
    else:
        # Zenodo data has no per-paper IDs, so we always compute A/B disruption
        # Reconstruct referenced_works from JSON if needed
        if "referenced_works_json" in works_df.columns and "referenced_works" not in works_df.columns:
            works_df["referenced_works"] = works_df["referenced_works_json"].apply(
                lambda x: json.loads(x) if pd.notna(x) else [])
        works_df = fallback_compute_disruption_openalex(client, works_df, max_papers=5000)
        works_df.to_csv(scored_cache_path, index=False)

    # ── Step 5: Author-level disruption metrics ──
    print("=" * 70, flush=True)
    print("Step 5/11: Compute author-level disruption metrics", flush=True)
    print("=" * 70, flush=True)

    author_metrics = compute_author_disruption_metrics(works_df, TOP_PERCENTILE)
    if author_metrics.empty:
        print("FATAL: No author-level metrics computed. Exiting.", flush=True)
        sys.exit(1)

    author_metrics.to_csv(OUTPUT_DIR / "author_disruption_metrics.csv", index=False)
    print(f"  Authors with metrics: {len(author_metrics):,}", flush=True)
    print(f"  Disruptors (any top-5%): {author_metrics['any_top5'].sum():,} "
          f"({100*author_metrics['any_top5'].mean():.1f}%)", flush=True)

    # ── Step 6: Distribution analysis ──
    print("=" * 70, flush=True)
    print("Step 6/11: Distribution analysis", flush=True)
    print("=" * 70, flush=True)

    dist_results = run_distribution_analysis(author_metrics, OUTPUT_DIR)

    # ── Step 7: Demographic profiling ──
    print("=" * 70, flush=True)
    print("Step 7/11: Demographic profiling", flush=True)
    print("=" * 70, flush=True)

    demo_results = run_demographic_analysis(author_metrics, author_demo, OUTPUT_DIR)

    # ── Step 8: Field heterogeneity ──
    print("=" * 70, flush=True)
    print("Step 8/11: Field-level heterogeneity", flush=True)
    print("=" * 70, flush=True)

    field_results = field_heterogeneity(author_metrics, OUTPUT_DIR)

    # ── Step 9: Temporal trends ──
    print("=" * 70, flush=True)
    print("Step 9/11: Temporal trends", flush=True)
    print("=" * 70, flush=True)

    temporal_results = temporal_trends(works_df, OUTPUT_DIR)

    # ── Step 10: Validation ──
    print("=" * 70, flush=True)
    print("Step 10/11: Validation (distribution comparison with Zenodo benchmark)", flush=True)
    print("=" * 70, flush=True)

    validation_results = run_validation_benchmark(works_df, zenodo_benchmark, OUTPUT_DIR)

    # ── Step 11: Report ──
    print("=" * 70, flush=True)
    print("Step 11/11: Writing report and figures", flush=True)
    print("=" * 70, flush=True)

    # Save all results as JSON
    all_results = {
        "distribution": dist_results,
        "demographics": demo_results,
        "field_heterogeneity": field_results,
        "temporal": temporal_results,
        "validation": validation_results,
    }
    (OUTPUT_DIR / "all_results.json").write_text(json.dumps(all_results, indent=2, default=str))

    write_report(
        author_demo, author_metrics, works_df,
        dist_results, demo_results, field_results,
        temporal_results, validation_results,
        OUTPUT_DIR,
    )

    print("=" * 70, flush=True)
    print("DONE. All outputs in:", OUTPUT_DIR, flush=True)
    print("=" * 70, flush=True)


if __name__ == "__main__":
    main()
