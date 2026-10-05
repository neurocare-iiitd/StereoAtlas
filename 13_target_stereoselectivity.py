#!/usr/bin/env python3
"""
12_target_stereoselectivity.py
========================================
Target-Centric Stereoselectivity Analysis (TSI Scoring)
StereoAtlas Pipeline — Script 12
Steps 1–11 as specified in the user objective.

INPUT FILES (existing pipeline outputs, not regenerated):
  - outputs/true_isomers/TRUE_ISOMER_PAIRS.csv
  - outputs/stereo_selectivity/family_views/FAMILY_MINMAX_FOLD.csv
  - outputs/stereo_selectivity/family_views/FAMILY_REFERENCE_FOLD.csv
  - outputs/literature_enrichment/TRUE_ISOMER_LITERATURE_CLEAN.csv

OUTPUT FILES (new directory: outputs/target_stereoselectivity/):
  - TARGET_STEREOSELECTIVITY_SUMMARY.csv
  - TARGET_PRIORITY_RANKING.csv
  - TARGET_CONSISTENCY_ANALYSIS.csv
  - TOP_STEREOSELECTIVE_TARGETS.csv  (top 25 by TSI)
  - TARGET_STEREOSELECTIVITY_REPORT.md

TSI FORMULA (documented):
  raw_score = W_FAMILIES   * norm(log10(n_families + 1))
            + W_PAIRS      * norm(log10(n_pairs + 1))
            + W_MEDIAN     * norm(log10(median_fold))
            + W_STRONG_PROP* proportion_strong_extreme (0-100)
            - W_OUTLIER_PEN* outlier_penalty (0 or 100)

  Weights: W_FAMILIES=0.20, W_PAIRS=0.10, W_MEDIAN=0.30,
           W_STRONG_PROP=0.25, W_OUTLIER_PEN=0.15

  outlier_penalty = 1 (scaled to 100) when:
    - At least 1 Extreme-class family exists
    - AND the number of families with fold>=10 equals n_Extreme
    - AND total families > 3
    (meaning ONLY extreme families cross the fold>=10 threshold -> outlier-driven)

  TSI = normalize_0_100(raw_score)
"""

import os

import os
# ── PATH CONFIGURATION ────────────────────────────────────────────────────────
# Set BASE_DIR to the root of your pipeline output directory
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
# ─────────────────────────────────────────────────────────────────────────────
import sys
import warnings
import logging
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

# ── Logging ──────────────────────────────────────────────────────────────────
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)-8s  %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)

# ── Paths ─────────────────────────────────────────────────────────────────────
BASE = Path(__file__).resolve().parent
OUT_DIR = BASE / "outputs" / "target_stereoselectivity"
OUT_DIR.mkdir(parents=True, exist_ok=True)

PAIRS_PATH  = BASE / "outputs" / "true_isomers"      / "TRUE_ISOMER_PAIRS.csv"
MINMAX_PATH = BASE / "outputs" / "stereo_selectivity" / "family_views" / "FAMILY_MINMAX_FOLD.csv"
REFOLD_PATH = BASE / "outputs" / "stereo_selectivity" / "family_views" / "FAMILY_REFERENCE_FOLD.csv"
LIT_PATH    = BASE / "outputs" / "literature_enrichment" / "TRUE_ISOMER_LITERATURE_CLEAN.csv"

# ── TSI weight constants ──────────────────────────────────────────────────────
W_FAMILIES    = 0.20
W_PAIRS       = 0.10
W_MEDIAN      = 0.30
W_STRONG_PROP = 0.25
W_OUTLIER_PEN = 0.15


# ═══════════════════════════════════════════════════════════════════════════════
# Utility helpers
# ═══════════════════════════════════════════════════════════════════════════════

def safe_read(path: Path, **kwargs) -> pd.DataFrame:
    log.info(f"Reading {path.name}  ({path.stat().st_size / 1e6:.1f} MB)")
    return pd.read_csv(path, low_memory=False, **kwargs)


def normalize_0_100(series: pd.Series) -> pd.Series:
    mn, mx = series.min(), series.max()
    if mx == mn:
        return pd.Series(50.0, index=series.index)
    return 100.0 * (series - mn) / (mx - mn)


def pct(s, q: float) -> float:
    return float(np.nanpercentile(s, q)) if len(s) else np.nan


# ═══════════════════════════════════════════════════════════════════════════════
# STEP 1 — Load data
# ═══════════════════════════════════════════════════════════════════════════════

def load_data():
    log.info("=" * 60)
    log.info("STEP 1 — Loading pipeline outputs")
    log.info("=" * 60)

    pairs  = safe_read(PAIRS_PATH)
    minmax = safe_read(MINMAX_PATH)
    refold = safe_read(REFOLD_PATH)
    lit    = safe_read(LIT_PATH)

    log.info(f"TRUE_ISOMER_PAIRS   : {pairs.shape[0]:>8,} rows | {pairs['target_chembl_id'].nunique()} targets")
    log.info(f"FAMILY_MINMAX_FOLD  : {minmax.shape[0]:>8,} rows")
    log.info(f"FAMILY_REFERENCE_FOLD:{refold.shape[0]:>8,} rows")
    log.info(f"LIT_CLEAN           : {lit.shape[0]:>8,} rows")

    minmax["fold_difference"] = pd.to_numeric(minmax["fold_difference"], errors="coerce")
    refold["fold_difference"] = pd.to_numeric(refold["fold_difference"], errors="coerce")

    minmax = minmax[np.isfinite(minmax["fold_difference"]) & (minmax["fold_difference"] > 0)].copy()
    refold = refold[np.isfinite(refold["fold_difference"]) & (refold["fold_difference"] > 0)].copy()

    return pairs, minmax, refold, lit


# ═══════════════════════════════════════════════════════════════════════════════
# STEP 1 — Build target summary
# ═══════════════════════════════════════════════════════════════════════════════

def build_target_summary(pairs: pd.DataFrame) -> pd.DataFrame:
    log.info("Building per-target compound/family/pair counts ...")

    grp = pairs.groupby(["target_chembl_id", "target_name"])

    summary = grp.agg(
        number_of_unique_compounds=("compound_chembl_id", "nunique"),
        number_of_stereochemical_families=("stereochemical_family_id", "nunique"),
        number_of_organisms=("organism", "nunique"),
    ).reset_index()

    # number of true isomer pairs = sum of C(n,2) per family
    pair_counts = (
        pairs.groupby(["target_chembl_id", "stereochemical_family_id"])["compound_chembl_id"]
        .count()
        .reset_index(name="members")
    )
    pair_counts["n_pairs"] = pair_counts["members"].apply(lambda n: max(0, n * (n - 1) // 2))
    pair_total = (
        pair_counts.groupby("target_chembl_id")["n_pairs"]
        .sum()
        .reset_index(name="number_of_true_isomer_pairs")
    )
    summary = summary.merge(pair_total, on="target_chembl_id", how="left")

    act_types = (
        pairs.groupby("target_chembl_id")["activity_type"]
        .apply(lambda x: "|".join(sorted(x.dropna().unique())))
        .reset_index(name="activity_types_present")
    )
    summary = summary.merge(act_types, on="target_chembl_id", how="left")

    log.info(f"  -> {len(summary):,} targets in summary")
    return summary


# ═══════════════════════════════════════════════════════════════════════════════
# STEP 2 — Fold difference profile per target
# ═══════════════════════════════════════════════════════════════════════════════

def build_fold_profile(minmax: pd.DataFrame) -> pd.DataFrame:
    log.info("=" * 60)
    log.info("STEP 2 — Fold difference profile per target")
    log.info("=" * 60)

    def target_profile(g):
        fd = g["fold_difference"].dropna()
        if len(fd) == 0:
            return pd.Series(dtype=float)
        return pd.Series({
            "max_fold"       : fd.max(),
            "min_fold"       : fd.min(),
            "median_fold"    : fd.median(),
            "mean_fold"      : fd.mean(),
            "std_fold"       : fd.std(),
            "iqr_fold"       : fd.quantile(0.75) - fd.quantile(0.25),
            "p90_fold"       : pct(fd, 90),
            "p95_fold"       : pct(fd, 95),
            "families_total" : len(fd),
            "families_ge2"   : int((fd >= 2).sum()),
            "families_ge5"   : int((fd >= 5).sum()),
            "families_ge10"  : int((fd >= 10).sum()),
            "families_ge20"  : int((fd >= 20).sum()),
            "families_ge50"  : int((fd >= 50).sum()),
            "families_ge100" : int((fd >= 100).sum()),
            "n_Mild"         : int((g["classification"] == "Mild").sum()),
            "n_Moderate"     : int((g["classification"] == "Moderate").sum()),
            "n_Strong"       : int((g["classification"] == "Strong").sum()),
            "n_Extreme"      : int((g["classification"] == "Extreme").sum()),
        })

    profile = (
        minmax.groupby(["target_chembl_id", "target_name"])
        .apply(target_profile)
        .reset_index()
    )
    log.info(f"  -> Fold profile computed for {len(profile):,} targets")
    return profile


# ═══════════════════════════════════════════════════════════════════════════════
# STEP 3 — Target Stereoselectivity Index (TSI)
# ═══════════════════════════════════════════════════════════════════════════════

def compute_tsi(summary: pd.DataFrame, profile: pd.DataFrame) -> pd.DataFrame:
    log.info("=" * 60)
    log.info("STEP 3 — Computing Target Stereoselectivity Index (TSI)")
    log.info("=" * 60)

    profile_cols = [
        "target_chembl_id",
        "max_fold", "min_fold", "median_fold", "mean_fold",
        "std_fold", "iqr_fold", "p90_fold", "p95_fold",
        "families_total",
        "families_ge2", "families_ge5", "families_ge10",
        "families_ge20", "families_ge50", "families_ge100",
        "n_Mild", "n_Moderate", "n_Strong", "n_Extreme",
    ]
    df = summary.merge(profile[profile_cols], on="target_chembl_id", how="left")

    df["families_total"] = df["families_total"].fillna(0)
    df["median_fold"]    = df["median_fold"].fillna(1.0)
    df["n_Strong"]       = df["n_Strong"].fillna(0)
    df["n_Extreme"]      = df["n_Extreme"].fillna(0)
    df["families_ge10"]  = df["families_ge10"].fillna(0)

    # Component A: family breadth (log-scaled)
    comp_A = normalize_0_100(np.log10(df["number_of_stereochemical_families"].clip(lower=1)))

    # Component B: pair count (log-scaled)
    comp_B = normalize_0_100(np.log10(df["number_of_true_isomer_pairs"].clip(lower=1)))

    # Component C: median fold (log-scaled)
    comp_C = normalize_0_100(np.log10(df["median_fold"].clip(lower=1)))

    # Component D: proportion Strong+Extreme (0-100)
    df["prop_strong_extreme"] = (
        (df["n_Strong"] + df["n_Extreme"])
        / df["families_total"].replace(0, np.nan)
    ).fillna(0)
    comp_D = df["prop_strong_extreme"] * 100

    # Component E: outlier penalty (100 = fully penalised)
    df["is_outlier_driven"] = (
        (df["n_Extreme"] >= 1)
        & (df["families_ge10"] == df["n_Extreme"])
        & (df["families_total"] > 3)
    ).astype(float)
    penalty = df["is_outlier_driven"] * 100

    # Raw composite score
    df["_raw_tsi"] = (
        W_FAMILIES    * comp_A
        + W_PAIRS     * comp_B
        + W_MEDIAN    * comp_C
        + W_STRONG_PROP * comp_D
        - W_OUTLIER_PEN * penalty
    )

    df["Target_Stereoselectivity_Index"] = normalize_0_100(df["_raw_tsi"]).round(2)
    df.drop(columns=["_raw_tsi"], inplace=True)

    log.info(f"  TSI range: {df['Target_Stereoselectivity_Index'].min():.1f} - {df['Target_Stereoselectivity_Index'].max():.1f}")
    log.info(f"  TSI mean={df['Target_Stereoselectivity_Index'].mean():.1f}  median={df['Target_Stereoselectivity_Index'].median():.1f}")
    return df


# ═══════════════════════════════════════════════════════════════════════════════
# STEP 4 — Classify targets by TSI
# ═══════════════════════════════════════════════════════════════════════════════

def classify_targets(df: pd.DataFrame) -> pd.DataFrame:
    log.info("STEP 4 — Classifying targets by TSI ...")
    tsi = df["Target_Stereoselectivity_Index"]
    conditions = [
        tsi >= 80,
        (tsi >= 60) & (tsi < 80),
        (tsi >= 40) & (tsi < 60),
        (tsi >= 20) & (tsi < 40),
    ]
    choices = [
        "Very High stereoselectivity",
        "High stereoselectivity",
        "Moderate stereoselectivity",
        "Low stereoselectivity",
    ]
    df["TSI_Class"] = np.select(conditions, choices, default="Minimal stereoselectivity")
    for cls, cnt in df["TSI_Class"].value_counts().items():
        log.info(f"    {cls:<38s}: {cnt:>5,}")
    return df


# ═══════════════════════════════════════════════════════════════════════════════
# STEP 5 — Consistency analysis
# ═══════════════════════════════════════════════════════════════════════════════

def consistency_analysis(df: pd.DataFrame) -> pd.DataFrame:
    log.info("STEP 5 — Consistency analysis ...")

    df["consistency_ratio"] = (
        df["families_ge10"].fillna(0) / df["families_total"].replace(0, np.nan)
    ).fillna(0).round(4)

    def consistency_label(row):
        ratio  = row["consistency_ratio"]
        n_fam  = row["families_total"]
        n_ext  = row.get("n_Extreme", 0) or 0
        n_ge10 = row.get("families_ge10", 0) or 0

        if n_fam < 2:
            return "Insufficient data"
        if ratio >= 0.75:
            return "Highly Consistent"
        if ratio >= 0.50:
            return "Consistent"
        if ratio >= 0.25:
            return "Partially Consistent"
        if n_ext >= 1 and n_ge10 == n_ext and n_fam > 3:
            return "Outlier-Driven"
        return "Inconsistent / Low"

    df["consistency_label"] = df.apply(consistency_label, axis=1)
    for lbl, cnt in df["consistency_label"].value_counts().items():
        log.info(f"    {lbl:<30s}: {cnt:>5,}")
    return df


# ═══════════════════════════════════════════════════════════════════════════════
# STEP 6 — Literature support
# ═══════════════════════════════════════════════════════════════════════════════

def literature_support(df: pd.DataFrame, lit: pd.DataFrame,
                       minmax: pd.DataFrame) -> pd.DataFrame:
    log.info("STEP 6 — Literature support ...")

    lit_pubs = lit[lit["publication_count"] > 0].copy()
    lit_pubs["publication_count"] = pd.to_numeric(lit_pubs["publication_count"], errors="coerce").fillna(0)

    lit_target = (
        lit_pubs.groupby("target_chembl_id")
        .agg(
            lit_validated_families=("stereochemical_family_id", "nunique"),
            lit_total_publications=("publication_count", "sum"),
        )
        .reset_index()
    )

    fam_with_lit = lit_pubs[["stereochemical_family_id", "target_chembl_id"]].drop_duplicates()
    minmax_lit   = minmax.merge(fam_with_lit,
                                on=["stereochemical_family_id", "target_chembl_id"],
                                how="inner")
    lit_fold = (
        minmax_lit.groupby("target_chembl_id")["fold_difference"]
        .agg(lit_max_fold="max", lit_median_fold="median")
        .reset_index()
    )

    df = df.merge(lit_target, on="target_chembl_id", how="left")
    df = df.merge(lit_fold,   on="target_chembl_id", how="left")
    df["lit_validated_families"] = df["lit_validated_families"].fillna(0).astype(int)
    df["lit_total_publications"]  = df["lit_total_publications"].fillna(0).astype(int)

    def evidence_level(row):
        fam  = row["lit_validated_families"]
        pubs = row["lit_total_publications"]
        if fam >= 5 or pubs >= 10:
            return "High"
        if fam >= 2 or pubs >= 3:
            return "Medium"
        if fam >= 1 or pubs >= 1:
            return "Low"
        return "None"

    df["literature_evidence_level"] = df.apply(evidence_level, axis=1)
    for lev, cnt in df["literature_evidence_level"].value_counts().items():
        log.info(f"    Evidence {lev:<8s}: {cnt:>5,} targets")
    return df


# ═══════════════════════════════════════════════════════════════════════════════
# STEP 7 — Biological insights per target
# ═══════════════════════════════════════════════════════════════════════════════

def biological_insights(df: pd.DataFrame, minmax: pd.DataFrame) -> pd.DataFrame:
    log.info("STEP 7 — Biological insights per target ...")

    best_fam = (
        minmax.sort_values("fold_difference", ascending=False)
        .groupby("target_chembl_id")
        .first()
        .reset_index()[["target_chembl_id", "stereochemical_family_id",
                        "compound_chembl_id", "fold_difference"]]
        .rename(columns={
            "stereochemical_family_id": "most_selective_family",
            "compound_chembl_id"      : "most_selective_compound",
            "fold_difference"         : "best_family_fold",
        })
    )

    strong_count = (
        minmax[minmax["classification"] == "Strong"]
        .groupby("target_chembl_id").size()
        .reset_index(name="n_strong_families")
    )
    extreme_count = (
        minmax[minmax["classification"] == "Extreme"]
        .groupby("target_chembl_id").size()
        .reset_index(name="n_extreme_families")
    )

    df = df.merge(best_fam,      on="target_chembl_id", how="left")
    df = df.merge(strong_count,  on="target_chembl_id", how="left")
    df = df.merge(extreme_count, on="target_chembl_id", how="left")
    df["n_strong_families"]  = df["n_strong_families"].fillna(0).astype(int)
    df["n_extreme_families"] = df["n_extreme_families"].fillna(0).astype(int)
    return df


# ═══════════════════════════════════════════════════════════════════════════════
# STEP 8 — Priority rankings
# ═══════════════════════════════════════════════════════════════════════════════

def priority_ranking(df: pd.DataFrame) -> pd.DataFrame:
    log.info("STEP 8 — Generating priority rankings ...")
    df = df.copy()
    # fillna(0) before ranking to handle targets with missing fold data
    df["rank_max_fold"]    = df["max_fold"].fillna(0).rank(ascending=False, method="min").astype(int)
    df["rank_median_fold"] = df["median_fold"].fillna(0).rank(ascending=False, method="min").astype(int)
    df["rank_TSI"]         = df["Target_Stereoselectivity_Index"].fillna(0).rank(ascending=False, method="min").astype(int)
    return df


# ═══════════════════════════════════════════════════════════════════════════════
# STEP 9 — Write output CSVs
# ═══════════════════════════════════════════════════════════════════════════════

SUMMARY_COLS = [
    "target_chembl_id", "target_name",
    "number_of_unique_compounds", "number_of_stereochemical_families",
    "number_of_true_isomer_pairs", "number_of_organisms",
    "activity_types_present",
    "max_fold", "min_fold", "median_fold", "mean_fold",
    "std_fold", "iqr_fold", "p90_fold", "p95_fold",
    "families_total", "families_ge2", "families_ge5",
    "families_ge10", "families_ge20", "families_ge50", "families_ge100",
    "n_Mild", "n_Moderate", "n_Strong", "n_Extreme",
    "prop_strong_extreme",
    "Target_Stereoselectivity_Index", "TSI_Class",
    "is_outlier_driven",
    "consistency_ratio", "consistency_label",
    "lit_validated_families", "lit_total_publications",
    "lit_max_fold", "lit_median_fold", "literature_evidence_level",
    "most_selective_family", "most_selective_compound", "best_family_fold",
    "n_strong_families", "n_extreme_families",
    "rank_max_fold", "rank_median_fold", "rank_TSI",
]


def write_outputs(df: pd.DataFrame):
    log.info("STEP 9 — Writing output CSVs ...")

    for col in SUMMARY_COLS:
        if col not in df.columns:
            df[col] = np.nan

    summary = df[SUMMARY_COLS].copy()

    # TARGET_STEREOSELECTIVITY_SUMMARY.csv
    p = OUT_DIR / "TARGET_STEREOSELECTIVITY_SUMMARY.csv"
    summary.to_csv(p, index=False)
    log.info(f"  Saved {p.name}  ({p.stat().st_size / 1e6:.2f} MB)")

    # TARGET_PRIORITY_RANKING.csv (sorted by TSI rank)
    ranking = summary.sort_values("rank_TSI").reset_index(drop=True)
    p = OUT_DIR / "TARGET_PRIORITY_RANKING.csv"
    ranking.to_csv(p, index=False)
    log.info(f"  Saved {p.name}")

    # TARGET_CONSISTENCY_ANALYSIS.csv
    cons_cols = [
        "target_chembl_id", "target_name",
        "number_of_stereochemical_families", "families_total",
        "families_ge10", "n_Extreme",
        "consistency_ratio", "consistency_label",
        "is_outlier_driven",
        "max_fold", "median_fold",
        "Target_Stereoselectivity_Index", "TSI_Class",
        "rank_TSI",
    ]
    for col in cons_cols:
        if col not in df.columns:
            df[col] = np.nan
    cons = df[cons_cols].sort_values("consistency_ratio", ascending=False).reset_index(drop=True)
    p = OUT_DIR / "TARGET_CONSISTENCY_ANALYSIS.csv"
    cons.to_csv(p, index=False)
    log.info(f"  Saved {p.name}")

    # TOP_STEREOSELECTIVE_TARGETS.csv (top 25 by TSI)
    top25 = summary.sort_values("rank_TSI").head(25).reset_index(drop=True)
    p = OUT_DIR / "TOP_STEREOSELECTIVE_TARGETS.csv"
    top25.to_csv(p, index=False)
    log.info(f"  Saved {p.name}  ({len(top25)} targets)")

    return summary


# ═══════════════════════════════════════════════════════════════════════════════
# STEP 10 — Additional insights
# ═══════════════════════════════════════════════════════════════════════════════

def additional_insights(df: pd.DataFrame, minmax: pd.DataFrame,
                        pairs: pd.DataFrame) -> dict:
    log.info("STEP 10 — Generating additional insights ...")
    insights = {}

    # A) Extreme-enriched targets
    insights["extreme_enriched"] = (
        df[df["n_Extreme"] >= 3]
        .sort_values("n_Extreme", ascending=False)
        [["target_chembl_id", "target_name", "n_Extreme",
          "number_of_stereochemical_families", "median_fold",
          "Target_Stereoselectivity_Index", "TSI_Class"]]
        .head(20)
    )

    # B) Many compounds, consistently low stereoselectivity
    insights["many_low"] = (
        df[(df["number_of_unique_compounds"] >= 30)
           & (df["median_fold"] < 2)
           & (df["families_total"] >= 5)]
        .sort_values("number_of_unique_compounds", ascending=False)
        [["target_chembl_id", "target_name",
          "number_of_unique_compounds", "number_of_stereochemical_families",
          "median_fold", "consistency_ratio", "TSI_Class"]]
        .head(20)
    )

    # C) Few compounds, high stereoselectivity
    insights["few_high"] = (
        df[(df["number_of_unique_compounds"] < 10)
           & (df["median_fold"] >= 5)]
        .sort_values("median_fold", ascending=False)
        [["target_chembl_id", "target_name",
          "number_of_unique_compounds", "number_of_stereochemical_families",
          "median_fold", "max_fold",
          "Target_Stereoselectivity_Index", "TSI_Class"]]
        .head(20)
    )

    # D) Activity-type specific: IC50 vs Ki per target
    act_stereo = (
        minmax.groupby(["target_chembl_id", "target_name", "activity_type"])
        .agg(median_fold_act=("fold_difference", "median"))
        .reset_index()
    )
    ic50 = act_stereo[act_stereo["activity_type"] == "IC50"].set_index("target_chembl_id")["median_fold_act"]
    ki   = act_stereo[act_stereo["activity_type"] == "Ki"].set_index("target_chembl_id")["median_fold_act"]
    shared = ic50.index.intersection(ki.index)
    if len(shared):
        act_compare = pd.DataFrame({
            "target_chembl_id"  : shared,
            "IC50_median_fold"  : ic50.loc[shared].values,
            "Ki_median_fold"    : ki.loc[shared].values,
        })
        act_compare["ratio_IC50_to_Ki"] = (
            act_compare["IC50_median_fold"] / act_compare["Ki_median_fold"].replace(0, np.nan)
        )
        act_compare = act_compare.merge(
            df[["target_chembl_id", "target_name"]], on="target_chembl_id", how="left"
        )
        insights["activity_type_specific"] = act_compare.sort_values(
            "ratio_IC50_to_Ki", ascending=False).head(20)
    else:
        insights["activity_type_specific"] = pd.DataFrame()

    # E) Species-specific variability
    sp_stereo = (
        minmax.groupby(["target_chembl_id", "target_name", "organism"])
        .agg(median_fold_sp=("fold_difference", "median"))
        .reset_index()
    )
    multi_sp = sp_stereo.groupby("target_chembl_id").filter(
        lambda g: g["organism"].nunique() > 1
    )
    if len(multi_sp):
        sp_range = (
            multi_sp.groupby(["target_chembl_id", "target_name"])
            .agg(
                max_sp_median=("median_fold_sp", "max"),
                min_sp_median=("median_fold_sp", "min"),
                n_species=("organism", "nunique"),
            )
            .reset_index()
        )
        sp_range["species_variability_ratio"] = (
            sp_range["max_sp_median"]
            / sp_range["min_sp_median"].replace(0, np.nan)
        )
        insights["species_specific"] = sp_range.sort_values(
            "species_variability_ratio", ascending=False).head(20)
    else:
        insights["species_specific"] = pd.DataFrame()

    # F) Correlation: n_families vs median fold
    corr_df = df[["number_of_stereochemical_families", "median_fold"]].dropna()
    if len(corr_df) > 2:
        insights["corr_families_vs_median_fold"] = round(
            corr_df["number_of_stereochemical_families"].corr(corr_df["median_fold"]), 4)
    else:
        insights["corr_families_vs_median_fold"] = None

    # G) Correlation: lit families vs TSI
    lit_df = df[["lit_validated_families", "Target_Stereoselectivity_Index"]].dropna()
    if len(lit_df) > 2:
        insights["corr_lit_vs_TSI"] = round(
            lit_df["lit_validated_families"].corr(
                lit_df["Target_Stereoselectivity_Index"]), 4)
    else:
        insights["corr_lit_vs_TSI"] = None

    log.info(f"  Corr(n_families, median_fold) = {insights['corr_families_vs_median_fold']}")
    log.info(f"  Corr(lit_families, TSI)        = {insights['corr_lit_vs_TSI']}")
    return insights


# ═══════════════════════════════════════════════════════════════════════════════
# STEP 11 — Generate markdown report
# ═══════════════════════════════════════════════════════════════════════════════

def fmt_table(sub: pd.DataFrame, cols: list, col_names: list = None) -> str:
    if sub is None or sub.empty:
        return "_No data available._\n"
    sub = sub[cols].copy()
    if col_names:
        sub.columns = col_names
    lines = ["| " + " | ".join(str(c) for c in sub.columns) + " |"]
    lines.append("|" + "|".join(["---"] * len(sub.columns)) + "|")
    for _, row in sub.iterrows():
        def fmt_val(v):
            if isinstance(v, float):
                if np.isnan(v):
                    return "N/A"
                if v > 1e6:
                    return f"{v:.2e}"
                return f"{v:.2f}"
            return str(v)
        lines.append("| " + " | ".join(fmt_val(v) for v in row) + " |")
    return "\n".join(lines) + "\n"


def generate_report(df: pd.DataFrame, insights: dict):
    log.info("STEP 11 — Generating TARGET_STEREOSELECTIVITY_REPORT.md ...")

    now   = datetime.now().strftime("%Y-%m-%d %H:%M")
    top10 = df.sort_values("rank_TSI").head(10)

    top_consistent = df[df["consistency_label"] == "Highly Consistent"].sort_values(
        "consistency_ratio", ascending=False).head(10)

    top_outlier = df[df["consistency_label"] == "Outlier-Driven"].sort_values(
        "max_fold", ascending=False).head(10)

    tsi_class_counts  = df["TSI_Class"].value_counts()
    cons_label_counts = df["consistency_label"].value_counts()

    docking_targets = df[
        df["TSI_Class"].isin(["Very High stereoselectivity", "High stereoselectivity"])
        & (df["consistency_ratio"] >= 0.50)
    ].sort_values("Target_Stereoselectivity_Index", ascending=False).head(10)

    ml_targets = df[
        (df["number_of_unique_compounds"] >= 30)
        & (df["number_of_stereochemical_families"] >= 5)
    ].sort_values("Target_Stereoselectivity_Index", ascending=False).head(10)

    atlas_targets = df[
        df["literature_evidence_level"].isin(["High", "Medium"])
        & (df["Target_Stereoselectivity_Index"] >= 60)
    ].sort_values("Target_Stereoselectivity_Index", ascending=False).head(10)

    tsi_table = "\n".join(
        f"| {cls} | {cnt:,} |" for cls, cnt in tsi_class_counts.items()
    )
    cons_table = "\n".join(
        f"| {lbl} | {cnt:,} |" for lbl, cnt in cons_label_counts.items()
    )

    corr_fam  = insights.get("corr_families_vs_median_fold", "N/A")
    corr_lit  = insights.get("corr_lit_vs_TSI", "N/A")

    report = f"""# Target-Centric Stereoselectivity Analysis Report
## All_species_v3 Pipeline · Generated {now}

---

## Executive Summary

This report presents a **target-centric stereoselectivity analysis** across
**{len(df):,} biological targets** extracted from the ChEMBL activity database.

The analysis uses three existing pipeline outputs without regenerating ChEMBL data
or recomputing fold differences:
- `TRUE_ISOMER_PAIRS.csv`
- `FAMILY_MINMAX_FOLD.csv`
- `FAMILY_REFERENCE_FOLD.csv`
- `TRUE_ISOMER_LITERATURE_CLEAN.csv` (literature enrichment)

The key deliverable is the **Target Stereoselectivity Index (TSI)** — a composite
0–100 score rewarding breadth, depth, and consistency of stereochemical discrimination.

---

## Dataset Overview

| Metric | Value |
|---|---|
| Total biological targets | {len(df):,} |
| Total stereochemical families (sum across targets) | {int(df["number_of_stereochemical_families"].sum()):,} |
| Total unique compounds | {int(df["number_of_unique_compounds"].sum()):,} |
| Total true isomer pairs | {int(df["number_of_true_isomer_pairs"].sum()):,} |
| Targets with ≥1 Extreme family | {int((df["n_Extreme"] > 0).sum()):,} |
| Targets with ≥1 Strong family | {int((df["n_Strong"] > 0).sum()):,} |
| Global median fold difference | {df["median_fold"].median():.2f} |
| Global max fold difference | {df["max_fold"].max():.2e} |

---

## TSI Class Distribution

| TSI Class | Target Count |
|---|---|
{tsi_table}

## Consistency Label Distribution

| Consistency Label | Target Count |
|---|---|
{cons_table}

---

## TSI Formula (Step 3 — Documented)

```
raw_score = W_FAMILIES    * norm(log10(n_families + 1))
          + W_PAIRS       * norm(log10(n_pairs + 1))
          + W_MEDIAN      * norm(log10(median_fold))
          + W_STRONG_PROP * proportion_strong_extreme   [0-100]
          - W_OUTLIER_PEN * outlier_penalty             [0 or 100]

Weights (sum = 1.00):
  W_FAMILIES    = 0.20  (breadth: more families -> higher score)
  W_PAIRS       = 0.10  (validated pair volume)
  W_MEDIAN      = 0.30  (central tendency: log-scaled median fold)
  W_STRONG_PROP = 0.25  (proportion of Strong+Extreme families)
  W_OUTLIER_PEN = 0.15  (penalty subtracted for outlier-driven targets)

outlier_penalty = 1 (scaled to 100) when:
  - n_Extreme >= 1
  - AND families_ge10 == n_Extreme       (ONLY extreme families pass fold>=10)
  - AND families_total > 3               (at least moderate family count)

TSI = normalize_0_100(raw_score)   [final range: 0 to 100]
```

TSI Classification thresholds:
- **Very High stereoselectivity** : TSI >= 80
- **High stereoselectivity**      : TSI 60–79
- **Moderate stereoselectivity**  : TSI 40–59
- **Low stereoselectivity**       : TSI 20–39
- **Minimal stereoselectivity**   : TSI < 20

---

## Top 10 Stereoselective Targets (by TSI)

{fmt_table(
    top10,
    ["target_chembl_id", "target_name", "number_of_stereochemical_families",
     "median_fold", "max_fold", "prop_strong_extreme",
     "Target_Stereoselectivity_Index", "TSI_Class", "consistency_label"],
    ["ChEMBL ID", "Target Name", "Families", "Median Fold",
     "Max Fold", "Prop S+E", "TSI", "TSI Class", "Consistency"]
)}

---

## Top 10 Most Consistent Targets (Highly Consistent)

{fmt_table(
    top_consistent,
    ["target_chembl_id", "target_name",
     "number_of_stereochemical_families", "families_ge10",
     "consistency_ratio", "median_fold",
     "Target_Stereoselectivity_Index", "TSI_Class"],
    ["ChEMBL ID", "Target Name", "Total Fam", "Fam >=10",
     "Consistency Ratio", "Median Fold", "TSI", "TSI Class"]
)}

---

## Top Outlier-Driven Targets

{fmt_table(
    top_outlier,
    ["target_chembl_id", "target_name", "number_of_stereochemical_families",
     "n_Extreme", "max_fold", "median_fold",
     "Target_Stereoselectivity_Index", "TSI_Class"],
    ["ChEMBL ID", "Target Name", "Families", "N Extreme",
     "Max Fold", "Median Fold", "TSI", "TSI Class"]
) if len(top_outlier) > 0 else "_No strongly outlier-driven targets identified._"}

---

## Ranking A — Highest Maximum Fold Difference

{fmt_table(
    df.sort_values("rank_max_fold").head(10),
    ["rank_max_fold", "target_chembl_id", "target_name",
     "max_fold", "most_selective_family",
     "most_selective_compound", "TSI_Class"],
    ["Rank A", "ChEMBL ID", "Target", "Max Fold",
     "Best Family", "Best Compound", "TSI Class"]
)}

---

## Ranking B — Highest Median Fold Difference

{fmt_table(
    df.sort_values("rank_median_fold").head(10),
    ["rank_median_fold", "target_chembl_id", "target_name",
     "median_fold", "number_of_stereochemical_families",
     "consistency_ratio", "TSI_Class"],
    ["Rank B", "ChEMBL ID", "Target", "Median Fold",
     "Families", "Consistency Ratio", "TSI Class"]
)}

---

## Ranking C — Highest Target Stereoselectivity Index

{fmt_table(
    df.sort_values("rank_TSI").head(10),
    ["rank_TSI", "target_chembl_id", "target_name",
     "Target_Stereoselectivity_Index", "TSI_Class",
     "median_fold", "consistency_ratio", "literature_evidence_level"],
    ["Rank C", "ChEMBL ID", "Target", "TSI", "Class",
     "Median Fold", "Consistency Ratio", "Lit Evidence"]
)}

---

## Step 10 — Additional Insights

### A) Targets Enriched for Extreme Stereoselectivity (n_Extreme >= 3)

{fmt_table(
    insights.get("extreme_enriched", pd.DataFrame()),
    ["target_chembl_id", "target_name", "n_Extreme",
     "number_of_stereochemical_families", "median_fold",
     "Target_Stereoselectivity_Index", "TSI_Class"],
    ["ChEMBL ID", "Target", "N Extreme", "Total Families",
     "Median Fold", "TSI", "TSI Class"]
)}

### B) Targets with Many Compounds but Consistently Low Stereoselectivity (>=30 cpds, median fold <2)

{fmt_table(
    insights.get("many_low", pd.DataFrame()),
    ["target_chembl_id", "target_name",
     "number_of_unique_compounds", "number_of_stereochemical_families",
     "median_fold", "consistency_ratio", "TSI_Class"],
    ["ChEMBL ID", "Target", "Compounds", "Families",
     "Median Fold", "Consistency Ratio", "TSI Class"]
)}

### C) Targets with Few Compounds but High Stereoselectivity (<10 cpds, median fold >=5)

{fmt_table(
    insights.get("few_high", pd.DataFrame()),
    ["target_chembl_id", "target_name",
     "number_of_unique_compounds", "median_fold", "max_fold",
     "Target_Stereoselectivity_Index", "TSI_Class"],
    ["ChEMBL ID", "Target", "Compounds", "Median Fold",
     "Max Fold", "TSI", "TSI Class"]
)}

### D) Activity-Type Specific Stereoselectivity (IC50 vs Ki)

{fmt_table(
    insights.get("activity_type_specific", pd.DataFrame()),
    ["target_chembl_id", "target_name",
     "IC50_median_fold", "Ki_median_fold", "ratio_IC50_to_Ki"],
    ["ChEMBL ID", "Target", "IC50 Median Fold", "Ki Median Fold", "IC50/Ki Ratio"]
) if not insights.get("activity_type_specific", pd.DataFrame()).empty else "_No targets with both IC50 and Ki data found._"}

> Targets with IC50/Ki ratio >> 1 may reflect assay-context effects
> (functional vs binding stereoselectivity divergence).

### E) Species-Specific Stereoselectivity Variability (multi-species targets)

{fmt_table(
    insights.get("species_specific", pd.DataFrame()),
    ["target_chembl_id", "target_name", "n_species",
     "max_sp_median", "min_sp_median", "species_variability_ratio"],
    ["ChEMBL ID", "Target", "N Species",
     "Max Sp. Median Fold", "Min Sp. Median Fold", "Variability Ratio"]
) if not insights.get("species_specific", pd.DataFrame()).empty else "_No multi-species targets found._"}

### F) Relationship: Number of Stereochemical Families vs Median Fold Difference

**Pearson r = {corr_fam}**

> A positive r indicates that targets with more stereochemical families
> tend to show higher median fold differences — implying that stereochemically
> diverse targets are also more stereoselective on average.
> A near-zero r would imply no systematic relationship between target depth
> and stereochemical potency discrimination.

### G) Relationship: Literature Validated Families vs Target Stereoselectivity Index

**Pearson r = {corr_lit}**

> A strong positive r validates the TSI as correlated with independent literature
> evidence for stereoselectivity. A weak r may indicate that the literature
> underrepresents newly identified stereoselective targets in ChEMBL.

---

## Biological Interpretation

### Targets Recommended for Docking / Molecular Dynamics

High TSI + high consistency ratio = systematic stereochemical discrimination
suitable for structure-activity stereospecific docking:

{fmt_table(
    docking_targets,
    ["target_chembl_id", "target_name", "median_fold",
     "consistency_ratio", "Target_Stereoselectivity_Index", "TSI_Class"],
    ["ChEMBL ID", "Target", "Median Fold",
     "Consistency", "TSI", "TSI Class"]
)}

### Targets Recommended for Machine Learning / Descriptor Analysis

Large compound sets (>=30) + multiple families (>=5) offer sufficient
chemical diversity for stereospecific QSAR or ML models:

{fmt_table(
    ml_targets,
    ["target_chembl_id", "target_name",
     "number_of_unique_compounds", "number_of_stereochemical_families",
     "median_fold", "Target_Stereoselectivity_Index", "TSI_Class"],
    ["ChEMBL ID", "Target", "Compounds", "Families",
     "Median Fold", "TSI", "TSI Class"]
)}

### Targets Recommended as StereoAtlas Case Studies

High literature evidence + high TSI = publication-ready stereoselective targets:

{fmt_table(
    atlas_targets,
    ["target_chembl_id", "target_name", "literature_evidence_level",
     "lit_validated_families", "lit_total_publications",
     "median_fold", "Target_Stereoselectivity_Index", "TSI_Class"],
    ["ChEMBL ID", "Target", "Lit Evidence", "Lit Families",
     "Lit Pubs", "Median Fold", "TSI", "TSI Class"]
)}

---

## Recommendations for Future Work

1. **Docking/MD studies**: Prioritise targets with TSI >= 60 and consistency_ratio >= 0.75.
   These exhibit systematic, reproducible stereochemical discrimination ideal
   for computational stereospecific binding predictions.

2. **Machine learning**: Targets with >=30 unique compounds and >=5 stereochemical
   families provide adequate training data for stereospecific activity models.

3. **Publication figures**: The top 10 by TSI and top 10 by consistency ratio
   offer the most compelling biological stereoselectivity narratives.

4. **StereoAtlas case studies**: High-literature-evidence targets (lit_evidence = High)
   with TSI >= 80 are the gold standard for StereoAtlas manuscript inclusion.

5. **Outlier investigation**: Targets classified as "Outlier-Driven" warrant
   structural scrutiny — single extreme-fold families may reflect unusual
   binding modes, data artefacts, or stereocenter assignment errors.

6. **Species divergence**: Multi-species targets with high species_variability_ratio
   are candidates for comparative pharmacology and cross-species selectivity studies.

7. **Activity-type discordance**: Targets where IC50/Ki ratio >> 1 may indicate
   functional vs binding assay context effects and merit mechanistic investigation.

---

## Output Files Generated

| File | Description |
|---|---|
| `TARGET_STEREOSELECTIVITY_SUMMARY.csv` | One row per target, all metrics |
| `TARGET_PRIORITY_RANKING.csv` | Sorted by TSI (rank_TSI) |
| `TARGET_CONSISTENCY_ANALYSIS.csv` | Consistency ratio + label per target |
| `TOP_STEREOSELECTIVE_TARGETS.csv` | Top 25 by TSI |
| `TARGET_STEREOSELECTIVITY_REPORT.md` | This report |

---

*Report generated by `12_target_stereoselectivity.py`*
*Pipeline: All_species_v3 | Analysis date: {now}*
"""

    report_path = OUT_DIR / "TARGET_STEREOSELECTIVITY_REPORT.md"
    report_path.write_text(report, encoding="utf-8")
    log.info(f"  Saved {report_path.name}  ({report_path.stat().st_size / 1e3:.1f} kB)")


# ═══════════════════════════════════════════════════════════════════════════════
# Main
# ═══════════════════════════════════════════════════════════════════════════════

def main():
    t0 = datetime.now()
    log.info("=" * 60)
    log.info("12_target_stereoselectivity.py -- START")
    log.info("=" * 60)

    pairs, minmax, refold, lit = load_data()

    target_summary = build_target_summary(pairs)
    fold_profile   = build_fold_profile(minmax)

    df = compute_tsi(target_summary, fold_profile)
    df = classify_targets(df)
    df = consistency_analysis(df)
    df = literature_support(df, lit, minmax)
    df = biological_insights(df, minmax)
    df = priority_ranking(df)

    write_outputs(df)

    insights = additional_insights(df, minmax, pairs)

    generate_report(df, insights)

    elapsed = (datetime.now() - t0).total_seconds()

    # Final console summary
    print("\n" + "=" * 60)
    print("TARGET STEREOSELECTIVITY ANALYSIS -- COMPLETE")
    print("=" * 60)
    print(f"  Elapsed time                 : {elapsed:.1f}s")
    print(f"  Targets analysed             : {len(df):>6,}")
    print(f"  Stereochemical families (sum): {int(df['number_of_stereochemical_families'].sum()):>6,}")
    print(f"  True isomer pairs (sum)      : {int(df['number_of_true_isomer_pairs'].sum()):>6,}")
    print(f"  Unique compounds (sum)       : {int(df['number_of_unique_compounds'].sum()):>6,}")
    print()
    print("  TSI Class distribution:")
    for cls, cnt in df["TSI_Class"].value_counts().items():
        print(f"    {cls:<40s}: {cnt:>5,}")
    print()
    print(f"  Mean TSI   : {df['Target_Stereoselectivity_Index'].mean():.1f}")
    print(f"  Median TSI : {df['Target_Stereoselectivity_Index'].median():.1f}")
    print()
    print("  Top 10 by TSI:")
    for _, row in df.sort_values("rank_TSI").head(10).iterrows():
        name = str(row["target_name"])[:48]
        print(f"    {row['rank_TSI']:>3}. {name:<48s}  TSI={row['Target_Stereoselectivity_Index']:.1f}")
    print()
    print("  Consistency breakdown:")
    for lbl, cnt in df["consistency_label"].value_counts().items():
        print(f"    {lbl:<35s}: {cnt:>5,}")
    print()
    print(f"  Output directory: {OUT_DIR}")
    print("=" * 60)


if __name__ == "__main__":
    main()
