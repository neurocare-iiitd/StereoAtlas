"""
All-Species v3 -- Step 10: TRUE Isomer Fold-Difference Analysis
===============================================================
Source: outputs/true_isomers/TRUE_ISOMER_PAIRS.csv  (already filtered)
        No DB access. No re-extraction. No re-aggregation.

Strategies retained in this release:
  S1 (MinMax): For each (family × target × organism × metric × unit) group
               with ≥2 distinct stereoforms, fold_difference = max_value / min_value.
               One representative row per group; the most biologically informative
               contrast captured by a single number.

  S2 (Summary): Aggregate group-level fold statistics into target-, species-, and
                metric-wise summary tables and a top-100 ranking.

Note: Full pairwise enumeration (all N×(N-1)/2 combinations per group) and S3
(all-vs-reference) have been removed from this manuscript version.

Outputs (outputs/stereo_selectivity/):
  FAMILY_MINMAX_FOLD.csv          -- one row per group (S1)
  target_wise_summary.csv
  species_wise_summary.csv
  metric_wise_summary.csv
  TOP100_STRONGEST_STEREOSELECTIVE_PAIRS.csv
  STEREOSELECTIVITY_SUMMARY.csv
"""

import pandas as pd
import numpy as np
import os, time

V3_DIR   = os.path.dirname(os.path.abspath(__file__))
SRC_FILE = os.path.join(V3_DIR, "outputs", "true_isomers", "TRUE_ISOMER_PAIRS.csv")
OUT_DIR  = os.path.join(V3_DIR, "outputs", "stereo_selectivity")
REP_DIR  = os.path.join(V3_DIR, "reports")
os.makedirs(OUT_DIR, exist_ok=True)

PAIR_KEYS = [
    "stereochemical_family_id",
    "target_name", "target_chembl_id",
    "organism",
    "activity_type",
    "activity_units",
]

FUNCTIONAL_METRICS = {
    "Potency","Emax","Efficacy","Activity","Inhibition",
    "Response","AC50","GI50","ED50","Imax","DC50","pD2"
}

# (upper_threshold, label) -- None means no upper bound
TIERS = [
    (2,    "No stereoselectivity (<2x)"),
    (5,    "Mild stereoselectivity (2-5x)"),
    (10,   "Moderate stereoselectivity (5-10x)"),
    (50,   "Strong stereoselectivity (10-50x)"),
    (None, "Extreme stereoselectivity (>50x)"),
]

def classify_fold(fold):
    if pd.isna(fold) or fold <= 0:
        return "Undetermined"
    for threshold, label in TIERS:
        if threshold is None or fold < threshold:
            return label
    return "Extreme stereoselectivity (>50x)"


# ---- Step 1: Load --------------------------------------------------------
def load_data():
    print(f"  Loading {os.path.basename(SRC_FILE)} ...")
    ti = pd.read_csv(SRC_FILE, low_memory=False)
    print(f"  {len(ti):,} rows")
    ti["mean_activity_value"] = pd.to_numeric(ti["mean_activity_value"], errors="coerce")
    if "target_chembl_id" not in ti.columns:
        ti["target_chembl_id"] = ""
    return ti


# ---- Step 2: Aggregate one representative value per (group+stereo_form) -
def aggregate_per_form(ti):
    agg = {}
    for col, func in [
        ("mean_activity_value", "median"),
        ("compound_chembl_id",  "first"),
        ("compound_name",       "first"),
        ("smiles_nostereo",     "first"),
        ("n_stereocenters",     "first"),
        ("n_measurements",      "sum"),
    ]:
        if col in ti.columns:
            agg[col] = pd.NamedAgg(column=col, aggfunc=func)

    pf = (
        ti.groupby(PAIR_KEYS + ["stereo_form", "stereochemical_variant"], dropna=False)
        .agg(**agg)
        .reset_index()
        .rename(columns={"mean_activity_value": "activity_value"})
    )
    return pf


# ---- Step 3: MinMax fold per group (S1) ---------------------------------
def compute_minmax(pf):
    """
    For each (group + stereo_form) slot, identify the variant with the
    minimum and maximum activity value. fold_difference = max / min.
    Returns one representative row per group.
    """
    rows = []
    for keys, grp in pf.groupby(PAIR_KEYS):
        valid = grp.dropna(subset=["activity_value"])
        valid = valid[valid["activity_value"] > 0]
        if len(valid) < 2:
            continue
        idx_min = valid["activity_value"].idxmin()
        idx_max = valid["activity_value"].idxmax()
        row_min = valid.loc[idx_min]
        row_max = valid.loc[idx_max]
        min_val = float(row_min["activity_value"])
        max_val = float(row_max["activity_value"])
        if min_val <= 0:
            continue
        fold = max_val / min_val
        log2f = np.log2(fold)
        key_dict = dict(zip(PAIR_KEYS, keys))
        rows.append({
            **key_dict,
            "smiles_nostereo_A":       row_min.get("smiles_nostereo", ""),
            "n_stereocenters_A":       row_min.get("n_stereocenters", np.nan),
            "compound_chembl_id_A":    row_min.get("compound_chembl_id", ""),
            "compound_name_A":         row_min.get("compound_name", ""),
            "compound_chembl_id_B":    row_max.get("compound_chembl_id", ""),
            "compound_name_B":         row_max.get("compound_name", ""),
            "stereo_form_A":           row_min.get("stereo_form", ""),
            "stereochemical_variant_A": row_min.get(VARIANT_COL, ""),
            "value_A":                 min_val,
            "stereo_form_B":           row_max.get("stereo_form", ""),
            "stereochemical_variant_B": row_max.get(VARIANT_COL, ""),
            "value_B":                 max_val,
            "fold_difference":         fold,
            "log2_fold":               log2f,
            "more_potent_variant":     row_min.get(VARIANT_COL, ""),
            "stereoselectivity_class": classify_fold(fold),
        })
    return pd.DataFrame(rows)


# (Fold computation is now inlined in compute_minmax above)


# ---- Step 5: Column ordering for S1 output ------------------------------
OUT_COLS = [
    "stereochemical_family_id",
    "smiles_nostereo_A", "n_stereocenters_A",
    "compound_chembl_id_A", "compound_name_A",
    "compound_chembl_id_B", "compound_name_B",
    "target_chembl_id", "target_name",
    "organism",
    "activity_type", "activity_units",
    "stereo_form_A", "stereochemical_variant_A", "value_A",
    "stereo_form_B", "stereochemical_variant_B", "value_B",
    "fold_difference", "log2_fold",
    "more_potent_variant",
    "stereoselectivity_class",
]

def tidy(pairs):
    cols = [c for c in OUT_COLS if c in pairs.columns]
    return pairs[cols].copy()


# ---- (Per-metric pair files removed — S1 master file covers all metrics) --


# ---- Step 7: Summary tables ---------------------------------------------
def build_summaries(pairs):
    v = pairs.dropna(subset=["fold_difference"])

    # Target-wise
    tgt = (
        v.groupby(["target_name","target_chembl_id"])
        .agg(
            pair_count=("fold_difference","count"),
            compound_count=("compound_chembl_id_A","nunique"),
            organism_count=("organism","nunique"),
            mean_fold=("fold_difference","mean"),
            median_fold=("fold_difference","median"),
            max_fold=("fold_difference","max"),
        )
        .reset_index().sort_values("max_fold", ascending=False)
    )
    tgt.to_csv(os.path.join(OUT_DIR, "target_wise_summary.csv"), index=False)
    print(f"    target_wise_summary.csv  -- {len(tgt):,} targets")

    # Species-wise
    spe = (
        v.groupby("organism")
        .agg(
            pair_count=("fold_difference","count"),
            compound_count=("compound_chembl_id_A","nunique"),
            target_count=("target_name","nunique"),
            mean_fold=("fold_difference","mean"),
            median_fold=("fold_difference","median"),
            max_fold=("fold_difference","max"),
        )
        .reset_index().sort_values("pair_count", ascending=False)
    )
    spe.to_csv(os.path.join(OUT_DIR, "species_wise_summary.csv"), index=False)
    print(f"    species_wise_summary.csv -- {len(spe):,} organisms")

    # Metric-wise
    met = (
        v.groupby("activity_type")
        .agg(
            pair_count=("fold_difference","count"),
            compound_count=("compound_chembl_id_A","nunique"),
            target_count=("target_name","nunique"),
            organism_count=("organism","nunique"),
            mean_fold=("fold_difference","mean"),
            median_fold=("fold_difference","median"),
            max_fold=("fold_difference","max"),
        )
        .reset_index().sort_values("pair_count", ascending=False)
    )
    met.to_csv(os.path.join(OUT_DIR, "metric_wise_summary.csv"), index=False)
    print(f"    metric_wise_summary.csv  -- {len(met):,} activity types")
    return tgt, spe, met


# ---- Step 8: Top-100 ranking --------------------------------------------
TOP_COLS = [
    "stereochemical_family_id",
    "compound_chembl_id_A","compound_name_A",
    "compound_chembl_id_B","compound_name_B",
    "target_name","target_chembl_id",
    "organism",
    "activity_type","activity_units",
    "stereochemical_variant_A","stereochemical_variant_B",
    "value_A","value_B",
    "fold_difference","log2_fold",
    "more_potent_variant",
    "stereoselectivity_class",
]

def build_top100(pairs):
    v = pairs.dropna(subset=["fold_difference"])
    cols = [c for c in TOP_COLS if c in v.columns]
    top = v.nlargest(100, "fold_difference")[cols]
    top.to_csv(os.path.join(OUT_DIR, "TOP100_STRONGEST_STEREOSELECTIVE_PAIRS.csv"), index=False)
    print(f"    TOP100_STRONGEST_STEREOSELECTIVE_PAIRS.csv saved")
    return top


# ---- Step 9: Global summary ---------------------------------------------
def build_global_summary(pairs):
    v = pairs.dropna(subset=["fold_difference"])
    n_total = len(pairs)
    n_valid = len(v)
    cdist = v["stereoselectivity_class"].value_counts()
    cpct  = (cdist / n_valid * 100).round(2)

    rows = []

    def make_row(label, sub_all, sub_valid):
        cd = sub_valid["stereoselectivity_class"].value_counts()
        cp = (cd / len(sub_valid) * 100).round(2) if len(sub_valid) else pd.Series(dtype=float)
        return {
            "metric":                   label,
            "total_pairs":              len(sub_all),
            "valid_pairs":              len(sub_valid),
            "compound_count":           sub_valid["compound_chembl_id_A"].nunique(),
            "target_count":             sub_valid["target_name"].nunique(),
            "organism_count":           sub_valid["organism"].nunique(),
            "mean_fold_difference":     round(sub_valid["fold_difference"].mean(), 2) if len(sub_valid) else None,
            "median_fold_difference":   round(sub_valid["fold_difference"].median(), 2) if len(sub_valid) else None,
            "max_fold_difference":      round(sub_valid["fold_difference"].max(), 2) if len(sub_valid) else None,
            "no_selectivity_pct":       cp.get("No stereoselectivity (<2x)", 0),
            "mild_selectivity_pct":     cp.get("Mild stereoselectivity (2-5x)", 0),
            "moderate_selectivity_pct": cp.get("Moderate stereoselectivity (5-10x)", 0),
            "strong_selectivity_pct":   cp.get("Strong stereoselectivity (10-50x)", 0),
            "extreme_selectivity_pct":  cp.get("Extreme stereoselectivity (>50x)", 0),
        }

    rows.append(make_row("ALL", pairs, v))

    for metric in ["IC50","Ki","Kd","EC50","Km","Potency","Emax","Activity","Inhibition"]:
        m_all   = pairs[pairs["activity_type"].str.strip().str.upper() == metric.upper()]
        m_valid = v[v["activity_type"].str.strip().str.upper() == metric.upper()]
        if len(m_valid) == 0:
            continue
        rows.append(make_row(metric, m_all, m_valid))

    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(OUT_DIR, "STEREOSELECTIVITY_SUMMARY.csv"), index=False)
    return df


# ---- Main ---------------------------------------------------------------
def main():
    print("=" * 65)
    print("Step 10: TRUE Isomer Fold-Difference Analysis (S1 + S2)")
    print("  Source: outputs/true_isomers/TRUE_ISOMER_PAIRS.csv")
    print("  Strategies: S1=MinMax, S2=Summary tables")
    print("=" * 65)
    t0 = time.time()

    # ---- [1/4] Load -------------------------------------------------------
    print("\n[1/4] Loading data...")
    ti = load_data()

    # ---- [2/4] Aggregate per (group + stereo_form) ------------------------
    print("[2/4] Aggregating per (group + stereo_form)...")
    pf = aggregate_per_form(ti)
    # Rename activity_value column for compatibility with compute_minmax
    if "activity_value" not in pf.columns and "mean_activity_value" in pf.columns:
        pf = pf.rename(columns={"mean_activity_value": "activity_value"})
    # Also ensure stereo_form column present
    if "stereo_form" not in pf.columns and "stereochemical_variant" in pf.columns:
        pf["stereo_form"] = pf["stereochemical_variant"]
    VARIANT_COL = "stereochemical_variant"
    print(f"  {len(pf):,} unique (context + stereo_form) slots")

    # ---- [3/4] S1 MinMax fold per group -----------------------------------
    print("[3/4] S1 — Computing MinMax fold difference per group...")
    pairs = compute_minmax(pf)
    pairs = tidy(pairs)
    pairs = pairs.sort_values("fold_difference", ascending=False)

    valid_n = pairs["fold_difference"].notna().sum()
    print(f"  Groups processed (S1)          : {len(pairs):,}")
    print(f"  Valid (both values > 0)        : {valid_n:,}")
    print(f"  Undetermined (zero/neg/NaN)    : {len(pairs) - valid_n:,}")

    cdist = pairs["stereoselectivity_class"].value_counts()
    for cls, cnt in cdist.items():
        pct = 100 * cnt / max(len(pairs), 1)
        print(f"    {cls:<46s}: {cnt:>8,} ({pct:.1f}%)")

    # Save S1 master file
    master_path = os.path.join(OUT_DIR, "FAMILY_MINMAX_FOLD.csv")
    pairs.to_csv(master_path, index=False)
    sz = os.path.getsize(master_path) / 1e6
    print(f"  Saved FAMILY_MINMAX_FOLD.csv -- {len(pairs):,} rows ({sz:.1f} MB)")

    # ---- [4/4] S2 Summary tables -----------------------------------------
    print("\n[4/4] S2 — Building summary tables...")
    build_summaries(pairs)
    build_top100(pairs)
    build_global_summary(pairs)

    # Key results
    v = pairs.dropna(subset=["fold_difference"])
    print(f"\n{'─'*65}")
    print("KEY RESULTS")
    print(f"{'─'*65}")
    print(f"  S1 groups processed            : {len(pairs):,}")
    print(f"  Valid (computable)             : {valid_n:,}")
    print(f"  Unique targets                 : {v['target_name'].nunique():,}")
    print(f"  Unique organisms               : {v['organism'].nunique():,}")
    print(f"  Mean fold difference           : {v['fold_difference'].mean():.1f}x")
    print(f"  Median fold difference         : {v['fold_difference'].median():.1f}x")
    print(f"  Max fold difference            : {v['fold_difference'].max():.0f}x")

    print("\n  Top 10 most stereoselective groups (S1):")
    t10_cols = [c for c in [
        "stereochemical_variant_A", "stereochemical_variant_B",
        "target_name", "organism", "activity_type", "activity_units",
        "value_A", "value_B", "fold_difference", "stereoselectivity_class"
    ] if c in v.columns]
    top10 = v.nlargest(10, "fold_difference")[t10_cols]
    for _, r in top10.iterrows():
        tgt  = str(r["target_name"])[:28]
        vA   = f"{r['value_A']:.3g}"
        vB   = f"{r['value_B']:.3g}"
        fold = f"{r['fold_difference']:.0f}x"
        print(f"    {str(r.get('stereochemical_variant_A','')):<5s} vs "
              f"{str(r.get('stereochemical_variant_B','')):<8s}"
              f"| {tgt:<30s} | {r['activity_type']:<6s}"
              f"| {vA} vs {vB}  => {fold}")

    elapsed = time.time() - t0
    print(f"\n  Time: {elapsed:.1f}s")
    print(f"  Outputs: {OUT_DIR}/")
    print(f"\n  Files generated:")
    for f in sorted(os.listdir(OUT_DIR)):
        if os.path.isfile(os.path.join(OUT_DIR, f)):
            sz = os.path.getsize(os.path.join(OUT_DIR, f)) / 1e6
            print(f"    {f:<55s} {sz:.1f} MB")
    print(f"\n[OK] Step 10 complete.")

if __name__ == "__main__":
    main()
