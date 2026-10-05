"""
All-Species v3 -- Step 11: Family-Level Stereoselectivity Analysis
==================================================================
Script:  11_family_view.py
Source:  outputs/true_isomers/TRUE_ISOMER_PAIRS.csv
Outputs: outputs/stereo_selectivity/family_views/

Strategies:
  S1: MinMax Fold       -- min-activity vs max-activity variant per family-context
                          One representative row per group; the most biologically
                          informative single contrast.
  S2: Reference Variant -- all variants compared to a single chosen reference,
                          useful for multi-stereoform families (n>=3).

Note: S3 (full pairwise enumeration of all N×(N-1)/2 combinations) has been
removed from this manuscript version. S1 and S2 together cover all relevant
analytical perspectives included in the final publication.

Additional outputs:
  TOP100_STEREOSELECTIVE_TARGETS.csv
  reports/11_family_view_report.txt

No DB access. No re-aggregation. Read-only source.
"""

import pandas as pd
import numpy as np
import os, time, math, sys
from itertools import combinations
from tqdm import tqdm
from datetime import datetime

# ---- Paths ---------------------------------------------------------------
V3_DIR    = os.path.dirname(os.path.abspath(__file__))
SRC_FILE  = os.path.join(V3_DIR, "outputs", "true_isomers", "TRUE_ISOMER_PAIRS.csv")
OUT_DIR   = os.path.join(V3_DIR, "outputs", "stereo_selectivity", "family_views")
METRIC_DIR= os.path.join(OUT_DIR, "by_metric")
SEL_DIR   = os.path.join(V3_DIR, "outputs", "stereo_selectivity")
REP_DIR   = os.path.join(V3_DIR, "reports")
os.makedirs(OUT_DIR, exist_ok=True)
os.makedirs(METRIC_DIR, exist_ok=True)

# ---- Group key -----------------------------------------------------------
PAIR_KEYS = [
    "stereochemical_family_id",
    "target_name", "target_chembl_id",
    "organism",
    "activity_type",
    "activity_units",
]
VARIANT_COL = "stereochemical_variant"
VALUE_COL   = "mean_activity_value"

# ---- Metrics for per-metric files ----------------------------------------
METRIC_FILES = {
    "IC50":       "IC50",
    "Ki":         "Ki",
    "Kd":         "Kd",
    "EC50":       "EC50",
    "Km":         "Km",
    "Vmax":       "Vmax",
    "Kcat":       "Kcat",
    "Kinact":     "Kinact",
    "Potency":    "Potency",
    "Emax":       "Emax",
    "Activity":   "Activity",
    "Inhibition": "Inhibition",
}

# ---- Classification ------------------------------------------------------
def classify(fold):
    """Classify fold difference into selectivity tier."""
    if pd.isna(fold) or fold <= 0:
        return "Undetermined"
    if fold < 2:   return "None"
    if fold < 5:   return "Mild"
    if fold < 10:  return "Moderate"
    if fold < 50:  return "Strong"
    return "Extreme"

def safe_fold(vA, vB):
    """Compute fold difference and log10. Returns (fold, log10_fold) or (NaN, NaN)."""
    try:
        vA, vB = float(vA), float(vB)
        if vA > 0 and vB > 0:
            fold = max(vA, vB) / min(vA, vB)
            return fold, math.log10(fold)
    except Exception:
        pass
    return float("nan"), float("nan")

# ---- Save helper ---------------------------------------------------------
def save(df, path):
    if len(df) == 0:
        return
    if os.path.exists(path):
        print(f"    [skip] exists: {os.path.basename(path)}")
        return
    df.to_csv(path, index=False)

def save_by_metric(df, suffix):
    """Save per-metric files into by_metric/."""
    for metric, label in METRIC_FILES.items():
        mask = df["activity_type"].str.strip().str.upper() == metric.upper()
        sub = df[mask]
        if len(sub):
            fpath = os.path.join(METRIC_DIR, f"{label}_{suffix}.csv")
            save(sub, fpath)

# ==========================================================================
# LOAD + PRE-AGGREGATE
# ==========================================================================
def load_and_aggregate():
    """
    Load TRUE_ISOMER_PAIRS.csv, aggregate one representative value per
    (PAIR_KEYS + VARIANT_COL) slot.
    Returns per_form DataFrame and raw stats.
    """
    print(f"  Source : {SRC_FILE}")
    ti = pd.read_csv(SRC_FILE, low_memory=False)
    ti[VALUE_COL] = pd.to_numeric(ti[VALUE_COL], errors="coerce")
    if "target_chembl_id" not in ti.columns:
        ti["target_chembl_id"] = ""
    if "n_measurements" not in ti.columns:
        ti["n_measurements"] = 1

    stats = {
        "total_rows":      len(ti),
        "unique_compounds": ti["compound_chembl_id"].nunique(),
        "unique_targets":   ti["target_name"].nunique(),
        "unique_organisms": ti["organism"].nunique(),
        "unique_families":  ti["stereochemical_family_id"].nunique(),
    }
    print(f"  Total rows         : {stats['total_rows']:,}")
    print(f"  Unique compounds   : {stats['unique_compounds']:,}")
    print(f"  Unique targets     : {stats['unique_targets']:,}")
    print(f"  Unique organisms   : {stats['unique_organisms']:,}")
    print(f"  Unique families    : {stats['unique_families']:,}")

    # Aggregate: one row per (group_key + variant)
    agg_fns = {
        VALUE_COL:            ("median"),
        "n_measurements":     ("sum"),
        "compound_chembl_id": ("first"),
        "compound_name":      ("first"),
        "n_stereocenters":    ("first"),
        "smiles_nostereo":    ("first"),
    }
    pf = (
        ti.groupby(PAIR_KEYS + [VARIANT_COL], dropna=False)
        .agg({k: v for k, v in agg_fns.items() if k in ti.columns})
        .reset_index()
    )
    pf = pf.rename(columns={VALUE_COL: "activity_value"})

    # Count n_stereoforms per group
    n_forms = (
        pf.groupby(PAIR_KEYS)[VARIANT_COL]
        .nunique()
        .reset_index()
        .rename(columns={VARIANT_COL: "n_stereoforms"})
    )
    pf = pf.merge(n_forms, on=PAIR_KEYS, how="left")

    # Only keep groups with >= 2 distinct variants
    pf = pf[pf["n_stereoforms"] >= 2].copy()

    total_groups = pf.groupby(PAIR_KEYS).ngroups
    print(f"  Groups with >=2 variants: {total_groups:,}")
    print(f"  Per-form slots          : {len(pf):,}")
    return pf, stats, total_groups

# ==========================================================================
# STRATEGY 1: MIN/MAX FOLD
# ==========================================================================
def strategy1_minmax(pf):
    """
    For each family group: identify variant with min value and variant with
    max value. Compute fold_difference = max/min.
    """
    rows = []
    groups = list(pf.groupby(PAIR_KEYS))

    for keys, grp in tqdm(groups, desc="  S1 MinMax", unit="grp", leave=False):
        grp_valid = grp.dropna(subset=["activity_value"])
        grp_valid = grp_valid[grp_valid["activity_value"] > 0]
        if len(grp_valid) < 2:
            continue
        try:
            n_forms    = int(grp["n_stereoforms"].iloc[0])
            idx_min    = grp_valid["activity_value"].idxmin()
            idx_max    = grp_valid["activity_value"].idxmax()
            row_min    = grp_valid.loc[idx_min]
            row_max    = grp_valid.loc[idx_max]
            min_val    = float(row_min["activity_value"])
            max_val    = float(row_max["activity_value"])
            if min_val <= 0:
                continue
            fold, log10f = safe_fold(min_val, max_val)
            key_dict = dict(zip(PAIR_KEYS, keys))
            rows.append({
                **key_dict,
                "compound_chembl_id": row_min["compound_chembl_id"],
                "compound_name":      row_min.get("compound_name", ""),
                "n_stereoforms":      n_forms,
                "min_variant":        row_min[VARIANT_COL],
                "min_value":          min_val,
                "max_variant":        row_max[VARIANT_COL],
                "max_value":          max_val,
                "fold_difference":    fold,
                "log10_fold_difference": log10f,
                "classification":     classify(fold),
            })
        except Exception as e:
            continue

    df = pd.DataFrame(rows).sort_values("fold_difference", ascending=False)
    return df

# ==========================================================================
# STRATEGY 2: REFERENCE VARIANT
# ==========================================================================
def pick_reference(grp):
    """
    Select the reference variant for a group.
    Priority: 1) highest n_measurements, 2) lowest activity_value, 3) alphabetical
    """
    # Step 1: highest assay count
    if "n_measurements" in grp.columns:
        max_count = grp["n_measurements"].max()
        candidates = grp[grp["n_measurements"] == max_count]
        if len(candidates) == 1:
            return candidates.iloc[0], "highest_count"
    else:
        candidates = grp

    # Step 2: lowest (most potent) value
    valid = candidates.dropna(subset=["activity_value"])
    valid = valid[valid["activity_value"] > 0]
    if len(valid) > 0:
        min_val = valid["activity_value"].min()
        potent = valid[valid["activity_value"] == min_val]
        if len(potent) == 1:
            return potent.iloc[0], "most_potent"
        candidates = potent

    # Step 3: alphabetical on variant
    if len(candidates) > 0:
        candidates = candidates.sort_values(VARIANT_COL)
        return candidates.iloc[0], "alphabetical"

    return grp.iloc[0], "alphabetical"

def strategy2_reference(pf):
    """
    For each group, pick a reference variant, then compare all others to it.
    """
    rows = []
    groups = list(pf.groupby(PAIR_KEYS))
    ref_reason_counts = {"highest_count": 0, "most_potent": 0, "alphabetical": 0}

    for keys, grp in tqdm(groups, desc="  S2 Reference", unit="grp", leave=False):
        try:
            n_forms = int(grp["n_stereoforms"].iloc[0])
            ref_row, reason = pick_reference(grp)
            ref_reason_counts[reason] = ref_reason_counts.get(reason, 0) + 1
            ref_variant = ref_row[VARIANT_COL]
            ref_value   = float(ref_row["activity_value"]) if pd.notna(ref_row["activity_value"]) else float("nan")

            key_dict = dict(zip(PAIR_KEYS, keys))
            comparisons = grp[grp[VARIANT_COL] != ref_variant]

            for _, comp_row in comparisons.iterrows():
                comp_variant = comp_row[VARIANT_COL]
                comp_value   = float(comp_row["activity_value"]) if pd.notna(comp_row["activity_value"]) else float("nan")
                fold, log10f = safe_fold(ref_value, comp_value)
                rows.append({
                    **key_dict,
                    "compound_chembl_id":    ref_row["compound_chembl_id"],
                    "compound_name":         ref_row.get("compound_name", ""),
                    "n_stereoforms":         n_forms,
                    "reference_variant":     ref_variant,
                    "reference_value":       ref_value,
                    "reference_selection_reason": reason,
                    "comparison_variant":    comp_variant,
                    "comparison_value":      comp_value,
                    "fold_difference":       fold,
                    "log10_fold_difference": log10f,
                    "classification":        classify(fold),
                })
        except Exception:
            continue

    df = pd.DataFrame(rows).sort_values("fold_difference", ascending=False)
    return df, ref_reason_counts


# ==========================================================================
# TARGET-LEVEL SUMMARY
# ==========================================================================
def build_target_top100(minmax_df):
    """Compute target-level fold stats from Strategy 1 output."""
    valid = minmax_df.dropna(subset=["fold_difference"])
    tgt = (
        valid.groupby(["target_name", "target_chembl_id"])
        .agg(
            n_families=("stereochemical_family_id", "nunique"),
            mean_fold_difference=("fold_difference", "mean"),
            median_fold_difference=("fold_difference", "median"),
            max_fold_difference=("fold_difference", "max"),
            n_extreme_pairs=("classification", lambda x: (x == "Extreme").sum()),
        )
        .reset_index()
        .sort_values("mean_fold_difference", ascending=False)
        .head(100)
    )
    return tgt

# ==========================================================================
# REPORT WRITER
# ==========================================================================
def write_report(lines, stats, total_groups, s1_df, ref_reasons, elapsed):
    lines_out = []

    def p(s=""):
        lines_out.append(s)

    p("=" * 70)
    p("  All-Species v3 -- Family View Stereoselectivity Report (S1 + S2)")
    p(f"  Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    p("=" * 70)
    p()
    p("SOURCE")
    p(f"  File     : {SRC_FILE}")
    p(f"  Rows     : {stats['total_rows']:,}")
    p(f"  Compounds: {stats['unique_compounds']:,}")
    p(f"  Targets  : {stats['unique_targets']:,}")
    p(f"  Families : {stats['unique_families']:,}")
    p(f"  Groups (>=2 variants): {total_groups:,}")

    # Strategy 1
    p()
    p("=" * 70)
    p("STRATEGY 1: MinMax Fold")
    p("=" * 70)
    p(f"  Groups processed: {len(s1_df):,}")
    valid_s1 = s1_df.dropna(subset=["fold_difference"])
    p(f"  Valid pairs     : {len(valid_s1):,}")
    p()
    p("  Classification breakdown:")
    for cls in ["None","Mild","Moderate","Strong","Extreme","Undetermined"]:
        n = (s1_df["classification"] == cls).sum()
        pct = 100 * n / max(len(s1_df), 1)
        p(f"    {cls:<15s}: {n:>8,}  ({pct:.1f}%)")
    p()
    p("  Top 5 highest fold differences:")
    top5 = valid_s1.head(5)
    for _, r in top5.iterrows():
        p(f"    {str(r['stereochemical_family_id']):<20s}"
          f"| {str(r.get('compound_name',''))[:20]:<22s}"
          f"| {str(r['target_name'])[:28]:<30s}"
          f"| {r['activity_type']:<6s}"
          f"| fold={r['fold_difference']:.2g}")

    # Strategy 2
    p()
    p("=" * 70)
    p("STRATEGY 2: Reference Variant")
    p("=" * 70)
    p(f"  Total comparisons generated: {len(lines):,}")
    p()
    p("  Reference selection breakdown:")
    for reason, count in ref_reasons.items():
        p(f"    {reason:<20s}: {count:,}")

    report_path = os.path.join(REP_DIR, "11_family_view_report.txt")
    with open(report_path, "w") as f:
        f.write("\n".join(lines_out))
    print(f"\n  Report saved: {report_path}")
    return lines_out

# ==========================================================================
# MAIN
# ==========================================================================
def main():
    print("=" * 70)
    print("11_family_view.py -- Family-Level Stereoselectivity Analysis")
    print(f"  Output : {OUT_DIR}")
    print("=" * 70)
    t0 = time.time()

    # ---- Load & Aggregate ------------------------------------------------
    print("\n[Load] Reading and pre-aggregating TRUE_ISOMER_PAIRS.csv...")
    try:
        pf, stats, total_groups = load_and_aggregate()
    except Exception as e:
        print(f"  ERROR during load: {e}")
        sys.exit(1)

    # ---- Strategy 1: MinMax ----------------------------------------------
    print("\n[S1] MinMax Fold Difference...")
    try:
        s1_df = strategy1_minmax(pf)
        print(f"  Groups processed : {len(s1_df):,}")
        valid_s1 = s1_df.dropna(subset=["fold_difference"])
        print(f"  Valid (computable): {len(valid_s1):,}")
        print(f"  Median fold      : {valid_s1['fold_difference'].median():.2f}x")
        cdist = s1_df["classification"].value_counts()
        for cls, n in cdist.items():
            print(f"    {cls:<15s}: {n:,}")

        save(s1_df, os.path.join(OUT_DIR, "FAMILY_MINMAX_FOLD.csv"))
        save_by_metric(s1_df, "MINMAX_FOLD")
        print(f"  Saved FAMILY_MINMAX_FOLD.csv ({len(s1_df):,} rows)")
    except Exception as e:
        print(f"  ERROR in Strategy 1: {e}")
        s1_df = pd.DataFrame()

    # ---- Strategy 2: Reference -------------------------------------------
    print("\n[S2] Reference Variant Analysis...")
    try:
        s2_df, ref_reasons = strategy2_reference(pf)
        print(f"  Comparisons generated: {len(s2_df):,}")
        print(f"  Reference selection reasons:")
        for reason, count in ref_reasons.items():
            print(f"    {reason:<20s}: {count:,}")

        save(s2_df, os.path.join(OUT_DIR, "FAMILY_REFERENCE_FOLD.csv"))
        save_by_metric(s2_df, "REFERENCE_FOLD")
        print(f"  Saved FAMILY_REFERENCE_FOLD.csv ({len(s2_df):,} rows)")
    except Exception as e:
        print(f"  ERROR in Strategy 2: {e}")
        s2_df, ref_reasons = pd.DataFrame(), {}

    #

    # ---- Target Top-100 --------------------------------------------------
    print("\n[Targets] Building TOP100_STEREOSELECTIVE_TARGETS.csv...")
    try:
        top100_tgt = build_target_top100(s1_df)
        tgt_path = os.path.join(SEL_DIR, "TOP100_STEREOSELECTIVE_TARGETS.csv")
        save(top100_tgt, tgt_path)
        print(f"  Saved TOP100_STEREOSELECTIVE_TARGETS.csv ({len(top100_tgt):,} targets)")
        print(f"\n  Top 5 targets by mean fold difference:")
        for _, r in top100_tgt.head(5).iterrows():
            print(f"    {str(r['target_name'])[:45]:<47s}"
                  f"mean={r['mean_fold_difference']:.2g}x  "
                  f"max={r['max_fold_difference']:.2g}x  "
                  f"extreme_n={int(r['n_extreme_pairs'])}")
    except Exception as e:
        print(f"  ERROR building target table: {e}")
        top100_tgt = pd.DataFrame()

    # ---- Manifest --------------------------------------------------------
    print(f"\n[Files] Output manifest (family_views/):")
    for f in sorted(os.listdir(OUT_DIR)):
        if os.path.isfile(os.path.join(OUT_DIR, f)):
            sz = os.path.getsize(os.path.join(OUT_DIR, f)) / 1e6
            print(f"  {f:<50s} {sz:.1f} MB")
    by_m_files = [f for f in os.listdir(METRIC_DIR) if f.endswith(".csv")]
    print(f"  by_metric/: {len(by_m_files)} files")

    # ---- Report ----------------------------------------------------------
    elapsed = time.time() - t0
    print(f"\n[Report] Writing 11_family_view_report.txt...")
    try:
        write_report(
            lines        = [] if s2_df.empty else list(range(len(s2_df))),
            stats        = stats,
            total_groups = total_groups,
            s1_df        = s1_df,
            ref_reasons  = ref_reasons,
            elapsed      = elapsed,
        )
    except Exception as e:
        print(f"  ERROR writing report: {e}")

    print(f"\n{'='*70}")
    print(f"  Total runtime : {elapsed:.1f}s ({elapsed/60:.1f} min)")
    print(f"  Output dir    : {OUT_DIR}")
    print(f"  [OK] 11_family_view.py complete.")
    print(f"{'='*70}")

if __name__ == "__main__":
    main()
