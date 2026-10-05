"""
All-Species v3 — Post-processing: Stereochemical Classification
================================================================
Reads master aggregated CSV from v3 (no DB re-extraction).
Classifies every record into four mutually-exclusive-for-single-center
(and overlapping-for-multi) categories:

  TRUE_ISOMERS  — group (smiles_nostereo + target + activity_type + organism)
                  has BOTH an "R" variant AND an "S" variant observed.
                  Also applies to multi-center families with ≥2 stereoforms.

  ONLY_R        — single-center compound, only "R" variant in its group.
                  Not in any TRUE_ISOMER group.

  ONLY_S        — single-center compound, only "S" variant in its group.
                  Not in any TRUE_ISOMER group.

  MULTI_CENTER  — n_stereocenters > 1 (R,R / R,S / S,R / S,S / R,R,S …).
                  Note: multi-center compounds can also appear in TRUE_ISOMERS
                  if their scaffold has multiple observed stereoforms.

Outputs:
  outputs/core_metrics/<METRIC>/
      TRUE_ISOMERS.csv
      ONLY_R.csv
      ONLY_S.csv
      MULTI_CENTER.csv

  outputs/species/<organism_name>/
      TRUE_ISOMERS.csv
      ONLY_R.csv
      ONLY_S.csv
      MULTI_CENTER.csv

  reports/
      metric_stereo_classification_summary.csv
      species_stereo_classification_summary.csv
"""

import pandas as pd
import re, os, sys, time
from pathlib import Path

# ── Paths ──────────────────────────────────────────────────────────────────
V3_DIR   = os.path.dirname(os.path.abspath(__file__))
AGG_FILE = os.path.join(V3_DIR, "aggregated", "ALL_SPECIES_ALL_METRICS_AGGREGATED.csv")
CORE_OUT = os.path.join(V3_DIR, "outputs", "core_metrics")
SPEC_OUT = os.path.join(V3_DIR, "outputs", "species")
REP_OUT  = os.path.join(V3_DIR, "reports")

# Core metrics to process (dedicated subfolders)
CORE_METRICS = ["IC50", "Ki", "Kd", "EC50", "Km", "Vmax", "Kcat", "Kinact",
                "Potency", "Emax", "Efficacy", "Activity", "Inhibition",
                "AC50", "GI50", "ED50"]

# Grouping key for TRUE ISOMER detection (same as original workflow)
GROUP_COLS = ["smiles_nostereo", "target_name", "activity_type", "organism"]

# ── Utilities ──────────────────────────────────────────────────────────────
def safe_name(name: str) -> str:
    """Filesystem-safe name."""
    s = re.sub(r"[^\w\s-]", "", str(name))
    s = re.sub(r"[\s/]+", "_", s.strip())
    return s[:80]

def save(df: pd.DataFrame, path: str):
    """Save DataFrame, create parent dirs, skip if empty."""
    if len(df) == 0:
        return
    os.makedirs(os.path.dirname(path), exist_ok=True)
    df.to_csv(path, index=False)

# ── Classification Logic ───────────────────────────────────────────────────
def classify(df: pd.DataFrame) -> pd.DataFrame:
    """
    Add a 'stereo_class' column to df.
    Logic (applied in priority order):
      MULTI_CENTER  : n_stereocenters > 1  (comma in variant)
      TRUE_ISOMERS  : single-center group with BOTH R and S present
                      OR multi-center family with ≥2 distinct stereo_form
      ONLY_R        : single-center, variant == 'R', no S counterpart
      ONLY_S        : single-center, variant == 'S', no R counterpart
    """
    df = df.copy()
    df["stereo_class"] = "UNCATEGORIZED"

    # ── Flag multi-center rows ─────────────────────────────────────────
    is_multi = df["n_stereocenters"] > 1
    df.loc[is_multi, "stereo_class"] = "MULTI_CENTER"

    # ── Identify TRUE_ISOMER groups (single-center) ────────────────────
    # A group is a TRUE pair if it has variant "R" AND variant "S"
    single_center = df[~is_multi].copy()

    if len(single_center) > 0:
        # For each group, collect the set of variants present
        group_variants = (
            single_center
            .groupby(GROUP_COLS)["stereochemical_variant"]
            .apply(lambda x: frozenset(str(v).strip() for v in x))
        )
        # TRUE if both "R" and "S" in the group
        true_groups = group_variants[
            group_variants.apply(lambda s: "R" in s and "S" in s)
        ].index

        # Build a set of (smiles_nostereo, target, activity, organism) tuples
        true_set = set(true_groups.tolist())

        def is_true(row):
            key = (row["smiles_nostereo"], row["target_name"],
                   row["activity_type"], row["organism"])
            return key in true_set

        true_mask_single = single_center.apply(is_true, axis=1)
        true_idx = single_center[true_mask_single].index
        df.loc[true_idx, "stereo_class"] = "TRUE_ISOMERS"

        # Remaining single-center that are not TRUE
        remaining = single_center[~true_mask_single]
        only_r_idx = remaining[
            remaining["stereochemical_variant"].str.strip() == "R"
        ].index
        only_s_idx = remaining[
            remaining["stereochemical_variant"].str.strip() == "S"
        ].index
        df.loc[only_r_idx, "stereo_class"] = "ONLY_R"
        df.loc[only_s_idx, "stereo_class"] = "ONLY_S"

    # ── TRUE ISOMERS for multi-center: families with ≥2 stereo_forms ──
    # (multi-center compounds whose family has both Stereo_1 and Stereo_2)
    multi_center_df = df[is_multi].copy()
    if len(multi_center_df) > 0:
        fam_forms = (
            multi_center_df
            .groupby("stereochemical_family_id")["stereo_form"]
            .nunique()
        )
        multi_true_families = fam_forms[fam_forms >= 2].index
        multi_true_idx = multi_center_df[
            multi_center_df["stereochemical_family_id"].isin(multi_true_families)
        ].index
        # Mark as both MULTI_CENTER and TRUE — keep a combined label
        df.loc[multi_true_idx, "stereo_class"] = "MULTI_CENTER_TRUE"

    return df

# ── Output writers ─────────────────────────────────────────────────────────
CATEGORY_FILES = {
    "TRUE_ISOMERS":      "TRUE_ISOMERS.csv",
    "ONLY_R":            "ONLY_R.csv",
    "ONLY_S":            "ONLY_S.csv",
    "MULTI_CENTER":      "MULTI_CENTER.csv",
    "MULTI_CENTER_TRUE": "MULTI_CENTER.csv",  # multi-center true → MULTI_CENTER.csv
}

def write_splits(classified: pd.DataFrame, out_dir: str) -> dict:
    """Write TRUE_ISOMERS/ONLY_R/ONLY_S/MULTI_CENTER.csv. Return counts."""
    counts = {"TRUE_ISOMERS": 0, "ONLY_R": 0, "ONLY_S": 0, "MULTI_CENTER": 0}
    os.makedirs(out_dir, exist_ok=True)

    for category, fname in [
        ("TRUE_ISOMERS",       "TRUE_ISOMERS.csv"),
        ("ONLY_R",             "ONLY_R.csv"),
        ("ONLY_S",             "ONLY_S.csv"),
    ]:
        sub = classified[classified["stereo_class"] == category]
        save(sub, os.path.join(out_dir, fname))
        counts[category] = len(sub)

    # MULTI_CENTER = MULTI_CENTER + MULTI_CENTER_TRUE combined
    multi = classified[classified["stereo_class"].isin(["MULTI_CENTER", "MULTI_CENTER_TRUE"])]
    save(multi, os.path.join(out_dir, "MULTI_CENTER.csv"))
    counts["MULTI_CENTER"] = len(multi)

    # TRUE_ISOMERS supplement: also include MULTI_CENTER_TRUE in TRUE_ISOMERS.csv
    true_all = classified[classified["stereo_class"].isin(["TRUE_ISOMERS", "MULTI_CENTER_TRUE"])]
    if len(true_all) > len(classified[classified["stereo_class"] == "TRUE_ISOMERS"]):
        save(true_all, os.path.join(out_dir, "TRUE_ISOMERS.csv"))
        counts["TRUE_ISOMERS"] = len(true_all)

    return counts

# ── Main ───────────────────────────────────────────────────────────────────
def main():
    print("=" * 65)
    print("Post-processing: Stereochemical Classification")
    print("  TRUE_ISOMERS | ONLY_R | ONLY_S | MULTI_CENTER")
    print("=" * 65)

    t0 = time.time()
    print(f"\nLoading: {AGG_FILE}")
    df = pd.read_csv(AGG_FILE, low_memory=False)
    print(f"  {len(df):,} rows, {df['organism'].nunique()} organisms")

    # ── Global classification ──────────────────────────────────────────
    print(f"\nClassifying all {len(df):,} records...")
    classified = classify(df)

    class_counts = classified["stereo_class"].value_counts()
    print(f"\n  Classification results:")
    for cls, cnt in class_counts.items():
        print(f"    {cls:<22s}: {cnt:>10,}")

    # ══════════════════════════════════════════════════════════════════
    # PART A: Per-metric subfolders in outputs/core_metrics/<METRIC>/
    # ══════════════════════════════════════════════════════════════════
    print(f"\n{'─'*65}")
    print(f"PART A: Per-metric classification files")
    print(f"{'─'*65}")

    metric_summary_rows = []

    for metric in CORE_METRICS:
        mask = classified["activity_type"].str.strip().str.upper() == metric.upper()
        sub = classified[mask]
        if len(sub) == 0:
            continue

        out_dir = os.path.join(CORE_OUT, metric)
        counts = write_splits(sub, out_dir)

        row = {
            "metric":              metric,
            "total_records":       len(sub),
            "compound_count":      sub["compound_chembl_id"].nunique(),
            "target_count":        sub["target_name"].nunique(),
            "organism_count":      sub["organism"].nunique(),
            "TRUE_ISOMER_count":   counts["TRUE_ISOMERS"],
            "ONLY_R_count":        counts["ONLY_R"],
            "ONLY_S_count":        counts["ONLY_S"],
            "MULTI_CENTER_count":  counts["MULTI_CENTER"],
        }
        metric_summary_rows.append(row)
        print(f"  {metric:<12s}: TRUE={counts['TRUE_ISOMERS']:>8,}  "
              f"R={counts['ONLY_R']:>8,}  S={counts['ONLY_S']:>8,}  "
              f"MULTI={counts['MULTI_CENTER']:>8,}")

    # Also process ALL remaining activity types not in CORE_METRICS
    covered_upper = {m.upper() for m in CORE_METRICS}
    others = classified[~classified["activity_type"].str.strip().str.upper().isin(covered_upper)]
    if len(others) > 0:
        out_dir = os.path.join(CORE_OUT, "OTHER")
        counts = write_splits(others, out_dir)
        row = {
            "metric": "OTHER",
            "total_records": len(others),
            "compound_count": others["compound_chembl_id"].nunique(),
            "target_count": others["target_name"].nunique(),
            "organism_count": others["organism"].nunique(),
            "TRUE_ISOMER_count": counts["TRUE_ISOMERS"],
            "ONLY_R_count": counts["ONLY_R"],
            "ONLY_S_count": counts["ONLY_S"],
            "MULTI_CENTER_count": counts["MULTI_CENTER"],
        }
        metric_summary_rows.append(row)
        print(f"  {'OTHER':<12s}: TRUE={counts['TRUE_ISOMERS']:>8,}  "
              f"R={counts['ONLY_R']:>8,}  S={counts['ONLY_S']:>8,}  "
              f"MULTI={counts['MULTI_CENTER']:>8,}")

    # Save metric summary
    metric_summary = pd.DataFrame(metric_summary_rows)
    metric_sum_path = os.path.join(REP_OUT, "metric_stereo_classification_summary.csv")
    metric_summary.to_csv(metric_sum_path, index=False)
    print(f"\n  ✓ Saved metric_stereo_classification_summary.csv")

    # ══════════════════════════════════════════════════════════════════
    # PART B: Per-species subfolders in outputs/species/<organism>/
    # ══════════════════════════════════════════════════════════════════
    print(f"\n{'─'*65}")
    print(f"PART B: Per-species classification files")
    print(f"{'─'*65}")

    organisms = sorted(classified["organism"].dropna().unique())
    print(f"  Generating classification files for {len(organisms)} organisms...")

    species_summary_rows = []
    for i, org in enumerate(organisms):
        sub = classified[classified["organism"] == org]
        org_dir = os.path.join(SPEC_OUT, safe_name(org))
        counts = write_splits(sub, org_dir)

        row = {
            "organism":            org,
            "total_records":       len(sub),
            "compound_count":      sub["compound_chembl_id"].nunique(),
            "target_count":        sub["target_name"].nunique(),
            "activity_types":      sub["activity_type"].nunique(),
            "TRUE_ISOMER_count":   counts["TRUE_ISOMERS"],
            "ONLY_R_count":        counts["ONLY_R"],
            "ONLY_S_count":        counts["ONLY_S"],
            "MULTI_CENTER_count":  counts["MULTI_CENTER"],
        }
        species_summary_rows.append(row)

        if (i + 1) % 50 == 0 or (i + 1) == len(organisms):
            elapsed = time.time() - t0
            print(f"  {i+1}/{len(organisms)} organisms done ({elapsed:.0f}s elapsed)", end="\r")

    print(f"\n  All {len(organisms)} organism folders written.")

    # Save species summary
    species_summary = (
        pd.DataFrame(species_summary_rows)
        .sort_values("total_records", ascending=False)
    )
    spec_sum_path = os.path.join(REP_OUT, "species_stereo_classification_summary.csv")
    species_summary.to_csv(spec_sum_path, index=False)
    print(f"  ✓ Saved species_stereo_classification_summary.csv")

    # ── Print top species ──────────────────────────────────────────────
    print(f"\n  Top 10 organisms by TRUE_ISOMER count:")
    top = species_summary.nlargest(10, "TRUE_ISOMER_count")[
        ["organism", "total_records", "TRUE_ISOMER_count", "ONLY_R_count",
         "ONLY_S_count", "MULTI_CENTER_count"]
    ]
    for _, r in top.iterrows():
        print(f"    {str(r['organism'])[:40]:<42s}: "
              f"TRUE={int(r['TRUE_ISOMER_count']):>7,}  "
              f"R={int(r['ONLY_R_count']):>7,}  "
              f"S={int(r['ONLY_S_count']):>7,}  "
              f"MULTI={int(r['MULTI_CENTER_count']):>7,}")

    # ── Grand summary ──────────────────────────────────────────────────
    total = time.time() - t0
    print(f"\n{'─'*65}")
    print(f"SUMMARY")
    print(f"{'─'*65}")
    print(f"  Core metric subfolders  : {len(CORE_METRICS)+1} (+ OTHER)")
    print(f"  Species subfolders      : {len(organisms)}")
    print(f"  Total TRUE_ISOMERS rows : {(classified['stereo_class'].isin(['TRUE_ISOMERS','MULTI_CENTER_TRUE'])).sum():,}")
    print(f"  Total ONLY_R rows       : {(classified['stereo_class'] == 'ONLY_R').sum():,}")
    print(f"  Total ONLY_S rows       : {(classified['stereo_class'] == 'ONLY_S').sum():,}")
    print(f"  Total MULTI_CENTER rows : {(classified['stereo_class'].isin(['MULTI_CENTER','MULTI_CENTER_TRUE'])).sum():,}")
    print(f"  Total time              : {total:.1f}s ({total/60:.1f} min)")
    print(f"\n✅ Post-processing complete.")

if __name__ == "__main__":
    main()
