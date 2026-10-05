"""
All-Species v3 — Step 5: True Isomer Datasets
==============================================
Replicates & extends original workflow's TRUE_ISOMERS / ONLY_R / ONLY_S split.
Groups by (smiles_nostereo, target_name, activity_type, organism).

Outputs in outputs/true_isomers/:
  TRUE_ISOMER_PAIRS.csv    — families with BOTH Stereo_1 AND Stereo_2 present
  ONLY_VARIANT_1.csv       — families with ONLY Stereo_1
  ONLY_VARIANT_2.csv       — families with ONLY Stereo_2
  MULTI_CENTER_STEREO.csv  — compounds with n_stereocenters > 2 in true-isomer families
"""
import pandas as pd
import os, sys
sys.path.insert(0, os.path.dirname(__file__))
from config import *

GROUP = ["smiles_nostereo", "target_name", "activity_type", "organism"]

def split_isomers(df):
    true_idx, s1_idx, s2_idx = [], [], []
    for _, grp in df.groupby(GROUP):
        forms = set(grp["stereo_form"].dropna().unique())
        if "Stereo_1" in forms and "Stereo_2" in forms:
            true_idx.extend(grp.index)
        elif forms == {"Stereo_1"}:
            s1_idx.extend(grp.index)
        elif forms == {"Stereo_2"}:
            s2_idx.extend(grp.index)
    return df.loc[true_idx], df.loc[s1_idx], df.loc[s2_idx]

def main():
    print("=" * 60)
    print("Step 5: True Isomer Datasets")
    print("=" * 60)

    print(f"\nLoading master aggregated CSV...")
    df = pd.read_csv(SRC_AGGREGATED, low_memory=False)
    print(f"  {len(df):,} rows")

    cols = [c for c in BASE_COLS if c in df.columns]

    print(f"  Splitting into TRUE / ONLY_1 / ONLY_2...")
    true_df, s1_df, s2_df = split_isomers(df)

    # ── TRUE_ISOMER_PAIRS ──────────────────────────────────────────────
    true_df[cols].to_csv(os.path.join(ISOMER_DIR, "TRUE_ISOMER_PAIRS.csv"), index=False)
    print(f"\n  TRUE_ISOMER_PAIRS.csv      : {len(true_df):,} rows")
    print(f"    Unique scaffolds           : {true_df['smiles_nostereo'].nunique():,}")
    print(f"    Unique targets             : {true_df['target_name'].nunique():,}")
    print(f"    Unique organisms           : {true_df['organism'].nunique():,}")
    print(f"    Activity types represented : {true_df['activity_type'].nunique():,}")

    # ── ONLY_VARIANT_1 ─────────────────────────────────────────────────
    s1_df[cols].to_csv(os.path.join(ISOMER_DIR, "ONLY_VARIANT_1.csv"), index=False)
    print(f"\n  ONLY_VARIANT_1.csv         : {len(s1_df):,} rows")

    # ── ONLY_VARIANT_2 ─────────────────────────────────────────────────
    s2_df[cols].to_csv(os.path.join(ISOMER_DIR, "ONLY_VARIANT_2.csv"), index=False)
    print(f"  ONLY_VARIANT_2.csv         : {len(s2_df):,} rows")

    # ── MULTI_CENTER_STEREO — true pairs with n_stereocenters > 2 ──────
    multi = true_df[true_df["n_stereocenters"] > 2][cols]
    multi.to_csv(os.path.join(ISOMER_DIR, "MULTI_CENTER_STEREO.csv"), index=False)
    print(f"\n  MULTI_CENTER_STEREO.csv    : {len(multi):,} rows (n_centers > 2)")
    if "n_stereocenters" in multi.columns:
        print(f"    Stereocenter distribution:")
        for n, cnt in multi["n_stereocenters"].value_counts().sort_index().head(10).items():
            print(f"      {n} centers: {cnt:,}")

    # ── Per-metric TRUE isomer files ───────────────────────────────────
    print(f"\n  Generating per-metric TRUE isomer files...")
    for metric in STEREO_METRICS:
        sub = true_df[
            true_df["activity_type"].str.strip().str.upper() == metric.upper()
        ][cols]
        if len(sub) == 0:
            continue
        sub.to_csv(os.path.join(ISOMER_DIR, f"{metric}_TRUE_ISOMERS.csv"), index=False)
        print(f"    {metric}_TRUE_ISOMERS.csv — {len(sub):,} rows, "
              f"{sub['organism'].nunique()} organisms")

    # ── Summary table ──────────────────────────────────────────────────
    isomer_summary = pd.DataFrame({
        "category":     ["TRUE_ISOMER_PAIRS", "ONLY_VARIANT_1", "ONLY_VARIANT_2", "MULTI_CENTER"],
        "records":      [len(true_df), len(s1_df), len(s2_df), len(multi)],
        "scaffolds":    [true_df["smiles_nostereo"].nunique(),
                         s1_df["smiles_nostereo"].nunique(),
                         s2_df["smiles_nostereo"].nunique(),
                         multi["smiles_nostereo"].nunique()],
        "organisms":    [true_df["organism"].nunique(), s1_df["organism"].nunique(),
                         s2_df["organism"].nunique(), multi["organism"].nunique()],
        "targets":      [true_df["target_name"].nunique(), s1_df["target_name"].nunique(),
                         s2_df["target_name"].nunique(), multi["target_name"].nunique()],
    })
    isomer_summary.to_csv(os.path.join(REPORTS_DIR, "true_isomer_summary.csv"), index=False)
    print(f"\n  Saved reports/true_isomer_summary.csv")

    print(f"\n✅ Step 5 complete.")

if __name__ == "__main__":
    main()
