"""
All-Species v3 — Step 4: Stereoisomer Files per Metric
=======================================================
For each metric in STEREO_METRICS:
  → outputs/stereoisomers/{METRIC}_stereoisomers.csv
     Only compounds that belong to stereochemical families
     (i.e., stereochemical_family_id has ≥2 distinct compounds).

Also generates:
  → outputs/stereoisomers/ALL_METRICS_STEREOISOMERS.csv (master)
"""
import pandas as pd
import os, sys
sys.path.insert(0, os.path.dirname(__file__))
from config import *

def filter_to_families(df: pd.DataFrame) -> pd.DataFrame:
    """Keep only rows belonging to families with ≥2 unique stereo_forms."""
    family_forms = (
        df.groupby("stereochemical_family_id")["stereo_form"]
        .nunique()
    )
    multi_families = family_forms[family_forms >= 2].index
    return df[df["stereochemical_family_id"].isin(multi_families)].copy()

def main():
    print("=" * 60)
    print("Step 4: Stereoisomer files per metric")
    print("=" * 60)

    print(f"\nLoading master aggregated CSV...")
    df = pd.read_csv(SRC_AGGREGATED, low_memory=False)
    print(f"  {len(df):,} rows")

    cols = [c for c in STEREO_COLS if c in df.columns]

    # ── Build master stereoisomer set (families with ≥2 stereoforms) ──
    print(f"  Identifying multi-stereoform families...")
    family_forms = (
        df.groupby("stereochemical_family_id")["stereo_form"].nunique()
    )
    multi_families = family_forms[family_forms >= 2].index
    stereo_df = df[df["stereochemical_family_id"].isin(multi_families)].copy()
    print(f"  Families with ≥2 stereoforms: {len(multi_families):,}")
    print(f"  Records in those families: {len(stereo_df):,}")

    # Save master
    stereo_df[cols].to_csv(
        os.path.join(STEREO_DIR, "ALL_METRICS_STEREOISOMERS.csv"), index=False)
    print(f"  Saved ALL_METRICS_STEREOISOMERS.csv")

    # ── Per-metric stereoisomer files ──────────────────────────────────
    print(f"\n  Generating per-metric stereoisomer files...")
    for metric in STEREO_METRICS:
        mask = stereo_df["activity_type"].str.strip().str.upper() == metric.upper()
        sub = stereo_df[mask][cols]
        if len(sub) == 0:
            continue
        fname = f"{metric}_stereoisomers.csv"
        sub.to_csv(os.path.join(STEREO_DIR, fname), index=False)
        n_fam = sub["stereochemical_family_id"].nunique()
        n_org = sub["organism"].nunique()
        print(f"    {fname} — {len(sub):,} rows, {n_fam:,} families, {n_org} organisms")

    # ── Also generate per-organism stereoisomer summary ────────────────
    print(f"\n  Generating organism stereoisomer summary...")
    org_stereo = (
        stereo_df.groupby("organism")
        .agg(
            stereo_records=("compound_chembl_id", "count"),
            stereo_families=("stereochemical_family_id", "nunique"),
            stereo_compounds=("compound_chembl_id", "nunique"),
            stereo_targets=("target_name", "nunique"),
        )
        .reset_index()
        .sort_values("stereo_records", ascending=False)
    )
    org_stereo.to_csv(
        os.path.join(REPORTS_DIR, "organism_stereo_summary.csv"), index=False)
    print(f"  Saved reports/organism_stereo_summary.csv ({len(org_stereo):,} organisms)")

    print(f"\n✅ Step 4 complete.")

if __name__ == "__main__":
    main()
