"""
All-Species v3 — Step 6: Stereochemistry Summary Report
========================================================
Generates reports/stereochemistry_summary.csv:
  Per-family summary: family_id, core_smiles, n_stereocenters,
  observed_variants, family_size, target_count, organism_count, activity_types

Reuses STEREOCHEMICAL_FAMILIES.csv from v2 if present,
otherwise rebuilds from master aggregated.
"""
import pandas as pd
import os, sys
sys.path.insert(0, os.path.dirname(__file__))
from config import *

V2_FAMILIES = os.path.join(V2_DIR, "STEREOCHEMICAL_FAMILIES.csv")
V2_SUMMARY  = os.path.join(V2_DIR, "STEREOCHEMISTRY_SUMMARY.csv")

def main():
    print("=" * 60)
    print("Step 6: Stereochemistry Summary")
    print("=" * 60)

    # Prefer v2 output if available (saves time)
    if os.path.exists(V2_FAMILIES):
        print(f"\n  Reusing STEREOCHEMICAL_FAMILIES.csv from v2...")
        fam = pd.read_csv(V2_FAMILIES, low_memory=False)
        import shutil
        shutil.copy2(V2_FAMILIES, os.path.join(REPORTS_DIR, "STEREOCHEMICAL_FAMILIES.csv"))
        if os.path.exists(V2_SUMMARY):
            shutil.copy2(V2_SUMMARY, os.path.join(REPORTS_DIR, "STEREOCHEMISTRY_SUMMARY.csv"))
        print(f"  Copied {len(fam):,} families → reports/")
    else:
        print(f"\n  Building from master aggregated CSV...")
        df = pd.read_csv(SRC_AGGREGATED, low_memory=False)
        fam = (
            df.groupby(["stereochemical_family_id", "smiles_nostereo", "n_stereocenters"])
            .agg(
                family_size=("compound_chembl_id", "nunique"),
                observed_variants=("stereochemical_variant",
                                   lambda x: "|".join(sorted(set(str(v) for v in x if pd.notna(v))))),
                n_stereoforms=("stereo_form", "nunique"),
                target_count=("target_name", "nunique"),
                organism_count=("organism", "nunique"),
                activity_type_count=("activity_type", "nunique"),
                total_records=("mean_activity_value", "count"),
            )
            .reset_index()
            .sort_values("family_size", ascending=False)
        )
        fam.insert(0, "family_id",
                   ["Family_{:05d}".format(i+1) for i in range(len(fam))])
        fam.to_csv(os.path.join(REPORTS_DIR, "STEREOCHEMICAL_FAMILIES.csv"), index=False)
        print(f"  Saved {len(fam):,} families → reports/STEREOCHEMICAL_FAMILIES.csv")

    # ── Stereochemistry summary (per variant per family) ──────────────
    print(f"\n  Loading master aggregated for per-variant summary...")
    df = pd.read_csv(SRC_AGGREGATED, low_memory=False)
    summary = (
        df.groupby([
            "stereochemical_family_id", "smiles_nostereo",
            "n_stereocenters", "stereochemical_variant", "stereo_form"
        ])
        .agg(
            family_size=("compound_chembl_id", "nunique"),
            target_count=("target_name", "nunique"),
            organism_count=("organism", "nunique"),
            activity_type_count=("activity_type", "nunique"),
            n_records=("mean_activity_value", "count"),
        )
        .reset_index()
        .sort_values(["stereochemical_family_id", "stereo_form"])
    )
    out = os.path.join(REPORTS_DIR, "stereochemistry_summary.csv")
    summary.to_csv(out, index=False)
    print(f"  Saved reports/stereochemistry_summary.csv ({len(summary):,} rows)")

    # ── Stats ──────────────────────────────────────────────────────────
    print(f"\n  Stereocenter distribution (families):")
    if "n_stereocenters" in fam.columns:
        dist = fam["n_stereocenters"].value_counts().sort_index()
        for n, cnt in dist.head(10).items():
            print(f"    {n} centers: {cnt:,} families")

    print(f"\n✅ Step 6 complete.")

if __name__ == "__main__":
    main()
