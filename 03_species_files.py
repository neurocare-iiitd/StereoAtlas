"""
All-Species v3 — Step 3: Species Files
=======================================
Generates one CSV per organism in outputs/species/.
Also generates species_summary.csv.
"""
import pandas as pd
import re, os, sys
sys.path.insert(0, os.path.dirname(__file__))
from config import *

def safe_filename(name: str) -> str:
    """Convert organism name to a safe filename."""
    safe = re.sub(r"[^\w\s-]", "", str(name))
    safe = re.sub(r"[\s/]+", "_", safe.strip())
    return safe[:80]  # cap length

def main():
    print("=" * 60)
    print("Step 3: Species Files")
    print("=" * 60)

    print(f"\nLoading master aggregated CSV...")
    df = pd.read_csv(SRC_AGGREGATED, low_memory=False)
    print(f"  {len(df):,} rows, {df['organism'].nunique():,} organisms")

    cols = [c for c in BASE_COLS if c in df.columns]

    organisms = df["organism"].dropna().unique()
    print(f"\n  Generating {len(organisms):,} species files...")

    summary_rows = []
    for i, org in enumerate(sorted(organisms)):
        sub = df[df["organism"] == org][cols].copy()
        fname = safe_filename(org) + ".csv"
        sub.to_csv(os.path.join(SPECIES_DIR, fname), index=False)

        summary_rows.append({
            "organism":          org,
            "filename":          fname,
            "activity_records":  len(sub),
            "compound_count":    sub["compound_chembl_id"].nunique(),
            "target_count":      sub["target_name"].nunique(),
            "activity_types":    sub["activity_type"].nunique(),
            "stereo_family_count": sub["stereochemical_family_id"].nunique(),
            "n_measurements_total": sub["n_measurements"].sum(),
        })

        if (i + 1) % 50 == 0 or (i + 1) == len(organisms):
            print(f"  {i+1}/{len(organisms)} species processed...", end="\r")

    print(f"\n  All {len(organisms):,} species files written.")

    # ── species_summary.csv ───────────────────────────────────────────
    summary = pd.DataFrame(summary_rows).sort_values("activity_records", ascending=False)
    out_path = os.path.join(REPORTS_DIR, "species_summary.csv")
    summary.to_csv(out_path, index=False)
    print(f"\n  Saved reports/species_summary.csv ({len(summary):,} rows)")
    print(f"\n  Top 10 organisms by record count:")
    print(summary[["organism", "activity_records", "compound_count", "target_count"]].head(10).to_string(index=False))

    print(f"\n✅ Step 3 complete.")

if __name__ == "__main__":
    main()
