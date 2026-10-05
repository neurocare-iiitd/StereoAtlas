"""
All-Species v3 — Step 1: Prepare Base Data
===========================================
Reuses ALL_SPECIES_ALL_METRICS_AGGREGATED.csv from v2.
Copies raw and aggregated CSVs into v3 folder structure.
No DB re-extraction needed.
"""
import pandas as pd
import shutil, os, sys
sys.path.insert(0, os.path.dirname(__file__))
from config import *

def main():
    print("=" * 60)
    print("Step 1: Prepare base data (reuse v2 outputs)")
    print("=" * 60)

    # ── Verify sources exist ───────────────────────────────────────────
    for src in [SRC_AGGREGATED, SRC_RAW_STEREO]:
        if not os.path.exists(src):
            print(f"❌ Missing: {src}")
            sys.exit(1)
        print(f"  ✓ Found: {os.path.basename(src)}")

    # ── Copy raw stereo file → raw/ ────────────────────────────────────
    raw_dest = os.path.join(RAW_DIR, "ALL_SPECIES_ALL_METRICS_RAW.csv")
    if not os.path.exists(raw_dest):
        print(f"\n  Copying raw stereo → raw/...")
        shutil.copy2(SRC_RAW_STEREO, raw_dest)
    else:
        print(f"\n  raw/ALL_SPECIES_ALL_METRICS_RAW.csv already exists, skipping copy.")

    # ── Copy aggregated → aggregated/ ─────────────────────────────────
    agg_dest = os.path.join(AGG_DIR, "ALL_SPECIES_ALL_METRICS_AGGREGATED.csv")
    if not os.path.exists(agg_dest):
        print(f"  Copying aggregated → aggregated/...")
        shutil.copy2(SRC_AGGREGATED, agg_dest)
    else:
        print(f"  aggregated/ALL_SPECIES_ALL_METRICS_AGGREGATED.csv already exists, skipping copy.")

    # ── Quick stats ────────────────────────────────────────────────────
    print(f"\n  Loading aggregated data for verification...")
    df = pd.read_csv(SRC_AGGREGATED, low_memory=False)
    print(f"  Rows           : {len(df):,}")
    print(f"  Organisms      : {df['organism'].nunique():,}")
    print(f"  Compounds      : {df['compound_chembl_id'].nunique():,}")
    print(f"  Targets        : {df['target_name'].nunique():,}")
    print(f"  Activity types : {df['activity_type'].nunique():,}")
    print(f"  Stereo families: {df['stereochemical_family_id'].nunique():,}")
    print(f"\n✅ Step 1 complete.")
    return df

if __name__ == "__main__":
    main()
