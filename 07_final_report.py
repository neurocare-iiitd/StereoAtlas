"""
All-Species v3 — Step 7: Final Report
======================================
Generates reports/FINAL_REPORT.txt with:
  - Overall counts (records, compounds, targets, organisms, families)
  - Counts per metric category
  - Counts per organism (top 50)
  - Counts per target (top 50)
  - Stereochemical class breakdown
  - Output file manifest with sizes
"""
import pandas as pd
import os, sys
from datetime import datetime
sys.path.insert(0, os.path.dirname(__file__))
from config import *

def fmt(n): return f"{n:,}" if isinstance(n, (int, float)) and not pd.isna(n) else str(n)

def dir_manifest(directory: str, label: str):
    """Return list of (label, filename, size_MB) for all CSVs in a directory."""
    rows = []
    if not os.path.isdir(directory):
        return rows
    for f in sorted(os.listdir(directory)):
        if f.endswith(".csv"):
            path = os.path.join(directory, f)
            rows.append((label, f, os.path.getsize(path) / 1e6))
    return rows

def main():
    print("=" * 60)
    print("Step 7: Final Report")
    print("=" * 60)

    lines = []
    def p(s=""):
        lines.append(s)
        print(s)

    p("=" * 72)
    p("  ChEMBL All-Species Stereochemistry Atlas v3 — FINAL REPORT")
    p(f"  Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    p("=" * 72)

    # Load master
    print(f"\nLoading master aggregated CSV for report...")
    df = pd.read_csv(SRC_AGGREGATED, low_memory=False)

    # ── 1. OVERALL ────────────────────────────────────────────────────
    p()
    p("── 1. OVERALL COUNTS ──────────────────────────────────────────────────")
    p(f"  Source               : ChEMBL (chembl_latest, all SINGLE PROTEIN targets)")
    p(f"  Aggregated records   : {fmt(len(df))}")
    p(f"  Unique compounds     : {fmt(df['compound_chembl_id'].nunique())}")
    p(f"  Unique targets       : {fmt(df['target_name'].nunique())}")
    p(f"  Unique organisms     : {fmt(df['organism'].nunique())}")
    p(f"  Activity types       : {fmt(df['activity_type'].nunique())}")
    p(f"  Stereo families      : {fmt(df['stereochemical_family_id'].nunique())}")
    p(f"  n_measurements total : {fmt(int(df['n_measurements'].sum()))}")

    # ── 2. STEREOCENTERS ──────────────────────────────────────────────
    p()
    p("── 2. STEREOCENTER DISTRIBUTION ──────────────────────────────────────")
    dist = df["n_stereocenters"].value_counts().sort_index()
    for n, cnt in dist.head(15).items():
        bar = "█" * min(int(cnt / dist.max() * 30), 30)
        p(f"  {int(n):>3} center(s): {fmt(cnt):>10}  {bar}")

    # ── 3. METRIC COUNTS ──────────────────────────────────────────────
    p()
    p("── 3. COUNTS PER METRIC (top 40) ─────────────────────────────────────")
    act_dist = df["activity_type"].value_counts().head(40)
    for atype, cnt in act_dist.items():
        p(f"  {str(atype):<35s}: {fmt(cnt)}")

    # ── 4. ORGANISM COUNTS ────────────────────────────────────────────
    p()
    p("── 4. COUNTS PER ORGANISM (top 50) ───────────────────────────────────")
    org_dist = df["organism"].value_counts().head(50)
    for org, cnt in org_dist.items():
        p(f"  {str(org):<55s}: {fmt(cnt)}")

    # ── 5. TARGET COUNTS ──────────────────────────────────────────────
    p()
    p("── 5. COUNTS PER TARGET (top 50) ─────────────────────────────────────")
    tgt_dist = df["target_name"].value_counts().head(50)
    for tgt, cnt in tgt_dist.items():
        p(f"  {str(tgt)[:52]:<54s}: {fmt(cnt)}")

    # ── 6. STEREOCHEMICAL CLASS ───────────────────────────────────────
    p()
    p("── 6. STEREOCHEMICAL VARIANT DISTRIBUTION (top 20) ───────────────────")
    var_dist = df["stereochemical_variant"].value_counts().head(20)
    for var, cnt in var_dist.items():
        p(f"  {str(var):<25s}: {fmt(cnt)}")

    p()
    p("── 7. STEREO FORM SUMMARY ────────────────────────────────────────────")
    form_dist = df["stereo_form"].value_counts()
    for form, cnt in form_dist.items():
        p(f"  {str(form):<15s}: {fmt(cnt)}")

    # True isomers
    ti_path = os.path.join(ISOMER_DIR, "TRUE_ISOMER_PAIRS.csv")
    if os.path.exists(ti_path):
        ti = pd.read_csv(ti_path, low_memory=False, usecols=["smiles_nostereo","organism","target_name"])
        p()
        p("── 8. TRUE ISOMER SUMMARY ────────────────────────────────────────────")
        p(f"  TRUE_ISOMER_PAIRS records : {fmt(len(ti))}")
        p(f"  Unique scaffolds           : {fmt(ti['smiles_nostereo'].nunique())}")
        p(f"  Organisms represented      : {fmt(ti['organism'].nunique())}")
        p(f"  Targets represented        : {fmt(ti['target_name'].nunique())}")

    # Species summary
    ss_path = os.path.join(REPORTS_DIR, "species_summary.csv")
    if os.path.exists(ss_path):
        ss = pd.read_csv(ss_path)
        p()
        p("── 9. SPECIES STATISTICS ─────────────────────────────────────────────")
        p(f"  Total organism files generated : {fmt(len(ss))}")
        p(f"  Largest organism (records)     : {ss.iloc[0]['organism']} ({fmt(int(ss.iloc[0]['activity_records']))})")
        p(f"  Median records per organism    : {ss['activity_records'].median():.0f}")

    # ── 10. OUTPUT MANIFEST ───────────────────────────────────────────
    p()
    p("── 10. OUTPUT FILE MANIFEST ──────────────────────────────────────────")
    manifest_dirs = [
        (RAW_DIR,    "raw/"),
        (AGG_DIR,    "aggregated/"),
        (CORE_DIR,   "outputs/core_metrics/"),
        (FUNC_DIR,   "outputs/functional_metrics/"),
        (ASSAY_DIR,  "outputs/assay_metrics/"),
        (ISOMER_DIR, "outputs/true_isomers/"),
        (STEREO_DIR, "outputs/stereoisomers/"),
        (REPORTS_DIR,"reports/"),
    ]
    total_size = 0
    for d, label in manifest_dirs:
        rows = dir_manifest(d, label)
        if rows:
            p(f"\n  [{label}]  ({len(rows)} files)")
            for _, fname, size_mb in rows:
                p(f"    ✓ {fname:<55s} {size_mb:.1f} MB")
                total_size += size_mb

    # Species dir separately (many files)
    if os.path.isdir(SPECIES_DIR):
        sp_files = [f for f in os.listdir(SPECIES_DIR) if f.endswith(".csv")]
        sp_size = sum(os.path.getsize(os.path.join(SPECIES_DIR, f)) for f in sp_files) / 1e6
        p(f"\n  [outputs/species/]  ({len(sp_files)} organism files, {sp_size:.0f} MB total)")
        total_size += sp_size

    p()
    p(f"  Total disk usage: {total_size:.0f} MB")
    p()
    p("=" * 72)
    p("  ✅ All-Species v3 Atlas Complete")
    p("=" * 72)

    out_path = os.path.join(REPORTS_DIR, "FINAL_REPORT.txt")
    with open(out_path, "w") as f:
        f.write("\n".join(lines))
    print(f"\nReport written → reports/FINAL_REPORT.txt")
    print(f"\n✅ Step 7 complete.")

if __name__ == "__main__":
    main()
