"""
All-Species v3 — Step 2: Core & Functional & Assay Metric Files
================================================================
Reads master aggregated CSV.
Outputs per-metric CSVs into:
  outputs/core_metrics/      (IC50, Ki, Kd, EC50, Km, Vmax, Kcat, Kinact)
  outputs/functional_metrics/ (Potency, Emax, Efficacy, GI50, ED50, ...)
  outputs/assay_metrics/      (% Control, Ratio_IC50, Selectivity_Index, ...)
Also dumps ALL remaining metrics into outputs/core_metrics/OTHER_METRICS.csv
"""
import pandas as pd
import os, sys
sys.path.insert(0, os.path.dirname(__file__))
from config import *

def save_metric(df: pd.DataFrame, metric_label: str, out_dir: str, filename: str):
    """Filter to a metric (case-insensitive), save, return row count."""
    mask = df["activity_type"].str.strip().str.upper() == metric_label.upper()
    sub = df[mask].copy()
    if len(sub) == 0:
        return 0
    # Keep only BASE_COLS that exist
    cols = [c for c in BASE_COLS if c in sub.columns]
    sub = sub[cols]
    path = os.path.join(out_dir, filename)
    sub.to_csv(path, index=False)
    return len(sub)

def save_metric_multi(df: pd.DataFrame, labels: list, out_dir: str, filename: str):
    """Filter to multiple metric names, save."""
    upper_labels = [l.upper() for l in labels]
    mask = df["activity_type"].str.strip().str.upper().isin(upper_labels)
    sub = df[mask].copy()
    if len(sub) == 0:
        return 0
    cols = [c for c in BASE_COLS if c in sub.columns]
    sub = sub[cols]
    path = os.path.join(out_dir, filename)
    sub.to_csv(path, index=False)
    return len(sub)

def main():
    print("=" * 60)
    print("Step 2: Metric-organized output files")
    print("=" * 60)

    print(f"\nLoading master aggregated CSV...")
    df = pd.read_csv(SRC_AGGREGATED, low_memory=False)
    print(f"  {len(df):,} rows loaded")

    # Track which activity types have been claimed
    covered = set()
    results = {}

    # ── 1. CORE METRICS ───────────────────────────────────────────────
    print(f"\n  [core_metrics/]")
    for metric, fname in CORE_METRICS.items():
        n = save_metric(df, metric, CORE_DIR, f"{fname}.csv")
        print(f"    {fname}.csv — {n:,} rows")
        if n > 0:
            covered.add(metric.upper())
        results[f"core/{fname}"] = n

    # ── 2. FUNCTIONAL METRICS ─────────────────────────────────────────
    print(f"\n  [functional_metrics/]")

    # Efficacy — combine Efficacy + Intrinsic_activity
    eff_variants = ["Efficacy", "Intrinsic activity", "Intrinsic_activity", "Alpha max", "Alpha_max"]
    n = save_metric_multi(df, eff_variants, FUNC_DIR, "Efficacy.csv")
    print(f"    Efficacy.csv — {n:,} rows")
    covered.update([v.upper() for v in eff_variants if v])

    # Inhibition — save separately from assay
    n = save_metric(df, "Inhibition", FUNC_DIR, "Inhibition.csv")
    print(f"    Inhibition.csv — {n:,} rows")
    covered.add("INHIBITION")

    # Other functional metrics
    for metric, fname in FUNCTIONAL_METRICS.items():
        if metric.upper() in covered:
            continue
        n = save_metric(df, metric, FUNC_DIR, f"{fname}.csv")
        print(f"    {fname}.csv — {n:,} rows")
        if n > 0:
            covered.add(metric.upper())
        results[f"func/{fname}"] = n

    # Emax already done in core, but also save in functional for easy access
    n_emax = save_metric(df, "Emax", FUNC_DIR, "Emax.csv")
    print(f"    Emax.csv — {n_emax:,} rows")

    # ── 3. ASSAY METRICS ──────────────────────────────────────────────
    print(f"\n  [assay_metrics/]")
    # Map: exact activity_type value → filename
    assay_map = {
        "% Control":          "Percent_Control.csv",
        "% Activity":         "Percent_Activity.csv",
        "% Inhibition":       "Percent_Inhibition.csv",
        "Percent Effect":     "Percent_Effect.csv",
        "% Enzyme Activity":  "Percent_Enzyme_Activity.csv",
        "Inhibition":         "Inhibition.csv",
        "Ratio IC50":         "Ratio_IC50.csv",
        "Ratio_IC50":         "Ratio_IC50.csv",
        "Ratio":              "Ratio.csv",
        "FC":                 "FC.csv",
        "Selectivity Index":  "Selectivity_Index.csv",
        "Selectivity_Index":  "Selectivity_Index.csv",
        "RBA":                "RBA.csv",
        "INH":                "INH.csv",
        "Activity":           "Activity.csv",
        "% of control":       "Percent_of_Control.csv",
        "% of activity":      "Percent_of_Activity.csv",
        "% of inhibition":    "Percent_of_Inhibition.csv",
        "Ratio IC50/EC50":    "Ratio_IC50_EC50.csv",
        "Ratio Ki":           "Ratio_Ki.csv",
    }
    seen_assay = {}  # fname → set of labels
    for atype, fname in assay_map.items():
        if fname not in seen_assay:
            seen_assay[fname] = []
        seen_assay[fname].append(atype)

    for fname, labels in seen_assay.items():
        n = save_metric_multi(df, labels, ASSAY_DIR, fname)
        print(f"    {fname} — {n:,} rows")
        covered.update([l.upper() for l in labels if n > 0])
        results[f"assay/{fname}"] = n

    # ── 4. DUMP REMAINING / RARE METRICS ──────────────────────────────
    rare_mask = ~df["activity_type"].str.strip().str.upper().isin(covered)
    rare_df = df[rare_mask].copy()
    print(f"\n  Remaining / rare metrics: {rare_df['activity_type'].nunique()} types, {len(rare_df):,} rows")
    if len(rare_df) > 0:
        cols = [c for c in BASE_COLS if c in rare_df.columns]
        rare_df[cols].to_csv(os.path.join(CORE_DIR, "OTHER_METRICS.csv"), index=False)
        print(f"    Saved → core_metrics/OTHER_METRICS.csv")
        # Also save per rare metric into assay folder
        for atype, grp in rare_df.groupby("activity_type"):
            safe = str(atype).replace(" ", "_").replace("/", "_").replace("%", "pct")[:50]
            grp[[c for c in BASE_COLS if c in grp.columns]].to_csv(
                os.path.join(ASSAY_DIR, f"{safe}.csv"), index=False)
        print(f"    Also saved {rare_df['activity_type'].nunique()} individual rare-metric files → assay_metrics/")

    print(f"\n✅ Step 2 complete.")

if __name__ == "__main__":
    main()
