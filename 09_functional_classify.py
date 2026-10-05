"""
All-Species v3 — Post-processing Extension: Functional Metric Classification
=============================================================================
Reads ONLY existing flat CSVs from outputs/functional_metrics/*.csv.
No DB access. No re-aggregation.

For each functional metric creates:
    outputs/functional_metrics/<METRIC>/
        TRUE_ISOMERS.csv
        ONLY_R.csv
        ONLY_S.csv
        MULTI_CENTER.csv

Classification rules (same as 08_postprocess_classify.py):
    MULTI_CENTER   : n_stereocenters > 1
    TRUE_ISOMERS   : group (smiles_nostereo + target_name + activity_type + organism)
                     has both "R" and "S" variant observed.
                     Multi-center families with >=2 distinct stereo_forms also qualify.
    ONLY_R         : single-center, only "R" in group, no paired "S"
    ONLY_S         : single-center, only "S" in group, no paired "R"

Also updates each species folder with:
    outputs/species/<organism>/functional_metrics_stereo_summary.csv
        (per-functional-metric counts within that organism)

Generates:
    reports/functional_metric_stereo_summary.csv
"""

import pandas as pd
import os, re, time
from pathlib import Path

# ── Paths ──────────────────────────────────────────────────────────────────
V3_DIR   = os.path.dirname(os.path.abspath(__file__))
FUNC_DIR = os.path.join(V3_DIR, "outputs", "functional_metrics")
SPEC_DIR = os.path.join(V3_DIR, "outputs", "species")
REP_DIR  = os.path.join(V3_DIR, "reports")

# All flat functional metric CSVs to process
FUNCTIONAL_METRICS = [
    "Potency", "Emax", "Efficacy", "GI50", "ED50",
    "Activity", "Response", "Inhibition",
    "AC50", "Imax", "DC50", "pD2",
]

GROUP_COLS = ["smiles_nostereo", "target_name", "activity_type", "organism"]

# ── Utilities ──────────────────────────────────────────────────────────────
def safe_name(s: str) -> str:
    n = re.sub(r"[^\w\s-]", "", str(s))
    return re.sub(r"[\s/]+", "_", n.strip())[:80]

def save(df: pd.DataFrame, path: str):
    if len(df) == 0:
        return
    os.makedirs(os.path.dirname(path), exist_ok=True)
    df.to_csv(path, index=False)

# ── Classification ─────────────────────────────────────────────────────────
def classify(df: pd.DataFrame) -> pd.DataFrame:
    """Add stereo_class column. Identical logic to 08_postprocess_classify.py."""
    df = df.copy()
    df["stereo_class"] = "UNCATEGORIZED"

    # Multi-center rows
    is_multi = df["n_stereocenters"] > 1
    df.loc[is_multi, "stereo_class"] = "MULTI_CENTER"

    # Single-center rows
    single = df[~is_multi].copy()
    if len(single) > 0:
        group_variants = (
            single.groupby(GROUP_COLS)["stereochemical_variant"]
            .apply(lambda x: frozenset(str(v).strip() for v in x))
        )
        true_groups = set(
            group_variants[
                group_variants.apply(lambda s: "R" in s and "S" in s)
            ].index.tolist()
        )

        def row_key(r):
            return (r["smiles_nostereo"], r["target_name"],
                    r["activity_type"], r["organism"])

        true_mask = single.apply(lambda r: row_key(r) in true_groups, axis=1)
        df.loc[single[true_mask].index, "stereo_class"] = "TRUE_ISOMERS"

        remaining = single[~true_mask]
        df.loc[remaining[remaining["stereochemical_variant"].str.strip() == "R"].index, "stereo_class"] = "ONLY_R"
        df.loc[remaining[remaining["stereochemical_variant"].str.strip() == "S"].index, "stereo_class"] = "ONLY_S"

    # Multi-center with >=2 stereoforms → MULTI_CENTER_TRUE
    multi_df = df[is_multi]
    if len(multi_df) > 0:
        fam_forms = multi_df.groupby("stereochemical_family_id")["stereo_form"].nunique()
        true_fams = set(fam_forms[fam_forms >= 2].index)
        mct_idx = multi_df[multi_df["stereochemical_family_id"].isin(true_fams)].index
        df.loc[mct_idx, "stereo_class"] = "MULTI_CENTER_TRUE"

    return df

def write_splits(classified: pd.DataFrame, out_dir: str) -> dict:
    """Write 4 CSV files, return counts dict."""
    os.makedirs(out_dir, exist_ok=True)
    counts = {}

    # TRUE_ISOMERS = TRUE_ISOMERS + MULTI_CENTER_TRUE
    true_df = classified[classified["stereo_class"].isin(["TRUE_ISOMERS", "MULTI_CENTER_TRUE"])]
    save(true_df, os.path.join(out_dir, "TRUE_ISOMERS.csv"))
    counts["TRUE_ISOMERS"] = len(true_df)

    only_r = classified[classified["stereo_class"] == "ONLY_R"]
    save(only_r, os.path.join(out_dir, "ONLY_R.csv"))
    counts["ONLY_R"] = len(only_r)

    only_s = classified[classified["stereo_class"] == "ONLY_S"]
    save(only_s, os.path.join(out_dir, "ONLY_S.csv"))
    counts["ONLY_S"] = len(only_s)

    multi = classified[classified["stereo_class"].isin(["MULTI_CENTER", "MULTI_CENTER_TRUE"])]
    save(multi, os.path.join(out_dir, "MULTI_CENTER.csv"))
    counts["MULTI_CENTER"] = len(multi)

    return counts

# ── Main ───────────────────────────────────────────────────────────────────
def main():
    print("=" * 65)
    print("Post-processing Extension: Functional Metric Classification")
    print("  Source: outputs/functional_metrics/*.csv (no DB access)")
    print("=" * 65)
    t_start = time.time()

    # ── PART A: Per-metric subfolders ──────────────────────────────────
    print(f"\n{'─'*65}")
    print("PART A: Classifying functional metric files")
    print(f"{'─'*65}")

    summary_rows = []

    # Build a combined DataFrame across all metrics for species work later
    all_classified_parts = []

    for metric in FUNCTIONAL_METRICS:
        src = os.path.join(FUNC_DIR, f"{metric}.csv")
        if not os.path.exists(src):
            print(f"  {metric:<14s}: file not found, skipping.")
            continue

        df = pd.read_csv(src, low_memory=False)
        if len(df) == 0:
            print(f"  {metric:<14s}: empty file, skipping.")
            continue

        # Ensure required columns exist
        missing = [c for c in ["n_stereocenters","stereochemical_variant",
                                "stereo_form","stereochemical_family_id",
                                "smiles_nostereo","target_name",
                                "activity_type","organism"]
                   if c not in df.columns]
        if missing:
            print(f"  {metric:<14s}: missing columns {missing}, skipping.")
            continue

        classified = classify(df)
        out_dir = os.path.join(FUNC_DIR, metric)
        counts = write_splits(classified, out_dir)

        row = {
            "metric":             metric,
            "total_records":      len(df),
            "compound_count":     df["compound_chembl_id"].nunique(),
            "target_count":       df["target_name"].nunique(),
            "organism_count":     df["organism"].nunique(),
            "TRUE_ISOMER_count":  counts["TRUE_ISOMERS"],
            "ONLY_R_count":       counts["ONLY_R"],
            "ONLY_S_count":       counts["ONLY_S"],
            "MULTI_CENTER_count": counts["MULTI_CENTER"],
        }
        summary_rows.append(row)

        print(f"  {metric:<14s}: "
              f"TRUE={counts['TRUE_ISOMERS']:>7,}  "
              f"R={counts['ONLY_R']:>7,}  "
              f"S={counts['ONLY_S']:>7,}  "
              f"MULTI={counts['MULTI_CENTER']:>7,}  "
              f"[{len(df):,} total]")

        # Keep classified slice for species work
        classified["_metric"] = metric
        all_classified_parts.append(classified)

    # Save functional metric summary
    func_summary = pd.DataFrame(summary_rows)
    out_path = os.path.join(REP_DIR, "functional_metric_stereo_summary.csv")
    func_summary.to_csv(out_path, index=False)
    print(f"\n  ✓ Saved reports/functional_metric_stereo_summary.csv")

    # ── PART B: Update species folders ────────────────────────────────
    print(f"\n{'─'*65}")
    print("PART B: Writing per-species functional metric stereo summaries")
    print(f"{'─'*65}")

    if not all_classified_parts:
        print("  No classified data available, skipping species update.")
    else:
        combined = pd.concat(all_classified_parts, ignore_index=True)
        print(f"  Combined functional records: {len(combined):,}")

        organisms = sorted(combined["organism"].dropna().unique())
        print(f"  Organisms with functional metric data: {len(organisms)}")

        species_func_rows = []
        for i, org in enumerate(organisms):
            org_sub = combined[combined["organism"] == org]
            org_dir = os.path.join(SPEC_DIR, safe_name(org))

            # Build per-metric breakdown for this organism
            org_metric_rows = []
            for metric_name, grp in org_sub.groupby("_metric"):
                true_n  = grp["stereo_class"].isin(["TRUE_ISOMERS","MULTI_CENTER_TRUE"]).sum()
                only_r  = (grp["stereo_class"] == "ONLY_R").sum()
                only_s  = (grp["stereo_class"] == "ONLY_S").sum()
                multi_n = grp["stereo_class"].isin(["MULTI_CENTER","MULTI_CENTER_TRUE"]).sum()
                org_metric_rows.append({
                    "organism":           org,
                    "metric":             metric_name,
                    "total_records":      len(grp),
                    "compound_count":     grp["compound_chembl_id"].nunique(),
                    "target_count":       grp["target_name"].nunique(),
                    "TRUE_ISOMER_count":  int(true_n),
                    "ONLY_R_count":       int(only_r),
                    "ONLY_S_count":       int(only_s),
                    "MULTI_CENTER_count": int(multi_n),
                })

            if org_metric_rows:
                os.makedirs(org_dir, exist_ok=True)
                org_df = pd.DataFrame(org_metric_rows)
                org_df.to_csv(
                    os.path.join(org_dir, "functional_metrics_stereo_summary.csv"),
                    index=False
                )
                # Aggregate row for the master species table
                species_func_rows.append({
                    "organism":           org,
                    "functional_records": len(org_sub),
                    "func_compound_count": org_sub["compound_chembl_id"].nunique(),
                    "func_target_count":  org_sub["target_name"].nunique(),
                    "func_metric_count":  org_sub["_metric"].nunique(),
                    "func_TRUE_ISOMERS":  int(org_sub["stereo_class"].isin(["TRUE_ISOMERS","MULTI_CENTER_TRUE"]).sum()),
                    "func_ONLY_R":        int((org_sub["stereo_class"] == "ONLY_R").sum()),
                    "func_ONLY_S":        int((org_sub["stereo_class"] == "ONLY_S").sum()),
                    "func_MULTI_CENTER":  int(org_sub["stereo_class"].isin(["MULTI_CENTER","MULTI_CENTER_TRUE"]).sum()),
                })

            if (i + 1) % 50 == 0 or (i + 1) == len(organisms):
                print(f"  {i+1}/{len(organisms)} species updated...", end="\r")

        print(f"\n  All {len(organisms)} species folders updated with functional_metrics_stereo_summary.csv")

        # Save master species-level functional summary
        sp_func_df = (
            pd.DataFrame(species_func_rows)
            .sort_values("functional_records", ascending=False)
        )
        sp_func_path = os.path.join(REP_DIR, "species_functional_stereo_summary.csv")
        sp_func_df.to_csv(sp_func_path, index=False)
        print(f"  ✓ Saved reports/species_functional_stereo_summary.csv")

        # Print top species
        print(f"\n  Top 10 species by functional TRUE_ISOMERS:")
        top = sp_func_df.nlargest(10, "func_TRUE_ISOMERS")
        for _, r in top.iterrows():
            print(f"    {str(r['organism'])[:42]:<44s}: "
                  f"TRUE={int(r['func_TRUE_ISOMERS']):>7,}  "
                  f"R={int(r['func_ONLY_R']):>7,}  "
                  f"S={int(r['func_ONLY_S']):>7,}  "
                  f"MULTI={int(r['func_MULTI_CENTER']):>7,}")

    # ── Grand totals ───────────────────────────────────────────────────
    elapsed = time.time() - t_start
    print(f"\n{'─'*65}")
    print("SUMMARY")
    print(f"{'─'*65}")
    if summary_rows:
        total_true  = sum(r["TRUE_ISOMER_count"]  for r in summary_rows)
        total_r     = sum(r["ONLY_R_count"]        for r in summary_rows)
        total_s     = sum(r["ONLY_S_count"]        for r in summary_rows)
        total_multi = sum(r["MULTI_CENTER_count"]  for r in summary_rows)
        print(f"  Functional metrics processed : {len(summary_rows)}")
        print(f"  Species folders updated      : {len(organisms) if all_classified_parts else 0}")
        print(f"  Total TRUE_ISOMERS           : {total_true:,}")
        print(f"  Total ONLY_R                 : {total_r:,}")
        print(f"  Total ONLY_S                 : {total_s:,}")
        print(f"  Total MULTI_CENTER           : {total_multi:,}")
    print(f"  Time                         : {elapsed:.1f}s")
    print(f"\n✅ Functional metric classification complete.")

if __name__ == "__main__":
    main()
