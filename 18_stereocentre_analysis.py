#!/usr/bin/env python3
"""
32d_stereocentre_analysis_v4.py
═══════════════════════════════════════════════════════════════════════════════
Stereocentre Analysis v4 — All_species_v3
Builds directly on 32c_stereocentre_analysis_v3.py (DO NOT modify v3)

Adds TWO things only:
  1. Grouped main table (STEREOCENTRE_COUNTS_v4.xlsx — Sheet 1)
       family_rank, is_family_start, family_member_index,
       stripped_smiles_verification added; sorted scaffold → binary
       Blue rows = first member of each family; alternating white/grey within
  2. Family verification sheet (Sheet 2 in same workbook)
       One row per family; green=full coverage, amber=partial, red=single

Run with:
  /opt/homebrew/Caskroom/miniforge/base/envs/descriptor_env/bin/python \
      32d_stereocentre_analysis_v4.py

READ-ONLY: No existing v3 files are overwritten.
"""

import csv
import itertools
import json
import os
import re
import sys
import numpy as np
import pandas as pd

try:
    import xlsxwriter
except ImportError:
    print("ERROR: xlsxwriter not found. Run with descriptor_env.")
    sys.exit(1)

# ─────────────────────────────────────────────────────────────────────────────
# PATHS
# ─────────────────────────────────────────────────────────────────────────────
BASE_DIR    = "BASE_DIR  # <-- set BASE_DIR to your pipeline root"
OUT_DIR     = os.path.join(BASE_DIR, "outputs", "stereocentre_analysis")
OUTPUTS_DIR = os.path.join(BASE_DIR, "outputs")
os.makedirs(OUT_DIR, exist_ok=True)

# ─────────────────────────────────────────────────────────────────────────────
# COLOURS
# ─────────────────────────────────────────────────────────────────────────────
COL_HEADER    = '#1F4E79'   # dark navy — header bg
COL_FAM_START = '#D6EAF8'   # light blue — first row of new family
COL_ODD       = '#FFFFFF'   # white — odd rows within family
COL_EVEN      = '#F8F9FA'   # very light grey — even rows within family
COL_GREEN_ROW = '#D5F5E3'   # green  — full coverage family
COL_AMBER_ROW = '#FEF9E7'   # amber  — partial coverage family
COL_RED_ROW   = '#FDEDEC'   # light red — single-entry family

# ─────────────────────────────────────────────────────────────────────────────
# STRING COLUMNS — always written with write_string() to preserve leading zeros
# ─────────────────────────────────────────────────────────────────────────────
FORCE_TEXT_COLS = {
    "stereo_binary_string", "stereo_label", "bracket_notation",
    "canonical_smiles", "scaffold_smiles", "stripped_smiles_verification",
    "verification_flag", "molecule_chembl_id", "molecule_name",
    "present_combinations", "absent_combinations", "present_labels", "absent_labels",
    "missing_binary_string", "missing_label", "missing_bracket",
    "compounds_present", "stereoselectivity_class",
    # family verification columns
    "members_chembl_ids", "all_stereo_labels", "all_binary_strings",
    "all_bracket_notations", "missing_labels",
}


# ═════════════════════════════════════════════════════════════════════════════
# FUNCTIONS FROM v3 — UNCHANGED
# ═════════════════════════════════════════════════════════════════════════════

def count_stereocentres(smiles: str):
    """Count R (@) and S (@@) stereocentres. @@ counted BEFORE @."""
    if not isinstance(smiles, str) or smiles.strip() == "" or smiles == "nan":
        return 0, 0, 0
    count_S   = len(re.findall(r'@@', smiles))
    remaining = smiles.replace('@@', '')
    count_R   = len(re.findall(r'@', remaining))
    return count_R, count_S, count_R + count_S


def encode_stereo_binary(smiles: str):
    """Walk SMILES left-to-right. @@ (S)→1, @ (R)→0. Returns (binary, label, bracket)."""
    if not isinstance(smiles, str) or smiles.strip() == "" or smiles == "nan":
        return "", "", ""
    bits = []
    i = 0
    while i < len(smiles):
        if smiles[i] == '@':
            if i + 1 < len(smiles) and smiles[i + 1] == '@':
                bits.append('1')
                i += 2
            else:
                bits.append('0')
                i += 1
        else:
            i += 1
    binary_str       = "".join(bits)
    label            = binary_str.replace('0', 'R').replace('1', 'S')
    bracket_notation = "".join('[@]' if b == '0' else '[@@]' for b in binary_str)
    return binary_str, label, bracket_notation


def verify_row(row) -> str:
    """Cross-check count and encoding; returns '' if OK, else error description."""
    n_R, n_S, n_tot = row["n_R_centres"], row["n_S_centres"], row["n_stereocentres"]
    bits = row["stereo_binary_string"]
    if not isinstance(bits, str) or bits == "":
        return "" if n_tot == 0 else "binary_string_empty"
    errors = []
    if len(bits) != n_R + n_S:
        errors.append(f"len(binary)={len(bits)} != n_R+n_S={n_R+n_S}")
    if bits.count('0') != n_R:
        errors.append(f"R_bits={bits.count('0')} != n_R={n_R}")
    if bits.count('1') != n_S:
        errors.append(f"S_bits={bits.count('1')} != n_S={n_S}")
    if n_R + n_S != n_tot:
        errors.append(f"n_R+n_S={n_R+n_S} != n_stereocentres={n_tot}")
    return "; ".join(errors)


def load_input_data():
    candidates = [
        os.path.join(OUTPUTS_DIR, "STEREOATLAS_MASTER_LITE.csv"),
        os.path.join(OUTPUTS_DIR, "STEREOATLAS_MASTER.csv"),
        os.path.join(OUTPUTS_DIR, "true_isomers", "TRUE_ISOMER_PAIRS.csv"),
        os.path.join(OUTPUTS_DIR, "true_isomers", "MULTI_CENTER_STEREO.csv"),
        os.path.join(OUTPUTS_DIR, "stereoisomers", "ALL_METRICS_STEREOISOMERS.csv"),
    ]
    df, loaded_path = None, None
    for path in candidates:
        if os.path.exists(path):
            print(f"  Loading: {path}")
            df = pd.read_csv(path, low_memory=False)
            loaded_path = path
            break
    if df is None:
        raise FileNotFoundError("No usable input file found.")
    print(f"  Loaded {len(df):,} rows from {os.path.basename(loaded_path)}")

    col_map = {}
    for c in ["canonical_smiles", "smiles", "standard_smiles"]:
        if c in df.columns: col_map["canonical_smiles"] = c; break
    for c in ["smiles_nostereo", "scaffold_smiles", "family_id",
              "scaffold_stripped", "stereochemical_family_id"]:
        if c in df.columns: col_map["scaffold_smiles"] = c; break
    for c in ["compound_chembl_id", "molecule_chembl_id", "chembl_id"]:
        if c in df.columns: col_map["molecule_chembl_id"] = c; break
    for c in ["compound_name", "molecule_name", "pref_name"]:
        if c in df.columns: col_map["molecule_name"] = c; break

    if "canonical_smiles" not in col_map:
        raise ValueError(f"No SMILES col found. Available: {df.columns.tolist()}")
    if "scaffold_smiles" not in col_map:
        raise ValueError(f"No scaffold col found. Available: {df.columns.tolist()}")

    print(f"  SMILES col  : '{col_map['canonical_smiles']}'")
    print(f"  Scaffold col: '{col_map['scaffold_smiles']}'")

    work = pd.DataFrame({
        "canonical_smiles"  : df[col_map["canonical_smiles"]].astype(str),
        "scaffold_smiles"   : df[col_map["scaffold_smiles"]].astype(str),
        "molecule_chembl_id": (df[col_map["molecule_chembl_id"]].astype(str)
                               if "molecule_chembl_id" in col_map else "UNKNOWN"),
        "molecule_name"     : (df[col_map["molecule_name"]]
                               if "molecule_name" in col_map else np.nan),
    })
    work = work[work["canonical_smiles"].str.strip().str.len() > 0].copy()
    work = work[work["canonical_smiles"] != "nan"].copy()
    work = work.drop_duplicates(subset=["molecule_chembl_id"]).reset_index(drop=True)
    print(f"  Unique compounds after dedup: {len(work):,}")
    return work


def load_fold_data():
    path = os.path.join(OUTPUTS_DIR, "stereo_selectivity", "family_views",
                        "FAMILY_MINMAX_FOLD.csv")
    if not os.path.exists(path):
        print("  FAMILY_MINMAX_FOLD.csv not found — fold columns will be blank.")
        return None
    df = pd.read_csv(path, low_memory=False)
    print(f"  Loaded fold data: {len(df):,} rows.")
    return df


def all_combinations(n: int):
    return ["".join(b) for b in itertools.product("01", repeat=n)]


def binary_to_label(b: str) -> str:
    return b.replace('0', 'R').replace('1', 'S')


def binary_to_bracket(b: str) -> str:
    return "".join('[@]' if c == '0' else '[@@]' for c in b)


def save_json(df: pd.DataFrame, path: str):
    records = []
    for _, row in df.iterrows():
        rec = {}
        for col in df.columns:
            v = row[col]
            if isinstance(v, (np.integer,)): v = int(v)
            elif isinstance(v, (np.floating,)): v = None if np.isnan(v) else float(v)
            elif isinstance(v, (np.bool_,)): v = bool(v)
            elif isinstance(v, float) and np.isnan(v): v = None
            else: v = str(v) if v is not None else None
            rec[col] = v
        records.append(rec)
    with open(path, "w", encoding="utf-8") as f:
        json.dump(records, f, indent=2, ensure_ascii=False)
    print(f"  Saved JSON : {os.path.basename(path)}  ({len(df):,} records)")


def save_csv_quoted(df: pd.DataFrame, path: str):
    df.to_csv(path, index=False, quoting=csv.QUOTE_ALL)
    print(f"  Saved CSV  : {os.path.basename(path)}  ({len(df):,} rows)")


# ═════════════════════════════════════════════════════════════════════════════
# NEW FUNCTIONS — v4 ADDITIONS
# ═════════════════════════════════════════════════════════════════════════════

def add_family_grouping_columns(df: pd.DataFrame) -> pd.DataFrame:
    """
    Sort by scaffold_smiles (alpha) then stereo_binary_string.
    Add: family_rank, family_member_index, is_family_start,
         stripped_smiles_verification.
    """
    # Sort scaffold alphabetically, then by binary string within family
    df = df.sort_values(
        ["scaffold_smiles", "stereo_binary_string"],
        key=lambda col: col.astype(str)
    ).reset_index(drop=True)

    # Assign family_rank: 1-indexed, ordered by scaffold_smiles alphabetically
    unique_scaffolds = sorted(df["scaffold_smiles"].unique(), key=str)
    rank_map = {s: i + 1 for i, s in enumerate(unique_scaffolds)}
    df["family_rank"] = df["scaffold_smiles"].map(rank_map).astype(int)

    # family_member_index: 1, 2, 3... within each family (resets per family)
    df["family_member_index"] = (
        df.groupby("scaffold_smiles").cumcount() + 1
    ).astype(int)

    # is_family_start: True only for first member
    df["is_family_start"] = df["family_member_index"] == 1

    # stripped_smiles_verification: explicit copy of scaffold_smiles
    df["stripped_smiles_verification"] = df["scaffold_smiles"]

    return df


def compute_family_coverage_lookup(stereo_df: pd.DataFrame) -> dict:
    """
    Compute family combination coverage for every scaffold.
    Returns dict keyed on scaffold_smiles string.
    """
    print(f"  Computing family coverage for {stereo_df['scaffold_smiles'].nunique():,} families ...")
    lookup = {}
    for scaffold, grp in stereo_df.groupby("scaffold_smiles"):
        present_combos = set(
            s for s in grp["stereo_binary_string"].unique()
            if isinstance(s, str) and s != ""
        )
        max_n      = int(grp["n_stereocentres"].max())
        n_for_enum = min(max_n, 5)
        if n_for_enum == 0:
            continue
        theoretical   = set(all_combinations(n_for_enum))
        present_at_n  = set(c for c in present_combos if len(c) == n_for_enum)
        absent_combos = theoretical - present_at_n
        max_comb      = 2 ** n_for_enum
        coverage_pct  = len(present_at_n) / max_comb * 100.0 if max_comb > 0 else 0.0
        has_full      = (absent_combos == set())
        missing_labels= ",".join(binary_to_label(c) for c in sorted(absent_combos))

        lookup[str(scaffold)] = {
            "coverage_pct"    : round(coverage_pct, 2),
            "has_full_coverage": has_full,
            "missing_labels"  : missing_labels,
            "max_combinations": max_comb,
            "present_count"   : len(present_at_n),
            "absent_count"    : len(absent_combos),
            "present_combos"  : ",".join(sorted(present_at_n)),
            "absent_combos"   : ",".join(sorted(absent_combos)),
            "present_labels"  : ",".join(binary_to_label(c) for c in sorted(present_at_n)),
            "n_for_enum"      : n_for_enum,
        }
    return lookup


def build_family_verification_df(
        grouped_df: pd.DataFrame,
        coverage_lookup: dict) -> pd.DataFrame:
    """
    Build one-row-per-family verification DataFrame.
    Members sorted by stereo_binary_string within each family.

    verification_pass uses PER-COMPOUND n_stereocentres (not family max)
    so that families with mixed stereocentre counts are not falsely flagged.
    A family is valid when each compound's binary string length == that
    compound's own n_stereocentres — regardless of other members.
    """
    rows = []
    for scaffold, grp in sorted(
            grouped_df.groupby("scaffold_smiles"),
            key=lambda x: str(x[0])):

        grp_s = grp.sort_values("stereo_binary_string", key=lambda c: c.astype(str))

        family_rank = int(grp_s["family_rank"].iloc[0])
        n_members   = len(grp_s)
        n_stereo    = int(grp_s["n_stereocentres"].max())  # family max (for display)

        chembl_ids     = ",".join(str(v) for v in grp_s["molecule_chembl_id"].tolist())
        stereo_labels  = ",".join(str(v) for v in grp_s["stereo_label"].tolist())
        binary_strings = ",".join(str(v) for v in grp_s["stereo_binary_string"].tolist())
        brackets       = ",".join(str(v) for v in grp_s["bracket_notation"].tolist())
        flag_list      = [str(v) for v in grp_s["verification_flag"].tolist()]

        # Per-compound n_stereocentres list (order matches binary_strings)
        per_compound_n = grp_s["n_stereocentres"].tolist()
        n_stereo_str   = ",".join(str(n) for n in per_compound_n)
        has_mixed_n    = len(set(per_compound_n)) > 1   # True if compounds differ in count

        cov            = coverage_lookup.get(str(scaffold), {})
        coverage_pct   = cov.get("coverage_pct",     0.0)
        has_full       = cov.get("has_full_coverage", False)
        missing_labels = cov.get("missing_labels",    "")

        # ── Verification pass — per-compound check (FIXED from v4.0) ─────────
        # Each binary string is checked against its OWN compound's n_stereocentres,
        # NOT the family maximum. Mixed-n families (valid data) now pass correctly.
        binary_list = [b for b in binary_strings.split(",") if b]
        labels_list = [l for l in stereo_labels.split(",") if l]

        count_ok = (len(labels_list) == n_members and len(binary_list) == n_members)
        flag_ok  = all(f == "" for f in flag_list)
        # Per-compound length check: zip each binary string with its compound's n
        len_ok   = all(
            len(b) == int(n)
            for b, n in zip(binary_list, per_compound_n)
        ) if count_ok else False
        vpass    = count_ok and flag_ok and len_ok
        # ─────────────────────────────────────────────────────────────────────

        rows.append({
            "family_rank"              : family_rank,
            "scaffold_smiles"          : str(scaffold),
            "n_members"                : n_members,
            "n_stereocentres_max"      : n_stereo,
            "n_stereocentres_per_member": n_stereo_str,
            "has_mixed_stereocentre_count": has_mixed_n,
            "members_chembl_ids"       : chembl_ids,
            "all_stereo_labels"        : stereo_labels,
            "all_binary_strings"       : binary_strings,
            "all_bracket_notations"    : brackets,
            "coverage_pct"             : coverage_pct,
            "missing_labels"           : missing_labels,
            "has_full_coverage"        : has_full,
            "verification_pass"        : vpass,
        })

    return pd.DataFrame(rows)


# ─────────────────────────────────────────────────────────────────────────────
# XLSXWRITER FORMAT FACTORY
# ─────────────────────────────────────────────────────────────────────────────

def _fmt(wb, bg, bold=False, italic=False, align='left', valign='vcenter',
         font_color='#000000', num_format=None, border=0, border_color='#CCCCCC'):
    """Create and return an xlsxwriter format object."""
    d = {
        'bg_color'  : bg,
        'align'     : align,
        'valign'    : valign,
        'font_color': font_color,
        'bold'      : bold,
        'italic'    : italic,
        'border'    : border,
        'border_color': border_color,
    }
    if num_format:
        d['num_format'] = num_format
    return wb.add_format(d)


# ─────────────────────────────────────────────────────────────────────────────
# TWO-SHEET XLSX WRITER
# ─────────────────────────────────────────────────────────────────────────────

def save_two_sheet_xlsx(df_main: pd.DataFrame,
                        df_verify: pd.DataFrame,
                        path: str):
    """
    Write STEREOCENTRE_COUNTS_v4.xlsx with two sheets:
      Sheet 1 "Stereocentre_Counts" — grouped, colour-coded rows
      Sheet 2 "Family_Verification" — one row per family, coverage colours

    Binary/label/notation columns always use write_string() (v3 rule preserved).
    """
    wb = xlsxwriter.Workbook(path, {'strings_to_numbers': False})

    # ── Shared header format ─────────────────────────────────────────────────
    hdr = wb.add_format({
        'bold': True, 'bg_color': COL_HEADER, 'font_color': '#FFFFFF',
        'align': 'center', 'valign': 'vcenter', 'border': 1,
        'text_wrap': True,
    })

    # ── Sheet-1 formats (keyed by bg colour) ─────────────────────────────────
    def s1_fmts(bg):
        """Return dict of formats for a given row background."""
        return {
            'base'    : _fmt(wb, bg),
            'bold_ctr': _fmt(wb, bg, bold=True,  align='center'),
            'str'     : _fmt(wb, bg, num_format='@'),
            'num'     : _fmt(wb, bg),
            'green_txt': _fmt(wb, bg, font_color='#006400', bold=True),   # True
            'pct'     : _fmt(wb, bg, num_format='0.00'),
        }

    FMT_START = s1_fmts(COL_FAM_START)   # blue  — first member of family
    FMT_ODD   = s1_fmts(COL_ODD)         # white — odd rows within family
    FMT_EVEN  = s1_fmts(COL_EVEN)        # grey  — even rows within family

    # ── Sheet-2 formats ───────────────────────────────────────────────────────
    def s2_fmts(bg):
        return {
            'base': _fmt(wb, bg),
            'str' : _fmt(wb, bg, num_format='@'),
            'bold': _fmt(wb, bg, bold=True),
            'pct' : _fmt(wb, bg, num_format='0.00'),
            'ctr' : _fmt(wb, bg, align='center'),
        }

    FMT_GREEN = s2_fmts(COL_GREEN_ROW)
    FMT_AMBER = s2_fmts(COL_AMBER_ROW)
    FMT_RED   = s2_fmts(COL_RED_ROW)
    FMT_VFAIL = wb.add_format({   # verification_pass = False cell
        'bold': True, 'font_color': '#FF0000', 'bg_color': '#FFFF00',
        'align': 'center', 'valign': 'vcenter',
    })

    # ─────────────────────────────────────────────────────────────────────────
    # SHEET 1 — Stereocentre_Counts
    # ─────────────────────────────────────────────────────────────────────────
    ws1 = wb.add_worksheet("Stereocentre_Counts")
    ws1.freeze_panes(1, 0)
    ws1.set_zoom(90)

    COLS1 = [
        "family_rank",
        "is_family_start",
        "family_member_index",
        "scaffold_smiles",
        "n_stereocentres",
        "molecule_chembl_id",
        "canonical_smiles",
        "stripped_smiles_verification",
        "n_R_centres",
        "n_S_centres",
        "stereo_binary_string",
        "stereo_label",
        "bracket_notation",
        "verification_flag",
    ]
    WIDTHS1 = {
        "family_rank": 12, "is_family_start": 14, "family_member_index": 20,
        "scaffold_smiles": 45, "n_stereocentres": 14, "molecule_chembl_id": 16,
        "canonical_smiles": 55, "stripped_smiles_verification": 45,
        "n_R_centres": 12, "n_S_centres": 12,
        "stereo_binary_string": 20, "stereo_label": 14,
        "bracket_notation": 28, "verification_flag": 22,
    }
    for ci, col in enumerate(COLS1):
        ws1.write(0, ci, col, hdr)
        ws1.set_column(ci, ci, WIDTHS1.get(col, 14))

    # Pre-build column-index lookup for fast access
    df_col_idx = {col: df_main.columns.get_loc(col) for col in COLS1 if col in df_main.columns}

    for ri, row_vals in enumerate(df_main.values, start=1):
        is_start   = bool(row_vals[df_col_idx["is_family_start"]])
        mem_idx    = int(row_vals[df_col_idx["family_member_index"]])

        if is_start:
            F = FMT_START
        elif mem_idx % 2 == 1:
            F = FMT_ODD
        else:
            F = FMT_EVEN

        for ci, col in enumerate(COLS1):
            if col not in df_col_idx:
                ws1.write_blank(ri, ci, None, F['base'])
                continue

            val = row_vals[df_col_idx[col]]

            # Handle NA
            if val is None or (isinstance(val, float) and np.isnan(val)):
                ws1.write_blank(ri, ci, None, F['base'])
                continue

            # ── family_rank: bold, centred ────────────────────────────────
            if col == "family_rank":
                ws1.write_number(ri, ci, int(val), F['bold_ctr'])

            # ── is_family_start: "True" green text, "" blank ──────────────
            elif col == "is_family_start":
                if bool(val):
                    ws1.write_string(ri, ci, "True", FMT_START['green_txt'])
                else:
                    ws1.write_string(ri, ci, "", F['base'])

            # ── family_member_index: centred number ───────────────────────
            elif col == "family_member_index":
                ws1.write_number(ri, ci, int(val), F['bold_ctr'])

            # ── boolean → string ─────────────────────────────────────────
            elif isinstance(val, (bool, np.bool_)):
                ws1.write_string(ri, ci, str(bool(val)), F['str'])

            # ── force-text columns ────────────────────────────────────────
            elif col in FORCE_TEXT_COLS:
                ws1.write_string(ri, ci, str(val), F['str'])

            # ── integer ───────────────────────────────────────────────────
            elif isinstance(val, (int, np.integer)):
                ws1.write_number(ri, ci, int(val), F['num'])

            # ── float ─────────────────────────────────────────────────────
            elif isinstance(val, (float, np.floating)):
                ws1.write_number(ri, ci, float(val), F['pct'])

            # ── fallback string ───────────────────────────────────────────
            else:
                ws1.write_string(ri, ci, str(val), F['str'])

    # ─────────────────────────────────────────────────────────────────────────
    # SHEET 2 — Family_Verification
    # ─────────────────────────────────────────────────────────────────────────
    ws2 = wb.add_worksheet("Family_Verification")
    ws2.freeze_panes(1, 0)
    ws2.set_zoom(90)

    COLS2 = [
        "family_rank",
        "scaffold_smiles",
        "n_members",
        "n_stereocentres_max",
        "n_stereocentres_per_member",
        "has_mixed_stereocentre_count",
        "members_chembl_ids",
        "all_stereo_labels",
        "all_binary_strings",
        "all_bracket_notations",
        "coverage_pct",
        "missing_labels",
        "has_full_coverage",
        "verification_pass",
    ]
    WIDTHS2 = {
        "family_rank": 12, "scaffold_smiles": 45, "n_members": 10,
        "n_stereocentres_max": 18, "n_stereocentres_per_member": 28,
        "has_mixed_stereocentre_count": 26,
        "members_chembl_ids": 50, "all_stereo_labels": 30,
        "all_binary_strings": 30, "all_bracket_notations": 45,
        "coverage_pct": 13, "missing_labels": 20,
        "has_full_coverage": 16, "verification_pass": 16,
    }
    STR_COLS2 = {
        "scaffold_smiles", "members_chembl_ids", "all_stereo_labels",
        "all_binary_strings", "all_bracket_notations", "missing_labels",
        "n_stereocentres_per_member",
    }

    for ci, col in enumerate(COLS2):
        ws2.write(0, ci, col, hdr)
        ws2.set_column(ci, ci, WIDTHS2.get(col, 14))

    dv_col_idx = {col: df_verify.columns.get_loc(col)
                  for col in COLS2 if col in df_verify.columns}

    for ri, row_vals in enumerate(df_verify.values, start=1):
        n_members = int(row_vals[dv_col_idx["n_members"]])
        has_full  = bool(row_vals[dv_col_idx["has_full_coverage"]])
        vpass     = bool(row_vals[dv_col_idx["verification_pass"]])

        # Row colour
        if has_full:
            F2 = FMT_GREEN
        elif n_members == 1:
            F2 = FMT_RED
        else:
            F2 = FMT_AMBER

        for ci, col in enumerate(COLS2):
            if col not in dv_col_idx:
                ws2.write_blank(ri, ci, None, F2['base'])
                continue

            val = row_vals[dv_col_idx[col]]

            # Handle NA
            if val is None or (isinstance(val, float) and np.isnan(val)):
                ws2.write_blank(ri, ci, None, F2['base'])
                continue

            # ── coverage_pct ──────────────────────────────────────────────
            if col == "coverage_pct":
                try:
                    ws2.write_number(ri, ci, float(val), F2['pct'])
                except (ValueError, TypeError):
                    ws2.write_string(ri, ci, str(val), F2['str'])

            # ── has_full_coverage ─────────────────────────────────────────
            elif col == "has_full_coverage":
                ws2.write_string(ri, ci, str(bool(val)), F2['str'])

            # ── verification_pass: red+yellow cell if False ────────────────
            elif col == "verification_pass":
                if not bool(val):
                    ws2.write_string(ri, ci, "False", FMT_VFAIL)
                else:
                    ws2.write_string(ri, ci, "True",  F2['str'])

            # ── family_rank: centred ──────────────────────────────────────
            elif col == "family_rank":
                ws2.write_number(ri, ci, int(val), F2['ctr'])

            # ── boolean ───────────────────────────────────────────────────
            elif isinstance(val, (bool, np.bool_)):
                ws2.write_string(ri, ci, str(bool(val)), F2['str'])

            # ── string columns ────────────────────────────────────────────
            elif col in STR_COLS2:
                ws2.write_string(ri, ci, str(val), F2['str'])

            # ── integer ───────────────────────────────────────────────────
            elif isinstance(val, (int, np.integer)):
                ws2.write_number(ri, ci, int(val), F2['base'])

            # ── float ─────────────────────────────────────────────────────
            elif isinstance(val, (float, np.floating)):
                ws2.write_number(ri, ci, float(val), F2['pct'])

            # ── fallback ──────────────────────────────────────────────────
            else:
                ws2.write_string(ri, ci, str(val), F2['str'])

    wb.close()
    print(f"  Saved 2-sheet xlsx: {os.path.basename(path)}")
    print(f"    Sheet 1 'Stereocentre_Counts'  : {len(df_main):,} rows")
    print(f"    Sheet 2 'Family_Verification'  : {len(df_verify):,} families")


# ─────────────────────────────────────────────────────────────────────────────
# SINGLE-SHEET XLSX WRITER (for FAMILY_COMBINATION_COVERAGE_v4 etc.)
# ─────────────────────────────────────────────────────────────────────────────

def save_xlsx_single(df: pd.DataFrame, path: str, sheet_name: str = "Sheet1"):
    """Single-sheet xlsx — same logic as v3 save_xlsx()."""
    wb = xlsxwriter.Workbook(path, {'strings_to_numbers': False})
    ws = wb.add_worksheet(sheet_name)

    hdr = wb.add_format({
        'bold': True, 'bg_color': COL_HEADER, 'font_color': '#FFFFFF',
        'align': 'center', 'valign': 'vcenter', 'border': 1,
    })
    str_fmt = wb.add_format({'num_format': '@', 'align': 'left'})
    num_fmt = wb.add_format({'align': 'right'})
    pct_fmt = wb.add_format({'num_format': '0.00', 'align': 'right'})
    bool_t  = wb.add_format({'font_color': '#006400', 'bold': True})
    bool_f  = wb.add_format({'font_color': '#8B0000'})

    for ci, col in enumerate(df.columns):
        ws.write(0, ci, col, hdr)
        ws.set_column(ci, ci, max(len(col) + 2, 12))

    ws.freeze_panes(1, 0)

    for ri, row_vals in enumerate(df.values, start=1):
        for ci, col in enumerate(df.columns):
            val = row_vals[ci]
            if val is None or (isinstance(val, float) and np.isnan(val)):
                ws.write_blank(ri, ci, None)
                continue
            if isinstance(val, (bool, np.bool_)):
                fmt = bool_t if val else bool_f
                ws.write_string(ri, ci, str(val), fmt)
            elif col in FORCE_TEXT_COLS:
                ws.write_string(ri, ci, str(val), str_fmt)
            elif col == "coverage_pct":
                try:
                    ws.write_number(ri, ci, float(val), pct_fmt)
                except (ValueError, TypeError):
                    ws.write_string(ri, ci, str(val), str_fmt)
            elif isinstance(val, (int, np.integer)):
                ws.write_number(ri, ci, int(val), num_fmt)
            elif isinstance(val, (float, np.floating)):
                ws.write_number(ri, ci, float(val), pct_fmt)
            else:
                ws.write_string(ri, ci, str(val), str_fmt)

    wb.close()
    print(f"  Saved xlsx : {os.path.basename(path)}  ({len(df):,} rows)")


# ═════════════════════════════════════════════════════════════════════════════
# MAIN
# ═════════════════════════════════════════════════════════════════════════════

def main():
    print()
    print("=" * 67)
    print("  STEREOCENTRE ANALYSIS v4 — All_species_v3")
    print(f"  Family grouping + verification sheet  [xlsxwriter {xlsxwriter.__version__}]")
    print("=" * 67)

    # ── Step 0: Load ──────────────────────────────────────────────────────────
    print("\n[Step 0] Loading input data ...")
    df      = load_input_data()
    fold_df = load_fold_data()

    # ── Step 1: Count ─────────────────────────────────────────────────────────
    print("\n[Step 1] Counting stereocentres ...")
    cts = df["canonical_smiles"].apply(count_stereocentres)
    df["n_R_centres"]     = [c[0] for c in cts]
    df["n_S_centres"]     = [c[1] for c in cts]
    df["n_stereocentres"] = [c[2] for c in cts]

    # ── Step 2: Encode ────────────────────────────────────────────────────────
    print("[Step 2] Building binary encodings ...")
    enc = df["canonical_smiles"].apply(encode_stereo_binary)
    df["stereo_binary_string"] = [e[0] for e in enc]
    df["stereo_label"]         = [e[1] for e in enc]
    df["bracket_notation"]     = [e[2] for e in enc]

    # ── Step 2b: Verify ───────────────────────────────────────────────────────
    print("[Step 2b] Verifying count-encoding consistency ...")
    df["verification_flag"] = df.apply(verify_row, axis=1)
    n_flag = int((df["verification_flag"] != "").sum())
    if n_flag == 0:
        print(f"  PASS: all {len(df):,} rows consistent.")
    else:
        print(f"  FAIL: {n_flag:,} rows have inconsistencies!")

    # ── Step 3: Add family grouping columns (NEW in v4) ───────────────────────
    print("\n[Step 3] Adding family grouping columns ...")
    df = add_family_grouping_columns(df)
    n_families  = df["family_rank"].max()
    n_fam_start = int(df["is_family_start"].sum())
    print(f"  Total unique families: {n_families:,}")
    print(f"  is_family_start=True:  {n_fam_start:,}  (should equal families)")

    # ── Step 4: Compute family coverage ───────────────────────────────────────
    print("\n[Step 4] Computing family combination coverage ...")
    coverage_lookup = compute_family_coverage_lookup(df)

    # ── Step 5: Build Family_Verification DataFrame ───────────────────────────
    print("\n[Step 5] Building Family_Verification table ...")
    df_verify = build_family_verification_df(df, coverage_lookup)
    n_vpass  = int(df_verify["verification_pass"].sum())
    n_vfail  = len(df_verify) - n_vpass
    n_full   = int(df_verify["has_full_coverage"].sum())
    n_single = int((df_verify["n_members"] == 1).sum())
    print(f"  Family_Verification rows:        {len(df_verify):,}")
    print(f"  verification_pass = True:        {n_vpass:,}  ({n_vpass/len(df_verify)*100:.1f}%)")
    print(f"  verification_pass = False:       {n_vfail:,}  ← {'INVESTIGATE' if n_vfail > 0 else 'OK'}")
    print(f"  Families with full coverage:     {n_full:,}  ({n_full/len(df_verify)*100:.1f}%)")
    print(f"  Families with 1 member only:     {n_single:,}  ({n_single/len(df_verify)*100:.1f}%)")

    # ── Step 6: Build final column order for main sheet ───────────────────────
    COLS_MAIN = [
        "family_rank", "is_family_start", "family_member_index",
        "scaffold_smiles", "n_stereocentres", "molecule_chembl_id",
        "canonical_smiles", "stripped_smiles_verification",
        "n_R_centres", "n_S_centres",
        "stereo_binary_string", "stereo_label", "bracket_notation",
        "verification_flag",
    ]
    df_main = df[[c for c in COLS_MAIN if c in df.columns]].copy()

    # ── Step 7: Save STEREOCENTRE_COUNTS_v4.xlsx (2 sheets) ──────────────────
    print("\n[File 1] STEREOCENTRE_COUNTS_v4.xlsx  (2 sheets) ...")
    xlsx4_path = os.path.join(OUT_DIR, "STEREOCENTRE_COUNTS_v4.xlsx")
    save_two_sheet_xlsx(df_main, df_verify, xlsx4_path)

    # ── Step 8: Build FAMILY_COMBINATION_COVERAGE_v4 ─────────────────────────
    print("\n[File 2] FAMILY_COMBINATION_COVERAGE_v4.xlsx ...")
    # Build from coverage_lookup, sorted by scaffold_smiles (= family_rank order)
    fam_rows = []
    for scaffold in sorted(coverage_lookup.keys(), key=str):
        grp   = df[df["scaffold_smiles"] == scaffold]
        rank  = int(grp["family_rank"].iloc[0]) if len(grp) > 0 else 0
        cov   = coverage_lookup[scaffold]
        n_max = int(grp["n_stereocentres"].max()) if len(grp) > 0 else 0
        fam_rows.append({
            "family_rank"                 : rank,
            "scaffold_smiles"             : scaffold,
            "n_compounds_in_family"       : len(grp),
            "max_stereocentres_in_family" : n_max,
            "capped_at_5"                : n_max > 5,
            "max_combinations_theoretical": cov["max_combinations"],
            "combinations_present_count"  : cov["present_count"],
            "combinations_absent_count"   : cov["absent_count"],
            "coverage_pct"               : cov["coverage_pct"],
            "present_combinations"        : cov["present_combos"],
            "absent_combinations"         : cov["absent_combos"],
            "present_labels"             : cov["present_labels"],
            "absent_labels"              : cov["missing_labels"],
            "has_full_coverage"          : cov["has_full_coverage"],
        })

    fam_df = pd.DataFrame(fam_rows).sort_values("family_rank").reset_index(drop=True)
    fam_path = os.path.join(OUT_DIR, "FAMILY_COMBINATION_COVERAGE_v4.xlsx")
    save_xlsx_single(fam_df, fam_path, "Family_Coverage")
    save_json(fam_df, os.path.join(OUT_DIR, "FAMILY_COMBINATION_COVERAGE_v4.json"))

    # ── Step 9: COMBINATION_SUMMARY_BY_N_v4 ───────────────────────────────────
    print("\n[File 3] COMBINATION_SUMMARY_BY_N_v4.xlsx ...")
    summary_rows = []
    for n in list(range(1, 6)) + ["6+"]:
        if n == "6+":
            sub   = fam_df[fam_df["max_stereocentres_in_family"] > 5]
            max_c = "32+ (capped)"
        else:
            sub   = fam_df[fam_df["max_stereocentres_in_family"] == n]
            max_c = 2 ** n
        if len(sub) == 0:
            continue
        full_cov = int((sub["has_full_coverage"] == True).sum())
        partial  = int(((sub["combinations_present_count"] > 1) &
                        (sub["has_full_coverage"] == False)).sum())
        single   = int((sub["combinations_present_count"] == 1).sum())
        summary_rows.append({
            "n_stereocentres"               : str(n),
            "n_families"                    : len(sub),
            "max_combinations_possible"     : str(max_c),
            "families_with_full_coverage"   : full_cov,
            "families_with_partial_coverage": partial,
            "families_with_single_entry"    : single,
            "mean_coverage_pct"             : round(float(sub["coverage_pct"].mean()), 2),
            "median_coverage_pct"           : round(float(sub["coverage_pct"].median()), 2),
        })
    summary_df = pd.DataFrame(summary_rows)
    sum_path = os.path.join(OUT_DIR, "COMBINATION_SUMMARY_BY_N_v4.xlsx")
    save_xlsx_single(summary_df, sum_path, "Summary_By_N")

    # ─────────────────────────────────────────────────────────────────────────
    # TERMINAL SUMMARY
    # ─────────────────────────────────────────────────────────────────────────
    total = len(df)
    cs = df["n_stereocentres"].value_counts().sort_index()
    gc = lambda n: int(cs.get(n, 0))
    gt5 = int(df[df["n_stereocentres"] > 5].shape[0])

    print()
    print("=" * 67)
    print("  STEREOCENTRE ANALYSIS v4 — COMPLETE")
    print("=" * 67)
    print(f"\n  Compound distribution:")
    print(f"  Total compounds analysed:          {total:>8,}")
    print(f"  Compounds with 1 stereocentre:     {gc(1):>8,}")
    print(f"  Compounds with 2 stereocentres:    {gc(2):>8,}")
    print(f"  Compounds with 3 stereocentres:    {gc(3):>8,}")
    print(f"  Compounds with 4 stereocentres:    {gc(4):>8,}")
    print(f"  Compounds with 5 stereocentres:    {gc(5):>8,}")
    print(f"  Compounds with >5 stereocentres:   {gt5:>8,}")

    print(f"\n  *** VERIFICATION (from v3 — unchanged) ***")
    print(f"  Count-encoding mismatches:         {n_flag:>8,}  ({'PASS' if n_flag == 0 else 'FAIL'})")

    print()
    print("  " + "─" * 45)
    print("  v4 ADDITIONS")
    print("  " + "─" * 45)
    print(f"  Total families in grouped table:   {int(n_families):>8,}")
    print(f"  Families with is_family_start=True:{n_fam_start:>8,}  "
          f"({'= families ✓' if n_fam_start == n_families else 'MISMATCH!'})")
    print(f"  Family_Verification sheet rows:    {len(df_verify):>8,}  "
          f"({'= families ✓' if len(df_verify) == n_families else 'MISMATCH!'})")
    print(f"  verification_pass = True:          {n_vpass:>8,}  ({n_vpass/len(df_verify)*100:.1f}%)")
    print(f"  verification_pass = False:         {n_vfail:>8,}  "
          f"({'OK' if n_vfail == 0 else '← INVESTIGATE'})")
    print(f"  Families with full coverage:       {n_full:>8,}  ({n_full/len(df_verify)*100:.1f}%)")
    print(f"  Families with 1 member only:       {n_single:>8,}  ({n_single/len(df_verify)*100:.1f}%)")
    print()
    print("  HOW TO VERIFY IN NUMBERS:")
    print("  1. Open STEREOCENTRE_COUNTS_v4.xlsx")
    print("  2. Filter column 'is_family_start' = 'True'")
    print("     → One row per family; check stripped_smiles_verification")
    print("       matches scaffold_smiles exactly ✓")
    print("  3. Switch to 'Family_Verification' sheet")
    print("     → Green rows = full coverage")
    print("     → Amber rows = partial coverage")
    print("     → Red rows   = single isomer (no pairing)")
    print("     → Yellow cell in verification_pass = investigate that family")
    print("  " + "─" * 45)
    print(f"  v4 files saved (outputs/stereocentre_analysis/):")
    print(f"    STEREOCENTRE_COUNTS_v4.xlsx          (2 sheets: {len(df_main):,} rows + {len(df_verify):,} families)")
    print(f"    FAMILY_COMBINATION_COVERAGE_v4.xlsx  ({len(fam_df):,} families, family_rank sorted)")
    print(f"    FAMILY_COMBINATION_COVERAGE_v4.json")
    print(f"    COMBINATION_SUMMARY_BY_N_v4.xlsx     ({len(summary_df)} rows)")
    print("  " + "─" * 45)
    print("  v3 files: NOT OVERWRITTEN ✓")
    print("=" * 67)
    print()


if __name__ == "__main__":
    main()
