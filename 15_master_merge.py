"""
All-Species v3 -- Literature Enrichment for TRUE Stereoisomer Pairs
====================================================================
Script:  14_master_merge.py
Source:  outputs/true_isomers/TRUE_ISOMER_PAIRS.csv
Output:  outputs/literature_enrichment/

Enrichment strategy:
  1. Load TRUE_ISOMER_PAIRS.csv (compound_chembl_id, target_chembl_id, activity_type present)
  2. Batch-query ChEMBL REST API (official, no scraping) for document_chembl_id per
     (molecule_chembl_id, target_chembl_id, standard_type) activity records.
     Batch size: 50 compound IDs per request.
  3. Fetch document metadata (pubmed_id, doi, journal, year, title, authors)
     for each unique document_chembl_id found.
  4. Aggregate publication metadata per (compound x target x activity_type).
  5. Classify reference quality: HIGH (>=3 pubs) / MEDIUM (2) / LOW (1) / NO_REFERENCE (0).
  6. Generate PubMed URLs for web-server use.

All API results are cached in:
  outputs/literature_enrichment/.cache/
  So repeat runs do NOT re-fetch already retrieved data.

No DB access. No PDF scraping. No unofficial sources.
Uses only: https://www.ebi.ac.uk/chembl/api/data/
"""

import pandas as pd
import numpy as np
import os

# ── PATH CONFIGURATION ──────────────────────────────────────────────
# Set BASE_DIR to the root of your pipeline output directory
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
# ────────────────────────────────────────────────────────────────────, sys, time, json, re
import urllib.request, urllib.parse
from datetime import datetime
from tqdm import tqdm

# ---- Paths ---------------------------------------------------------------
V3_DIR   = os.path.dirname(os.path.abspath(__file__))
SRC_FILE = os.path.join(V3_DIR, "outputs", "true_isomers", "TRUE_ISOMER_PAIRS.csv")
OUT_DIR  = os.path.join(V3_DIR, "outputs", "literature_enrichment")
CACHE_DIR= os.path.join(OUT_DIR, ".cache")
REP_DIR  = os.path.join(V3_DIR, "reports")
os.makedirs(OUT_DIR,  exist_ok=True)
os.makedirs(CACHE_DIR,exist_ok=True)

# ---- ChEMBL API ----------------------------------------------------------
CHEMBL_BASE  = "https://www.ebi.ac.uk/chembl/api/data"
BATCH_SIZE   = 50        # compound IDs per activity request
ACTIVITY_LIMIT = 1000    # max activities to retrieve per batch
REQUEST_DELAY  = 0.35    # seconds between requests (rate-limiting)
MAX_RETRIES    = 3

# ---- Cache helpers -------------------------------------------------------
def cache_path(key_type, key_val):
    safe = re.sub(r"[^\w]", "_", str(key_val))[:80]
    return os.path.join(CACHE_DIR, f"{key_type}__{safe}.json")

def cache_load(key_type, key_val):
    path = cache_path(key_type, key_val)
    if os.path.exists(path):
        with open(path) as f:
            return json.load(f)
    return None

def cache_save(key_type, key_val, data):
    with open(cache_path(key_type, key_val), "w") as f:
        json.dump(data, f)

# ---- API request helper --------------------------------------------------
def api_get(url, retries=MAX_RETRIES):
    """GET a ChEMBL REST URL, return parsed JSON or None on failure."""
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=15) as resp:
                return json.loads(resp.read())
        except Exception as e:
            if attempt < retries - 1:
                time.sleep(1.5 * (attempt + 1))
            else:
                return None
    return None

# ==========================================================================
# STEP 1: Load source data
# ==========================================================================
def load_source():
    print(f"  Source: {SRC_FILE}")
    ti = pd.read_csv(SRC_FILE, low_memory=False)
    print(f"  Rows              : {len(ti):,}")
    print(f"  Unique compounds  : {ti['compound_chembl_id'].nunique():,}")
    print(f"  Unique targets    : {ti['target_chembl_id'].nunique():,}")
    print(f"  Unique organisms  : {ti['organism'].nunique():,}")
    return ti

# ==========================================================================
# STEP 2: Batch-fetch activity -> document_chembl_id mapping
# ==========================================================================
def fetch_activity_docs_batch(compound_ids):
    """
    Query ChEMBL activity API for a batch of compound IDs.
    Returns list of activity dicts with document info.
    Handles pagination automatically.
    """
    results = []
    # Comma-delimited for ChEMBL __in filter
    mol_filter = ",".join(compound_ids)
    fields = "molecule_chembl_id,target_chembl_id,standard_type,document_chembl_id,document_journal,document_year"
    offset = 0
    while True:
        url = (
            f"{CHEMBL_BASE}/activity.json?"
            f"molecule_chembl_id__in={urllib.parse.quote(mol_filter)}"
            f"&limit={ACTIVITY_LIMIT}&offset={offset}"
            f"&fields={fields}"
        )
        data = api_get(url)
        if not data or "activities" not in data:
            break
        batch_acts = data["activities"]
        for act in batch_acts:
            if act.get("document_chembl_id"):
                results.append({
                    "molecule_chembl_id": act.get("molecule_chembl_id",""),
                    "target_chembl_id":   act.get("target_chembl_id",""),
                    "standard_type":      act.get("standard_type",""),
                    "document_chembl_id": act.get("document_chembl_id",""),
                    "document_journal":   act.get("document_journal",""),
                    "document_year":      act.get("document_year",""),
                })
        total_count = data.get("page_meta",{}).get("total_count", 0)
        offset += ACTIVITY_LIMIT
        if offset >= total_count:
            break
        time.sleep(REQUEST_DELAY)
    return results

def build_activity_doc_map(ti):
    """
    Batch compounds in groups of BATCH_SIZE, fetch all activity records
    (no target/type filter — let ChEMBL return all docs for those compounds).
    Filter client-side to only keep entries for compounds/targets/types in our dataset.
    Returns DataFrame: molecule_chembl_id, target_chembl_id, standard_type,
                       document_chembl_id, document_journal, document_year
    """
    print("  Building activity -> document mapping via ChEMBL REST API...")
    print(f"  Batch size: {BATCH_SIZE} compounds per request")

    unique_mols = ti["compound_chembl_id"].dropna().unique().tolist()
    print(f"  Unique compounds to fetch: {len(unique_mols):,}")

    # Build allowed set for client-side filter
    allowed = set(
        zip(ti["compound_chembl_id"], ti["target_chembl_id"], ti["activity_type"])
    )

    all_rows = []
    cache_hits = 0
    api_calls  = 0
    batches = [unique_mols[i:i+BATCH_SIZE] for i in range(0, len(unique_mols), BATCH_SIZE)]

    for batch_idx, batch in enumerate(tqdm(batches,
                                           desc="  Fetching activity docs",
                                           unit="batch")):
        cache_key = f"batch_{batch_idx}"
        cached = cache_load("actdoc", cache_key)
        if cached is not None:
            all_rows.extend(cached)
            cache_hits += 1
        else:
            time.sleep(REQUEST_DELAY)
            results = fetch_activity_docs_batch(batch)
            cache_save("actdoc", cache_key, results)
            all_rows.extend(results)
            api_calls += 1

    print(f"  API calls made       : {api_calls:,}")
    print(f"  Cache hits           : {cache_hits:,}")
    print(f"  Raw activity records : {len(all_rows):,}")

    if not all_rows:
        return pd.DataFrame(columns=[
            "molecule_chembl_id","target_chembl_id","standard_type",
            "document_chembl_id","document_journal","document_year"
        ])

    df = pd.DataFrame(all_rows).drop_duplicates()
    # Client-side filter: keep only (compound, target, activity_type) combos in our dataset
    df["_key"] = list(zip(df["molecule_chembl_id"], df["target_chembl_id"], df["standard_type"]))
    df = df[df["_key"].isin(allowed)].drop(columns=["_key"])
    print(f"  After filtering to dataset: {len(df):,} records")
    return df

# ==========================================================================
# STEP 3: Fetch document metadata (pubmed_id, doi, journal, year, title)
# ==========================================================================
def fetch_document_metadata(document_chembl_ids):
    """
    For each unique document_chembl_id, fetch full metadata from ChEMBL.
    Returns DataFrame with one row per document_chembl_id.
    """
    unique_docs = list(set(d for d in document_chembl_ids if d and str(d) != "nan"))
    print(f"  Fetching metadata for {len(unique_docs):,} unique documents...")

    doc_rows = []
    cache_hits = 0
    api_calls  = 0

    for doc_id in tqdm(unique_docs, desc="  Fetching doc metadata", unit="doc"):
        cached = cache_load("doc", doc_id)
        if cached is not None:
            doc_rows.append(cached)
            cache_hits += 1
            continue

        url = f"{CHEMBL_BASE}/document/{doc_id}.json"
        time.sleep(REQUEST_DELAY)
        data = api_get(url)
        api_calls += 1

        if data:
            rec = {
                "document_chembl_id": doc_id,
                "pubmed_id":   str(data.get("pubmed_id","") or ""),
                "doi":         str(data.get("doi","") or ""),
                "journal":     str(data.get("journal","") or ""),
                "year":        str(data.get("year","") or ""),
                "title":       str(data.get("title","") or ""),
                "authors":     str(data.get("authors","") or ""),
                "doc_type":    str(data.get("doc_type","") or ""),
            }
        else:
            rec = {
                "document_chembl_id": doc_id,
                "pubmed_id":"","doi":"","journal":"","year":"",
                "title":"","authors":"","doc_type":"",
            }
        cache_save("doc", doc_id, rec)
        doc_rows.append(rec)

    print(f"  API calls made: {api_calls:,}  |  Cache hits: {cache_hits:,}")
    return pd.DataFrame(doc_rows) if doc_rows else pd.DataFrame()

# ==========================================================================
# STEP 4: Aggregate publication metadata per (compound × target × activity_type)
# ==========================================================================
def aggregate_pubs(ti, act_doc_df, doc_meta_df):
    """
    Join activity->doc->metadata back to TRUE_ISOMERS.
    For each row in ti, aggregate all linked documents into semicolon-separated strings.
    """
    if act_doc_df.empty or doc_meta_df.empty:
        print("  WARNING: No document data available from API. Filling with NO_REFERENCE.")
        ti = ti.copy()
        for col in ["publication_count","pubmed_ids","dois","document_chembl_ids",
                    "journals","publication_years","reference_links","reference_quality"]:
            ti[col] = "NO_REFERENCE" if col == "reference_quality" else ("" if "count" not in col else 0)
        return ti

    # Merge act_doc with doc_meta
    enriched_docs = act_doc_df.merge(
        doc_meta_df,
        on="document_chembl_id",
        how="left"
    )

    # Aggregate per (molecule_chembl_id, target_chembl_id, standard_type)
    def join_unique(series):
        vals = [str(v) for v in series.dropna() if str(v).strip() and str(v) != "nan"]
        return ";".join(sorted(set(vals)))

    pub_agg = (
        enriched_docs
        .groupby(["molecule_chembl_id","target_chembl_id","standard_type"])
        .agg(
            publication_count    =("document_chembl_id","nunique"),
            document_chembl_ids  =("document_chembl_id", join_unique),
            pubmed_ids           =("pubmed_id",          join_unique),
            dois                 =("doi",                join_unique),
            journals             =("journal",            join_unique),
            publication_years    =("year",               join_unique),
        )
        .reset_index()
        .rename(columns={
            "molecule_chembl_id": "compound_chembl_id",
            "standard_type":      "activity_type",
        })
    )

    # Join back to TRUE_ISOMERS
    ti_enriched = ti.merge(
        pub_agg,
        on=["compound_chembl_id","target_chembl_id","activity_type"],
        how="left"
    )

    # Fill missing
    ti_enriched["publication_count"] = ti_enriched["publication_count"].fillna(0).astype(int)
    for col in ["document_chembl_ids","pubmed_ids","dois","journals","publication_years"]:
        ti_enriched[col] = ti_enriched[col].fillna("")

    # PubMed URLs
    def build_links(pmid_str):
        if not pmid_str or pmid_str == "nan":
            return ""
        pmids = [p.strip() for p in pmid_str.split(";") if p.strip() and p.strip() != "nan"]
        return ";".join(f"https://pubmed.ncbi.nlm.nih.gov/{p}/" for p in pmids if p.isdigit())

    ti_enriched["reference_links"] = ti_enriched["pubmed_ids"].apply(build_links)

    # Reference quality classification
    def ref_quality(n):
        if n >= 3: return "HIGH"
        if n == 2: return "MEDIUM"
        if n == 1: return "LOW"
        return "NO_REFERENCE"

    ti_enriched["reference_quality"] = ti_enriched["publication_count"].apply(ref_quality)

    return ti_enriched

# ==========================================================================
# STEP 5: Build output tables
# ==========================================================================
def save_main_output(ti_enriched):
    """TRUE_ISOMER_LITERATURE.csv"""
    cols = [
        "stereochemical_family_id",
        "compound_chembl_id","compound_name",
        "target_chembl_id","target_name",
        "organism",
        "activity_type","activity_units",
        "stereo_form","stereochemical_variant",
        "mean_activity_value",
        "publication_count",
        "pubmed_ids","dois",
        "document_chembl_ids",
        "journals","publication_years",
        "reference_quality","reference_links",
    ]
    out_cols = [c for c in cols if c in ti_enriched.columns]
    out = ti_enriched[out_cols].copy()
    path = os.path.join(OUT_DIR, "TRUE_ISOMER_LITERATURE.csv")
    out.to_csv(path, index=False)
    print(f"  Saved TRUE_ISOMER_LITERATURE.csv ({len(out):,} rows)")
    return out

def save_reference_summary(ti_enriched):
    """TRUE_ISOMER_REFERENCE_SUMMARY.csv — one row per compound"""
    summary = (
        ti_enriched.groupby(["compound_chembl_id","compound_name"])
        .agg(
            total_publications=("publication_count","max"),
            total_targets     =("target_name","nunique"),
            total_species     =("organism","nunique"),
            total_metrics     =("activity_type","nunique"),
        )
        .reset_index()
        .sort_values("total_publications", ascending=False)
    )
    path = os.path.join(OUT_DIR, "TRUE_ISOMER_REFERENCE_SUMMARY.csv")
    summary.to_csv(path, index=False)
    print(f"  Saved TRUE_ISOMER_REFERENCE_SUMMARY.csv ({len(summary):,} rows)")
    return summary

def save_top100_referenced(ti_enriched):
    """TOP100_MOST_REFERENCED_STEREOISOMERS.csv"""
    top = (
        ti_enriched.groupby(["compound_chembl_id","compound_name"])
        .agg(
            publication_count=("publication_count","max"),
            target_count     =("target_name","nunique"),
            organism_count   =("organism","nunique"),
            metric_count     =("activity_type","nunique"),
        )
        .reset_index()
        .sort_values("publication_count", ascending=False)
        .head(100)
    )
    path = os.path.join(OUT_DIR, "TOP100_MOST_REFERENCED_STEREOISOMERS.csv")
    top.to_csv(path, index=False)
    print(f"  Saved TOP100_MOST_REFERENCED_STEREOISOMERS.csv ({len(top):,} rows)")
    return top

def save_quality_breakdown(ti_enriched):
    """REFERENCE_QUALITY_BREAKDOWN.csv"""
    qb = (
        ti_enriched.groupby(["reference_quality","activity_type"])
        .agg(
            pair_count   =("compound_chembl_id","count"),
            compound_count=("compound_chembl_id","nunique"),
            target_count =("target_name","nunique"),
        )
        .reset_index()
    )
    path = os.path.join(OUT_DIR, "REFERENCE_QUALITY_BREAKDOWN.csv")
    qb.to_csv(path, index=False)
    print(f"  Saved REFERENCE_QUALITY_BREAKDOWN.csv ({len(qb):,} rows)")
    return qb

# ==========================================================================
# STEP 6: Report
# ==========================================================================
def write_report(ti_enriched, top100, elapsed):
    lines = []
    def p(s=""): lines.append(s)

    p("=" * 70)
    p("  All-Species v3 -- Literature Enrichment Report")
    p(f"  Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    p("=" * 70)
    p()
    p("SOURCE")
    p(f"  File: {SRC_FILE}")
    p(f"  Total TRUE isomer pairs processed: {len(ti_enriched):,}")

    # Reference quality distribution
    p()
    p("REFERENCE QUALITY DISTRIBUTION")
    qd = ti_enriched["reference_quality"].value_counts()
    for q, n in qd.items():
        pct = 100*n/len(ti_enriched)
        p(f"  {q:<15s}: {n:>8,} pairs ({pct:.1f}%)")

    # Pairs with/without references
    with_ref = (ti_enriched["publication_count"] > 0).sum()
    without  = (ti_enriched["publication_count"] == 0).sum()
    p()
    p("COVERAGE")
    p(f"  Pairs with references   : {with_ref:,}")
    p(f"  Pairs without references: {without:,}")
    p(f"  Coverage %              : {100*with_ref/max(len(ti_enriched),1):.1f}%")

    # PubMed and DOI counts
    has_pmid = ti_enriched["pubmed_ids"].fillna("").str.len() > 0
    has_doi  = ti_enriched["dois"].fillna("").str.len() > 0
    total_pmids = ti_enriched["pubmed_ids"].fillna("").str.split(";").apply(
        lambda x: len([v for v in x if v.strip() and v.strip().isdigit()])
    ).sum()
    p()
    p("PUBLICATION IDENTIFIERS")
    p(f"  Rows with PubMed IDs    : {has_pmid.sum():,}")
    p(f"  Rows with DOIs          : {has_doi.sum():,}")
    p(f"  Total PubMed ID links   : {int(total_pmids):,}")

    # Top referenced compounds
    p()
    p("TOP 10 MOST-REFERENCED STEREOISOMERS")
    for _, r in top100.head(10).iterrows():
        name = str(r["compound_name"])[:35] if pd.notna(r["compound_name"]) else str(r["compound_chembl_id"])
        p(f"  {name:<37s}: {int(r['publication_count'])} pubs  "
          f"| {int(r['target_count'])} targets  "
          f"| {int(r['organism_count'])} organisms")

    # Top targets
    p()
    p("TOP 10 MOST-REFERENCED TARGETS")
    tgt_summary = (
        ti_enriched.groupby("target_name")["publication_count"]
        .agg(["sum","mean","max"])
        .sort_values("max", ascending=False)
        .head(10)
    )
    for tgt, row in tgt_summary.iterrows():
        p(f"  {str(tgt)[:45]:<47s}: max={int(row['max'])} pubs  mean={row['mean']:.1f}")

    p()
    p(f"  Total runtime: {elapsed:.1f}s")
    p("=" * 70)

    report_path = os.path.join(REP_DIR, "12_literature_enrichment_report.txt")
    with open(report_path, "w") as f:
        f.write("\n".join(lines))
    print(f"\n  Report saved: {report_path}")

# ==========================================================================
# MAIN
# ==========================================================================
def main():
    print("=" * 70)
    print("14_master_merge.py -- ChEMBL Literature Enrichment")
    print(f"  API  : {CHEMBL_BASE}")
    print(f"  Cache: {CACHE_DIR}")
    print("=" * 70)
    t0 = time.time()

    # ---- 1. Load ---------------------------------------------------------
    print("\n[1/6] Loading TRUE_ISOMER_PAIRS.csv...")
    try:
        ti = load_source()
    except Exception as e:
        print(f"  ERROR: {e}"); sys.exit(1)

    # ---- 2. Activity -> document mapping ---------------------------------
    print("\n[2/6] Fetching activity -> document_chembl_id (batched, cached)...")
    try:
        act_doc_df = build_activity_doc_map(ti)
    except Exception as e:
        print(f"  ERROR in activity doc fetch: {e}")
        act_doc_df = pd.DataFrame()

    # ---- 3. Document metadata --------------------------------------------
    print("\n[3/6] Fetching document metadata (pubmed_id, doi, journal, year)...")
    try:
        if not act_doc_df.empty and "document_chembl_id" in act_doc_df.columns:
            doc_meta_df = fetch_document_metadata(act_doc_df["document_chembl_id"].dropna())
        else:
            doc_meta_df = pd.DataFrame()
            print("  No documents to fetch.")
    except Exception as e:
        print(f"  ERROR in doc metadata fetch: {e}")
        doc_meta_df = pd.DataFrame()

    # Save intermediate caches as CSV for inspection
    if not act_doc_df.empty:
        act_doc_df.to_csv(os.path.join(OUT_DIR, "activity_document_map.csv"), index=False)
        print(f"  Saved activity_document_map.csv ({len(act_doc_df):,} rows)")
    if not doc_meta_df.empty:
        doc_meta_df.to_csv(os.path.join(OUT_DIR, "document_metadata.csv"), index=False)
        print(f"  Saved document_metadata.csv ({len(doc_meta_df):,} rows)")

    # ---- 4. Aggregate publications per compound/target/metric -----------
    print("\n[4/6] Aggregating publication metadata per compound-target-activity...")
    try:
        ti_enriched = aggregate_pubs(ti, act_doc_df, doc_meta_df)
        q_dist = ti_enriched["reference_quality"].value_counts()
        for q, n in q_dist.items():
            print(f"  {q:<15s}: {n:,} pairs")
    except Exception as e:
        print(f"  ERROR in aggregation: {e}")
        ti_enriched = ti.copy()
        ti_enriched["publication_count"] = 0
        for col in ["pubmed_ids","dois","document_chembl_ids","journals",
                    "publication_years","reference_links"]:
            ti_enriched[col] = ""
        ti_enriched["reference_quality"] = "NO_REFERENCE"

    # ---- 5. Save outputs -------------------------------------------------
    print("\n[5/6] Saving output files...")
    try:
        save_main_output(ti_enriched)
        save_reference_summary(ti_enriched)
        save_top100_referenced(ti_enriched)
        save_quality_breakdown(ti_enriched)
    except Exception as e:
        print(f"  ERROR saving outputs: {e}")

    # ---- 6. Report -------------------------------------------------------
    print("\n[6/6] Writing report...")
    elapsed = time.time() - t0
    try:
        top100 = pd.read_csv(os.path.join(OUT_DIR, "TOP100_MOST_REFERENCED_STEREOISOMERS.csv"))
        write_report(ti_enriched, top100, elapsed)
    except Exception as e:
        print(f"  ERROR writing report: {e}")

    # ---- Manifest --------------------------------------------------------
    print(f"\n  Output files:")
    for f in sorted(os.listdir(OUT_DIR)):
        if f.endswith(".csv"):
            sz = os.path.getsize(os.path.join(OUT_DIR, f)) / 1e6
            n  = sum(1 for _ in open(os.path.join(OUT_DIR, f))) - 1
            print(f"  {f:<55s} {n:>8,} rows  {sz:.1f} MB")
    cache_files = len([f for f in os.listdir(CACHE_DIR) if f.endswith(".json")])
    print(f"  .cache/: {cache_files:,} cached API responses")

    print(f"\n  Total runtime: {elapsed:.1f}s ({elapsed/60:.1f} min)")
    print(f"\n[OK] 14_master_merge.py complete.")

if __name__ == "__main__":
    main()
