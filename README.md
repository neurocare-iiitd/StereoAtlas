# StereoAtlas Pipeline

Automated 15-script Python pipeline generating the StereoAtlas dataset of
experimentally measured bioactivity differences between stereoisomers,
derived from ChEMBL v37 (1,252,051 records · 879 activity types · 155,526 true isomer pairs).

**Manuscript:** Saisaran Thangaraj & N. Arul Murugan. *StereoAtlas: A Large-Scale Dataset of Chiral Compounds with Distinct Bioactivities Curated from ChEMBL37.*
Scientific Data (2026). DOI: [ZENODO_DOI_PENDING — update before final submission]

---

## Repository Structure

```
StereoAtlas_GitHub_Submission/
├── 00_descriptor_environment_setup.sh   # Environment setup
├── 01_prepare_base.py                   # Core pipeline (Steps 1–15)
├── 02_metric_files.py
├── ...
├── 15_master_merge.py
├── 16_atropisomer_detection.py          # Supplementary analyses
├── 17_atropisomer_verification.py
├── 18_stereocentre_analysis.py
├── requirements.txt
└── README.md
```

---

## Core Pipeline Scripts (01–15)

| Script | Step | Purpose |
|--------|------|---------|
| `01_prepare_base.py`         | 1  | ChEMBL v37 extraction and base data preparation (1,252,051 records) |
| `02_metric_files.py`         | 2  | Per-metric file generation (IC₅₀, Kᵢ, Kd, EC₅₀, Kₘ, Eₘₐₓ, 873 additional) |
| `03_species_files.py`        | 3  | Per-species file generation (471 organisms) |
| `04_stereoisomer_files.py`   | 4  | Stereochemical family construction via RDKit SMILES stripping (377,232 families) |
| `05_true_isomers.py`         | 5  | True isomer pair identification — same target, organism, metric (155,526 pairs) |
| `06_stereo_summary.py`       | 6  | Summary statistics across families, targets, and organisms |
| `07_final_report.py`         | 7  | Pipeline execution report |
| `08_postprocess_classify.py` | 8  | Five-tier stereoselectivity classification (None / Mild / Moderate / Strong / Extreme) |
| `09_functional_classify.py`  | 9  | Functional metric classification (Potency, Emax, Efficacy, etc.) |
| `10_fold_difference.py`      | 10 | Fold difference engine — Strategy 1 MinMax and Strategy 2 Reference-Form |
| `11_family_view.py`          | 11 | Family-level fold views (MinMax and Reference-Form strategies only) |
| `12_literature_enrichment.py`| 12 | Literature cross-reference against PubMed; generates TRUE_ISOMER_LITERATURE_CLEAN.csv |
| `13_target_stereoselectivity.py` | 13 | Target Stereoselectivity Index (TSI, 0–100) across 2,732 targets |
| `14_dataset_coverage_score.py`   | 14 | Dataset Coverage Score (DCS, 0–100) across all targets |
| `15_master_merge.py`         | 15 | Master merge of all pipeline outputs into STEREOATLAS_MASTER.csv (60.4 MB) |

---

## Supplementary Analysis Scripts (16–18)

| Script | Purpose |
|--------|---------|
| `16_atropisomer_detection.py`    | SMARTS-based atropisomer scaffold flagging (388 families detected) |
| `17_atropisomer_verification.py` | Atropisomer verification via optical rotation naming and known-drug cross-reference |
| `18_stereocentre_analysis.py`    | Stereocentre enumeration and configurational coverage (76,576 compound–family entries) |

---

## Setup

**Step 1 — Install dependencies**

```bash
bash 00_descriptor_environment_setup.sh
```

Or manually:

```bash
pip install rdkit==2023.09.6 pandas==3.0.2 numpy==2.4.4 openpyxl xlsxwriter tqdm
```

**Step 2 — Place ChEMBL v37 SQLite database**

Download `chembl_37.db` from https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/chembl_37/
and place it in the pipeline root directory.

**Step 3 — Run the pipeline**

```bash
python 01_prepare_base.py
python 02_metric_files.py
python 03_species_files.py
python 04_stereoisomer_files.py
python 05_true_isomers.py
python 06_stereo_summary.py
python 07_final_report.py
python 08_postprocess_classify.py
python 09_functional_classify.py
python 10_fold_difference.py
python 11_family_view.py
python 12_literature_enrichment.py
python 13_target_stereoselectivity.py
python 14_dataset_coverage_score.py
python 15_master_merge.py
```

Each script reads only from `outputs/` of the previous step — no database
access after Script 01. Scripts 16–18 can be run independently after Script 05.

**Smoke-test mode** (validates outputs without a full run):

```bash
python 01_prepare_base.py --smoke-test
```

---

## Requirements

Tested on Ubuntu 20.04 · CPU-only · reproducible on any standard workstation.

```
python==3.14.3
rdkit==2023.09
pandas==3.0.2
numpy==2.4.4
openpyxl>=3.1
xlsxwriter>=3.1
tqdm>=4.65
```

---

## Key Outputs

| File | Size | Description |
|------|------|-------------|
| `FAMILY_MINMAX_FOLD.csv`           | 10.7 MB  | One max fold difference per family–target–organism–metric group. Primary analysis file. |
| `FAMILY_REFERENCE_FOLD.csv`        | 13.3 MB  | Fold differences relative to a designated reference stereoform per family. |
| `TRUE_ISOMER_LITERATURE_CLEAN.csv` | 32.4 MB  | Literature-validated pairs with PubMed IDs and DOIs (Script 12). |
| `STEREOATLAS_MASTER.csv`           | 60.4 MB  | Full merged master file combining all pipeline outputs (Script 15). |
| `STEREOATLAS_MASTER_LITE.csv`      | 31.2 MB  | Filtered master: true isomer pairs with defined stereo class and fold ≥ 2× (Script 15). |
| `by_metric/`                       | ~30 MB   | Per-metric breakdowns (IC₅₀, Kᵢ, Kd, EC₅₀, and 875 additional types). |
| `STEREOCENTRE_COUNTS_v4.xlsx`      | 8.1 MB   | Per-compound stereocentre counts with configurational coverage (Script 18). |
| `MISSING_COMBINATIONS.csv`         | ~6 MB    | Absent stereoisomer configurations — candidate synthesis targets (Script 18). |
| `TARGET_STEREOSELECTIVITY_SUMMARY.csv` | 1.1 MB | TSI scores across 2,732 targets (Script 13). |
| `TARGET_PRIORITY_RANKING.csv`      | 1.1 MB   | Targets ranked by TSI descending (Script 13). |
| `TARGET_DATASET_COVERAGE.csv`      | 713 KB   | DCS scores and coverage classifications (Script 14). |
| `ATROPISOMER_FOLD_SUMMARY.csv`     | 119 KB   | Candidate atropisomeric families with fold differences (Scripts 16–17). |
| `protein_class_extreme_summary.csv`| Variable | Protein class breakdown of Extreme (≥50×) stereoselectivity pairs. |

---

## Quick Start — Python

```python
import pandas as pd

# Load primary analysis file
df = pd.read_csv('FAMILY_MINMAX_FOLD.csv', low_memory=False)

# Filter to Extreme stereoselectivity pairs (fold ≥ 50×)
extreme = df[df['stereoselectivity_class'] == 'Extreme']
print(f"Extreme pairs: {len(extreme):,}")

# Top targets by pair count
top = df.groupby('target_chembl_id').size().sort_values(ascending=False).head(10)
print(top)
```

## Quick Start — R

```r
library(data.table)
df <- fread('FAMILY_MINMAX_FOLD.csv')
extreme <- df[stereoselectivity_class == 'Extreme']
cat("Extreme pairs:", nrow(extreme), "\n")
```

---

## Stereoselectivity Tiers

| Tier | Fold Difference | Biological interpretation |
|------|----------------|--------------------------|
| None | < 2× | Within typical assay variability |
| Mild | 2–5× | Weak but reproducible stereoselectivity |
| Moderate | 5–10× | Notable stereoselectivity |
| Strong | 10–50× | Clinically relevant regime |
| Extreme | > 50× | Chiral-switch justification threshold |

---

## Dataset Statistics

| Metric | Value |
|--------|-------|
| Raw ChEMBL records | 1,252,051 |
| Records after QC | 995,043 (79.5% retention) |
| Unique compounds | 420,721 |
| Stereochemical families | 377,232 |
| True isomer pairs | 155,526 |
| IC₅₀ pairs | 38,446 |
| Biological targets | 4,519 |
| Organisms | 471 |
| Activity types | 879 |
| Maximum fold difference | 195,000,000× (oxytocin receptor, Kd) |
| Homo sapiens fraction | 74.6% |

---

## Data Availability

Full dataset deposited at Zenodo (CC BY 4.0):
**[ZENODO_DOI_PENDING — update before final submission]**

---

## Citation

```bibtex
@article{thangaraj2026stereoatlas,
  title   = {StereoAtlas: A Large-Scale Dataset of Chiral Compounds with
             Distinct Bioactivities Curated from ChEMBL37},
  author  = {Thangaraj, Saisaran and Murugan, N. Arul},
  journal = {Scientific Data},
  year    = {2026},
  doi     = {ZENODO_DOI_PENDING}
}
```

---

## Licence

**Code** (this repository): MIT License — see `LICENSE` file.

**Data** (Zenodo deposit): Creative Commons CC BY 4.0.
Data may be used, shared, and built upon with appropriate attribution.

## New in This Version: Single-form Family Files

Three additional files describe families where only one stereoisomer was tested:

| File | Description |
|------|-------------|
| `ONLY_R_FAMILIES.csv` | 99,467 families where only the all-R configuration was tested (n=1–5) |
| `ONLY_S_FAMILIES.csv` | 88,356 families where only the all-S configuration was tested (n=1–5) |
| `MIXED_SINGLE_FAMILIES.csv` | 114,573 families (n≥2) where only one mixed R/S stereoform was tested |

These files correct the previously inflated coverage estimates that arose from
considering paired families alone. The corrected mean configurational coverage
across all 377,204 stereochemical families is:

| n stereocentres | Theoretical combos | Total families | Corrected coverage |
|---|---|---|---|
| 1 | 2 | 138,781 | 55.1% |
| 2 | 4 | 106,034 | 28.1% |
| 3 | 8 | 46,550 | 14.0% |
| 4 | 16 | 28,420 | 6.9% |
| 5 | 32 | 15,572 | 3.4% |

## Output Files

All data files are deposited on Zenodo (DOI: [10.5281/zenodo.XXXXXXX](https://doi.org/10.5281/zenodo.XXXXXXX)).

| File | Size | Description |
|------|------|-------------|
| `STEREOATLAS_MASTER.csv` | 60.4 MB | Full merged dataset: all stereoisomer pairs with bioactivity, fold differences, and classifications |
| `STEREOATLAS_MASTER_LITE.csv` | 31.2 MB | Lightweight version with key columns only |
| `FAMILY_MINMAX_FOLD.csv` | 10.7 MB | S1 strategy: one representative MinMax fold row per (family × target × metric) |
| `FAMILY_REFERENCE_FOLD.csv` | 13.3 MB | S2 strategy: fold vs reference (least potent form) per group |
| `TRUE_ISOMER_LITERATURE_CLEAN.csv` | 32.4 MB | Literature-enriched dataset with PubMed references |
| `by_metric/` | ~30 MB | Per-metric breakdowns (IC50, Ki, EC50, Kd, potency) |
| `STEREOCENTRE_COUNTS_v4.xlsx` | 8.1 MB | Per-compound stereocentre enumeration and family verification |
| `MISSING_COMBINATIONS.csv` | 28.4 MB | Untested stereoform combinations per family |
| `ONLY_R_FAMILIES.csv` | 39.4 MB | Families with only the all-R form tested |
| `ONLY_S_FAMILIES.csv` | 34.5 MB | Families with only the all-S form tested |
| `MIXED_SINGLE_FAMILIES.csv` | 94.8 MB | Families with only one mixed-config form tested |
| `TARGET_STEREOSELECTIVITY_SUMMARY.csv` | 1.1 MB | Per-target stereoselectivity summary |
| `TARGET_PRIORITY_RANKING.csv` | 1.1 MB | Targets ranked by stereoselectivity evidence |
| `TARGET_DATASET_COVERAGE.csv` | 713 KB | Coverage score per target |
| `ATROPISOMER_FOLD_SUMMARY.csv` | 119 KB | Atropisomer pair fold differences |
| `PROBABLE_ATROPISOMERS_BY_SCAFFOLD.csv` | 290 KB | Predicted atropisomers by scaffold |
