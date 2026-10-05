# StereoAtlas Pipeline

Automated 15-script Python pipeline generating the StereoAtlas dataset of
experimentally measured bioactivity differences between stereoisomers,
derived from ChEMBL v37 (1,252,051 records · 879 activity types · 155,526 true isomer pairs).

**Manuscript:** Saisaran Thangaraj & N. Arul Murugan. *StereoAtlas: A Large-Scale Dataset of Chiral Compounds with Distinct Bioactivities Curated from ChEMBL37.*
Scientific Data (2026).

**Data deposit (Zenodo):** https://doi.org/10.5281/zenodo.23158759

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

## Data

All output files are deposited at Zenodo under CC BY 4.0:

**https://doi.org/10.5281/zenodo.23158759**

The Zenodo deposit contains the full dataset (STEREOATLAS_MASTER.csv, FAMILY_MINMAX_FOLD.csv,
FAMILY_REFERENCE_FOLD.csv, by_metric.zip, and all supporting files) with detailed
descriptions of every file and column.

---

## Citation

```bibtex
@article{thangaraj2026stereoatlas,
  title   = {StereoAtlas: A Large-Scale Dataset of Chiral Compounds with
             Distinct Bioactivities Curated from ChEMBL37},
  author  = {Thangaraj, Saisaran and Murugan, N. Arul},
  year    = {2026},
  doi     = {10.5281/zenodo.23158759}
}
```

---

## Licence

**Code** (this repository): MIT License — see `LICENSE` file.

**Data** (Zenodo deposit): Creative Commons CC BY 4.0.
Data may be used, shared, and built upon with appropriate attribution.
