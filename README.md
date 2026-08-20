# Computational Drug Discovery

End-to-end template for building a machine-learning model of **human acetylcholinesterase (AChE, CHEMBL220)** inhibitors from [ChEMBL](https://www.ebi.ac.uk/chembl/) bioactivity data.

The workflow follows the classic Computational Drug Discovery series (data collection → EDA → fingerprints → models), but the logic lives in a Python package so notebooks stay thin and each step is rerunnable from the command line.

## Pipeline

| Part | Notebook | Package entry | Status |
|------|----------|---------------|--------|
| 1. Bioactivity data | `notebooks/01_bioactivity_data.ipynb` | `python -m drug_discovery part1` | Implemented |
| 2. Exploratory data analysis | `notebooks/02_exploratory_data_analysis.ipynb` | `python -m drug_discovery part2` | Template |
| 3. Descriptor calculation | `notebooks/03_descriptor_calculation.ipynb` | `python -m drug_discovery part3` | Template |
| 4. Regression models | `notebooks/04_regression_models.ipynb` | `python -m drug_discovery part4` | Template |
| 5. Compare algorithms | `notebooks/05_compare_algorithms.ipynb` | `python -m drug_discovery part5` | Template |

## Project layout

```
drug-discovery/
├── notebooks/          # narrative + figures
├── src/drug_discovery/ # reusable pipeline code
├── tests/              # unit tests for cleaning/labeling
├── data/raw/           # ChEMBL download
├── data/processed/     # curated tables for later parts
├── data/descriptors/   # fingerprint matrices (Part 3)
├── models/             # fitted models (Parts 4–5)
├── reports/figures/    # EDA and model plots
└── archive/            # original Colab notebook
```

## Setup

Python 3.10+. From the repo root:

```bash
python3.10 -m venv .venv
source .venv/bin/activate
pip install -e ".[notebooks,dev]"
```

Later parts also need:

```bash
pip install -e ".[eda,ml]"
```

Or install everything listed in `requirements.txt`.

## Run Part 1

Downloads IC50 records for CHEMBL220, drops missing/duplicate structures, and labels compounds:

- **active** — IC50 ≤ 1000 nM
- **inactive** — IC50 ≥ 10,000 nM
- **intermediate** — in between

```bash
python -m drug_discovery part1
```

Reuse an existing download:

```bash
python -m drug_discovery part1 --skip-fetch
```

Outputs:

- `data/raw/acetylcholinesterase_01_bioactivity_data_raw.csv`
- `data/processed/acetylcholinesterase_02_bioactivity_data_preprocessed.csv`
- `data/processed/acetylcholinesterase_03_bioactivity_data_curated.csv`

To point the same pipeline at another protein, edit `src/drug_discovery/config.py`.

## Tests

```bash
pytest
```
