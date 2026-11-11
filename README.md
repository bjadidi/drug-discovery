# Computational Drug Discovery

End-to-end QSAR pipeline that predicts the potency (pIC50) of **human acetylcholinesterase (AChE, CHEMBL220)** inhibitors from [ChEMBL](https://www.ebi.ac.uk/chembl/) bioactivity data. AChE is the target of Alzheimer's drugs such as donepezil.

The steps follow the classic computational drug discovery workflow (data collection → EDA → fingerprints → models), but with two changes that matter in practice:

* **Everything is a tested Python package**, not a pile of notebooks: each step reruns from the command line, and notebooks only narrate.
* **Models are judged on a scaffold split** as well as a random split. Holding out whole Bemis–Murcko scaffolds asks "how well do we predict a *new* chemical series?", which a random split overstates because near-duplicate analogues land in both train and test.

RDKit Morgan fingerprints replace the PaDEL Java tool from the original tutorial, so there is no Java dependency.

## Pipeline

| Part | Notebook | Package entry | Status |
|------|----------|---------------|--------|
| 1. Bioactivity data | `notebooks/01_bioactivity_data.ipynb` | `python -m drug_discovery part1` | Implemented |
| 2. Exploratory data analysis: Lipinski descriptors, pIC50, Mann–Whitney tests | `notebooks/02_exploratory_data_analysis.ipynb` | `python -m drug_discovery part2` | Implemented |
| 3. Morgan fingerprints and Murcko scaffolds | `notebooks/03_descriptor_calculation.ipynb` | `python -m drug_discovery part3` | Implemented |
| 4. Random forest: random vs. scaffold split | `notebooks/04_regression_models.ipynb` | `python -m drug_discovery part4` | Implemented |
| 5. Compare 6 algorithms + mean baseline (scaffold split) | `notebooks/05_compare_algorithms.ipynb` | `python -m drug_discovery part5` | Implemented |

Run everything:

```bash
python -m drug_discovery all          # downloads ChEMBL data, then parts 2-5
```

Outputs: curated CSVs in `data/processed/`, fingerprints in `data/descriptors/`, figures in `reports/figures/`, and all metrics in `reports/metrics.json`.

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

The tests run parts 2–5 end to end on a small offline fixture (real drug-like SMILES with synthetic IC50 values), so CI needs no network access to ChEMBL.
