"""Project paths and Part 1 defaults.

Change TARGET_QUERY / TARGET_CHEMBL_ID here to point the same pipeline
at a different protein without rewriting notebooks.
"""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]

DATA_RAW = ROOT / "data" / "raw"
DATA_PROCESSED = ROOT / "data" / "processed"
DATA_DESCRIPTORS = ROOT / "data" / "descriptors"
MODELS_DIR = ROOT / "models"
FIGURES_DIR = ROOT / "reports" / "figures"

TARGET_QUERY = "acetylcholinesterase"
TARGET_CHEMBL_ID = "CHEMBL220"
TARGET_NAME = "Human Acetylcholinesterase"
STANDARD_TYPE = "IC50"
STANDARD_UNITS = "nM"

# IC50 thresholds in nM (standard ChEMBL-tutorial cutoffs)
ACTIVE_NM = 1000.0
INACTIVE_NM = 10000.0

RAW_CSV = DATA_RAW / "acetylcholinesterase_01_bioactivity_data_raw.csv"
PREPROCESSED_CSV = DATA_PROCESSED / "acetylcholinesterase_02_bioactivity_data_preprocessed.csv"
CURATED_CSV = DATA_PROCESSED / "acetylcholinesterase_03_bioactivity_data_curated.csv"

BIOACTIVITY_FIELDS = [
    "activity_id",
    "assay_chembl_id",
    "molecule_chembl_id",
    "canonical_smiles",
    "standard_type",
    "standard_relation",
    "standard_value",
    "standard_units",
    "pchembl_value",
    "document_year",
    "target_chembl_id",
    "target_organism",
    "target_pref_name",
]
