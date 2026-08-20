"""Part 1: download and curate ChEMBL bioactivity data."""

from __future__ import annotations

import pandas as pd

from .chembl import fetch_bioactivity, search_targets
from .config import (
    BIOACTIVITY_FIELDS,
    CURATED_CSV,
    DATA_PROCESSED,
    DATA_RAW,
    PREPROCESSED_CSV,
    RAW_CSV,
    STANDARD_TYPE,
    TARGET_CHEMBL_ID,
    TARGET_QUERY,
)
from .preprocess import (
    drop_duplicate_smiles,
    drop_missing,
    filter_standard_units,
    label_bioactivity,
    select_columns,
)


def run_part1(
    target_chembl_id: str = TARGET_CHEMBL_ID,
    skip_fetch: bool = False,
    search: bool = True,
) -> pd.DataFrame:
    """Fetch IC50 data, clean it, label activity class, and write CSVs."""
    DATA_RAW.mkdir(parents=True, exist_ok=True)
    DATA_PROCESSED.mkdir(parents=True, exist_ok=True)

    if skip_fetch and RAW_CSV.exists():
        raw = pd.read_csv(RAW_CSV)
        print(f"Loaded cached raw data: {RAW_CSV} ({len(raw):,} rows)")
    else:
        if search:
            print(f"Searching ChEMBL for {TARGET_QUERY!r}...")
            targets = search_targets(TARGET_QUERY)
            print(f"Found {len(targets)} target hits. Using {target_chembl_id}.")
            if "target_chembl_id" in targets.columns:
                match = targets[targets["target_chembl_id"] == target_chembl_id]
                if not match.empty:
                    pref = match.iloc[0].get("pref_name", "")
                    org = match.iloc[0].get("organism", "")
                    print(f"  {pref} ({org})")

        print(f"Fetching {STANDARD_TYPE} activities for {target_chembl_id}...")
        raw = fetch_bioactivity(
            target_chembl_id,
            standard_type=STANDARD_TYPE,
            fields=BIOACTIVITY_FIELDS,
        )
        raw.to_csv(RAW_CSV, index=False)
        print(f"Saved raw data -> {RAW_CSV} ({len(raw):,} rows)")

    cleaned = drop_missing(raw)
    cleaned = filter_standard_units(cleaned)
    cleaned = drop_duplicate_smiles(cleaned)

    preprocessed = select_columns(cleaned)
    preprocessed.to_csv(PREPROCESSED_CSV, index=False)
    print(f"Saved preprocessed -> {PREPROCESSED_CSV} ({len(preprocessed):,} rows)")

    curated = label_bioactivity(preprocessed)
    curated.to_csv(CURATED_CSV, index=False)
    print(f"Saved curated -> {CURATED_CSV} ({len(curated):,} rows)")
    print("\nBioactivity class counts:")
    print(curated["class"].value_counts().to_string())
    return curated
