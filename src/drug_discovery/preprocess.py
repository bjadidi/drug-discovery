"""Cleaning and labeling for ChEMBL IC50 tables."""

from __future__ import annotations

import pandas as pd

from .config import ACTIVE_NM, INACTIVE_NM, STANDARD_UNITS

KEEP_COLS = ["molecule_chembl_id", "canonical_smiles", "standard_value"]


def drop_missing(df: pd.DataFrame) -> pd.DataFrame:
    """Drop rows missing potency or structure."""
    return df[df["standard_value"].notna() & df["canonical_smiles"].notna()].copy()


def filter_standard_units(df: pd.DataFrame, units: str = STANDARD_UNITS) -> pd.DataFrame:
    """Keep a single potency unit so IC50 thresholds stay comparable."""
    if "standard_units" not in df.columns:
        return df.copy()
    return df[df["standard_units"] == units].copy()


def drop_duplicate_smiles(df: pd.DataFrame) -> pd.DataFrame:
    return df.drop_duplicates(subset=["canonical_smiles"]).copy()


def select_columns(df: pd.DataFrame, columns: list[str] | None = None) -> pd.DataFrame:
    columns = columns or KEEP_COLS
    missing = [c for c in columns if c not in df.columns]
    if missing:
        raise KeyError(f"Missing expected columns: {missing}")
    return df[columns].copy()


def label_bioactivity(
    df: pd.DataFrame,
    active_nm: float = ACTIVE_NM,
    inactive_nm: float = INACTIVE_NM,
) -> pd.DataFrame:
    """Label compounds active / inactive / intermediate from IC50 (nM)."""
    out = df.copy()
    values = out["standard_value"].astype(float)
    labels = pd.Series("intermediate", index=out.index, dtype="object")
    labels = labels.mask(values <= active_nm, "active")
    labels = labels.mask(values >= inactive_nm, "inactive")
    out["class"] = labels
    return out
