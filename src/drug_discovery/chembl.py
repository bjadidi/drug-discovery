"""ChEMBL helpers: target search and bioactivity download."""

from __future__ import annotations

import pandas as pd

from .config import BIOACTIVITY_FIELDS, STANDARD_TYPE


def _client():
    from chembl_webresource_client.new_client import new_client

    return new_client


def search_targets(query: str) -> pd.DataFrame:
    """Return ChEMBL target hits for a free-text query."""
    target = _client().target
    hits = target.search(query)
    return pd.DataFrame.from_dict(hits)


def fetch_bioactivity(
    target_chembl_id: str,
    standard_type: str = STANDARD_TYPE,
    fields: list[str] | None = None,
    progress_every: int = 1000,
) -> pd.DataFrame:
    """Download bioactivity records for a target, with optional column subset."""
    activity = _client().activity
    query = activity.filter(target_chembl_id=target_chembl_id).filter(
        standard_type=standard_type
    )
    if fields:
        query = query.only(fields)

    records: list[dict] = []
    for i, rec in enumerate(query, start=1):
        records.append(rec)
        if progress_every and i % progress_every == 0:
            print(f"  fetched {i:,} records...")

    print(f"  fetched {len(records):,} records total")
    if not records:
        return pd.DataFrame(columns=fields or BIOACTIVITY_FIELDS)
    return pd.DataFrame.from_records(records)
