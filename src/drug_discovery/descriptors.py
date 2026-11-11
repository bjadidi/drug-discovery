"""Part 3: molecular fingerprints and Bemis-Murcko scaffolds.

Uses RDKit Morgan (ECFP4-style) fingerprints instead of the PaDEL Java tool used in
the original tutorial: no Java dependency, and fast enough to recompute on demand.
"""

from __future__ import annotations

import numpy as np
import pandas as pd


def smiles_to_fingerprint(smiles: str, radius: int = 2, n_bits: int = 2048) -> np.ndarray | None:
    from rdkit import Chem
    from rdkit.Chem import rdFingerprintGenerator

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=radius, fpSize=n_bits)
    return gen.GetFingerprintAsNumPy(mol).astype(np.uint8)


def murcko_scaffold(smiles: str) -> str:
    from rdkit import Chem
    from rdkit.Chem.Scaffolds import MurckoScaffold

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return ""
    return MurckoScaffold.MurckoScaffoldSmiles(mol=mol, includeChirality=False)


def build_ml_dataset(df: pd.DataFrame, smiles_col: str = "canonical_smiles", target: str = "pIC50",
                     radius: int = 2, n_bits: int = 2048) -> tuple[np.ndarray, np.ndarray, list[str]]:
    """Return X (fingerprints), y (target) and scaffolds; rows with invalid SMILES are dropped."""
    X, y, scaffolds = [], [], []
    for smi, val in zip(df[smiles_col], df[target]):
        fp = smiles_to_fingerprint(smi, radius, n_bits)
        if fp is None:
            continue
        X.append(fp)
        y.append(float(val))
        scaffolds.append(murcko_scaffold(smi))
    return np.vstack(X), np.array(y), scaffolds
