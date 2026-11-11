"""Part 2: Lipinski descriptors, pIC50, and statistical comparison of actives vs. inactives."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

LIPINSKI_COLS = ["MW", "LogP", "NumHDonors", "NumHAcceptors"]


def lipinski(df: pd.DataFrame, smiles_col: str = "canonical_smiles") -> pd.DataFrame:
    """Append Lipinski rule-of-five descriptors computed with RDKit; drops unparsable SMILES."""
    from rdkit import Chem
    from rdkit.Chem import Crippen, Descriptors, Lipinski

    rows, keep = [], []
    for idx, smi in df[smiles_col].items():
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            continue
        keep.append(idx)
        rows.append([Descriptors.MolWt(mol), Crippen.MolLogP(mol),
                     Lipinski.NumHDonors(mol), Lipinski.NumHAcceptors(mol)])
    desc = pd.DataFrame(rows, columns=LIPINSKI_COLS, index=keep)
    return df.loc[keep].join(desc)


def ic50_to_pic50(df: pd.DataFrame, value_col: str = "standard_value", cap_nm: float = 1e8) -> pd.DataFrame:
    """pIC50 = -log10(IC50 in molar). Values above cap_nm are capped to avoid negative pIC50."""
    out = df.copy()
    nm = out[value_col].astype(float).clip(upper=cap_nm, lower=1e-3)
    out["pIC50"] = -np.log10(nm * 1e-9)
    return out.drop(columns=[value_col])


def drop_intermediate(df: pd.DataFrame) -> pd.DataFrame:
    return df[df["class"] != "intermediate"].copy()


def mannwhitney(df: pd.DataFrame, descriptor: str, alpha: float = 0.05) -> dict:
    """Two-sided Mann-Whitney U test of a descriptor between active and inactive compounds."""
    from scipy.stats import mannwhitneyu

    active = df.loc[df["class"] == "active", descriptor]
    inactive = df.loc[df["class"] == "inactive", descriptor]
    stat, p = mannwhitneyu(active, inactive, alternative="two-sided")
    return {"descriptor": descriptor, "U": float(stat), "p": float(p),
            "different_distribution": bool(p < alpha),
            "median_active": float(active.median()), "median_inactive": float(inactive.median())}


def save_eda_figures(df: pd.DataFrame, out_dir: Path) -> list[Path]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir.mkdir(parents=True, exist_ok=True)
    paths = []
    fig, axes = plt.subplots(1, 5, figsize=(17, 3.4))
    for ax, col in zip(axes, ["pIC50", *LIPINSKI_COLS]):
        data = [df.loc[df["class"] == c, col] for c in ("active", "inactive")]
        ax.boxplot(data, tick_labels=["active", "inactive"])
        ax.set_title(col)
        ax.grid(alpha=0.3)
    fig.tight_layout()
    p = out_dir / "eda_boxplots.png"
    fig.savefig(p, dpi=150)
    paths.append(p)
    fig, ax = plt.subplots(figsize=(5.5, 4.2))
    for c, color in (("active", "#2b6cb0"), ("inactive", "#c05621")):
        sub = df[df["class"] == c]
        ax.scatter(sub["MW"], sub["LogP"], s=12, alpha=0.6, label=c, color=color)
    ax.set_xlabel("Molecular weight")
    ax.set_ylabel("LogP")
    ax.set_title("Chemical space: MW vs LogP")
    ax.legend()
    fig.tight_layout()
    p = out_dir / "eda_mw_logp.png"
    fig.savefig(p, dpi=150)
    paths.append(p)
    plt.close("all")
    return paths
