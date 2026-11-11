"""Runs parts 2-5 on a small offline fixture (real drug-like SMILES, synthetic IC50 values)."""
import numpy as np
import pandas as pd
import pytest

pytest.importorskip("rdkit")

from drug_discovery import parts  # noqa: E402
from drug_discovery.descriptors import murcko_scaffold, smiles_to_fingerprint  # noqa: E402
from drug_discovery.eda import ic50_to_pic50  # noqa: E402
from drug_discovery.models import scaffold_split  # noqa: E402
from drug_discovery.preprocess import label_bioactivity  # noqa: E402

SMILES = [
    "COc1cc2c(cc1OC)C(=O)C(CC1CCN(Cc3ccccc3)CC1)C2",  # donepezil
    "CN(C)C(=O)Oc1cccc(C(C)N(C)C)c1", "CCN(C)C(=O)Oc1cccc(C(C)N(C)C)c1",  # carbamates
    "COc1ccc2c(c1)C13CCCC1N(C)CC3O2", "Nc1c2c(nc3ccccc13)CCCC2",  # tacrine-like
    "Nc1c2c(nc3cc(Cl)ccc13)CCCC2", "Nc1c2c(nc3ccc(Cl)cc13)CCCC2", "Nc1c2c(nc3ccccc13)CCC2",
    "c1ccc2c(c1)ccc1ccccc12", "CCOC(=O)c1ccccc1", "CC(=O)Oc1ccccc1C(=O)O", "CC(C)Cc1ccc(C(C)C(=O)O)cc1",
    "c1ccc(CN2CCC(CC2)c2ccccc2)cc1", "c1ccc(CN2CCC(CCc3ccccc3)CC2)cc1", "O=C(c1ccccc1)N1CCC(Cc2ccccc2)CC1",
    "COc1ccc(CCN2CCCC2)cc1", "CN1CCC(c2ccc(Cl)cc2)CC1", "OC(c1ccccc1)(c1ccccc1)C1CCNCC1",
    "Cc1ccccc1NC(=O)C1CCCN1", "CC(C)NCC(O)c1ccc(O)c(O)c1", "Oc1ccc2c(c1)CCN(C)C2", "COc1ccc2[nH]ccc2c1",
    "CN1CCN(c2ccccc2)CC1", "Clc1ccc(CN2CCOCC2)cc1", "c1ccc(-c2ccccn2)nc1", "O=C1CCCN1Cc1ccccc1",
    "CC1=CC(=O)c2ccccc2C1=O", "O=c1cc(-c2ccccc2)oc2ccccc12", "COc1cc(C=CC(=O)O)ccc1O", "CCCCN1CCCCC1C(=O)Nc1c(C)cccc1C",
]
IC50 = [6, 300, 450, 2000, 120, 40, 60, 250, 50000, 80000, 90000, 70000, 900, 500, 1500, 15000, 30000, 8000,
        60000, 40000, 12000, 25000, 20000, 35000, 45000, 18000, 22000, 9000, 50000, 3000]


@pytest.fixture()
def curated(tmp_path):
    df = pd.DataFrame({"molecule_chembl_id": [f"M{i}" for i in range(len(SMILES))],
                       "canonical_smiles": SMILES, "standard_value": IC50})
    path = tmp_path / "curated.csv"
    label_bioactivity(df).to_csv(path, index=False)
    return path


def test_pic50_conversion():
    out = ic50_to_pic50(pd.DataFrame({"standard_value": [1000.0, 1.0]}))
    assert np.allclose(out["pIC50"], [6.0, 9.0])


def test_fingerprint_and_scaffold():
    fp = smiles_to_fingerprint("CCOC(=O)c1ccccc1")
    assert fp.shape == (2048,) and fp.sum() > 0
    assert smiles_to_fingerprint("not a smiles") is None
    assert murcko_scaffold("CCOC(=O)c1ccccc1") == "c1ccccc1"


def test_scaffold_split_keeps_groups_together():
    scaff = ["a"] * 10 + ["b"] * 4 + ["c"] * 3 + ["d"] * 2 + ["e"]
    tr, te = scaffold_split(scaff, test_size=0.25, seed=1)
    assert set(tr).isdisjoint(te) and len(tr) + len(te) == len(scaff)
    assert {scaff[i] for i in tr}.isdisjoint({scaff[i] for i in te})


def test_parts_2_to_5_run_end_to_end(curated, tmp_path):
    pic = tmp_path / "pic50.csv"
    npz = tmp_path / "fp.npz"
    metrics = tmp_path / "metrics.json"
    df = parts.run_part2(curated, pic, tmp_path / "figs", metrics)
    assert "pIC50" in df and (tmp_path / "figs" / "eda_boxplots.png").exists()
    parts.run_part3(pic, npz)
    res = parts.run_part4(npz, tmp_path / "figs", metrics)
    assert set(res) == {"random split", "scaffold split"}
    rows = parts.run_part5(npz, metrics)
    assert any(r["model"].startswith("Baseline") for r in rows)
    assert metrics.exists()
