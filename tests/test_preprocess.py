import pandas as pd

from drug_discovery.preprocess import drop_missing, label_bioactivity


def test_label_bioactivity_thresholds():
    df = pd.DataFrame({"standard_value": [100, 1000, 5000, 10000, 20000]})
    out = label_bioactivity(df)
    assert list(out["class"]) == [
        "active",
        "active",
        "intermediate",
        "inactive",
        "inactive",
    ]


def test_drop_missing_requires_value_and_smiles():
    df = pd.DataFrame(
        {
            "standard_value": [1.0, None, 3.0],
            "canonical_smiles": ["CCO", "CC", None],
        }
    )
    cleaned = drop_missing(df)
    assert len(cleaned) == 1
    assert cleaned.iloc[0]["canonical_smiles"] == "CCO"
