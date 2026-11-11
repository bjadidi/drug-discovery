"""Parts 2-5 as runnable steps. Each reads the previous step's CSV and writes its own outputs."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from .config import CURATED_CSV, DATA_DESCRIPTORS, DATA_PROCESSED, FIGURES_DIR, ROOT

PIC50_CSV = DATA_PROCESSED / "acetylcholinesterase_04_bioactivity_data_pic50.csv"
FP_NPZ = DATA_DESCRIPTORS / "acetylcholinesterase_05_morgan_fp.npz"
METRICS_JSON = ROOT / "reports" / "metrics.json"


def _update_metrics(key: str, value, path: Path = METRICS_JSON) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    data = json.loads(path.read_text()) if path.exists() else {}
    data[key] = value
    path.write_text(json.dumps(data, indent=2))


def run_part2(curated_csv: Path = CURATED_CSV, out_csv: Path = PIC50_CSV, fig_dir: Path = FIGURES_DIR,
              metrics_path: Path = METRICS_JSON) -> pd.DataFrame:
    from .eda import LIPINSKI_COLS, drop_intermediate, ic50_to_pic50, lipinski, mannwhitney, save_eda_figures

    df = ic50_to_pic50(lipinski(pd.read_csv(curated_csv)))
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_csv, index=False)
    two_class = drop_intermediate(df)
    tests = [mannwhitney(two_class, d) for d in ["pIC50", *LIPINSKI_COLS]]
    save_eda_figures(two_class, fig_dir)
    _update_metrics("part2_mannwhitney", tests, metrics_path)
    _update_metrics("dataset", {"compounds": int(len(df)),
                                "class_counts": df["class"].value_counts().to_dict()}, metrics_path)
    for t in tests:
        print(f"{t['descriptor']:14s} p={t['p']:.2e}  different={t['different_distribution']}")
    return df


def run_part3(pic50_csv: Path = PIC50_CSV, out_npz: Path = FP_NPZ) -> Path:
    from .descriptors import build_ml_dataset

    X, y, scaffolds = build_ml_dataset(pd.read_csv(pic50_csv))
    out_npz.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(out_npz, X=X, y=y, scaffolds=np.array(scaffolds, dtype=object))
    print(f"fingerprints: {X.shape}, unique scaffolds: {len(set(scaffolds))} -> {out_npz}")
    return out_npz


def _load(npz: Path):
    d = np.load(npz, allow_pickle=True)
    return d["X"], d["y"], list(d["scaffolds"])


def run_part4(npz: Path = FP_NPZ, fig_dir: Path = FIGURES_DIR, metrics_path: Path = METRICS_JSON) -> dict:
    from sklearn.model_selection import train_test_split

    from .models import evaluate_regressor, random_vs_scaffold, train_random_forest

    X, y, scaffolds = _load(npz)
    res = random_vs_scaffold(X, y, scaffolds)
    _update_metrics("part4_random_forest", res, metrics_path)
    for k, v in res.items():
        print(f"{k:15s} R2={v['r2']:.3f} RMSE={v['rmse']:.3f} MAE={v['mae']:.3f}")
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    tr, te = train_test_split(np.arange(len(y)), test_size=0.2, random_state=0)
    model = train_random_forest(X[tr], y[tr])
    pred = model.predict(X[te])
    fig, ax = plt.subplots(figsize=(4.6, 4.4))
    ax.scatter(y[te], pred, s=12, alpha=0.6)
    lo, hi = min(y.min(), pred.min()), max(y.max(), pred.max())
    ax.plot([lo, hi], [lo, hi], "k--", lw=1)
    ax.set_xlabel("Experimental pIC50")
    ax.set_ylabel("Predicted pIC50")
    ax.set_title(f"Random forest (random split), R² = {evaluate_regressor(model, X[te], y[te])['r2']:.2f}")
    fig.tight_layout()
    fig_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(fig_dir / "rf_parity.png", dpi=150)
    plt.close(fig)
    return res


def run_part5(npz: Path = FP_NPZ, metrics_path: Path = METRICS_JSON) -> list[dict]:
    from .models import compare_algorithms

    X, y, scaffolds = _load(npz)
    rows = compare_algorithms(X, y, scaffolds)
    _update_metrics("part5_scaffold_split_comparison", rows, metrics_path)
    print("| Model | R² | RMSE | MAE |\n|---|---|---|---|")
    for r in rows:
        print(f"| {r['model']} | {r['r2']:.3f} | {r['rmse']:.3f} | {r['mae']:.3f} |")
    return rows
