"""Parts 4-5: QSAR regression on pIC50, with random vs. scaffold splits and an algorithm comparison.

A random split leaks near-duplicate chemistry between train and test and flatters the model.
A scaffold split (whole Bemis-Murcko scaffolds held out) is closer to the real question:
"how well do we predict a *new* chemical series?" Both are reported.
"""

from __future__ import annotations

from collections import defaultdict

import numpy as np
from sklearn.ensemble import ExtraTreesRegressor, HistGradientBoostingRegressor, RandomForestRegressor
from sklearn.linear_model import Ridge
from sklearn.metrics import mean_absolute_error, r2_score, root_mean_squared_error
from sklearn.model_selection import train_test_split
from sklearn.neighbors import KNeighborsRegressor
from sklearn.svm import SVR


def scaffold_split(scaffolds: list[str], test_size: float = 0.2, seed: int = 0) -> tuple[np.ndarray, np.ndarray]:
    """Assign whole scaffold groups to train or test (largest groups to train first, then shuffled rest)."""
    groups: dict[str, list[int]] = defaultdict(list)
    for i, s in enumerate(scaffolds):
        groups[s].append(i)
    rng = np.random.default_rng(seed)
    big = sorted([g for g in groups.values() if len(g) > len(scaffolds) * test_size / 2], key=len, reverse=True)
    small = [g for g in groups.values() if len(g) <= len(scaffolds) * test_size / 2]
    rng.shuffle(small)
    n_test = int(round(test_size * len(scaffolds)))
    train, test = [], []
    for g in big + small:
        (test if len(test) + len(g) <= n_test and g not in big else train).extend(g)
    return np.array(sorted(train)), np.array(sorted(test))


def evaluate_regressor(model, X_test, y_test) -> dict:
    pred = model.predict(X_test)
    return {"r2": float(r2_score(y_test, pred)),
            "rmse": float(root_mean_squared_error(y_test, pred)),
            "mae": float(mean_absolute_error(y_test, pred))}


def train_random_forest(X_train, y_train, seed: int = 0) -> RandomForestRegressor:
    return RandomForestRegressor(n_estimators=500, n_jobs=-1, random_state=seed).fit(X_train, y_train)


def random_vs_scaffold(X, y, scaffolds, seed: int = 0) -> dict:
    out = {}
    tr, te = train_test_split(np.arange(len(y)), test_size=0.2, random_state=seed)
    out["random split"] = evaluate_regressor(train_random_forest(X[tr], y[tr], seed), X[te], y[te])
    tr, te = scaffold_split(scaffolds, 0.2, seed)
    out["scaffold split"] = evaluate_regressor(train_random_forest(X[tr], y[tr], seed), X[te], y[te])
    return out


def candidate_models(seed: int = 0) -> dict:
    return {
        "Random forest": RandomForestRegressor(n_estimators=300, n_jobs=-1, random_state=seed),
        "Extra trees": ExtraTreesRegressor(n_estimators=300, n_jobs=-1, random_state=seed),
        "Gradient boosting": HistGradientBoostingRegressor(random_state=seed),
        "SVR (RBF)": SVR(C=5.0),
        "k-NN (k=5)": KNeighborsRegressor(n_neighbors=5, metric="jaccard"),
        "Ridge": Ridge(alpha=10.0),
        "Baseline (predict mean)": None,
    }


def compare_algorithms(X, y, scaffolds, seed: int = 0) -> list[dict]:
    """Fit each model on the same scaffold split; returns rows sorted by test R^2."""
    tr, te = scaffold_split(scaffolds, 0.2, seed)
    rows = []
    for name, model in candidate_models(seed).items():
        if model is None:
            pred = np.full(len(te), y[tr].mean())
            rows.append({"model": name, "r2": float(r2_score(y[te], pred)),
                         "rmse": float(root_mean_squared_error(y[te], pred)),
                         "mae": float(mean_absolute_error(y[te], pred))})
            continue
        Xtr = X[tr].astype(bool) if name.startswith("k-NN") else X[tr]
        Xte = X[te].astype(bool) if name.startswith("k-NN") else X[te]
        model.fit(Xtr, y[tr])
        rows.append({"model": name, **evaluate_regressor(model, Xte, y[te])})
    return sorted(rows, key=lambda r: r["r2"], reverse=True)
