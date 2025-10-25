from dataclasses import dataclass
import numpy as np, pandas as pd
from sklearn.ensemble import GradientBoostingClassifier, GradientBoostingRegressor
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split
from sklearn.pipeline import Pipeline
from .metrics import regression_summary, classification_summary
from .preprocessing import median_split

@dataclass
class GBConfig:
    mode: str = "classifier"   # "classifier" or "regressor"
    test_size: float = 0.4
    params: dict = None
    select_top_genes: int = 300

def run_gb(expr: pd.DataFrame, y: np.ndarray, cfg: GBConfig, standardize: bool=True):
    X = expr.values
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=cfg.test_size, random_state=42, stratify=(y if cfg.mode=="classifier" else None))
    steps = []
    if standardize: steps.append(("scaler", StandardScaler()))
    if cfg.mode == "classifier":
        model = GradientBoostingClassifier(random_state=42, **(cfg.params or {}))
    else:
        model = GradientBoostingRegressor(random_state=42, **(cfg.params or {}))
    steps.append(("gb", model))
    pipe = Pipeline(steps)
    pipe.fit(X_train, y_train)
    y_pred = pipe.predict(X_test)

    if cfg.mode == "classifier":
        metrics = classification_summary(y_test, y_pred)
        importances = pipe.named_steps["gb"].feature_importances_
    else:
        metrics = regression_summary(y_test, y_pred)
        importances = pipe.named_steps["gb"].feature_importances_

    order = np.argsort(importances)[::-1]
    order = order[importances[order] > 0]  # keep non-zero importances
    top = order[:cfg.select_top_genes]
    top_genes = expr.columns[top]
    return top_genes, importances, metrics
