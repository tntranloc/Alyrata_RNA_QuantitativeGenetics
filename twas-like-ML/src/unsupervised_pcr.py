from dataclasses import dataclass
import numpy as np, pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LinearRegression
from sklearn.pipeline import Pipeline
from sklearn.model_selection import RepeatedKFold, GridSearchCV, cross_val_score

@dataclass
class PCRConfig:
    leakage_safe: bool = False
    pca_variance_cutoff: float = 0.95
    cv_splits: int = 5
    cv_repeats: int = 100
    select_top_genes: int = 300
    grid_k: tuple = (10, 20, 50, 100, 200, 300)
    randomized_pca: bool = False

def run_pcr(expr: pd.DataFrame, y: np.ndarray, cfg: PCRConfig, standardize: bool=True):
    X = expr.values
    if not cfg.leakage_safe:
        # Your original approach: scale + PCA on full data, then CV on PC scores
        scaler = StandardScaler() if standardize else None
        Xs = scaler.fit_transform(X) if scaler else X

        pca = PCA(n_components=cfg.pca_variance_cutoff, svd_solver="full")
        if cfg.randomized_pca:
            pca = PCA(n_components=cfg.pca_variance_cutoff, svd_solver="randomized", iterated_power=2, random_state=1)

        pc_scores = pca.fit_transform(Xs)
        # optional: drop last PC (kept for fidelity with your script)
        if pc_scores.shape[1] > 1:
            pc_scores = pc_scores[:, :-1]

        # rank PCs with simple 1 - R^2 score (your “pseudo p-value” idea)
        def pseudo_pvals(scores, response):
            vals = []
            for i in range(scores.shape[1]):
                m = LinearRegression().fit(scores[:, [i]], response)
                preds = m.predict(scores[:, [i]])
                rss = ((response - preds)**2).sum()
                tss = ((response - response.mean())**2).sum()
                r2 = 1 - rss/tss
                vals.append(1 - r2)
            return np.array(vals)

        pvals = pseudo_pvals(pc_scores, y)
        order = np.argsort(pvals)

        # CV to choose how many ordered PCs
        rkf = RepeatedKFold(n_splits=cfg.cv_splits, n_repeats=cfg.cv_repeats, random_state=1)
        results = []
        for k in range(1, len(order)+1):
            Xk = pc_scores[:, order[:k]]
            r2 = cross_val_score(LinearRegression(), Xk, y, scoring="r2", cv=rkf, n_jobs=-1).mean()
            rmse = (-cross_val_score(LinearRegression(), Xk, y, scoring="neg_root_mean_squared_error", cv=rkf, n_jobs=-1)).mean()
            mae = (-cross_val_score(LinearRegression(), Xk, y, scoring="neg_mean_absolute_error", cv=rkf, n_jobs=-1)).mean()
            results.append((k, r2, rmse, mae))
        k_opt, r2_opt, rmse_opt, mae_opt = max(results, key=lambda t: t[1])
        X_opt = pc_scores[:, order[:k_opt]]
        final = LinearRegression().fit(X_opt, y)

        # back-transform to gene space
        beta_pc = final.coef_
        loadings = pca.components_[order[:k_opt], :]
        beta_genes_std = beta_pc @ loadings
        if scaler:
            beta_genes = beta_genes_std / scaler.scale_
        else:
            beta_genes = beta_genes_std

        abs_beta = np.abs(beta_genes)
        top_idx = np.argsort(abs_beta)[-cfg.select_top_genes:]
        top_genes = expr.columns[top_idx]

        metrics = dict(RSquared=r2_opt, RMSE=rmse_opt, MAE=mae_opt)
        return top_genes, beta_genes, metrics

    else:
        # Leakage-safe pipeline: scale + PCA + OLS inside CV with grid over k
        steps = []
        if standardize:
            steps.append(("scaler", StandardScaler()))
        solver = "randomized" if cfg.randomized_pca else "full"
        steps.append(("pca", PCA(svd_solver=solver, iterated_power=2 if solver=="randomized" else None, random_state=1)))
        steps.append(("ols", LinearRegression()))
        pipe = Pipeline(steps)

        rkf = RepeatedKFold(n_splits=cfg.cv_splits, n_repeats=cfg.cv_repeats, random_state=1)
        search = GridSearchCV(pipe, {"pca__n_components": list(cfg.grid_k)}, scoring="r2", cv=rkf, n_jobs=-1, refit=True)
        search.fit(X, y)
        best = search.best_estimator_

        pca_f = best.named_steps["pca"]
        beta_pc = best.named_steps["ols"].coef_
        if standardize:
            scaler_f = best.named_steps["scaler"]
            beta_genes_std = pca_f.components_.T @ beta_pc
            beta_genes = beta_genes_std / scaler_f.scale_
        else:
            beta_genes = pca_f.components_.T @ beta_pc

        abs_beta = np.abs(beta_genes)
        top_idx = np.argsort(abs_beta)[-cfg.select_top_genes:]
        top_genes = expr.columns[top_idx]

        metrics = dict(RSquared=float(search.best_score_))
        return top_genes, beta_genes, metrics
