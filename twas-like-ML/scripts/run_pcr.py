import yaml, argparse, pandas as pd, numpy as np
from pathlib import Path
from src.utils import set_seed, ensure_dir
from src.data_io import load_data
from src.preprocessing import relative_fitness, keep_top_variance
from src.unsupervised_pcr import run_pcr, PCRConfig

def main(cfg_path):
    cfg = yaml.safe_load(open(cfg_path))
    set_seed(cfg.get("random_seed", 42))
    outdir = Path(cfg["outputs_dir"]); ensure_dir(outdir)

    expr, pheno = load_data(
        cfg["data"]["expression_csv"],
        cfg["data"]["phenotypes_csv"],
        cfg["data"].get("sample_index_col", 0)
    )

    trait = cfg["data"]["trait"]
    y = pheno[trait].astype(float)
    if cfg["data"].get("use_relative_fitness", True):
        y = relative_fitness(y)
    y = y.values

    if cfg["preprocessing"]["filter_low_variance"]["enable"]:
        expr = keep_top_variance(expr, cfg["preprocessing"]["filter_low_variance"]["top_n"])

    pcr_cfg = PCRConfig(
        leakage_safe=cfg["unsupervised_pcr"]["leakage_safe"],
        pca_variance_cutoff=cfg["unsupervised_pcr"]["pca_variance_cutoff"],
        cv_splits=cfg["unsupervised_pcr"]["cv"]["n_splits"],
        cv_repeats=cfg["unsupervised_pcr"]["cv"]["n_repeats"],
        select_top_genes=cfg["unsupervised_pcr"]["select_top_genes"],
        grid_k=tuple(cfg["unsupervised_pcr"]["grid_k"]),
        randomized_pca=cfg["unsupervised_pcr"]["randomized_pca"],
    )

    top_genes, beta_genes, metrics = run_pcr(expr, y, pcr_cfg, standardize=cfg["preprocessing"]["standardize_features"])
    pd.DataFrame({"Top Genes": top_genes}).to_csv(outdir/"top_genes_pcr.csv", index=False)
    pd.Series(beta_genes, index=expr.columns, name="beta_gene").to_csv(outdir/"pcr_beta_gene.csv")
    pd.Series(metrics).to_csv(outdir/"pcr_metrics.csv")
    print("PCR metrics:", metrics)
    print("saved:", outdir/"top_genes_pcr.csv")

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="config/default.yaml")
    args = ap.parse_args()
    main(args.config)
