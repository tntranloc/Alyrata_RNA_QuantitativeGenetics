import yaml, argparse, pandas as pd, numpy as np
from pathlib import Path
from src.utils import set_seed, ensure_dir
from src.data_io import load_data
from src.preprocessing import relative_fitness, keep_top_variance, median_split
from src.supervised_gb import run_gb, GBConfig

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
    y_cont = pheno[trait].astype(float)
    if cfg["data"].get("use_relative_fitness", True):
        y_cont = relative_fitness(y_cont)

    # choose target per mode
    mode = cfg["supervised_gb"]["mode"]
    y = median_split(y_cont) if mode == "classifier" else y_cont

    if cfg["preprocessing"]["filter_low_variance"]["enable"]:
        expr = keep_top_variance(expr, cfg["preprocessing"]["filter_low_variance"]["top_n"])

    gb_cfg = GBConfig(mode=mode,
                      test_size=cfg["supervised_gb"]["test_size"],
                      params=cfg["supervised_gb"]["params"],
                      select_top_genes=cfg["supervised_gb"]["select_top_genes"])

    top_genes, importances, metrics = run_gb(expr, y.values, gb_cfg, standardize=cfg["preprocessing"]["standardize_features"])
    pd.DataFrame({"Gene": expr.columns, "Importance": importances}).to_csv(outdir/"gb_importances.csv", index=False)
    pd.DataFrame({"Top Genes": top_genes}).to_csv(outdir/"top_genes_gb.csv", index=False)
    pd.Series(metrics).to_csv(outdir/"gb_metrics.csv")
    print("GB metrics:", metrics)
    print("saved:", outdir/"top_genes_gb.csv")

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="config/default.yaml")
    args = ap.parse_args()
    main(args.config)
