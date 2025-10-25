import argparse, yaml, pandas as pd
from pathlib import Path
from src.utils import ensure_dir
from src.overlap import compute_overlap

def main(cfg_path):
    cfg = yaml.safe_load(open(cfg_path))
    outdir = Path(cfg["outputs_dir"]); ensure_dir(outdir)

    pcr = pd.read_csv(outdir/"top_genes_pcr.csv")["Top Genes"].tolist()
    gb  = pd.read_csv(outdir/"top_genes_gb.csv")["Top Genes"].tolist() if "Top Genes" in pd.read_csv(outdir/"top_genes_gb.csv").columns \
          else pd.read_csv(outdir/"top_genes_gb.csv")["Gene"].tolist()  # fallback

    # universe = all genes in expression (best if you pass it explicitly; here we infer from saved beta file)
    try:
        universe_size = pd.read_csv(outdir/"pcr_beta_gene.csv", index_col=0).shape[0]
    except Exception:
        # if not available, fall back to union size
        universe_size = len(set(pcr) | set(gb))

    inter, stats = compute_overlap(pcr, gb, universe_size)
    pd.DataFrame({"Overlap Genes": inter}).to_csv(cfg["overlap"]["out_csv"], index=False)
    pd.Series(stats).to_csv(outdir/"overlap_stats.csv")
    print("overlap:", stats, "saved:", cfg["overlap"]["out_csv"])

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="config/default.yaml")
    args = ap.parse_args()
    main(args.config)
