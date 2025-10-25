# TWAS-like ML for Expression→Trait

This repo provides a small, reproducible pipeline for linking transcriptome (gene counts) to a trait using:
1) **Unsupervised PCR** (PCA → regression → back-transform to gene coefficients),
2) **Supervised Gradient Boosting** (regressor or classifier),
3) **Overlap** of top genes across methods + enrichment-ready outputs.

## Quickstart
```bash
# 1) install
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt

# 2) put your data
# - data/expression.csv : rows=samples, cols=genes (numeric)
# - data/phenotypes.csv : rows=samples, includes column named 'rosetteSize' (default)

# 3) run
python scripts/run_pcr.py --config config/default.yaml
python scripts/run_gb.py  --config config/default.yaml
python scripts/run_overlap.py --config config/default.yaml
