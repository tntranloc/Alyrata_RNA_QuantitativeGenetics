### TWAS-like Machine Learning: Expression → Trait

This module implements two complementary analyses linking transcriptome (gene counts) to phenotype:
1) **Unsupervised PCR**: PCA on standardized expression → regress trait on PC scores → back-transform coefficients to per-gene effects (“selection gradients”) → rank genes.
2) **Supervised GB**: Gradient Boosting (regressor or classifier) on expression to predict the trait → rank genes by model importance.
3) It also computes **overlap** between top genes from the two approaches and reports enrichment statistics.

#### Data format
- `data/expression.csv`: rows = samples, columns = genes (numeric).  
- `data/phenotypes.csv`: rows = samples, includes a column named `rosetteSize` (default).  
Both files must share the same sample index (row names). Set `config/default.yaml` if your column names differ.

#### Quickstart
```bash
# 1) install
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt

# 2) put your data
# - data/expression.csv : rows=samples, cols=genes (numeric)
# - data/phenotypes.csv : rows=samples, includes column named according to the phenotypes (default)

# 3) run
# Run unsupervised PCR
python scripts/run_pcr.py --config config/default.yaml

# Run supervised GB (classifier by default; switch to regressor in YAML)
python scripts/run_gb.py --config config/default.yaml

# Overlap top genes between methods
python scripts/run_overlap.py --config config/default.yaml

