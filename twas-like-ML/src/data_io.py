import pandas as pd
from typing import Optional, Tuple

def load_data(expr_csv: str, pheno_csv: str, sample_index_col: Optional[int] = 0
             ) -> Tuple[pd.DataFrame, pd.DataFrame]:
    expr = pd.read_csv(expr_csv, index_col=sample_index_col)
    pheno = pd.read_csv(pheno_csv, index_col=sample_index_col)
    # align rows
    common = expr.index.intersection(pheno.index)
    expr = expr.loc[common].sort_index()
    pheno = pheno.loc[common].sort_index()
    return expr, pheno
