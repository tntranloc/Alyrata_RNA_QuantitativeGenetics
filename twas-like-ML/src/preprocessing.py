import numpy as np
import pandas as pd

def relative_fitness(y: pd.Series) -> pd.Series:
    return y / float(y.mean())

def keep_top_variance(expr: pd.DataFrame, top_n: int) -> pd.DataFrame:
    var = expr.var(axis=0, numeric_only=True)
    keep = var.nlargest(top_n).index
    return expr.loc[:, keep]

def median_split(y: pd.Series):
    median = y.median()
    return (y >= median).map({True: "large", False: "small"})
