import os, random, numpy as np

def set_seed(seed: int = 42):
    import numpy as np, random
    random.seed(seed); np.random.seed(seed)

def ensure_dir(p):
    os.makedirs(p, exist_ok=True)
