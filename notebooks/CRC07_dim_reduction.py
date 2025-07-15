from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc

print("Load data...", flush=True)
data_path = Path("/path/to/celltypist_folder/results_5000_1e2.h5ad")
adata = sc.read(data_path)

print("Perform PCA...", flush=True)
sc.pp.pca(adata)

print("Calculate neighbors...", flush=True)
# calculate neighbors
sc.pp.neighbors(adata)

print("Perform UMAP...", flush=True)
# dimensionality reduction
sc.tl.umap(adata)

print("Perform Leiden clustering...", flush=True)
# leiden clustering
sc.tl.leiden(adata)

# save results
out_path = data_path.parent / "results_5000_1e2_dimred.h5ad"
print(f"Save results to {out_path}.")
adata.write(out_path)