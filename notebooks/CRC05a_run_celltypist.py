#!/usr/bin/env python3
import time
from pathlib import Path

import celltypist
import scanpy as sc

n_cells = 5000
norm = "1e2"

celltypist_path = Path("/path/to/celltypist_folder/")
out_path = celltypist_path / "results"
out_path.mkdir(exist_ok=True) # create output directory

ad_path = celltypist_path / f"adata_rounded_pp{norm}.h5ad"
ref_path = celltypist_path / f"reference_{n_cells}_pp{norm}.h5ad"

# load data
adata = sc.read(ad_path)
reference = sc.read(ref_path)

# perform subsampling
print("Subsampling...")
sc.pp.sample(adata, n=1000)

# Add `check_expression = False` to bypass expression check with only a subset of genes.
t_start = time.time()
model = celltypist.train(reference, 'cell_type_middle', check_expression = False, n_jobs = -1, max_iter = None)
t_end = time.time()
print(f"Time elapsed: {(t_end - t_start)/60} minutes", flush=True)

# # Save the model.
model_path = celltypist_path / f'crc_atlas_model_max{n_cells}{norm}.pkl'
model.write(model_path)

# CellTypist prediction without over-clustering and majority-voting.
t_start = time.time()
predictions = celltypist.annotate(adata, model = str(model_path), majority_voting = True)
t_end = time.time()
print(f"Time elapsed: {(t_end - t_start)/60} minutes", flush=True)

# convert to adata
result = predictions.to_adata()

results_path = out_path / f'celltypist_rounded_{n_cells}_{norm}.h5ad'
print(f"Save results to {str(results_path)}", flush=True)
result.write(results_path)

print("Finished.")