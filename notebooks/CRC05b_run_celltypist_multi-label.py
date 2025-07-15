#!/usr/bin/env python3
import time
from pathlib import Path

import celltypist
import scanpy as sc

n_cells = 5000

celltypist_path = Path("/path/to/celltypist_folder/")
ad_path = celltypist_path / "adata_pp.h5ad"
ref_path = celltypist_path / f"reference_{n_cells}.h5ad"

# load data
adata = sc.read(ad_path)
reference = sc.read(ref_path)

# Add `check_expression = False` to bypass expression check with only a subset of genes.
t_start = time.time()
model = celltypist.train(reference, 'cell_type_middle', check_expression = False, n_jobs = 10,
                         #max_iter = 100
                        )
t_end = time.time()
print(f"Time elapsed: {(t_end - t_start)/60} minutes")

# # Save the model.
model_path = celltypist_path / f'crc_atlas_model_max{n_cells}.pkl'
model.write(model_path)

# CellTypist prediction
t_start = time.time()
predictions = celltypist.annotate(adata, model = str(model_path),
                                  majority_voting = True, mode = 'prob match', p_thres = 0.5)
t_end = time.time()
print(f"Time elapsed: {t_end - t_start} seconds")

# convert to adata
result = predictions.to_adata(insert_prob = True)

results_path = celltypist_path / f'results_{n_cells}_multi-label.h5ad'
print(f"Save results to {str(results_path)}")
result.write(results_path)