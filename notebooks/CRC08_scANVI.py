#!/usr/bin/env python3

# check GPU availability
import torch

print(f"GPU availability: {torch.cuda.is_available()}")

print("Import packages", flush=True)
from pathlib import Path

import scanpy as sc
import scvi

torch.multiprocessing.set_sharing_strategy('file_system')

# parameters
adata_name = "celltypist_rounded_5000_1e2"
scvi_batch_key = "uid"
scanvi_labels_key = "majority_voting_simple"
gene_likelihood = "nb"
dispersion = "gene-batch"
max_epochs_scvi = 100
max_epochs_scanvi = 50
out_name = f"{gene_likelihood}_{dispersion}_{max_epochs_scvi}_{max_epochs_scanvi}"
subsample = False
num_workers = 5

data_path = Path(f"/dss/dssfs02/lwp-dss-0001/pn57fo/pn57fo-dss-0000/projects/2301-CRC/celltypist/results/{adata_name}.h5ad")
output_path = Path(f"/dss/dssfs02/lwp-dss-0001/pn57fo/pn57fo-dss-0000/projects/2301-CRC/scANVI_batch_correction/results_gpu__{adata_name}__{scanvi_labels_key}__{out_name}")
output_path.mkdir(exist_ok=True) # create output directory

# setup model paths
scvi_path = output_path / f"scVI_model"
scanvi_path = output_path / f"scANVI_model"

# output path for result anndata
scvi_adata_out_path = output_path / f"scvi.h5ad"
adata_out_path = output_path / f"scanvi.h5ad"

print("Load data.", flush=True)
adata = sc.read(data_path)

print("Calculate n_samples_per_label")
min_counts = int(adata.obs[scanvi_labels_key].value_counts().min())
if min_counts < 100:
    n_samples_per_label = 100
else:
    n_samples_per_label = min_counts

print(f"scVI batch key: {scvi_batch_key}")
print(f"scANVI labels key: {scanvi_labels_key}")
print(f"scANVI n_samples_per_label: {n_samples_per_label}")
print(f"Parameters: {out_name}")
print(f"Data path: {data_path}")
print(f"Output path: {output_path}")
print(f"scVI model path: {scvi_path}")
print(f"scANVI model path: {scanvi_path}")
print(f"scVI anndata path: {scvi_adata_out_path}")
print(f"scANVI anndata path: {adata_out_path}")
print(f"Subsampling: {subsample}")
print(f"Number of workers: {num_workers}")

if subsample:
    # perform subsampling
    print("Subsampling...")
    sc.pp.sample(adata, n=1000)

print("Dimensionality reduction and clustering before scVI.", flush=True)
sc.pp.pca(adata)
sc.pp.neighbors(adata)
#sc.tl.leiden(adata, key_added="leiden")
sc.tl.umap(adata, min_dist=0.3)

# print("Plot UMAP before scVI.", flush=True)
# sc.pl.embedding(
#     adata,
#     basis="X_umap",
#     color=["uid", "majority_voting", "predicted_labels"],
#     ncols=1,
#     save=f"_{out_name}.png"
# )

print("Train scVI model.", flush=True)
scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key=scvi_batch_key)
scvi_model = scvi.model.SCVI(adata, n_layers=2, n_latent=30, gene_likelihood=gene_likelihood, dispersion=dispersion)
scvi_model.train(
    max_epochs=max_epochs_scvi,
    train_size=0.9,
    log_every_n_steps=1,
    check_val_every_n_epoch=1,
    #early_stopping=True, early_stopping_patience=5,
    datasplitter_kwargs={"num_workers": num_workers}
    )

print("Save scVI model.", flush=True)
scvi_model.save(scvi_path, overwrite=True)

# get latent space representation
SCVI_LATENT_KEY = "X_scVI"
adata.obsm[SCVI_LATENT_KEY] = scvi_model.get_latent_representation()

print("Dimensionality reduction and clustering of scVI results.", flush=True)
scvi_neighbor_key = "neighbors_scvi"
sc.pp.neighbors(adata, use_rep=SCVI_LATENT_KEY, key_added=scvi_neighbor_key)
#sc.tl.leiden(adata, neighbors_key=scvi_neighbor_key, key_added="leiden_scvi")
sc.tl.umap(adata, neighbors_key=scvi_neighbor_key, key_added="X_umap_scvi", min_dist=0.3)

# print("Plot UMAP of scVI results.")
# sc.pl.embedding(
#     adata,
#     basis="X_umap_scvi",
#     color=["uid", "majority_voting", "predicted_labels"],
#     ncols=1,
#     save=f"_{out_name}.png"
# )

print("Save scVI results.", flush=True)
adata.write(scvi_adata_out_path)
# adata = sc.read(out_path)

#scvi_path = data_path.parent / f"scVI_model_rounded_{out_name}"
print(f"Loading scVI model from {scvi_path}")
scvi_model = scvi.model.SCVI.load(scvi_path, adata)

print("Train scANVI model", flush=True)
scanvi_model = scvi.model.SCANVI.from_scvi_model(
    scvi_model,
    adata=adata,
    labels_key=scanvi_labels_key,
    unlabeled_category="Unknown",
)

scanvi_model.train(
    max_epochs=max_epochs_scanvi,
    datasplitter_kwargs={"num_workers": num_workers},
    train_size=0.9,
    log_every_n_steps=1,
    check_val_every_n_epoch=1,
    n_samples_per_label=n_samples_per_label,
)

print("Save scANVI model.", flush=True)
scanvi_model.save(scanvi_path, overwrite=True)

# get latent space
SCANVI_LATENT_KEY = "X_scANVI"
adata.obsm[SCANVI_LATENT_KEY] = scanvi_model.get_latent_representation(adata)

print("Dimensionality reduction on scANVI results.", flush=True)
scanvi_neighbors_key = "neighbors_scanvi"
sc.pp.neighbors(adata, use_rep=SCANVI_LATENT_KEY, key_added=scanvi_neighbors_key)
sc.tl.umap(adata, neighbors_key=scanvi_neighbors_key, key_added="X_umap_scanvi", min_dist=0.3)

print("Save scANVI results.", flush=True)
adata.write(adata_out_path)

# print("Plot UMAP of scANVI results.", flush=True)
# sc.pl.embedding(
#     adata,
#     basis="X_umap_scanvi",
#     color=["uid", "majority_voting", "predicted_labels"],
#     ncols=1,
#     save=f"_{out_name}.png"
# )

print("Finished.")