from pathlib import Path
import sys
path_root = Path(__file__).parents[2]
sys.path.append(str(path_root))
print(sys.path)

import scvi
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import os

# download the data
datadir = "/workspaces/cuda_devcon/data/peakvi"
scvi_dir_untrained = "/workspaces/cuda_devcon/saved_models/untrained/peakvi/pbmc5k"
scvi_dir_trained = "/workspaces/cuda_devcon/saved_models/trained/peakvi/pbmc5k"
fig_dir = '/workspaces/cuda_devcon/figures/peakvi/pbmc5k'

if not os.path.exists(datadir):
    print(f"creating {datadir}")
    os.makedirs(datadir)
if not os.path.exists(scvi_dir_trained):
    print(f"creating {scvi_dir_trained}")
    os.makedirs(scvi_dir_trained)
if not os.path.exists(scvi_dir_untrained):
    print(f"creating {scvi_dir_untrained}")
    os.makedirs(scvi_dir_untrained)
if not os.path.exists(fig_dir):
    print(f"creating {fig_dir}")
    os.makedirs(fig_dir)


sc.settings.figdir=fig_dir

os.system(f"wget https://cf.10xgenomics.com/samples/cell-atac/1.2.0/atac_pbmc_5k_nextgem/atac_pbmc_5k_nextgem_filtered_peak_bc_matrix.tar.gz -P {datadir}")
os.system(f"tar -xvf {datadir}/atac_pbmc_5k_nextgem_filtered_peak_bc_matrix.tar.gz --directory {datadir}")

# initialize scvi settings
scvi.settings.seed = 420
sc.set_figure_params(figsize=(4, 4), frameon=False)

# or use this methods to read 10x data directly
adata = scvi.data.read_10x_atac(f"{datadir}/filtered_peak_bc_matrix")

# analysis
print(adata.shape)
# compute the threshold: 5% of the cells
min_cells = int(adata.shape[0] * 0.05)
# in-place filtering of regions
sc.pp.filter_genes(adata, min_cells=min_cells)
print(adata.shape)

scvi.model.PEAKVI.setup_anndata(adata)
pvi = scvi.model.PEAKVI(adata)
pvi.save(dir_path=scvi_dir_untrained, overwrite=True)

# train the model
pvi.train()
pvi.save(dir_path=scvi_dir_trained, overwrite=True)

# analyze the latent space
latent = pvi.get_latent_representation()
adata.obsm["X_PeakVI"] = latent
print(latent.shape)

# compute the k-nearest-neighbor graph that is used in both clustering and umap algorithms
sc.pp.neighbors(adata, use_rep="X_PeakVI")
# compute the umap
sc.tl.umap(adata, min_dist=0.2)
# cluster the space (we use a lower resolution to get fewer clusters than the default)
sc.tl.leiden(adata, key_added="cluster_pvi", resolution=0.2)
sc.pl.umap(adata, color='cluster_pvi',save="_cluster")