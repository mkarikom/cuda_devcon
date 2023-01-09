from pathlib import Path

import scvi
import anndata
import numpy as np
import pandas as pd
import scanpy as sc
import episcanpy as epi
import matplotlib.pyplot as plt
import mypackage

# dirs
datadir = "/workspaces/cuda_devcon/scvi_work/data/multi_acvi" # the downloaded data
scvi_dir_untrained = "/workspaces/cuda_devcon/scvi_work/saved_models/untrained/gmax_pbmc"
scvi_dir_trained = "/workspaces/cuda_devcon/scvi_work/saved_models/trained/gmax_pbmc"
adata_dir = "/workspaces/cuda_devcon/scvi_work/saved_adata/gmax/test_pbmc"
fig_dir = '/workspaces/cuda_devcon/scvi_work/figures/gmax/test_pbmc'

# create and set figures directory
if not os.path.exists(fig_dir):
    print(f"creating fig dir: {fig_dir}")
    os.makedirs(fig_dir)
sc.settings.figdir = Path(fig_dir)

# download the data
if not os.path.exists(f"{datadir}/filtered_feature_bc_matrix"):
    print("downloading data to {datadir}/filtered_feature_bc_matrix")
    os.system(f"wget https://cf.10xgenomics.com/samples/cell-arc/2.0.0/pbmc_unsorted_10k/pbmc_unsorted_10k_filtered_feature_bc_matrix.tar.gz -P {datadir}")
    os.system(f"tar -xvf {datadir}/pbmc_unsorted_10k_filtered_feature_bc_matrix.tar.gz --directory {datadir}")
    os.system(f"gunzip -f {datadir}/filtered_feature_bc_matrix/*.gz")
else:
    print("found dir {datadir}/filtered_feature_bc_matrix, skipping download")

# initialize scvi settings
scvi.settings.seed = 420
sc.set_figure_params(figsize=(4, 4), frameon=False)

# create anndata
print(f"processing data")
adata = scvi.data.read_10x_multiome(f"{datadir}/filtered_feature_bc_matrix")
adata.var_names_make_unique()

# split a 10X dataset with 12012 cells into 3 datasets, taking either rna, atac, or both
n = 4004
adata_rna = adata[:n, adata.var.modality == "Gene Expression"].copy()
adata_paired = adata[n:2 * n].copy()
adata_atac = adata[2 * n:, adata.var.modality == "Peaks"].copy()

# inspect the rna data
sc.pl.highest_expr_genes(adata_paired,n_top=10,save=True)
adata_paired.var['mt'] = adata_paired.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(adata_paired, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

sc.pl.violin(adata_paired, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],jitter=0.4, multi_panel=True,save=True)
sc.pl.scatter(adata_paired, x='total_counts', y='pct_counts_mt',save="_pct-mt")
sc.pl.scatter(adata_paired, x='total_counts', y='n_genes_by_counts',save="_n_genes")

# do QC on the multimodal data using the RNA annotations as a guide
adata_paired = adata_paired[adata_paired.obs.n_genes_by_counts < 30000, :]
adata_paired = adata_paired[adata_paired.obs.pct_counts_mt < 20, :]

## get the ATAC clusters using PeakVI
# create and train the PeakVI model
adata_clust = adata_paired.copy()
scvi.model.PEAKVI.setup_anndata(adata_clust)
pvi = scvi.model.PEAKVI(adata_clust)
pvi.train()

if not os.path.exists(adata_dir):
    print(f"creating {adata_dir}")
    os.makedirs(adata_dir)
print(f"saving peakvi adata in {adata_dir}")
adata_clust.write(f"{adata_dir}/adata_peakvi.h5ad")
print(f"saving peakvi model in {scvi_dir_trained}")
pvi.save(prefix="trained_peakvi_",dir_path=scvi_dir_trained,overwrite=True)


## load the peakvi for metadata
adata_clust = anndata.read_h5ad(f"{adata_dir}/adata_peakvi.h5ad")
pvi = scvi.model.PEAKVI.load(prefix="trained_peakvi_",dir_path=scvi_dir_trained, adata=adata_clust)


# add the PeakVI latent representation to the observation matrices for adata_paired
latent = pvi.get_latent_representation()
adata_clust.obsm["X_PeakVI"] = latent
adata_paired.obsm["X_PeakVI"] = latent

# compute the k-nearest-neighbor graph that is used in both clustering and umap algorithms
sc.pp.neighbors(adata_clust, use_rep="X_PeakVI")
# compute the umap
sc.tl.umap(adata_clust, min_dist=0.2)

# cluster the space (we use a lower resolution to get fewer clusters than the default)
sc.tl.leiden(adata_clust, resolution=0.8, key_added='leiden80')
sc.tl.leiden(adata_clust, resolution=0.4, key_added='leiden40')
sc.tl.leiden(adata_clust, resolution=0.2, key_added='leiden20')
sc.tl.leiden(adata_clust, resolution=0.1, key_added='leiden10')
sc.tl.leiden(adata_clust, resolution=0.05, key_added='leiden05')

# plot the atac clusters
sc.pl.umap(adata_clust, color='leiden80',save=f"_atac_leiden80.pdf")
sc.pl.umap(adata_clust, color='leiden40',save=f"_atac_leiden40.pdf")
sc.pl.umap(adata_clust, color='leiden20',save=f"_atac_leiden20.pdf")
sc.pl.umap(adata_clust, color='leiden10',save=f"_atac_leiden10.pdf")
sc.pl.umap(adata_clust, color='leiden05',save=f"_atac_leiden05.pdf")

# copy the ATAC clustering to adata_paired, note that the columns 'leiden*' in adata_mvi.obs will be NaN unless the modality=='paired'
adata_paired.obs = adata_clust.obs

# organize multiome data and sort features
# this concatenates [paired,expression,accessibility]
adata_mvi = scvi.data.organize_multiome_anndatas(multi_anndata=adata_paired, rna_anndata=adata_rna, atac_anndata=adata_atac)
adata_mvi = adata_mvi[:, adata_mvi.var["modality"].argsort()].copy()

# add the 'Unknown' category for unpaired data clustering
for c in ['leiden80','leiden40','leiden20','leiden10','leiden05']:
    adata_mvi.obs[c] = adata_mvi.obs[c].cat.add_categories('Unknown')
    adata_mvi.obs[c] = adata_mvi.obs[c].fillna('Unknown')

# filter genes
print(adata_mvi.shape)
sc.pp.filter_genes(adata_mvi, min_cells=int(adata_mvi.shape[0] * 0.01))
print(adata_mvi.shape)

# adata_mvi.layers["counts"] will be stored as REGISTRY_KEYS.X_KEY during setup_anndata
adata_mvi.layers["counts"] = adata_mvi.X.copy() # preserve counts

# set up the model
# note that adata.layers["counts"] (must be explicitly saved as part of the scvi workflow) is not created in acvi following the mvi workflow
# todo make sure that the new acvi mulitmodal likelihood is appropriate for the 
mypackage.LinearMULTIACVI.setup_anndata(
    adata_mvi, 
    layer="counts",
    batch_key='modality',
    categorical_covariate_keys=['leiden80','leiden40','leiden20','leiden10','leiden05']
)

mvi = mypackage.LinearMULTIACVI(
    adata_mvi, 
    n_proteins=0,
    n_latent=16,
    n_genes=(adata_mvi.var['modality']=='Gene Expression').sum(),
    n_regions=(adata_mvi.var['modality']=='Peaks').sum(),
)
mvi.view_anndata_setup()

if not os.path.exists(adata_dir):
    print(f"creating {adata_dir}")
    os.makedirs(adata_dir)
print(f"saving processed adata in {adata_dir}")
adata_mvi.write(f"{adata_dir}/adata_acvi_multi.h5ad")
print(f"saving model in {scvi_dir_untrained}")
mvi.save(prefix="untrained_multi_",dir_path=scvi_dir_untrained,overwrite=True)

# load the saved model and data
adata_mvi = anndata.read_h5ad(f"{adata_dir}/adata_acvi_multi.h5ad")
mvi = mypackage.LinearMULTIACVI.load(prefix="untrained_multi_",dir_path=scvi_dir_untrained, adata=adata_mvi)

mvi.train()

# # save the trained model
# if not os.path.exists(scvi_dir_trained):
#     print(f"creating {scvi_dir_trained}")
#     os.makedirs(scvi_dir_trained)
# print(f"saving trained model")
# mvi.save(prefix="trained_multi_",dir_path=scvi_dir_trained,overwrite=True)

# # load the trained model
# adata_mvi = anndata.read_h5ad(f"{adata_dir}/adata_acvi_multi.h5ad")
# mvi = mypackage.LinearMULTIACVI.load(prefix="trained_multi_",dir_path=scvi_dir_trained, adata=adata_mvi)

# # inspect convergence metrics
# train_elbo = mvi.history['elbo_train'][1:]
# test_elbo = mvi.history['elbo_validation']
# epochs = np.arange(1,20,1)
# fig, ax = plt.subplots()
# ax.plot(train_elbo, label='train')
# ax.plot(test_elbo, label='test')
# ax.legend()
# ax.set_title('elbo metrics')
# print('saving elbo metrics plot')
# plt.savefig(f'{fig_dir}/elbo.pdf')

# # inspect the latent dimension
# Z_hat = mvi.get_latent_representation()
# for i, z in enumerate(Z_hat.T):
#     adata_mvi.obs[f'Z_{i}'] = z

# fig = plt.figure(figsize=(12, 8))

# for f in range(0, 9, 2):
#     plt.subplot(2, 3, int(f / 2) + 1)
#     plt.scatter(adata_mvi.obs[f'Z_{f}'], adata_mvi.obs[f'Z_{f + 1}'], marker='.', s=4, label='Cells')
#     plt.xlabel(f'Z_{f}')
#     plt.ylabel(f'Z_{f + 1}')

# plt.subplot(2, 3, 6)
# plt.scatter(adata_mvi.obs[f'Z_{f}'], adata_mvi.obs[f'Z_{f + 1}'], marker='.', label='Cells', s=4)
# plt.scatter(adata_mvi.obs[f'Z_{f}'], adata_mvi.obs[f'Z_{f + 1}'], c='w', label=None)
# plt.gca().set_frame_on(False)
# plt.gca().axis('off')

# lgd = plt.legend(scatterpoints=3, loc='upper left')
# for handle in lgd.legendHandles:
#     handle.set_sizes([200])


# plt.tight_layout()
# print('saving latent dim')
# plt.savefig(f'{fig_dir}/latent.pdf')

# # inspect the loadings
# loadings = mvi.get_loadings()
# loadings.head()

# # #############################################
# # # OPTION1: Use ALL GENES
# # #############################################

# # #### for genes x loadings matrix: M x N
# # ### get the correlation between all the genes: M x M
# # ## cluster the correlation matrix of the genes (using ONLY the weights across all of the atac clusters)
# # # this shows which genes are similar in their multi-scale linear relationship between accessibility and transcription
# # gene_trans = loadings.transpose()
# # gene_corr = gene_trans.corr()

# #############################################
# # OPTION2: Use the highly variable genes to save time
# #############################################

# scvi.data.poisson_gene_selection(adata_mvi)
# hvg_trans_ind = adata_mvi.var[(adata_mvi.var['highly_variable']==True) & (adata_mvi.var['modality']=='Gene Expression')].index

# #### for genes x loadings matrix: M x N
# ### get the correlation between all the genes: M x M
# ## cluster the correlation matrix of the genes (using ONLY the weights across all of the atac clusters)
# # this shows which genes are similar in their multi-scale linear relationship between accessibility and transcription
# gene_trans = loadings.loc[hvg_trans_ind,loadings.columns.str.startswith('K')].transpose()
# gene_corr = gene_trans.corr()

# import seaborn as sns
# print('computing and plotting the gene clustermap')
# sns.clustermap(gene_corr, figsize=(20,20))
# plt.savefig(f'{fig_dir}/gene_corr_clust.pdf', dpi=300)

# # ### get the correlation between latent dim loadings and atac loadings loadings 
# # ## cluster the correlation matrix of all the weights
# # # this shows whether certain transcriptional latent dimensions are related to accessibility at a particular scale
# weights_corr = loadings.corr()

# # subset the weights matrix so that 
# print('computing and plotting the z x weight cluster map')
# sns.clustermap(weights_corr.loc[weights_corr.index.str.startswith('Z'),weights_corr.columns.str.startswith('K')], figsize=(20,20))
# plt.savefig(f'{fig_dir}/latdim_atac_clust.pdf', dpi=300)
