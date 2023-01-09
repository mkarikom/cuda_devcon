import sys
from pathlib import Path
from importlib import reload  # Python 3.4+
import scvi
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt

sys.path.append("scvi_work/helper_functions")
from data_handling_new import load_inflam_new
from atac_clustering_newer import add_atac_clusters_newer
from backup_utils import makepath

import os
os.environ["TF_CPP_MIN_LOG_LEVEL"] = "0"

# dirs
datadir = "/workspaces/cuda_devcon/scvi_work/data/boukhaled_2022_processed_multi_geoquery" # the signac-processed multiome data
unzipdatadir = "/workspaces/cuda_devcon/scvi_work/data/boukhaled_2022_processed_multi_geoquery_unzipped"
metadir = "/workspaces/cuda_devcon/scvi_work/data/boukhaled_2022_geoquery/metadata" # the phenotype data and gene activity data
scvi_dir_untrained = "/workspaces/cuda_devcon/scvi_work/saved_models/untrained/multivi/inflam"
scvi_dir_trained = "/workspaces/cuda_devcon/scvi_work/saved_models/trained/multivi/inflam"
adata_dir = "/workspaces/cuda_devcon/scvi_work/saved_adata/multivi/inflam"
fig_dir = '/workspaces/cuda_devcon/scvi_work/figures/multivi/inflam'
serializedir = '/workspaces/cuda_devcon/scvi_work/serialized/multivi/inflam'

makepath(fig_dir)
sc.settings.figdir = Path(fig_dir)

makepath(scvi_dir_untrained)
data = load_inflam_new(datadir,metadir,serializedir)

# integrate and concatenate the data into a single anndata object
adata_list = list(data['adata'].values())
gact_list = list(data['gact'].values())

adata_concat_full = adata_list[0].concatenate(*adata_list[1:],index_unique=None)
gact_concat_full = gact_list[0].concatenate(*gact_list[1:],index_unique=None)

# only take the barcodes with both Multiome and GACT
adata_concat = adata_concat_full[adata_concat_full.obs.index.intersection(gact_concat_full.obs.index)]
gact_concat = gact_concat_full[gact_concat_full.obs.index.intersection(adata_concat_full.obs.index)]

# genes in adata_concat RNA but not in gact_concat (~18249)
del_adata_var = adata_concat[:,adata_concat.var.modality=="Gene Expression"].var.index.difference(gact_concat.var.index)
keep_adata_var = adata_concat.var.index.difference(del_adata_var)

# genes in gact_concat AND adata_concat RNA (~18352)
keep_gact_common = adata_concat[:,adata_concat.var.modality=="Gene Expression"].var.index.intersection(gact_concat.var.index)

# only take the adata_concat.var.modality="Gene Expression" that have a match in the GACT
adata_concat = adata_concat[:,keep_adata_var]
gact_concat = gact_concat[:,keep_gact_common]


# run the atac clustering
# initialize scvi settings
scvi.settings.seed = 420
sc.set_figure_params(figsize=(4, 4), frameon=False)

# inspect the rna data
sc.pl.highest_expr_genes(adata_concat,n_top=10,save=True)
adata_concat.var['mt'] = adata_concat.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(adata_concat, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

sc.pl.violin(adata_concat, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],jitter=0.4, multi_panel=True,save=True)
sc.pl.scatter(adata_concat, x='total_counts', y='pct_counts_mt',save="_pct-mt")
sc.pl.scatter(adata_concat, x='total_counts', y='n_genes_by_counts',save="_n_genes")

# do QC on the multimodal data using the RNA annotations as a guide
adata_concat = adata_concat[adata_concat.obs.n_genes_by_counts < 30000, :]
adata_concat = adata_concat[adata_concat.obs.pct_counts_mt < 20, :]

## get the ATAC clusters using PeakVI
makepath(adata_dir)
adata_concat = add_atac_clusters_newer(adata_concat,adata_dir,scvi_dir_trained,device="cuda:0")

# organize multiome data and sort features
# this concatenates [paired,expression,accessibility]
# adata_mvi = scvi.data.organize_multiome_anndatas(multi_anndata=adata_concat, rna_anndata=adata_rna, atac_anndata=adata_atac)
# adata_mvi = adata_mvi[:, adata_mvi.var["modality"].argsort()].copy()
adata_mvi = scvi.data.organize_multiome_anndatas(multi_anndata=adata_concat)

# # add the 'Unknown' category for unpaired data clustering
# for c in ['leiden80','leiden40','leiden20','leiden10','leiden05']:
#     adata_mvi.obs[c] = adata_mvi.obs[c].cat.add_categories('Unknown')
#     adata_mvi.obs[c] = adata_mvi.obs[c].fillna('Unknown')

# # filter genes
# print(adata_mvi.shape)
# sc.pp.filter_genes(adata_mvi, min_cells=int(adata_mvi.shape[0] * 0.01))
# print(adata_mvi.shape)

# adata_mvi.layers["counts"] will be stored as REGISTRY_KEYS.X_KEY during setup_anndata
adata_mvi.layers["counts"] = adata_mvi.X.copy() # preserve counts

# set up the model
# note that adata.layers["counts"] (must be explicitly saved as part of the scvi workflow) is not created in acvi following the mvi workflow
# todo make sure that the new acvi mulitmodal likelihood is appropriate for the 

scvi.model.MULTIVI.setup_anndata(adata_mvi, batch_key='modality')

mvi = scvi.model.MULTIVI(
    adata_mvi,
    n_latent=16,
    n_genes=(adata_mvi.var['modality']=='Gene Expression').sum(),
    n_regions=(adata_mvi.var['modality']=='Peaks').sum(),
)
mvi.view_anndata_setup()

print(f"saving processed adata in {adata_dir}")
adata_mvi.write(f"{adata_dir}/multi.h5ad")
print(f"saving model in {scvi_dir_untrained}")
mvi.save(prefix="untrained_multi_",dir_path=scvi_dir_untrained,overwrite=True)

# # load the saved model and data
# adata_mvi = anndata.read_h5ad(f"{adata_dir}/adata_acvi_multi.h5ad")
# mvi = mypackage.LinearMULTIACVI.load(prefix="untrained_multi_",dir_path=scvi_dir_untrained, adata=adata_mvi)

mvi.train()

# save the trained model
makepath(scvi_dir_trained)
print(f"saving trained model")
mvi.save(prefix="trained_multi_",dir_path=scvi_dir_trained,overwrite=True)

# # load the trained model
# adata_mvi = anndata.read_h5ad(f"{adata_dir}/adata_acvi_multi.h5ad")
# mvi = mypackage.LinearMULTIACVI.load(prefix="trained_multi_",dir_path=scvi_dir_trained, adata=adata_mvi)

# inspect convergence metrics
train_elbo = mvi.history['elbo_train'][1:]
test_elbo = mvi.history['elbo_validation']
epochs = np.arange(1,20,1)
fig, ax = plt.subplots()
ax.plot(train_elbo, label='train')
ax.plot(test_elbo, label='test')
ax.legend()
ax.set_title('elbo metrics')
print('saving elbo metrics plot')
plt.savefig(f'{fig_dir}/elbo.pdf')

# inspect the latent dimension
Z_hat = mvi.get_latent_representation()
for i, z in enumerate(Z_hat.T):
    adata_mvi.obs[f'Z_{i}'] = z

fig = plt.figure(figsize=(12, 8))

for f in range(0, 9, 2):
    plt.subplot(2, 3, int(f / 2) + 1)
    plt.scatter(adata_mvi.obs[f'Z_{f}'], adata_mvi.obs[f'Z_{f + 1}'], marker='.', s=4, label='Cells')
    plt.xlabel(f'Z_{f}')
    plt.ylabel(f'Z_{f + 1}')

plt.subplot(2, 3, 6)
plt.scatter(adata_mvi.obs[f'Z_{f}'], adata_mvi.obs[f'Z_{f + 1}'], marker='.', label='Cells', s=4)
plt.scatter(adata_mvi.obs[f'Z_{f}'], adata_mvi.obs[f'Z_{f + 1}'], c='w', label=None)
plt.gca().set_frame_on(False)
plt.gca().axis('off')

lgd = plt.legend(scatterpoints=3, loc='upper left')
for handle in lgd.legendHandles:
    handle.set_sizes([200])


plt.tight_layout()
print('saving latent dim')
plt.savefig(f'{fig_dir}/latent.pdf')

# inspect the loadings
loadings = mvi.get_loadings()
loadings.head()


# # get the highly variable genes to save time
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