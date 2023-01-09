import os
import scvi
import scanpy as sc
import anndata

def add_atac_clusters_newer(adata_concat,adata_dir,scvi_dir_trained,device="cuda:1"):
    if not (os.path.exists(f"{adata_dir}/adata_peakvi.h5ad") and os.path.exists(f"{scvi_dir_trained}/trained_peakvi_model.pt")):
        print(f"no atac clustering detected, running peakvi")
        adata_clust = adata_concat[:,adata_concat.var.modality == 'Peaks'].copy()
        scvi.model.PEAKVI.setup_anndata(adata_clust)
        pvi = scvi.model.PEAKVI(adata_clust)

        print(f"using device {device}")
        pvi.to_device(device)
        pvi.train(use_gpu=device)

        # add the PeakVI latent representation to the observation matrices for adata_concat
        latent = pvi.get_latent_representation()
        adata_clust.obsm["X_PeakVI"] = latent

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

        # copy the ATAC clustering to adata_concat, note that the columns 'leiden*' in adata_mvi.obs will be NaN unless the modality=='paired'
        adata_concat.obs = adata_clust.obs
        adata_concat.obsm["X_PeakVI"] = latent
        # # pre-run and coment this for interactive repl pasting
        if not os.path.exists(adata_dir):
            print(f"creating {adata_dir}")
            os.makedirs(adata_dir)
        print(f"saving peakvi adata in {adata_dir}")
        adata_clust.write(f"{adata_dir}/adata_peakvi.h5ad")
        print(f"saving peakvi model in {scvi_dir_trained}")
        pvi.save(prefix="trained_peakvi_",dir_path=scvi_dir_trained,overwrite=True)
    else:
        ## load the peakvi for metadata
        adata_clust = anndata.read_h5ad(f"{adata_dir}/adata_peakvi.h5ad")
        pvi = scvi.model.PEAKVI.load(prefix="trained_peakvi_",dir_path=scvi_dir_trained, adata=adata_clust)
    
        # copy the ATAC clustering to adata_concat, note that the columns 'leiden*' in adata_mvi.obs will be NaN unless the modality=='paired'
        adata_concat.obs = adata_clust.obs
        adata_concat.obsm["X_PeakVI"] = adata_clust.obsm["X_PeakVI"]

    return adata_concat