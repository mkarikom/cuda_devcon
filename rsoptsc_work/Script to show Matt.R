#Libraries to load

library(dplyr)
library(Seurat)
library(patchwork)
library(sctransform)
library(ggplot2)
library(rgl)
require(data.table)
library(tidyverse)
library(Matrix)
library(RSoptSC)
library(ggrepel)

setwd("rsoptsc_work")

#load seurat object 
pbmc <- readRDS("snRAWfiles/cv_73_11_02_22.rds")
pbmc #sample cv_73,  1,384 nuclei , has been log normalized and found 1000 variable features

SimilaritySeurat <- function(obj, minvar_gene=0,minvar_cell=0,...){
  
  # make sure there are no zero variance cells
  genes = names(apply(Seurat::GetAssayData(obj)[which(apply(Seurat::GetAssayData(obj)[Seurat::VariableFeatures(obj),], 1, var)>minvar_gene),], 1, var))
  cells = names(apply(Seurat::GetAssayData(obj)[which(apply(Seurat::GetAssayData(obj)[Seurat::VariableFeatures(obj),], 2, var)>minvar_cell),], 2, var))
  
  SimilarityM(data=Seurat::GetAssayData(obj)[genes,cells],...)
}

sim = 
  SimilaritySeurat(
    obj=pbmc[,1:500],
    lambda = 0.05, 
    dims = 3,
    pre_embed_method = 'tsne',
    perplexity = 20, 
  )

low_dim_mapping <- RepresentationMap(similarity_matrix = sim$W,
                                     flat_embedding_method = 'tsne',
                                     join_components = TRUE,
                                     perplexity = 35,
                                     theta = 0.5,
                                     normalize = FALSE,
                                     pca = TRUE,
                                     pca_center = TRUE,
                                     pca_scale = TRUE,
                                     dims = 2,
                                     initial_dims = 2)

clusters <- ClusterCells(similarityMatrix = sim$W, n_comp = 15, .options='p')
H <- clusters$H
labels <- clusters$labels
n_clusters <- length(unique(clusters$labels))


pdf("snPlots/eigs.pdf")
plot(c(1:20), 
     clusters$ensemble$eigs$val[1:20],
     xlab = NA,
     ylab = 'eigenvalues',
     main = 'Eigenvalues of the Graph Laplacian')
dev.off()

# define a scheme for the coloring.  This is a required parameter for plotting discrete (factor) data
colorscale <- ColorHue(n = length(unique(labels)))
colorscale <- colorscale$hex.1.n.
# plot clusters

pdf("snPlots/featureplot_clusters.pdf")
FeatureScatterPlot(flat_embedding = low_dim_mapping$flat_embedding,
                   feature = as.factor(labels),
                   title = "NMF Cluster Labeling",
                   subtitle = "t-SNE Embedding",
                   featurename = "Cluster ID",
                   colorscale = colorscale)
dev.off()

markers <- GetMarkerTable(counts_data = Seurat::GetAssayData(pbmc)[,1:500],
                          cluster_labels = labels,
                          H = H,
                          n_sorted = 25)

pdf("snPlots/featureplot_marker_heatmap.pdf")
PlotTopN_Grid(data = Seurat::GetAssayData(pbmc)[,1:500],
              cluster_labels = labels,
              markers = markers$all,
              n_features = 15)
dev.off()

