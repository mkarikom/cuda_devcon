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

#assigning data directory 
data_dir <- "snRAWfiles"
#setting working directory 
setwd(data_dir)

#load seurat object 
pbmc <- readRDS("~/snRAWfiles/cv_73_11_02_22.rds")
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
    pre_embed_method = 'umap',
    perplexity = 20, 
  )

#throws error message 
#Error: umap: missing arguments: n_components







