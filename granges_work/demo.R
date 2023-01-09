if (!require("BiocManager"))
    install.packages("BiocManager")
BiocManager::install("GenomicRanges")

library(GenomicRanges)