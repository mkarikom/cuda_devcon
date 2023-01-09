setRepositories(ind=1:3) # needed to automatically install Bioconductor dependencies
install.packages("RCurl") # needs libcurl4-openssl-dev
install.packages("GenomicRanges")
install.packages("Signac")
install.packages("future")