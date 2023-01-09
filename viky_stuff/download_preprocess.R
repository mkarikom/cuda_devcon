library(GEOquery)
library(readr)
library(stringr)
library(Seurat)

setwd("/workspaces/cuda_devcon/viky_stuff")
datadir = "data"
geoacc = "GSE159812"

# filter the list based on Viky 8 datasets
gsmkeep = c(
    "GSM4848457",
    "GSM4848459",
    "GSM4848460",
    "GSM4848461",
    "GSM4848451",
    "GSM4848452",
    "GSM4848453",
    "GSM4848455")

# gsmkeep = c(
#     "GSM4848457",
#     "GSM4848459",
#     "GSM4848460",
#     "GSM4848461")

gse <- getGEO(geoacc, GSEMatrix = TRUE)

dir.create(file.path(datadir,"raw"),recursive=TRUE)
filePaths = getGEOSuppFiles(geoacc,baseDir = datadir)
untar(rownames(filePaths),exdir = file.path(datadir,"raw"))

rawfiles = list.files(file.path(datadir,"raw"))
pdata = pData(gse[[1]])
seur_list = list()
for(geo in pdata$geo_accession){
    print(sprintf("detecting accession:%s",geo))
    if(geo %in% gsmkeep){
        dir.create(file.path(datadir,"zipped",geo),recursive=TRUE)
        detected = rawfiles[str_detect(rawfiles,geo)]
        for(dtfile in detected){
            file.copy(
                file.path(datadir,"raw",dtfile),
                file.path(datadir,"zipped",geo))

            file.rename(
                file.path(datadir,"zipped",geo,dtfile),
                file.path(datadir,"zipped",geo,tail(str_split(dtfile,"_")[[1]],1))
            )
        }
        gsmdir = file.path(datadir,"zipped",geo)
        print(sprintf("reading 10x from %s",gsmdir))
        expression_matrix <- Read10X(data.dir = gsmdir)
        seur = CreateSeuratObject(counts = expression_matrix)

        seur@meta.data["sex"] = rep(pdata[geo,"Sex:ch1"],dim(seur)[2])
        seur@meta.data["disease"] = rep(pdata[geo,"disease group:ch1"],dim(seur)[2])

        seur[["percent.mt"]] <- PercentageFeatureSet(seur, pattern = "^MT-")
        seur <- subset(seur, subset = percent.mt < 5)

        seur_list[[geo]] = seur
    }else{
        print(sprintf("%s is not in keep list, skipping",geo))
    }
}

# varialble features
seur_list_norm <- lapply(X = seur_list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# common variable features
features <- SelectIntegrationFeatures(object.list = seur_list_norm)

# integration features
anchors <- FindIntegrationAnchors(object.list = seur_list_norm, anchor.features = features)

# integrate data
combined <- IntegrateData(anchorset = anchors)

dir.create(file.path(datadir,"saved"))
saveRDS(combined,file.path(datadir,"saved","combined.RDS"))