# Delete environmental variables
rm(list=ls())

# Libraries
library(minfi)
library(minfiData)
library(wateRmelon)
library(limma)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylation450kmanifest)

# Load methylation data

metadata <- getGEO(GEO_ID, GSEMatrix=TRUE)
methyl_phenodata <- phenoData(data$`GSE103008-GPL13534_series_matrix.txt.gz`)

read.metharray.sheet(methyl_phenodata)

read.metharray.exp(targets=methyl_phenodata, extended = TRUE, force = TRUE,
                   base = NULL, recursive = FALSE, verbose = TRUE)

