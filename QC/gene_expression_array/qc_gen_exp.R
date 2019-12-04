# Delete environmental variables
rm(list=ls())

# Libraries
library(oligo)
library(GEOquery)

# Set working directory and GEO ID
setwd("/home/guille/Guille/MBHS/2Y/TFM/QC/gene_expression_array/")
GEO_ID = "GSE124026"

# Download raw data from GEO database in the current directory
getGEOSuppFiles(GEO_ID, baseDir=getwd())

# Untar the raw data
path_to_GEO_folder = paste(getwd(), "/", GEO_ID, sep = "")
path_to_GEO_raw_tar_data = paste(path_to_GEO_folder, "/", GEO_ID, "_RAW.tar", sep = "")
untar(path_to_GEO_raw_tar_data, exdir = path_to_GEO_folder)
#dir.create(path_to_download)

# Obtain a vector of CEL files names 
path_to_GEO_raw_data = paste(path_to_GEO_folder, "/", GEO_ID, "_RAW", sep= "")
celfiles_list = list.files(path = path_to_GEO_raw_data, pattern = ".CEL.*", all.files =
                        FALSE, full.names = FALSE, recursive = FALSE, ignore.case = FALSE)
celfiles_list

# Obtain Series Matrix File(s), i.e., the samples metadata

#getGEOSuppFiles(GEO_ID, makeDirectory = TRUE, baseDir = getwd(), fetch_files = TRUE, filter_regex = NULL)
#setwd("/home/guille/Guille/MBHS/2Y/TFM/QC/gene_expression_array/GSE124026/prueba")

data <- getGEO(GEO_ID, GSEMatrix=TRUE)
phenodata <- phenoData(data$GSE124026_series_matrix.txt.gz)

# Set working directory where CEL files are stored
setwd(path_to_GEO_raw_data)

# Read CEL files + add metadata from OLIGO package
cel_files <- read.celfiles(filenames=celfiles_list, phenoData=phenodata)
cel_files


prueba <- read.celfiles(filenames=celfiles_list)


## MA-plot

xl <- c(2.8, 4)
yl <- c(-1, 1)
MAplot(cel_files[,1:3], pairs=TRUE, ylim=yl, xlim=xl)


## Box-plot

boxplot(cel_files, which="all", transfo=log2, nsample=10000)




