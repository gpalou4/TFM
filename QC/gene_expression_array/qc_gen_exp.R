# Delete environmental variables
rm(list=ls())

# Libraries
library(oligo)
library(GEOquery)

### DOWNLOAD DATA ###

# Set working directory and GEO ID
current_directory = "/home/guille/Guille/MBHS/2Y/TFM/QC/gene_expression_array/"
setwd(current_directory)
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

### QUALITY CONTROL ###

# Set working directory where CEL files are stored
setwd(path_to_GEO_raw_data)

# Read CEL files + add metadata from OLIGO package
cel_files <- read.celfiles(filenames=celfiles_list, phenoData=phenodata)
cel_files

# imagenes

oligo::image(cel_files)

## MA-plot

# Important to take into account the X and Y scales, it should be more or less
# the same in all samples

# Save the MA plots in a pdf in the working directory
setwd(current_directory)
pdf("MA_plots")
MAplot(cel_files)
dev.off()

## Box-plot

# HACER PLOT POR BATCHES (1200 INDV)

# Level of summarization? probeset, core, full or extended?
oligo::boxplot(cel_files, target = "probeset", transfo=log2, nsample=10000)
showMethods(boxplot, class = "FeatureSet", includeDefs=T)

## Density plot

# hacer bloques de 20 para visualziar bien los colores y poder descartar el malo..
oligo::hist(cel_files, target = "probeset")

## Probe Level Models: RLE y NUSE
#si RLE da mal, en NUSE da peor, y descartar

fit_PLM <- fitProbeLevelModel(cel_files)

# Relative Log Expression (RLE) boxplot
RLE(fit_PLM)

# Normalized Unscaled Standard Errors (NUSE)
# Everything above 1 is bad quality samples

NUSE(fit_PLM)

### PRE-PROCESSING ###

# Background substraction, normalization and summarization by RMA method
## preprocessing step by spte or simultaneously with rma()
# me transforma a ExpressionSet

cel_files_norm <- rma(cel_files)








