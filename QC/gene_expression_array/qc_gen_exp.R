# Delete environmental variables
rm(list=ls())

# Libraries
library(oligo)
library(GEOquery)

# Set working directory and GEO ID
setwd("/home/guille/Guille/MBHS/2Y/TFM/QC/gene_expression_array/")
GEO_ID = "GSE124026"

# Download data from GEO database
getGEOSuppFiles(GEO_ID, baseDir=getwd())

?getGEO

# Untar the raw data
path_to_GEO_folder = paste(getwd(), "/", GEO_ID, sep = "")
path_to_GEO_raw_tar_data = paste(path_to_GEO_folder, "/", GEO_ID, "_RAW.tar", sep = "")
untar(path_to_GEO_raw_tar_data, exdir = path_to_GEO_folder)
#dir.create(path_to_download)

# Obtain CEL files names 
path_to_GEO_raw_data = paste(path_to_GEO_folder, "/", GEO_ID, "_RAW", sep= "")
celfiles = list.files(path = path_to_GEO_raw_data, pattern = ".CEL.*", all.files =
                        FALSE, full.names = FALSE, recursive = FALSE, ignore.case = FALSE)
celfiles

# Obtain Series Matrix File(s), i.e., the samples metadata

getGEOSuppFiles(GEO_ID, makeDirectory = TRUE, baseDir = getwd(), fetch_files = TRUE, filter_regex = NULL)

# Set working directory where CEL files are stored
setwd(celfiles_path)

# Obtain metadata from samples



# Read CEL files from OLIGO package
read.celfiles(filenames=celfiles)












