# Delete environmental variables
rm(list=ls())

# Libraries
library(oligo)
library(GEOquery)

# Set working directory
setwd("/home/guille/Guille/MBHS/2Y/TFM/github/TFM/QC/gene_expression_array/")

# Download data from GEO database

# 

celfiles_path = "/home/guille/Guille/MBHS/2Y/TFM/QC/gene_expression_array/HuEx-1_0-st-v2-colon-cancer-data-set-agcc-files"
celfiles = list.files(path = celfiles_path, pattern = ".CEL$", all.files =
                        FALSE, full.names = FALSE, recursive = FALSE, ignore.case = FALSE)
