rm(list=ls())

date()

cat("\n\n ######################## START LOADING LIBRARIES ######################## \n\n", file = stdout())

library(oligo)
cat("\n\n oligo LOADED \n\n", file = stdout())

cat("\n\n ######################## FINISHED LOADING LIBRARIES ######################## \n\n", file = stdout())

cat("\n\n ######################## 1.READ INPUT SCRIPT BEGINS ######################## \n\n", file = stdout())

cat("\n\n... SETTING WORKING DIRECTORY ...\n\n", file = stdout())
working.dir <- "/home/guille/Escritorio/Guille/TFM/QC/gene_expression_array"
setwd(working.dir)
getwd()

cat("\n\n... OBTAINING DIRECTORY OF RAW DATA ...\n\n", file = stdout())
##rawdata.dir <- "/projects/regicor/Guille/TFM/transcriptomics/raw_data/phe000002.v7.FHS_SABRe_project3_OFF_CCS.raw-data-cel.c1"
rawdata.dir <- "/home/guille/Escritorio/Guille/TFM/QC/gene_expression_array/mRNA_rawdata"
rawdata.dir

cat("\n\n... READING SAMPLES METADATA DATA ...\n\n", file = stdout())
# El plate es inventado
samples.metadata <- read.AnnotatedDataFrame(file = "phs000363.v17.pht002944.v4.p11.c1.Framingham_SABRe_Sample_Attributes.HMB-IRB-MDS.subset.txt", 
                             header = TRUE, sep ="\t")

cat("\n\n... SETTING WORKING DIRECTORY WHERE RAW DATA IS STORED ...\n\n", file = stdout())
setwd(rawdata.dir)
getwd()

list.files(rawdata.dir)

cat("\n\n... READING RAW DATA (.CEL FILES) AND METADATA ...\n\n", file = stdout())

# Read CEL files + add metadata from OLIGO package
#cel_files <- read.celfiles(filenames=celfiles_list, phenoData=phenodata)
EFSet <- read.celfiles(filenames=list.files(rawdata.dir), phenoData = samples.metadata)
cat("\n// ExonFeatureSet //\n", file = stdout())
getClass(EFSet)

cat("\n\n... SAVING ExonFeatureSet OBJECT ...\n\n", file = stdout())

save(EFSet, file="/home/guille/Escritorio/Guille/TFM/QC/gene_expression_array/1.read_input/ExonFeatureSet.RData")
##save(ExonFeatureSet, file="/projects/regicor/Guille/TFM/transcriptomics/QC/1.read_input/ExonFeatureSet.RData")

cat("\n\n ######################## 1.READ INPUT SCRIPT ENDS ######################## \n\n", file = stdout())

