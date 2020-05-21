rm(list=ls())
start_time <- Sys.time()
date()

cat("\n\n ######################## START LOADING LIBRARIES ######################## \n\n", file = stdout())

# library(ff)
# cat("\n\n ff LOADED \n\n", file = stdout())
library(oligo)
cat("\n\n oligo LOADED \n\n", file = stdout())


cat("\n\n ######################## FINISHED LOADING LIBRARIES ######################## \n\n", file = stdout())

cat("\n\n ######################## 1.READ INPUT SCRIPT BEGINS ######################## \n\n", file = stdout())

#Sys.setenv(R_THREADS=5)

# ldPath()
# ocSamples()
# ocSamples(50)## changing default to 50
# ocProbesets()
# ocProbesets(100)#


cat("\n\n... SETTING WORKING DIRECTORY ...\n\n", file = stdout())
working.dir <- "/projects/regicor/Guille/TFM/transcriptomics"
setwd(working.dir)
getwd()

cat("\n\n... SETTING RAW DATA DIRECTORY ...\n\n", file = stdout())
rawdata.dir <- "/projects/regicor/Guille/TFM/transcriptomics/raw_data/final_raw_data/"
rawdata.dir

cat("\n\n... READING SAMPLES METADATA DATA ...\n\n", file = stdout())

samples.metadata <- read.csv(file = "phenotype/QC/samples_metadata_merged_QC_1818.csv",
                             header=TRUE, stringsAsFactors = F, sep=",")

rownames(samples.metadata) <- substr(samples.metadata$cel_files,1,23)

## BORRAR
# SUBSET DE PRUEBA 1300 INDIVIDUOS
#samples.metadata <- samples.metadata[c(1:10,70:80,100:110,300:310,900:910),]
#samples.metadata <- samples.metadata[1201:1818,]
samples.metadata <- samples.metadata[51:100,]

# MOVE THE SELECTED 1300 SAMPLES TO FINAL_RAW_DATA FILE
# 
# library(filesstrings)
# 
# rawdata.dir.2 <- "raw_data/phe000002.v7.FHS_SABRe_project3_OFF_CCS.raw-data-cel.c1/"
# raw.data.file.names <- substr(list.files(rawdata.dir.2),1,27)
# 
# lapply(samples.metadata$cel_files, function(ID) {
# 
#   file.matched <- grep(ID, raw.data.file.names, value = TRUE)
# 
#   if (length(file.matched) == 1) {
# 
#     file.move(paste("raw_data/phe000002.v7.FHS_SABRe_project3_OFF_CCS.raw-data-cel.c1/",file.matched,".gz",sep ="")
#               , "raw_data/final_raw_data/")
#   }
# 
#  })

## BORRAR

samples.metadata <- AnnotatedDataFrame(data = samples.metadata)
cat("\nSamples metadata\n", file = stdout())
str(samples.metadata)
dim(samples.metadata)
head(pData(samples.metadata))

# In read.celfiles(filenames = list.rawfiles, phenoData = samples.metadata) :
#   'channel' automatically added to varMetadata in phenoData.
# SOLVE:
# x <- data.frame(x, channel = "_ALL_")
# varMetadata(pd) <- x

cat("\n\n... SETTING WORKING DIRECTORY WHERE RAW DATA IS STORED ...\n\n", file = stdout())
setwd(rawdata.dir)
getwd()

list.rawfiles <- list.celfiles(rawdata.dir, listGzipped = TRUE)
head(list.rawfiles)

cat("\n\n... READING RAW DATA (.CEL FILES) AND METADATA ...\n\n", file = stdout())
#cat("\n\n MEMORY USAGE BEFORE \n\n", file = stderr())
#system("cat /proc/meminfo 1>&2") # redirect stdout to stderr file
# Read CEL files + add metadata from OLIGO package
ExonFeatureSet <- read.celfiles(filenames = list.rawfiles, phenoData = samples.metadata)
#cat("\n\n MEMORY USAGE AFTER \n\n", file = stderr())
#system("cat /proc/meminfo 1>&2") # redirect stdout to stderr file
cat("\n// ExonFeatureSet //\n", file = stdout())
getClass(ExonFeatureSet)

cat("\n\n... SAVING ExonFeatureSet OBJECT ...\n\n", file = stdout())

setwd(working.dir)

save(ExonFeatureSet, file="QC/1.read_input/ExonFeatureSet_51_100.RData")

date()

cat("\nComputational time\n", file = stdout())
end_time <- Sys.time()
print(end_time - start_time)

cat("\n\n ######################## 1.READ INPUT SCRIPT ENDS ######################## \n\n", file = stdout())











    