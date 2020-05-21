rm(list=ls())

date()

cat("\n\n ######################## START LOADING LIBRARIES ######################## \n\n", file = stdout())

library(minfi)
cat("\n\n minfi LOADED \n\n", file = stdout())
library(minfiData)
cat("\n\n minfiData LOADED \n\n", file = stdout())
library(methylumi)
cat("\n\n methylumi LOADED \n\n", file = stdout())
library(lumi)
cat("\n\n lumi LOADED \n\n", file = stdout())
library(wateRmelon)
cat("\n\n wateRmelon LOADED \n\n", file = stdout())
library(snpStats)
cat("\n\n snpStats LOADED \n\n", file = stdout())
library(Biobase)
cat("\n\n BioBase LOADED \n\n", file = stdout())
library(limma)
cat("\n\n limma LOADED \n\n", file = stdout())
library(sva)
cat("\n\n sva LOADED \n\n", file = stdout())
library(data.table)
cat("\n\n data.table LOADED \n\n", file = stdout())
library(sandwich)
cat("\n\n sandwich LOADED \n\n", file = stdout())
library(lmtest)
cat("\n\n lmtest LOADED \n\n", file = stdout())
library(gap)
cat("\n\n gap LOADED \n\n", file = stdout())
library(ggplot2)
cat("\n\n ggplot2 LOADED \n\n", file = stdout())
library(gdata)
cat("\n\n gdata LOADED \n\n", file = stdout())
library(IlluminaHumanMethylation450kmanifest)
cat("\n\n IlluminaHumanMethylation450kmanifest LOADED \n\n", file = stdout())

cat("\n\n ######################## FINISHED LOADING LIBRARIES ######################## \n\n", file = stdout())

cat("\n\n ######################## 1.READ INPUT SCRIPT BEGINS ######################## \n\n", file = stdout())

cat("\n\n... SETTING WORKING DIRECTORY ...\n\n", file = stdout())
setwd("/projects/regicor/Guille/TFM/methylation")
getwd()

cat("\n\n... SETTING RAW DATA DIRECTORY ...\n\n", file = stdout())
rawdata.dir <- "/projects/regicor/Guille/TFM/methylation/final_raw_data/"
rawdata.dir
#list.files(rawdata.dir)

cat("\n\n... READING SAMPLES METADATA DATA ...\n\n", file = stdout())

#cat("\nReformatting and saving samples metadata\n", file = stdout())
# 2123 individuals (from 2618 initially, because of NA of CVD.incidence)
# samples.metadata <- read.csv(file = "phenotype/samples_metadata_merged_filt.csv",
#                     header=TRUE, stringsAsFactors = F, sep=",")

samples.metadata <- read.csv(file = "phenotype/samples_metadata/QC/samples_metadata_merged_QC_2620.csv",
                             header=TRUE, stringsAsFactors = F, sep=",")

# write.csv(samples.metadata.filt, file = "phenotype/samples_metadata_merged_filt.csv",
#           row.names = FALSE)

# IMPORTANT: csv must be on same folder as raw_data to be read properly
# write.csv(samples.metadata, file = "final_raw_data/samples_metadata_merged_filt.csv",
#            row.names = FALSE)

write.csv(samples.metadata, file = "final_raw_data/samples_metadata_merged_QC_2620.csv",
          row.names = FALSE)

cat("\nReading new samples metadata\n", file = stdout())
# Give directory where .csv file is located

samples.metadata <- read.metharray.sheet(base = "final_raw_data", pattern = "samples_metadata_merged_QC_2620.csv$", 
                                         ignore.case = TRUE, recursive = TRUE, verbose = TRUE)
cat("\nSamples metadata\n", file = stdout())
head(samples.metadata)
# A new column is added to the metadata file: "basename", necessary for read.metharray.exp function
# So let's save this new samples metadata file
write.csv(samples.metadata,  file = "final_raw_data/samples_metadata_merged_QC_2620_targets.csv",
          row.names = FALSE)

cat("\n\n... READING RAW DATA (.IDAT FILES) ...\n\n", file = stdout())

RGChannelSetExtended <- read.metharray.exp(base = rawdata.dir,
                             targets = samples.metadata, 
                             extended = TRUE, 
                             force= TRUE,
                             recursive = FALSE, 
                             verbose = TRUE)

cat("\n// RGChannelSetExtended //\n", file = stdout())
getClass(RGChannelSetExtended)

cat("\n\n... SAVING RGChannelSetExtended OBJECT ...\n\n", file = stdout())

save(RGChannelSetExtended, file="QC/1.read_input/RGChannelSetExtended_2620.RData")

cat("\n\n ######################## 1.READ INPUT SCRIPT ENDS ######################## \n\n", file = stdout())

