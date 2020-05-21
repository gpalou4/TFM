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
library("IlluminaHumanMethylation450kmanifest")
cat("\n\n IlluminaHumanMethylation450kmanifest LOADED \n\n", file = stdout())

cat("\n\n ######################## FINISHED LOADING LIBRARIES ######################## \n\n", file = stdout())

cat("\n\n ######################## 2.FILTERING-2 SCRIPT BEGINS ######################## \n\n", file = stdout())

cat("\n\n... SETTING WORKING DIRECTORY ...\n\n", file = stdout())
setwd("/projects/regicor/Guille/TFM/methylation/QC")
getwd()

cat("\n\n... LOADING Methylset OBJECT ...\n\n", file = stdout())

Methylset.name <- load("2.filtering/Methylset_2620_filt.RData")
Methylset <- get(Methylset.name)
cat("\n// Methylset.filt //\n", file = stdout())
getClass(Methylset)
cat("\n// Illumina 450k Manifest //\n", file = stdout())
getManifest(Methylset)

cat("\n\n... 2) REMOVING CROSS-REACTIVE PROBES, NO CPGS, AND FROM ILLUMINA MANIFEST ...\n\n", file = stdout())

## 450k cpgs_to_discard ##
cpgs.discard <- scan("2.filtering/cpgs_to_discard/cpgs_to_discard.txt", what = character(), quote = "")

cat("\nList of probes that should be removed\n", file = stdout())
length(cpgs.discard)

cat("\nNumber of probes to be removed from our data\n", file = stdout())
length(which(Methylset@NAMES%in%cpgs.discard))

Methylset.filt <- Methylset[-which(Methylset@NAMES%in%cpgs.discard),]
# o probes <- filter_probes[filter_probes@NAMES%nin%discard$probes_discard,] #%nin% del paquete Hmisc es para no quedarse (nin) con las del primer argumento que están en el segundo

cat("\n// MethylSet filtered //\n", file = stdout())
getClass(Methylset.filt)

cat("\n\n... SAVING FILTERED PROBES ...\n\n", file = stdout())
save(Methylset.filt, file="/projects/regicor/Guille/TFM/methylation/QC/2.filtering/Methylset_2620_filt2.RData")
#system("cat /proc/meminfo")

date()

cat("\n\n ######################## SCRIPT ENDS ######################## \n\n", file = stdout())
