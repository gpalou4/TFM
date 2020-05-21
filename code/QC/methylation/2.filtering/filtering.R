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

cat("\n\n ######################## 2.FILTERING SCRIPT BEGINS ######################## \n\n", file = stdout())

cat("\n\n... SETTING WORKING DIRECTORY ...\n\n", file = stdout())
setwd("/projects/regicor/Guille/TFM/methylation/QC")
getwd()

cat("\n\n... LOADING RGChannelSetExtended OBJECT ...\n\n", file = stdout())

RGCSext.name <- load("1.read_input/RGChannelSetExtended_2620.RData")
RGCSext <- get(RGCSext.name)
cat("\n// RGChannelSetExtended //\n", file = stdout())
getClass(RGCSext)
cat("\n// Illumina 450k Manifest //\n", file = stdout())
getManifest(RGCSext)

cat("\n\n... 1) REMOVING BAD QUALITY PROBES AND SAMPLES ...\n\n", file = stdout())

Methylset <- pfilter(RGCSext)

cat("\n// RGChannelSetExtended filtered --> MethylSet //\n", file = stdout())
getClass(Methylset)

cat("\n\n... 2) REMOVING CROSS-REACTIVE PROBES, NO CPGS, AND FROM ILLUMINA MANIFEST ...\n\n", file = stdout())

# ## 450k cpgs_to_discard ##
# cpgs.discard <- scan("2.filtering/cpgs_to_discard/cpgs_to_discard.txt", what = character(), quote = "")
# 
# cat("\nList of probes that should be removed\n", file = stdout())
# length(cpgs.discard)
# 
# cat("\nNumber of probes to be removed from our data\n", file = stdout())
# length(which(Methylset@NAMES%in%cpgs.discard))
# 
# Methylset.filt <- Methylset[-which(Methylset@NAMES%in%cpgs.discard),]
# # o probes <- filter_probes[filter_probes@NAMES%nin%discard$probes_discard,] #%nin% del paquete Hmisc es para no quedarse (nin) con las del primer argumento que están en el segundo
# 
# cat("\n// MethylSet filtered //\n", file = stdout())
# getClass(Methylset.filt)
# 
# cat("\n\n ... TESTING ... \n\n", file = stdout())

# prueba de que se han eliminado con un cpg que hay que eliminar
# which(Methylset.filt@NAMES=="cg00004260")
# [1] 208215
# which(probes@NAMES=="cg00004260")
# integer(0)

#Comprobación que las 4 muestras que no pasan el QC de pfilter() se han eliminado
# samples_rm<-c("202884930004_R01C01",
#               "202995740056_R04C01",
#               "203259060022_R08C01",
#               "203084910186_R08C01")
# length(filter_probes@colData$Sample_ID)
# # [1] 204
# length(which(filter_probes@colData$Sample_ID%in%samples_rm))
# # [1] 0
# # Las 4 muestras se eliminaron en el pfilter

cat("\n\n... SAVING FILTERED PROBES ...\n\n", file = stdout())
save(Methylset, file="/projects/regicor/Guille/TFM/methylation/QC/2.filtering/Methylset_2620_filt.RData")
#system("cat /proc/meminfo")

date()

cat("\n\n ######################## SCRIPT ENDS ######################## \n\n", file = stdout())
