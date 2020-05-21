rm(list=ls())

date()

cat("\n\n ######################## START LOADING LIBRARIES ######################## \n\n", file = stdout())

library(minfi)
cat("\n\n minfi LOADED \n\n", file = stdout())
library(wateRmelon)
cat("\n\n wateRmelon LOADED \n\n", file = stdout())
library(RnBeads)
cat("\n\n RnBeads LOADED \n\n", file = stdout())
library(lumi)
cat("\n\n lumi LOADED \n\n", file = stdout())
library(IlluminaHumanMethylation450kmanifest)
cat("\n\n IlluminaHumanMethylation450kmanifest LOADED \n\n", file = stdout())

cat("\n\n ######################## FINISHED LOADING LIBRARIES ######################## \n\n", file = stdout())


cat("\n\n ######################## 6.NORMALIZATION OF M VALUES SCRIPT BEGINS ######################## \n\n", file = stdout())

cat("\n\n... SETTING WORKING DIRECTORY ...\n\n", file = stdout())
setwd("/projects/regicor/Guille/TFM/methylation/QC")
getwd()

# cat("\n\n... LOADING Methylset.filt2 OBJECT ...\n\n", file = stdout())
# Methylset.filt.name <- load("4.sex/Methylset_2620_filt3.RData")
# Methylset.filt <- get(Methylset.filt.name)
# cat("\n// Methylset.filt3 //\n", file = stdout())
# Methylset.filt
# 
# cat("\n\n... NORMALIZING BETAS WITH DASEN ...\n\n", file = stdout())
# 
# Methylset.norm <- dasen(Methylset.filt)
# cat("\n// Methylset.norm //\n", file = stdout())
# getClass(Methylset.norm)
# 
# cat("\n// BETAS BEFORE DASEN //\n", file = stdout())
# betas <- (Methylset.filt@assays$data$Meth / (Methylset.filt@assays$data$Meth + Methylset.filt@assays$data$Unmeth + 100))
# betas[1:5,1:5]
# 
# cat("\n// BETAS AFTER DASEN //\n", file = stdout())
# betas.norm <- (Methylset.norm@assays$data$Meth / (Methylset.norm@assays$data$Meth + Methylset.norm@assays$data$Unmeth + 100))
# betas.norm[1:5,1:5]
# 
# cat("\n\n... SAVING NORMALIZED BETAS ...\n\n", file = stdout())
# #save(betas.norm, file="6.normalization/betas_norm.RData")
# 
# cat("\n\n... SAVING NORMALIZED M VALUES ...\n\n", file = stdout())
# 
# cat("\n// M VALUES BEFORE DASEN //\n", file = stdout())
# mvals <- log2((Methylset.filt@assays$data$Meth + 1) / (Methylset.filt@assays$data$Unmeth + 1))
# mvals[1:5,1:5]
# 
# cat("\n// METH VALUES BEFORE DASEN //\n", file = stdout())
# 
# Methylset.filt@assays$data$Meth[1:5,1:5]
# 
# cat("\n// UNMETH VALUES BEFORE DASEN //\n", file = stdout())
# 
# Methylset.filt@assays$data$Unmeth[1:5,1:5]
# 
# cat("\n// M VALUES WITH BETA2M TRANSFORMATION //\n", file = stdout())
# mvals.norm <- beta2m(betas.norm)
# mvals.norm[1:5,1:5]
# 
# cat("\n// M VALUES AFTER DASEN //\n", file = stdout())
# mvals.norm <- log2((Methylset.norm@assays$data$Meth + 1) / (Methylset.norm@assays$data$Unmeth + 1))
# mvals.norm[1:5,1:5]
# 
# cat("\n// METH VALUES AFTER DASEN //\n", file = stdout())
# 
# Methylset.filt@assays$data$Meth[1:5,1:5]
# 
# cat("\n// UNMETH VALUES AFTER DASEN //\n", file = stdout())
# 
# Methylset.filt@assays$data$Unmeth[1:5,1:5]

#save(mvals.norm, file="6.normalization/m_norm.RData")



###COMPROBACION MVALS ESTAN NORMALIZADOS###

cat("\n\n... LOADING betas norm ...\n\n", file = stdout())
betas.norm.name <- load("6.normalization/betas_norm.RData")
betas.norm <- get(betas.norm.name)
cat("\n// betas.norm //\n", file = stdout())
betas.norm[1:5,1:5]

cat("\n\n... LOADING M vals norm ...\n\n", file = stdout())
mvals.norm.name <- load("6.normalization/m_norm.RData")
mvals.norm <- get(mvals.norm.name)
cat("\n// mvals.norm //\n", file = stdout())
mvals.norm[1:5,1:5]

cat("\n// mvals.norm with beta2m transformation//\n", file = stdout())
mvals.norm.2 <- beta2m(betas.norm)
mvals.norm.2[1:5,1:5]

cat("\n// M VALUES BEFORE DASEN //\n", file = stdout())

cat("\n\n... LOADING Methylset.filt2 OBJECT ...\n\n", file = stdout())
Methylset.filt.name <- load("4.sex/Methylset_2620_filt3.RData")
Methylset.filt <- get(Methylset.filt.name)
cat("\n// mvals //\n", file = stdout())
mvals <- log2((Methylset.filt@assays$data$Meth + 1) / (Methylset.filt@assays$data$Unmeth + 1))
mvals[1:5,1:5]






cat("\n\n ######################## 6.NORMALIZATION OF M VALUES SCRIPT ENDS ######################## \n\n", file = stdout())

