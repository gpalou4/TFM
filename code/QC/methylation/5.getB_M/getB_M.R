rm(list=ls())

date()

cat("\n\n ######################## START LOADING LIBRARIES ######################## \n\n", file = stdout())

library(minfi)
cat("\n\n minfi LOADED \n\n", file = stdout())

cat("\n\n ######################## FINISHED LOADING LIBRARIES ######################## \n\n", file = stdout())

cat("\n\n ######################## 5.GET BETA AND M VALUES SCRIPT BEGINS ######################## \n\n", file = stdout())

cat("\n\n... SETTING WORKING DIRECTORY ...\n\n", file = stdout())
setwd("/projects/regicor/Guille/TFM/methylation/QC")
getwd()

cat("\n\n... LOADING Methylset.filt OBJECT ...\n\n", file = stdout())
Methylset.filt.name <- load("4.sex/Methylset_filt2.RData")
Methylset.filt <- get(Methylset.filt.name)
cat("\n// Methylset.filt //\n", file = stdout())
Methylset.filt

cat("\n\n... 1) BETA VALUES ...\n\n", file = stdout())

cat("\nCalculating Beta values\n", file = stdout())
betas <- (Methylset.filt@assays$data$Meth / (Methylset.filt@assays$data$Meth + Methylset.filt@assays$data$Unmeth + 100))
betas[1:5,1:5]
cat("\n\nSaving Beta values\n", file = stdout())
save(betas, file="5.getB_M/beta_values.RData")

cat("\n\n... 2) M VALUES ...\n\n", file = stdout())

cat("\nCalculating M values\n", file = stdout())
mvals <- log2((Methylset.filt@assays$data$Meth + 1) / (Methylset.filt@assays$data$Unmeth + 1))
mvals[1:5,1:5]
cat("\n\nSaving M values\n", file = stdout())
save(mvals, file="5.getB_M/m_values.RData")

date()

cat("\n\n ######################## 5.GET BETA AND M VALUES SCRIPT ENDS ######################## \n\n", file = stdout())

