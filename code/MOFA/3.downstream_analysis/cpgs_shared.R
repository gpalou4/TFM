rm(list=ls())
start_time <- Sys.time()
date()

cat("\n\n ######################## START LOADING LIBRARIES ######################## \n\n", file = stdout())

library(ggplot2)
cat("\n\n ggplot2 LOADED \n\n", file = stdout())
library(MOFA2)
cat("\n\n MOFA LOADED \n\n", file = stdout())


cat("\n\n ######################## FINISHED LOADING LIBRARIES ######################## \n\n", file = stdout())

cat("\n\n ######################## 3. MOFA DOWNSTREAM ANALYSIS: VARIANCE DECOMPOSITION SCRIPT BEGINS ######################## \n\n", file = stdout())

cat("\n\n... SETTING WORKING DIRECTORY ...\n\n", file = stdout())
setwd("/projects/regicor/Guille/TFM")
getwd()

cat("\n\n... LOADING MOFA TRAINED OBJECT ...\n\n", file = stdout())

groups <- "nogroups"
samples <- "914_2055"

MOFA.model.sd <- load_model(paste("MOFA/2.training_model/MOFA_trained_",groups,"_",samples,".hdf5",sep = ""))
MOFA.model.sd

groups <- "ewas_nogroups"
samples <- "914_2055"

MOFA.model.ewas <- load_model(paste("MOFA/2.training_model/MOFA_trained_",groups,"_",samples,".hdf5",sep = ""))
MOFA.model.ewas


cat("\n\n... SHARED CPGS ...\n\n", file = stdout())

table(rownames(MOFA.model.sd@data$meth$group1)%in%rownames(MOFA.model.ewas@data$meth$group1))













