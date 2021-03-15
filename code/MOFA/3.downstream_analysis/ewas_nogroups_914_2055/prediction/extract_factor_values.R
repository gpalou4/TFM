rm(list=ls())
start_time <- Sys.time()
date()

cat("\n\n ######################## START LOADING LIBRARIES ######################## \n\n", file = stdout())

library(ggplot2)
cat("\n\n ggplot2 LOADED \n\n", file = stdout())
library(MOFA2)
cat("\n\n MOFA LOADED \n\n", file = stdout())

cat("\n\n ######################## FINISHED LOADING LIBRARIES ######################## \n\n", file = stdout())

cat("\n\n ######################## 3. MOFA DOWNSTREAM ANALYSIS: EXTRACT FACTOR VALUES SCRIPT BEGINS ######################## \n\n", file = stdout())

cat("\n\n... SETTING WORKING DIRECTORY ...\n\n", file = stdout())
setwd("/projects/regicor/Guille/TFM")
getwd()

cat("\n\n... 1) LOADING DATA ...\n\n", file = stdout())

samples <- "914_2055"
groups <- "nogroups"

cat("\n\n... 1.1) LOADING SAMPLES METADATA ...\n\n", file = stdout())

MOFA.covariates <- read.csv(file = paste("MOFA/MOFA_covariates_order_",groups,"_",samples,".csv", sep = ""),
                            header=TRUE, stringsAsFactors = F, sep=",")
# Metadata has to contain the columns 'sample' (shareid) and 'group' (cvd) for MOFA to work
colnames(MOFA.covariates)[1] <- "sample"
# Change sample names from integers to characters
MOFA.covariates$sample <- as.character(MOFA.covariates$sample)
MOFA.covariates[1:3,]
dim(MOFA.covariates)
# 2080 24

samples <- "914_2055"
groups <- "ewas_nogroups"

cat("\n\n... 1.2) LOADING MOFA TRAINED OBJECT ...\n\n", file = stdout())

MOFA.model <- load_model(paste("MOFA/2.training_model/MOFA_trained_",groups,"_",samples,".hdf5", sep = ""))
MOFA.model

cat("\n\n... 2) ADDING SAMPLES METADATA (COVARIABLES) TO THE TRAINED MODEL OBJECT ...\n\n", file = stdout())

samples_metadata(MOFA.model) <- MOFA.covariates[-2]

cat("\n\n... 3) INSPECTING THE MOFA OBJECT ...\n\n", file = stdout())

cat("\nSlot names of the MOFA model\n", file = stdout())
slotNames(MOFA.model)
cat("\nMetadata stored in the MOFA model\n", file = stdout())
head(MOFA.model@samples_metadata)

cat("\n\n... 4) EXTRACT AND SAVE ALL FACTOR VALUES ...\n\n", file = stdout())

factors <- get_factors(MOFA.model,
                       groups = "all",
                       factors = "all"
)
lapply(factors,dim)
lapply(factors,head)

# cat("\n\nInterested factor\n\n", file = stdout())
# factor <- 2
# factor.values <- factors$group1[,factor]
# head(factor.values)

save(factors, file = paste("MOFA/3.downstream_analysis/",groups,"_",samples,"/prediction/factors_values.RData", sep = ""))

# cat("\n\n... 5) EXTRACT AND SAVE ALL TOP 20 HEATMAP FEATURES EXPRESSION VALUES ...\n\n", file = stdout())
# 
# meth.30features.heatmap.name <- load(file = paste("MOFA/3.downstream_analysis/",groups,"_",samples,"/variance_decomposition/meth_30features_heatmap.RData", sep = ""))
# meth.30features.heatmap <- get(meth.30features.heatmap.name)
# cat("\n// Top 20 features from Factor 12 heatmap for Methylation //\n", file = stdout())
# meth.30features.heatmap[[12]]
# 
# cat("\nM matrix subset for that specific features\n", file = stdout())
# MOFA.data <- get_data(MOFA.model)
# M.matrix <- t(MOFA.data$meth$group1[rownames(MOFA.data$meth$group1)%in%meth.30features.heatmap[[12]],])
# dim(M.matrix)
# M.matrix[1:5,1:5]
# save(M.matrix, file = paste("MOFA/3.downstream_analysis/",groups,"_",samples,"/prediction/features_expression.RData", sep = ""))
# 

## TOP 20 FEATURES WEIGHTS FOR A FACTOR
# features.weights <- get_weights(MOFA.model)
# 
# top20.features.weights <- features.weights$meth[meth.30features.heatmap[[12]],]
# top20.features.weights[1:5,1:5]
# 
# save(top20.features.weights, file = paste("MOFA/3.downstream_analysis/",groups,"_",samples,"/prediction/features_weights.RData", sep = ""))

date()

cat("\n\nComputational time\n\n", file = stdout())
end_time <- Sys.time()
end_time - start_time

cat("\n\n ######################## 3. MOFA DOWNSTREAM ANALYSIS: EXTRACT FACTOR VALUES SCRIPT ENDS ######################## \n\n", file = stdout())


