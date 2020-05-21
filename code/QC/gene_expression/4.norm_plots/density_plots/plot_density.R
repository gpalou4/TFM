rm(list=ls())

date()

cat("\n\n ######################## START LOADING LIBRARIES ######################## \n\n", file = stdout())

library(oligo)
cat("\n\n oligo LOADED \n\n", file = stdout())

cat("\n\n ######################## FINISHED LOADING LIBRARIES ######################## \n\n", file = stdout())

cat("\n\n ######################## 4.DENSITY NORMALIZED PLOTS SCRIPT BEGINS ######################## \n\n", file = stdout())

cat("\n\n... SETTING WORKING DIRECTORY ...\n\n", file = stdout())
setwd("/projects/regicor/Guille/TFM/transcriptomics")
getwd()

cat("\n\n... LOADING NORMALIZED ExpressionSet OBJECT ...\n\n", file = stdout())
ESet.name <- load("QC/3.normalization/ESet_norm.RData")
ESet <- get(ESet.name)
cat("\n\n // ExpressionSet //\n", file = stdout())
getClass(ESet)

cat("\n\n... DENSITY PLOTS BY BATCHES ...\n\n", file = stdout())

cat("\nBatches\n", file = stdout())
samples.batches <- sort(unique(ESet@phenoData@data$Batch))
samples.batches

# Use in-house function
source("QC/2.plots/plots_functions.R")

density.outlier.samples <- plotbyBatch(Set = ESet, batches = samples.batches, plot.type = "density",
            directory = "QC/4.norm_plots/density_plots")

cat("\nOutlier samples list by batches\n\n", file = stdout())
names(density.outlier.samples) <- samples.batches
density.outlier.samples

cat("\n\n ... SAVING OUTLIER SAMPLES LIST ...\n\n", file = stdout())

save(density.outlier.samples, file = "QC/4.norm_plots/density_plots/density_outlier_samples.RData")

date()


cat("\n\n ######################## 4.DENSITY NORMALIZED PLOTS SCRIPT ENDS ######################## \n\n", file = stdout())


