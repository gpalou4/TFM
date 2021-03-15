rm(list=ls())
start_time <- Sys.time()
date()

cat("\n\n ######################## START LOADING LIBRARIES ######################## \n\n", file = stdout())

library(oligo)
cat("\n\n oligo LOADED \n\n", file = stdout())

cat("\n\n ######################## FINISHED LOADING LIBRARIES ######################## \n\n", file = stdout())

cat("\n\n ######################## 2.PLOT BOXPLOTS SCRIPT BEGINS ######################## \n\n", file = stdout())

cat("\n\n... SETTING WORKING DIRECTORY ...\n\n", file = stdout())
setwd("/projects/regicor/Guille/TFM/transcriptomics")
getwd()

cat("\n\n... LOADING ExonFeatureSet OBJECT ...\n\n", file = stdout())
EFSet.name <- load("QC/1.read_input/ExonFeatureSet.RData")
EFSet <- get(EFSet.name)
cat("\n\n // ExonFeatureSet //\n", file = stdout())
getClass(EFSet)

cat("\n\n... BOX PLOTS BY BATCHES ...\n\n", file = stdout())

cat("\nBatches\n", file = stdout())
samples.batches <- sort(unique(EFSet@phenoData@data$Batch))
samples.batches

# Use in-house function
source("QC/2.plots/plots_functions.R")

boxplot.outlier.samples <- plotbyBatch(Set = EFSet, batches = samples.batches, plot.type = "boxplot",
                              directory = "QC/2.plots/boxplots")

cat("\nOutlier samples list by batches\n\n", file = stdout())
names(boxplot.outlier.samples) <- samples.batches
boxplot.outlier.samples

cat("\n\n ... SAVING OUTLIER SAMPLES LIST ...\n\n", file = stdout())

save(boxplot.outlier.samples, file = "QC/2.plots/boxplots/boxplot_outlier_samples.RData")

date()

cat("\nComputational time\n", file = stdout())
end_time <- Sys.time()
print(end_time - start_time)

cat("\n\n ######################## 2.PLOT BOXPLOTS SCRIPT ENDS ######################## \n\n", file = stdout())


