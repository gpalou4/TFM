rm(list=ls())
start_time <- Sys.time()
date()

cat("\n\n ######################## START LOADING LIBRARIES ######################## \n\n", file = stdout())

library(oligo)
cat("\n\n oligo LOADED \n\n", file = stdout())

cat("\n\n ######################## FINISHED LOADING LIBRARIES ######################## \n\n", file = stdout())

cat("\n\n ######################## 4.PLOT NORMALIZED BOXPLOTS SCRIPT BEGINS ######################## \n\n", file = stdout())

cat("\n\n... SETTING WORKING DIRECTORY ...\n\n", file = stdout())
setwd("/projects/regicor/Guille/TFM/transcriptomics")
getwd()

cat("\n\n... LOADING NORMALIZED ExpressionSet OBJECT ...\n\n", file = stdout())
ESet.name <- load("QC/3.normalization/ESet_norm.RData")
ESet <- get(ESet.name)
cat("\n\n // ExpressionSet //\n", file = stdout())
getClass(ESet)

cat("\n\n... BOX PLOTS BY BATCHES ...\n\n", file = stdout())

cat("\nBatches\n", file = stdout())
samples.batches <- sort(unique(ESet@phenoData@data$Batch))
samples.batches

# Use in-house function
source("QC/2.plots/plots_functions.R")

boxplot.outlier.samples <- plotbyBatch(Set = ESet, batches = samples.batches, plot.type = "boxplot",
                              directory = "QC/4.norm_plots/boxplots")

cat("\nOutlier samples list by batches\n\n", file = stdout())
names(boxplot.outlier.samples) <- samples.batches
boxplot.outlier.samples

cat("\n\n ... SAVING OUTLIER SAMPLES LIST ...\n\n", file = stdout())

save(boxplot.outlier.samples, file = "QC/4.norm_plots/boxplots/boxplot_outlier_samples.RData")

date()

cat("\nComputational time\n", file = stdout())
end_time <- Sys.time()
print(end_time - start_time)

cat("\n\n ######################## 4.PLOT NORMALIZED BOXPLOTS SCRIPT ENDS ######################## \n\n", file = stdout())


