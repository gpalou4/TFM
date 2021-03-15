rm(list=ls())

date()

cat("\n\n ######################## START LOADING LIBRARIES ######################## \n\n", file = stdout())

library(oligo)
cat("\n\n oligo LOADED \n\n", file = stdout())

cat("\n\n ######################## FINISHED LOADING LIBRARIES ######################## \n\n", file = stdout())

cat("\n\n ######################## 2.PLOTS SCRIPT BEGINS ######################## \n\n", file = stdout())

cat("\n\n... SETTING WORKING DIRECTORY ...\n\n", file = stdout())
setwd("/projects/regicor/Guille/TFM/transcriptomics")
getwd()

cat("\n\n... LOADING ExonFeatureSet OBJECT ...\n\n", file = stdout())
EFSet.name <- load("QC/1.read_input/ExonFeatureSet.RData")
EFSet <- get(EFSet.name)
getClass(EFSet)

# cat("\n\n... PLOTTING ARRAY IMAGES ...\n\n", file = stdout())
# pdf(paste(working.dir,"/2.plots/images_plots.pdf", sep = ""))
# oligo::image(EFSet, which = 1)
# plot(1)
# dev.off()

cat("\n\n... MA PLOT FOR EACH INDIVIDUAL ...\n\n", file = stdout())

pdf("QC/2.plots/MA_plots.pdf")
MAplot(EFSet)
dev.off()

cat("\n\n... BOX PLOTS BY BATCHES ...\n\n", file = stdout())

# CHANGE THIS: WE NEED A VECTOR OF BATCHES
samples.batches <- unique(EFSet@phenoData@data$Batch)

# Use in-house function
source("QC/2.plots/plots_functions.R")

plotbyBatch(ExonFeatureSet = EFSet, batches = samples.batches, plot.type = "boxplot",
            directory = "QC/2.plots")

cat("\n\n... DENSITY PLOTS BY BATCHES ...\n\n", file = stdout())

plotbyBatch(ExonFeatureSet = EFSet, batches = samples.batches, plot.type = "density",
            directory = "QC/2.plots")

cat("\n\n... RELATIVE LOG EXPRESSION PLOTS BY BATCHES ...\n\n", file = stdout())

plotbyBatch(ExonFeatureSet = EFSet, batches = samples.batches, plot.type = "RLE",
            directory = "QC/2.plots")

cat("\n\n... NORMALIZED UNSCALED STANDARD ERRORS PLOTS BY BATCHES ...\n\n", file = stdout())

plotbyBatch(ExonFeatureSet = EFSet, batches = samples.batches, plot.type = "NUSE",
            directory = "QC/2.plots")

cat("\n\n ######################## 2.PLOTS SCRIPT ENDS ######################## \n\n", file = stdout())


