rm(list=ls())

date()

cat("\n\n ######################## START LOADING LIBRARIES ######################## \n\n", file = stdout())

library(oligo)
cat("\n\n oligo LOADED \n\n", file = stdout())

cat("\n\n ######################## FINISHED LOADING LIBRARIES ######################## \n\n", file = stdout())

cat("\n\n ######################## 2.IMAGE PLOTS SCRIPT BEGINS ######################## \n\n", file = stdout())

cat("\n\n... SETTING WORKING DIRECTORY ...\n\n", file = stdout())
setwd("/projects/regicor/Guille/TFM/transcriptomics")
getwd()

cat("\n\n... LOADING ExonFeatureSet OBJECT ...\n\n", file = stdout())
EFSet.name <- load("QC/1.read_input/ExonFeatureSet.RData")
EFSet <- get(EFSet.name)
cat("\n\n // ExonFeatureSet //\n", file = stdout())
getClass(EFSet)

cat("\n\n... PLOTTING ARRAY IMAGES ...\n\n", file = stdout())

pdf("QC/2.plots/image_plots/images_plots.pdf")
oligo::image(EFSet)
plot(1)
dev.off()


cat("\n\n ######################## 2.IMAGE PLOTS SCRIPT ENDS ######################## \n\n", file = stdout())


