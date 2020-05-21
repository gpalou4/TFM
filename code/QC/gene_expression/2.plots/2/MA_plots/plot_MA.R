rm(list=ls())
start_time <- Sys.time()
date()

cat("\n\n ######################## START LOADING LIBRARIES ######################## \n\n", file = stdout())

library(oligo)
cat("\n\n oligo LOADED \n\n", file = stdout())

cat("\n\n ######################## FINISHED LOADING LIBRARIES ######################## \n\n", file = stdout())

cat("\n\n ######################## 2.MA PLOTS SCRIPT BEGINS ######################## \n\n", file = stdout())

cat("\n\n... SETTING WORKING DIRECTORY ...\n\n", file = stdout())
setwd("/projects/regicor/Guille/TFM/transcriptomics")
getwd()

cat("\n\n... LOADING ExonFeatureSet OBJECT ...\n\n", file = stdout())
EFSet.name <- load("QC/1.read_input/ExonFeatureSet.RData")
EFSet <- get(EFSet.name)
getClass(EFSet)

cat("\n\n... MA PLOT FOR EACH INDIVIDUAL ...\n\n", file = stdout())

pdf("QC/2.plots/MA_plots/MA_plots.pdf")
MAplot <- MAplot(EFSet)
dev.off()

cat("\nComputational time\n", file = stdout())
end_time <- Sys.time()
print(end_time - start_time)

cat("\n\n ######################## 2.MA PLOTS SCRIPT ENDS ######################## \n\n", file = stdout())


