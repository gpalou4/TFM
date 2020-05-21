rm(list=ls())

date()

cat("\n\n ######################## START LOADING LIBRARIES ######################## \n\n", file = stdout())

library(oligo)
cat("\n\n oligo LOADED \n\n", file = stdout())

cat("\n\n ######################## FINISHED LOADING LIBRARIES ######################## \n\n", file = stdout())

cat("\n\n ######################## 4.MA NORMALIZED PLOTS SCRIPT BEGINS ######################## \n\n", file = stdout())

cat("\n\n... SETTING WORKING DIRECTORY ...\n\n", file = stdout())
setwd("/projects/regicor/Guille/TFM/transcriptomics")
getwd()

cat("\n\n... LOADING NORMALIZED ExpressionSet OBJECT ...\n\n", file = stdout())
ESet.name <- load("QC/3.normalization/ESet_norm.RData")
ESet <- get(ESet.name)
cat("\n\n // ExpressionSet //\n", file = stdout())
getClass(ESet)

cat("\n\n... MA PLOT FOR EACH INDIVIDUAL ...\n\n", file = stdout())

pdf("QC/4.norm_plots/MA_plots/MA_plots.pdf")
oligo::MAplot(ESet)
dev.off()

cat("\n\n ######################## 4.MA NORMALIZED PLOTS SCRIPT ENDS ######################## \n\n", file = stdout())


