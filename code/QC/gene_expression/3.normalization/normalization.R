rm(list=ls())
start_time <- Sys.time()
date()

cat("\n\n ######################## START LOADING LIBRARIES ######################## \n\n", file = stdout())

library(oligo)
cat("\n\n oligo LOADED \n\n", file = stdout())

cat("\n\n ######################## FINISHED LOADING LIBRARIES ######################## \n\n", file = stdout())

cat("\n\n ######################## 3.NORMALIZATION SCRIPT BEGINS ######################## \n\n", file = stdout())

sessionInfo()

cat("\n\n... SETTING WORKING DIRECTORY ...\n\n", file = stdout())
setwd("/projects/regicor/Guille/TFM/transcriptomics")
getwd()

cat("\n\n... LOADING ExonFeatureSet OBJECT ...\n\n", file = stdout())
EFSet.name <- load("QC/1.read_input/ExonFeatureSet.RData")
EFSet <- get(EFSet.name)
cat("\n\n // ExonFeatureSet //\n", file = stdout())
getClass(EFSet)

cat("\n\n... SUBSTRACTION, NORMALIZATION AND SUMMARIZATION BY RMA METHOD ...\n\n", file = stdout())

ESet.norm <- rma(EFSet, background = TRUE, normalize = TRUE, target = "core")

cat("\n\n... CHECKING THE RESULTING GENE EXPRESSION MATRIX ...\n\n", file = stdout())

cat("\n // ExpressionSet //\n", file = stdout())
getClass(ESet.norm)
cat("\n //Gene expression matrix// \n", file = stdout())
assayDataElement(ESet.norm,'exprs')[1:5,1:5]

cat("\n\n ... SAVING NORMALIZED EXPRESSION SET ...\n\n", file = stdout())

save(ESet.norm, file = "QC/3.normalization/ESet_norm.RData")

date()

cat("\nComputational time\n", file = stdout())
end_time <- Sys.time()
print(end_time - start_time)

cat("\n\n ######################## 3.NORMALIZATION SCRIPT ENDS ######################## \n\n", file = stdout())