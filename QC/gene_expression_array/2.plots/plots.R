rm(list=ls())

date()

cat("\n\n ######################## START LOADING LIBRARIES ######################## \n\n", file = stdout())

library(oligo)
cat("\n\n oligo LOADED \n\n", file = stdout())

cat("\n\n ######################## FINISHED LOADING LIBRARIES ######################## \n\n", file = stdout())

cat("\n\n ######################## 2.FILTERING SCRIPT BEGINS ######################## \n\n", file = stdout())

cat("\n\n... SETTING WORKING DIRECTORY ...\n\n", file = stdout())
##rawdata.dir <- "/projects/regicor/Guille/TFM/transcriptomics/raw_data/phe000002.v7.FHS_SABRe_project3_OFF_CCS.raw-data-cel.c1"
working.dir <- "/home/guille/Escritorio/Guille/TFM/QC/gene_expression_array"
setwd(working.dir)
getwd()

cat("\n\n... LOADING ExonFeatureSet OBJECT ...\n\n", file = stdout())
EFSet.name <- load(paste(working.dir,"/1.read_input/ExonFeatureSet.RData",sep = ""))
EFSet <- get(EFSet.name)

# cat("\n\n... PLOTTING ARRAY IMAGES ...\n\n", file = stdout())
# pdf(paste(working.dir,"/2.plots/images_plots.pdf", sep = ""))
# oligo::image(EFSet, which = 1)
# plot(1)
# dev.off()

cat("\n\n... MA PLOTS ...\n\n", file = stdout())

pdf(paste(working.dir,"/2.plots/MA_plots.pdf", sep = ""))
MAplot(EFSet)
dev.off()

cat("\n\n... BOX PLOTS BY BATCHES ...\n\n", file = stdout())

# We must do plot by batches
# Obtain a vector of batches
samples.batches <- unique(EFSet@phenoData@data$PLATE)

# Use in-house function
source(paste(working.dir,"/functions/plots.R", sep = ""))
plotbyBatch(ExonFeatureSet = EFSet, batches = samples.batches, plot.type = "boxplot",
             directory = paste(working.dir, "/2.plots", sep = ""))

cat("\n\n... DENSITY PLOTS BY BATCHES ...\n\n", file = stdout())

plotbyBatch(ExonFeatureSet = EFSet, batches = samples.batches, plot.type = "density",
            directory = paste(working.dir, "/2.plots", sep = ""))

cat("\n\n... RELATIVE LOG EXPRESSION PLOTS BY BATCHES ...\n\n", file = stdout())

plotbyBatch(ExonFeatureSet = EFSet, batches = samples.batches, plot.type = "RLE",
            directory = paste(working.dir, "/2.plots", sep = ""))

cat("\n\n... NORMALIZED UNSCALED STANDARD ERRORS PLOTS BY BATCHES ...\n\n", file = stdout())

plotbyBatch(ExonFeatureSet = EFSet, batches = samples.batches, plot.type = "NUSE",
            directory = paste(working.dir, "/2.plots", sep = ""))

cat("\n\n... SUBSTRACTION, NORMALIZATION AND SUMMARIZATION BY RMA METHOD ...\n\n", file = stdout())

EFSet.norm <- rma(EFSet, background = TRUE, normalize = TRUE, target = "core")
# One of the following values:  ’core’, ’full’, ’extended’, ’probeset’.  
# Used only with Gene ST and Exon ST designs. difference?

#The ’rma’method allows for two targets:  ’probeset’ (target=’probeset’) and ’transcript’ 
#(target=’core’,  target=’full’, target=’extended’)

cat("\n\n... CHECKING THE RESULTING GENE EXPRESSION MATRIX ...\n\n", file = stdout())

cat("\n// ExpressionSet //\n", file = stdout())
getClass(EFSet.norm)
cat("\nGene expression matrix\n", file = stdout())
assayDataElement(EFSet.norm,'exprs')[1:5,1:5]

# alternatives
#exprs(cel_files_norm)[1:5,1:5]
#EFSet.norm@assayData$exprs[1:5,1:5]

## checking logCPM ###

library("edgeR")























