rm(list=ls())
start_time <- Sys.time()
date()

cat("\n\n ######################## START LOADING LIBRARIES ######################## \n\n", file = stdout())

library(limma)
cat("\n\n limma LOADED \n\n", file = stdout())
library(oligo)
cat("\n\n oligo LOADED \n\n", file = stdout())
library(sva)
cat("\n\n sva LOADED \n\n", file = stdout())

cat("\n\n ######################## FINISHED LOADING LIBRARIES ######################## \n\n", file = stdout())

sessionInfo()

cat("\n\n ######################## 6. BATCH EFFECT CONTROL SCRIPT BEGINS ######################## \n\n", file = stdout())

cat("\n\n... SETTING WORKING DIRECTORY ...\n\n", file = stdout())
setwd("/projects/regicor/Guille/TFM/transcriptomics")
getwd()

cat("\n\n... LOADING SAMPLES METADATA FILE THAT CONTAINS METADATA FOR THE MATCHED 1522 INDIVIDUALS...\n\n", file = stdout())

samples.metadata <- read.csv(file="phenotype/Analysis/samples_metadata_merged.csv",
                             header=TRUE, stringsAsFactors = FALSE, sep=",")
rownames(samples.metadata) <- substr(samples.metadata$cel_files,1,nchar(samples.metadata$cel_files)-4)

# 1522 individuals
dim(samples.metadata)
head(samples.metadata)
names(samples.metadata)

cat("\n\n... LOADING SAMPLES METADATA FILE THAT CONTAINS METADATA FOR ALL THE INDIVIDUALS (1818) ...\n\n", file = stdout())
mrna.samples.metadata.1818 <- read.csv(file = "phenotype/Analysis/samples_metadata_1818.csv",
                                       header=TRUE, stringsAsFactors = F, sep=",")
rownames(mrna.samples.metadata.1818) <- substr(mrna.samples.metadata.1818$cel_files,1,nchar(mrna.samples.metadata.1818$cel_files)-4)
# 1818 individuals
dim(mrna.samples.metadata.1818)
head(mrna.samples.metadata.1818)
names(mrna.samples.metadata.1818)

cat("\n\n... LOADING NORMALIZED AND FILTERED ExpressionSet OBJECT ...\n\n", file = stdout())
ESet.name <- load("QC/5.filtering/ESet_norm_filt.RData")
ESet <- get(ESet.name)
cat("\n\n // ExpressionSet //\n", file = stdout())
getClass(ESet)

# samples metadata subset
cat("\n\n... MATCHING SAMPLE NAMES FROM SAMPLES METADATA (1522) WITH ExpressionSet SAMPLE NAMES (1110) ...\n\n", file = stdout())
samples.metadata <- samples.metadata[rownames(samples.metadata)%in%colnames(ESet),]
cat("\n// Samples metadata subset //\n", file = stdout())
dim(samples.metadata) # 938 30

cat("\n\n... MATCHING SAMPLE NAMES FROM SAMPLES METADATA (1818) WITH ExpressionSet SAMPLE NAMES (1110) ...\n\n", file = stdout())
mrna.samples.metadata.1818 <- mrna.samples.metadata.1818[rownames(mrna.samples.metadata.1818)%in%colnames(ESet),]
cat("\n// Samples metadata subset //\n", file = stdout())
dim(mrna.samples.metadata.1818) # 1110 9

# ExpressionSet subset
cat("\n\n... MATCHING ExpressionSet SAMPLE NAMES (1110) WITH THE SAMPLES METADATA SUBSET (938)  ...\n\n", file = stdout())
ESet.subset <- ESet[,colnames(ESet)%in%rownames(samples.metadata)]
cat("\n\n // ExpressionSet subset//\n", file = stdout())
getClass(ESet.subset) # 938

cat("\nOrdering samples from ExpressionSet matched (938) and samples from samples metadata (938) \n", file = stdout())
sample.ids <- match(colnames(ESet.subset), rownames(samples.metadata))
samples.metadata <- samples.metadata[sample.ids,]
cat("\nCheck\n", file = stdout())
cat("\nSamples metadata\n", file = stdout())
rownames(samples.metadata)[1:20]
head(samples.metadata)
dim(samples.metadata)
cat("\nExpressionSet\n", file = stdout())
colnames(ESet.subset)[1:20]
dim(ESet.subset)

cat("\nOrdering samples from ExpressionSet and samples from samples metadata (1110)\n", file = stdout())
sample.ids <- match(colnames(ESet), rownames(mrna.samples.metadata.1818))
mrna.samples.metadata.1818 <- mrna.samples.metadata.1818[sample.ids,]
cat("\nCheck\n", file = stdout())
cat("\nSamples metadata\n", file = stdout())
rownames(mrna.samples.metadata.1818)[1:20]
head(mrna.samples.metadata.1818)
dim(mrna.samples.metadata.1818)
cat("\nExpressionSet\n", file = stdout())
colnames(ESet)[1:20]
dim(ESet)

cat("\n\n... CVD TABLE BY POSSIBLE BATCH EFFECT VARIABLES ...\n\n", file = stdout())

table(data.frame(CVD=samples.metadata$cvd,Batch=samples.metadata$Batch))
table(data.frame(CVD=samples.metadata$cvd,Row=samples.metadata$Row))
table(data.frame(CVD=samples.metadata$cvd,Column=samples.metadata$Column))
#table(data.frame(CVD=samples.metadata$cvd,Plate=samples.metadata$Plate_Location))

## QUITAR DATOS DE EXPRESION DE CHR X E Y?

cat("\n\n... 1) MDS PLOTS BY POSSIBLE BATCH EFFECT VARIABLES ...\n\n", file = stdout())

colfunc <- colorRampPalette(c("red","yellow","springgreen","royalblue"))
colfunc(20)

cat("\n\n... BATCH 1: Batch for 938 ...\n\n", file = stdout())

# batch <- as.numeric(as.character(factor(samples.metadata$Batch)))
# cat("\nBatches\n", file = stdout())
# sort(unique(batch))
# length(sort(unique(batch)))
# outcome <- factor(samples.metadata$cvd)
# batch.colors <- c()
# for (i in 1:length(batch)) {
#   #batch[i] <- colfunc(20)[batch[i]]
#   color <- colfunc(20)[batch[i]-6]
#   batch.colors <- c(batch.colors,color)
#   }
# 
# cat("\n\n... SAVING PLOT ...\n\n", file = stdout())
# 
# # layout(matrix(1:2,ncol=2), width = c(2,1),height = c(1,1))
# png("QC/6.batch_effect/batch_batch_FULL_cvd_938.png", height=3600, width=6000, res=600, units="px")
# plotMDS(exprs(ESet.subset), col = batch.colors, labels = outcome, cex = 0.7,cex.axis = 1, cex.lab = 1.6, main = "MDS for CVD colored by Batch")
# legend("bottomleft", paste("Batch", sort(unique(batch))), fill=colfunc(20), inset=0.005, cex = 0.6)
# #legend("bottomleft", paste("Batch", sort(unique(batch))), fill=sort(unique(batch)), inset=0.005, cex = 0.6)
# 
# # legend_image <- as.raster(matrix(colfunc(25), ncol=1))
# # plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'legend title')
# # text(x=1.5, y = seq(0,1,l=2), labels = seq(7,30,l=2))
# # rasterImage(legend_image, 0, 0, 1,1)
# legend_image <- as.raster(matrix(colfunc(25), ncol=1))
# text(x=1.30, y = seq(-1.40,-1.05,l=2), labels = seq(7,30,l=2))
# rasterImage(legend_image, 1.25, -1.45, 1.10, -1)
# 
# dev.off()

cat("\n\n... BATCH 1: Batch for 1110 ...\n\n", file = stdout())

dim(mrna.samples.metadata.1818)
batch <- as.numeric(as.character(factor(mrna.samples.metadata.1818$Batch)))
cat("\nBatches\n", file = stdout())
sort(unique(batch))
length(sort(unique(batch)))
outcome <- factor(mrna.samples.metadata.1818$cvd)
batch.colors <- c()
for (i in 1:length(batch)) {
  #batch[i] <- colfunc(20)[batch[i]]
  color <- colfunc(20)[batch[i]-6]
  batch.colors <- c(batch.colors,color)
}

cat("\n\n... SAVING PLOT ...\n\n", file = stdout())

dim(exprs(ESet))
length(batch.colors)

png("QC/6.batch_effect/batch_batch_FULL_cvd_1110.png", height=3600, width=6000, res=600, units="px")
plotMDS(exprs(ESet), col = batch.colors, labels = outcome, cex = 0.7,cex.axis = 1, cex.lab = 1.6, main = "MDS for CVD colored by Batch")
legend("bottomleft", paste("Batch", sort(unique(batch))), fill=colfunc(20), inset=0.005, cex = 0.6)

legend_image <- as.raster(matrix(colfunc(25), ncol=1))
text(x=1.30, y = seq(-1.40,-1.05,l=2), labels = seq(7,30,l=2))
rasterImage(legend_image, 1.25, -1.45, 1.10, -1)

dev.off()

# cat("\n\n... BATCH 2: ROW ...\n\n", file = stdout())
# 
# batch <- as.integer(factor(samples.metadata$Row))
# outcome <- factor(samples.metadata$cvd)
# 
# cat("\n\n... SAVING PLOT ...\n\n", file = stdout())
# 
# png("QC/6.batch_effect/batch_row_FULL_cvd.png", height=3600, width=6000, res=600, units="px")
# plotMDS(exprs(ESet), col = batch.colors, labels = outcome, cex = 0.7,cex.axis = 1, cex.lab = 1.6, main = "MDS for CVD colored by Row")
# legend("bottomleft", paste("Row", sort(unique(batch))), fill=colfunc(20), inset=0.005, cex = 0.6)
# dev.off()
# 
# cat("\n\n... BATCH 3: Column ...\n\n", file = stdout())
# 
# batch <- as.integer(factor(samples.metadata$Column))
# outcome <- factor(samples.metadata$cvd)
# 
# cat("\n\n... SAVING PLOT ...\n\n", file = stdout())
# 
# png("QC/6.batch_effect/batch_column_FULL_cvd.png", height=3600, width=6000, res=600, units="px")
# plotMDS(exprs(ESet), col = batch.colors, labels = outcome, cex = 0.7,cex.axis = 1, cex.lab = 1.6, main = "MDS for CVD colored by Column")
# legend("bottomleft", paste("Column", sort(unique(batch))), fill=colfunc(20), inset=0.005, cex = 0.6)
# dev.off()
# 
# cat("\n\n... BATCH 4: GENDER ...\n\n", file = stdout())
# 
# batch <- as.integer(factor(samples.metadata$sex))
# outcome <- factor(samples.metadata$cvd)
# 
# cat("\n\n... SAVING PLOT ...\n\n", file = stdout())
# 
# png("QC/6.batch_effect/batch_sex_FULL_cvd.png", height=3600, width=6000, res=600, units="px")
# plotMDS(exprs(ESet), col = batch.colors, labels = outcome, cex = 0.7,cex.axis = 1, cex.lab = 1.6, main = "MDS for CVD colored by Gender")
# legend("bottomleft", paste("sex", sort(unique(batch))), fill=colfunc(20), inset=0.005, cex = 0.6)
# dev.off()

cat("\n\n... 2) BATCH EFFECT REMOVAL USING COMBAT ...\n\n", file = stdout())

# For subset (938)
ESet.subset.batch <- ComBat(exprs(ESet.subset), batch = factor(samples.metadata$Batch), par.prior = TRUE)
# For all (1110)
ESet.batch <- ComBat(exprs(ESet), batch = factor(mrna.samples.metadata.1818$Batch), par.prior = TRUE)
cat("\n\n // ExpressionSet before batch effect removal //\n\n", file = stdout())
exprs(ESet.subset)[1:5,1:5]
dim(exprs(ESet.subset))
cat("\n\n // ExpressionSet after batch effect removal //\n\n", file = stdout())
ESet.subset.batch[1:5,1:5]
dim(ESet.subset.batch)

cat("\n\n... SAVING NEW EXPRESSION SET OBJECT ...\n\n", file = stdout())

save(ESet.subset.batch, file = "QC/6.batch_effect/ESet_norm_filt_batch_938.RData")
save(ESet.batch, file = "QC/6.batch_effect/ESet_norm_filt_batch_1110.RData")

cat("\n\n... SAVING NEW BATCH EFFECT PLOT for 938 ...\n\n", file = stdout())

# batch <- as.numeric(as.character(factor(samples.metadata$Batch)))
# cat("\nBatches\n", file = stdout())
# sort(unique(batch))
# length(sort(unique(batch)))
# outcome <- factor(samples.metadata$cvd)
# batch.colors <- c()
# for (i in 1:length(batch)) {
#   #batch[i] <- colfunc(20)[batch[i]]
#   color <- colfunc(20)[batch[i]-6]
#   batch.colors <- c(batch.colors,color)
# }
# 
# png("QC/6.batch_effect/batch_batch_FULL_cvd_938_combat.png", height=3600, width=6000, res=600, units="px")
# plotMDS(ESet.subset.batch, col = batch.colors, labels = outcome, cex = 0.7,cex.axis = 1, cex.lab = 1.6, main = "MDS for CVD colored by Batch")
# legend("bottomleft", paste("Batch", sort(unique(batch))), fill=colfunc(20), inset=0.005, cex = 0.6)
# 
# legend_image <- as.raster(matrix(colfunc(25), ncol=1))
# text(x=1.30, y = seq(-0.9,-0.55,l=2), labels = seq(7,30,l=2))
# rasterImage(legend_image, 1.25, -0.95, 1.10, -0.5)
# dev.off()

cat("\n\n... SAVING NEW BATCH EFFECT PLOT for 1110 ...\n\n", file = stdout())

batch <- as.numeric(as.character(factor(mrna.samples.metadata.1818$Batch)))
cat("\nBatches\n", file = stdout())
sort(unique(batch))
length(sort(unique(batch)))
outcome <- factor(mrna.samples.metadata.1818$cvd)
batch.colors <- c()
for (i in 1:length(batch)) {
  #batch[i] <- colfunc(20)[batch[i]]
  color <- colfunc(20)[batch[i]-6]
  batch.colors <- c(batch.colors,color)
}

cat("\n\n... SAVING PLOT ...\n\n", file = stdout())

png("QC/6.batch_effect/batch_batch_FULL_cvd_1110_combat.png", height=3600, width=6000, res=600, units="px")
plotMDS(exprs(ESet.batch), col = batch.colors, labels = outcome, cex = 0.7,cex.axis = 1, cex.lab = 1.6, main = "MDS for CVD colored by Batch")
legend("bottomleft", paste("Batch", sort(unique(batch))), fill=colfunc(20), inset=0.005, cex = 0.6)

legend_image <- as.raster(matrix(colfunc(25), ncol=1))
text(x=1.30, y = seq(-0.9,-0.55,l=2), labels = seq(7,30,l=2))
rasterImage(legend_image, 1.25, -0.95, 1.10, -0.5)

dev.off()

date()

cat("\nComputational time\n", file = stdout())
end_time <- Sys.time()
print(end_time - start_time)


cat("\n\n ######################## 6. BATCH EFFECT CONTROL SCRIPT ENDS ######################## \n\n", file = stdout())
