rm(list=ls())
start_time <- Sys.time()
date()

cat("\n\n ######################## START LOADING LIBRARIES ######################## \n\n", file = stdout())

library(limma)
cat("\n\n limma LOADED \n\n", file = stdout())
library(minfi)
cat("\n\n mifi LOADED \n\n", file = stdout())
library(IlluminaHumanMethylation450kmanifest)
cat("\n\n IlluminaHumanMethylation450kmanifest LOADED \n\n", file = stdout())

cat("\n\n ######################## FINISHED LOADING LIBRARIES ######################## \n\n", file = stdout())

sessionInfo()

cat("\n\n ######################## 9. BATCH EFFECT CONTROL SCRIPT BEGINS ######################## \n\n", file = stdout())

cat("\n\n... SETTING WORKING DIRECTORY ...\n\n", file = stdout())
setwd("/projects/regicor/Guille/TFM/methylation")
getwd()

cat("\n\n... LOADING SAMPLES METADATA FILE THAT CONTAINS ALL METADATA, I.E, NO NA'S IN CVD  (2080 individuals) ...\n\n", file = stdout())

samples.metadata <- read.csv(file="phenotype/samples_metadata/Analysis/samples_metadata_merged.csv",
                             header=TRUE, stringsAsFactors = FALSE, sep=",")

# 2080 individuals
dim(samples.metadata)
head(samples.metadata)
names(samples.metadata)

cat("\n\n... LOADING M values MATRIX FROM ALL 2596 INDIVIDUALS (NOT ONLY THE 1522 MATCHED WITH MRNA) ...\n\n", file = stdout())
mvals.name <- load("QC/6.normalization/m_norm_4SD.RData")
mvals <- get(mvals.name)
cat("\n// M values //\n", file = stdout())
mvals[1:5,1:5]
dim(mvals)

### BORRAR ###

# subset
#mvals <- mvals[50000:60000,]

### BORRAR ###

# Match M vals matrix sample names with sample names from samples metadata
cat("\n\n... MATCHING M VALS MATRIX SAMPLE NAMES (2596) WITH SAMPLE NAMES FROM SAMPLES METADATA (2080)  ...\n\n", file = stdout())
mvals.filt <- mvals[,colnames(mvals)%in%samples.metadata$LABID]
cat("\n// M values filtered //\n", file = stdout())
mvals.filt[1:5,1:5]
dim(mvals.filt)

cat("\n\n... CVD TABLE BY POSSIBLE BATCH EFFECT VARIABLES ...\n\n", file = stdout())

table(data.frame(CVD=samples.metadata$cvd,Sample_Well=samples.metadata$Sample.Well))
table(data.frame(CVD=samples.metadata$cvd,Plate=samples.metadata$Plate))
table(data.frame(CVD=samples.metadata$cvd,Slide=samples.metadata$Slide))
table(data.frame(CVD=samples.metadata$cvd,Array=samples.metadata$Array))

cat("\n\n... LOADING CPGS METADATA FROM ILLUMINA 450KS MANIFEST ...\n\n", file = stdout())
Illumina.450k.manifest <- read.csv("phenotype/phg000492.v1.FHS_DNAMethylation.marker-info.MULTI/HumanMethylation450_15017482_v.1.1.csv",
                                   skip = 7, stringsAsFactors=F, header=T, sep=",", quote="")
head(Illumina.450k.manifest)
dim(Illumina.450k.manifest)
cat("\nNumber of CpGs located in Chromosome Y\n", file = stdout())
#table(Illumina.450k.manifest$CHR=="X")
cpgs.chrY <- Illumina.450k.manifest[which(Illumina.450k.manifest$CHR == "Y"),"IlmnID"]
length(cpgs.chrY)
cat("\nNumber of CpGs located in Chromosome X\n", file = stdout())
#table(Illumina.450k.manifest$CHR=="Y")
cpgs.chrX <- Illumina.450k.manifest[which(Illumina.450k.manifest$CHR == "X"),"IlmnID"]
length(cpgs.chrX)

###Comprobar que los CpGs del manifiesto estén en el mismo orden que la DB
###Sino estarás borrando los CpGs de un orden del manifiesto que no corresponden
###con los de la matriz

cat("\n\n... REMOVING CPGS LOCATED IN BOTH X AND Y CHROMOSOMES ...\n\n", file = stdout())
# mvals.filt.Y <- mvals.filt[-which(Illumina.450k.manifest$CHR=="Y"),]
# mvals.filt.XY <- mvals.filt.Y[-which(Illumina.450k.manifest$CHR=="X"),]
mvals.filt.Y<- mvals.filt[!rownames(mvals.filt)%in%cpgs.chrY,]
mvals.filt.XY<- mvals.filt.Y[!rownames(mvals.filt.Y)%in%cpgs.chrX,]
cat("\n// M values XY filtered //\n", file = stdout())
mvals.filt.XY[1:5,1:5]
dim(mvals.filt.XY)

cat("\n\n... MDS PLOTS BY POSSIBLE BATCH EFFECT VARIABLES ...\n\n", file = stdout())

# batch <- as.integer(factor(samples.metadata$cvd))
# outcome <- factor(samples.metadata$Plate)

# se_gender <- se.filt4[,!is.na(se.filt4$gender)]
# dge_gender <- dge.filt4[,!is.na(se.filt4$gender)]
# batch <- as.integer(factor(se_gender$gender))
# outcome <- paste(factor(se_gender$type), substr(colnames(dge_gender), 9, 12), sep="-")

cat("\n\n... BATCH 1: ARRAY ...\n\n", file = stdout())

batch <- as.integer(factor(samples.metadata$Array))
outcome <- factor(samples.metadata$cvd)

cat("\n\n... SAVING PLOT ...\n\n", file = stdout())

png("QC/9.batch_effect/batch_array_noXY_FULL_cvd.png", height=3600, width=6000, res=600, units="px")
plotMDS(mvals.filt.XY, col = batch, labels = outcome, cex = 0.7,cex.axis = 1, cex.lab = 1.6, main = "MDS for CVD colored by Array")
legend("bottomleft", paste("Array", sort(unique(batch))), fill=sort(unique(batch)), inset=0.005, cex = 0.6)
dev.off()

cat("\n\n... BATCH 2: PLATE ...\n\n", file = stdout())

batch <- as.integer(factor(samples.metadata$Plate))
outcome <- factor(samples.metadata$cvd)

cat("\n\n... SAVING PLOT ...\n\n", file = stdout())

png("QC/9.batch_effect/batch_plate_noXY_FULL_cvd.png", height=3600, width=6000, res=600, units="px")
plotMDS(mvals.filt.XY, col = batch, labels = outcome, cex = 0.7,cex.axis = 1, cex.lab = 1.6, main = "MDS for CVD colored by Plate")
legend("left", paste("Plate", sort(unique(batch))), fill=sort(unique(batch)), inset=0.005, cex = 0.4)
dev.off()

cat("\n\n... BATCH 3: SAMPLE WELL ...\n\n", file = stdout())

batch <- as.integer(factor(samples.metadata$Sample.Well))
outcome <- factor(samples.metadata$cvd)

cat("\n\n... SAVING PLOT ...\n\n", file = stdout())

png("QC/9.batch_effect/batch_sampleWell_noXY_FULL_cvd.png", height=3600, width=6000, res=600, units="px")
plotMDS(mvals.filt.XY, col = batch, labels = outcome, cex = 0.7,cex.axis = 1, cex.lab = 1.6, main = "MDS for CVD colored by Sample Well")
legend("bottomleft", paste("Sample Well", sort(unique(batch))), fill=sort(unique(batch)), inset=0.005, cex = 0.15)
dev.off()

cat("\n\n... BATCH 4: SLIDE ...\n\n", file = stdout())

batch <- as.integer(factor(samples.metadata$Slide))
outcome <- factor(samples.metadata$cvd)

cat("\n\n... SAVING PLOT ...\n\n", file = stdout())

png("QC/9.batch_effect/batch_slide_noXY_FULL_cvd.png", height=3600, width=6000, res=600, units="px")
plotMDS(mvals.filt.XY, col = batch, labels = outcome, cex = 0.7,cex.axis = 1, cex.lab = 1.6, main = "MDS for CVD colored by Slide")
legend("bottomleft", paste("Slide", sort(unique(batch))), fill=sort(unique(batch)), inset=0.005, cex = 0.1)
dev.off()

cat("\n\n... BATCH 4: GENDER ...\n\n", file = stdout())

batch <- as.integer(factor(samples.metadata$sex))
outcome <- factor(samples.metadata$cvd)

cat("\n\n... SAVING PLOT ...\n\n", file = stdout())

png("QC/9.batch_effect/batch_sex_noXY_FULL_cvd.png", height=3600, width=6000, res=600, units="px")
plotMDS(mvals.filt.XY, col = batch, labels = outcome, cex = 0.7,cex.axis = 1, cex.lab = 1.6, main = "MDS for CVD colored by Gender")
legend("bottomleft", paste("Sex", sort(unique(batch))), fill=sort(unique(batch)), inset=0.005, cex = 1)
dev.off()

date()

cat("\nComputational time\n", file = stdout())
end_time <- Sys.time()
print(end_time - start_time)


cat("\n\n ######################## 9. BATCH EFFECT CONTROL SCRIPT ENDS ######################## \n\n", file = stdout())
