rm(list=ls())
start_time <- Sys.time()
date()

cat("\n\n ######################## START LOADING LIBRARIES ######################## \n\n", file = stdout())

library(MOFA2)
cat("\n\n MOFA LOADED \n\n", file = stdout())
library(Biobase)
cat("\n\n Biobase LOADED \n\n", file = stdout())
library(MultiAssayExperiment)
cat("\n\n MultiAssayExperiment LOADED \n\n", file = stdout())

sessionInfo()

cat("\n\n ######################## FINISHED LOADING LIBRARIES ######################## \n\n", file = stdout())

cat("\n\n ######################## 1.MOFA OBJECT SCRIPT BEGINS ######################## \n\n", file = stdout())

cat("\n\n... SETTING WORKING DIRECTORY ...\n\n", file = stdout())
setwd("/projects/regicor/Guille/TFM")
getwd()

groups <- "nogroups"
samples <- "914_2055"

cat("\n\n... 1) LOADING DATA ...\n\n", file = stdout())

cat("\n\n... 1.1) LOADING GENE EXPRESSION DATA: ExpressionSet (its matrix) ...\n\n", file = stdout())
ESet.name <- load("transcriptomics/QC/6.batch_effect/ESet_norm_filt_batch_1110.RData")
ESet <- get(ESet.name)
cat("\n\n // ExpressionSet //\n", file = stdout())
ESet[1:5,1:5]
dim(ESet)

cat("\n\n... 1.2) LOADING METHYLATION DATA: M matrix values ...\n\n", file = stdout())
M.matrix.name <- load("methylation/QC/6.normalization/m_norm_4SD.RData")
M.matrix <- get(M.matrix.name)
cat("\n\n // M values //\n", file = stdout())
dim(M.matrix)
M.matrix[1:5,1:5]

cat("\n\n... 2) LOADING SAMPLES METADATA ...\n\n", file = stdout())

cat("\n\n... 2.1) LOADING METHYLATION SAMPLES METADATA THAT CONTAINS SVA (2080) ...\n\n", file = stdout())
cpgs.samples.metadata.sva <- read.csv(file = "methylation/phenotype/samples_metadata/Analysis/samples_metadata_merged_sva.csv",
                                  header=TRUE, stringsAsFactors = F, sep=",")
head(cpgs.samples.metadata.sva)
dim(cpgs.samples.metadata.sva)
# 2080 30

cat("\n\n... 2.2) LOADING mRNA SAMPLES METADATA ...\n\n", file = stdout())
mrna.samples.metadata <- read.csv(file = "transcriptomics/phenotype/Analysis/samples_metadata_merged.csv",
                                      header=TRUE, stringsAsFactors = F, sep=",")
head(mrna.samples.metadata)
dim(mrna.samples.metadata)
# 1522 30

# # Samples metadata file with batch
# mrna.samples.metadata.batch <- read.csv(file = "transcriptomics/phe000002.v7.FHS_SABRe_project3.sample-info.MULTI/sample-info/NHLBI_FHS_DCC_files/MasterFile_Gene_OFF_2446.txt",
#                                    header=TRUE, stringsAsFactors = F, sep="\t")
# 
# # Filter only these fields: cel_files Batch Plate Row Column Plate_location
# mrna.samples.metadata.batch <- mrna.samples.metadata.batch[,c(1,3,5,6,7)]
# dim(mrna.samples.metadata.batch)
# # 2446 5 (OFFs + CCS)
# 
# # Samples metadata file with shareid
# mrna.samples.metadata.IDs <- read.csv(file = "transcriptomics/phenotype/phs000363.v17.pht002944.v4.p11.c1.Framingham_SABRe_Sample_Attributes.HMB-IRB-MDS_clean.txt",
#                                  header=TRUE, stringsAsFactors = F, sep="\t")
# 
# # Check if raw data files match these metadata file
# rawdata.dir <- "/projects/regicor/Guille/TFM/transcriptomics/raw_data/final_raw_data/"
# raw.data.file.names <- substr(list.files(rawdata.dir),1,27)
# length(raw.data.file.names)
# # 1893 OFFs
# table(mrna.samples.metadata.IDs$sampid%in%raw.data.file.names)
# # FALSE: 9058, TRUE: 1818
# # Remove samples from the metadata that do not match with raw data file names
# mrna.samples.metadata.IDs.filt <- mrna.samples.metadata.IDs[which(mrna.samples.metadata.IDs$sampid%in%raw.data.file.names),]
# dim(mrna.samples.metadata.IDs.filt)
# # 1818 5
# head(mrna.samples.metadata.IDs.filt)
# colnames(mrna.samples.metadata.IDs.filt)[3] <- "cel_files"
# 
# # Merge both dataframes based on the same named column of IDs (cel_files)
# mrna.merged.samples.metadata <- merge(mrna.samples.metadata.batch, mrna.samples.metadata.IDs.filt, by = "cel_files")
# head(mrna.merged.samples.metadata)
# dim(mrna.merged.samples.metadata)
# # 1818 9
# 
# # write.csv(mrna.merged.samples.metadata, file = "transcriptomics/phenotype/Analysis/samples_metadata_1768.csv",
# #           row.names = FALSE)

cat("\n\n... 3) ADDING MRNA SAMPLES METADATA TO METHYLATION SAMPLES METADATA  ...\n\n", file = stdout())

names(cpgs.samples.metadata.sva) <- c("ID_METH","shareid", "sex", "age", "tot_chol", "hdl_chol", "ldl_chol", "trig", "sbp", "dbp", "glucose", "weight", "height", "waist_u", "smoke", "CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran", "cvd", "prev_ex8_cvd_date", "SAMPID", "Sample.Well_METH", "Plate_METH", "Image.date_METH", "Slide_METH", "Array_METH", "sva1_METH")
#names(mrna.merged.samples.metadata) <- c("ID_MRNA", "Batch_MRNA", "Row_MRNA", "Column_MRNA", "Plate_location_MRNA", "dbGaP_Sample_ID", "shareid", "BODY_SITE", "ANALYTE_TYPE")
names(mrna.samples.metadata) <- c("shareid", "ID_MRNA", "Batch_MRNA", "Row_MRNA", "Column_MRNA", "Plate_location_MRNA", "dbGaP_Sample_ID", "BODY_SITE", "ANALYTE_TYPE", "sex", "age", "tot_chol", "hdl_chol", "ldl_chol", "trig", "sbp", "dbp", "glucose", "weight", "height", "waist_u", "smoke", "CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran", "cvd", "prev_ex8_cvd_date")

# Merge mRNA samples metadata and methylation samples metadata (containing Batch NA's in meth samples with no mRNA data, obviously)
common.columns <- c("shareid", "sex", "age", "tot_chol", "hdl_chol", "ldl_chol", "trig", "sbp", "dbp", "glucose", "weight", "height", "waist_u", "smoke", "CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran", "cvd", "prev_ex8_cvd_date")
MOFA.samples.metadata <- merge(mrna.samples.metadata, cpgs.samples.metadata.sva,
                                 by = common.columns, all = TRUE)
dim(MOFA.samples.metadata)
# 2080 38
head(MOFA.samples.metadata)

cat("\n\n... 4) REMOVE COLUMNS THAT ARE NOT COVARIATES FOR THE MOFA MODEL (KEEP THE OUTCOME OF INTEREST: 'CVD') ...\n\n", file = stdout())

covariates <- c("shareid","cvd","sex","age","tot_chol","hdl_chol","trig","sbp","dbp","glucose","weight","height","waist_u","smoke","CD8T","CD4T","NK","Bcell","Mono","Gran","Array_METH","sva1_METH","Batch_MRNA")
MOFA.covariates <- MOFA.samples.metadata[,covariates]
head(MOFA.covariates)
dim(MOFA.covariates)
# 2080 23

cat("\n\n... SAVING NEW SAMPLES METADATA MERGED AND COVARIATES THAT WILL BE USED ON MOFA ...\n\n", file = stdout())

write.csv(MOFA.samples.metadata, file = paste("MOFA/MOFA_samples_metadata_",groups,"_",samples,".csv", sep = ""),
          row.names = FALSE)
write.csv(MOFA.covariates, file = paste("MOFA/MOFA_covariates_",groups,"_",samples,".csv", sep = ""),
          row.names = FALSE)

cat("\n\n... 5) SUBSETTING SAMPLES METADATA TO MATCH THE SAME INDIVIDUALS OF THE DATA MATRICES ...\n\n", file = stdout())

cat("\n\n... 5.1) SUBSETTING METHYLATION SAMPLES WITH THE NEW SAMPLES METADATA  ...\n\n", file = stdout())

# NOTE: here we use MOFA.samples.metadata (and not MOFA.covariates) because it contains the ID's for both METH and mRNA
# necessary to do the subsetting

M.matrix.subset <- M.matrix[,colnames(M.matrix)%in%MOFA.samples.metadata$ID_METH]
cat("\n\n // M values matrix subset //\n", file = stdout())
dim(M.matrix.subset)
M.matrix.subset[1:5,1:5]

cat("\nChanging sample names by their corresponding shareid\n", file = stdout())
sample.names.index <- match(colnames(M.matrix.subset), MOFA.samples.metadata$ID_METH)
sample.names.shareid <- MOFA.samples.metadata$shareid[sample.names.index]
colnames(M.matrix.subset) <- sample.names.shareid

cat("\n\n // M values matrix subset, names changed to IDs //\n", file = stdout())
dim(M.matrix.subset)
M.matrix.subset[1:5,1:5]

cat("\nRemove samples from batch 15, because they seem to be outliers in TSNE plots\n", file = stdout())
tsne.batch.outliers <- na.omit(MOFA.samples.metadata[MOFA.samples.metadata$Batch_MRNA == 15, "shareid"])
tsne.batch.outliers <- tsne.batch.outliers[!is.na(tsne.batch.outliers)]
tsne.batch.outliers
length(tsne.batch.outliers)
# 24

M.matrix.subset <- M.matrix.subset[,!colnames(M.matrix.subset)%in%tsne.batch.outliers]
cat("\n\n // M values matrix subset batch 15 filtered //\n", file = stdout())
M.matrix.subset[1:5,1:5]
dim(M.matrix.subset)
# 420764 914

cat("\n\n... 5.2) SUBSETTING GENE EXPRESSION SAMPLES  ...\n\n", file = stdout())

ESet.subset <- ESet[,colnames(ESet)%in%substr(MOFA.samples.metadata$ID_MRNA,1,nchar(MOFA.samples.metadata$ID_MRNA)-4)]
cat("\n\n // ExpressionSet subset //\n", file = stdout())
ESet.subset[1:5,1:5]
dim(ESet.subset)

cat("\nChanging sample names by their corresponding shareid\n", file = stdout())

sample.names.index <- match(colnames(ESet.subset), substr(MOFA.samples.metadata$ID_MRNA,1,nchar(MOFA.samples.metadata$ID_MRNA)-4))
sample.names.shareid <- MOFA.samples.metadata$shareid[sample.names.index]
colnames(ESet.subset) <- sample.names.shareid
table(is.na(colnames(ESet.subset)))

cat("\n\n // ExpressionSet subset, names changed to IDs //\n", file = stdout())
ESet.subset[1:5,1:5]
dim(ESet.subset)

# cat("\nConfirm both datasets contain the same samples\n", file = stdout())
# table(colnames(ESet)%in%colnames(M.matrix.subset))
# 
# cat("\nRemove samples from batch 15, because they seem to be outliers in TSNE plots\n", file = stdout())
# tsne.batch.outliers <- na.omit(MOFA.samples.metadata[MOFA.samples.metadata$Batch_MRNA == 15, "shareid"])
# tsne.batch.outliers <- tsne.batch.outliers[!is.na(tsne.batch.outliers)]
# tsne.batch.outliers
# length(tsne.batch.outliers)
# # 24

ESet.subset <- ESet.subset[,!colnames(ESet.subset)%in%tsne.batch.outliers]
cat("\n\n // ExpressionSet subset batch 15 filtered //\n", file = stdout())
ESet.subset[1:5,1:5]
dim(ESet.subset)
# 20567 914

cat("\nAdding shareid to rownames in Covariates and ordering\n", file = stdout())

# Very important as.character: sample names must be of the same type of class in both samples metadata and data matrices.
rownames(MOFA.covariates) <- as.character(MOFA.samples.metadata$shareid)
MOFA.covariates.order <- MOFA.covariates[order(rownames(MOFA.covariates)),]
head(MOFA.covariates.order)
dim(MOFA.covariates.order)
# 2080 23

cat("\n\n... 6) SUBSETTING DATA VALUES  ...\n\n", file = stdout())

cat("\n\n... 6.1) REMOVING CpGs FROM CHR X AND Y ...\n\n", file = stdout())

cat("\n\n... LOADING CpGs METADATA FROM ILLUMINA 450KS MANIFEST ...\n\n", file = stdout())
Illumina.450k.manifest <- read.csv("methylation/phenotype/phg000492.v1.FHS_DNAMethylation.marker-info.MULTI/HumanMethylation450_15017482_v.1.1.csv",
                                   skip = 7, stringsAsFactors=F, header=T, sep=",", quote="")

head(Illumina.450k.manifest)
dim(Illumina.450k.manifest)
cat("\nNumber of CpGs located in Chromosome Y\n", file = stdout())
cpgs.chrY <- Illumina.450k.manifest[which(Illumina.450k.manifest$CHR == "Y"),"IlmnID"]
length(cpgs.chrY)
cat("\nNumber of CpGs located in Chromosome X\n", file = stdout())
cpgs.chrX <- Illumina.450k.manifest[which(Illumina.450k.manifest$CHR == "X"),"IlmnID"]
length(cpgs.chrX)

cat("\n\n... REMOVING CPGS LOCATED IN BOTH X AND Y CHROMOSOMES ...\n\n", file = stdout())
M.matrix.subset.filt.Y <- M.matrix.subset[!rownames(M.matrix.subset)%in%cpgs.chrY,]
M.matrix.subset.filt.XY <- M.matrix.subset.filt.Y[!rownames(M.matrix.subset.filt.Y)%in%cpgs.chrX,]
cat("\n// M values XY filtered //\n", file = stdout())
M.matrix.subset.filt.XY[1:5,1:5]
dim(M.matrix.subset.filt.XY)

cat("\n\n... 6.2) LOADING THE LIST WITH THE 20.000 MOST VARIABLE CPGS AND SUBSETTING THEM ...\n\n", file = stdout())

cpgs.sd.name <- load("methylation/QC/10.reducing_dimensions/SD/cpgs_SD_ranked.RData")
cpgs.sd <- get(cpgs.sd.name)
cat("\n\n // CpGs ranked by SD //\n", file = stdout())
length(cpgs.sd)
cat("\nTop 5\n", file = stdout())
cpgs.sd[1:5]
cat("\nLast 5\n", file = stdout())
cpgs.sd[19996:20000]
cat("\n// M values SD filtered //\n", file = stdout())
M.matrix.subset.filt.XY.sd <- M.matrix.subset.filt.XY[rownames(M.matrix.subset.filt.XY)%in%names(cpgs.sd),]
M.matrix.subset.filt.XY.sd[1:5,1:5]
dim(M.matrix.subset.filt.XY.sd)

cat("\n\n... 6.3) REMOVING genes FROM CHR X AND Y ...\n\n", file = stdout())

cat("\n\n... LOADING genes METADATA ...\n\n", file = stdout())
genes.metadata.batch1 <- read.csv(file = "transcriptomics/phenotype/phe000002.v7.FHS_SABRe_project3.marker-info.MULTI/GPL5175.txt",
                                  header=TRUE, stringsAsFactors = F, sep="\t", skip = 14)
head(genes.metadata.batch1)
dim(genes.metadata.batch1)
genes.metadata.batch2 <- read.csv(file = "transcriptomics/phenotype/phe000002.v7.FHS_SABRe_project3.marker-info.MULTI/GPL5188.txt",
                                  header=TRUE, stringsAsFactors = F, sep="\t", skip = 18)
head(genes.metadata.batch2)
dim(genes.metadata.batch2)
cat("\nNumber of genes located in Chromosome Y\n", file = stdout())
genes.chrY.1 <- genes.metadata.batch1[which(genes.metadata.batch1$seqname == "chrY"),"ID"]
length(genes.chrY.1)
genes.chrY.2 <- genes.metadata.batch2[which(genes.metadata.batch2$seqname == "chrY"),"ID"]
length(genes.chrY.2)
genes.chrY <- c(genes.chrY.1, genes.chrY.2)
length(genes.chrY)
cat("\nNumber of genes located in Chromosome X\n", file = stdout())
genes.chrX.1 <- genes.metadata.batch1[which(genes.metadata.batch1$seqname == "chrX"),"ID"]
length(genes.chrX.1)
genes.chrX.2 <- genes.metadata.batch2[which(genes.metadata.batch2$seqname == "chrX"),"ID"]
length(genes.chrX.2)
genes.chrX <- c(genes.chrX.1, genes.chrX.2)
length(genes.chrX)

cat("\n\n... REMOVING GENES LOCATED IN BOTH X AND Y CHROMOSOMES ...\n\n", file = stdout())
ESet.subset.filt.Y <- ESet.subset[!rownames(ESet.subset)%in%genes.chrY,]
ESet.subset.filt.XY <- ESet.subset.filt.Y[!rownames(ESet.subset.filt.Y)%in%genes.chrX,]
cat("\n// ExpressionSet subset XY filtered //\n", file = stdout())
ESet.subset.filt.XY[1:5,1:5]
dim(ESet.subset.filt.XY)

cat("\n\n... 7) CREATING A LIST WITH DATA MATRICES (MRNA + METHYLATION)  ...\n\n", file = stdout())

# Also order alphabetically the sample names for each data set, and converting to matrices instead of dataframes
data.matrices.list <- list(mRNA = as.matrix(ESet.subset.filt.XY[,order(colnames(ESet.subset.filt.XY))]), meth = as.matrix(M.matrix.subset.filt.XY.sd[,order(colnames(M.matrix.subset.filt.XY.sd))]))

cat("\nChecking dimensionalities, samples are columns, features are rows\n", file = stdout())

lapply(data.matrices.list, dim)

cat("\nAnd confirm both datasets are matrices and contain the same sample names as the Samples metadata\n", file = stdout())

lapply(data.matrices.list, function(matrix) {
    print(class(matrix))
    print(paste("Class of rows and columns: ",class(rownames(matrix))," and ",class(colnames(matrix)), sep = ""))
    print(matrix[1:5,1:5])
    table(colnames(matrix)%in%rownames(MOFA.covariates.order))
})

cat("\nOrder the data matrices samples IDs with the samples IDs from the samples metadata\n", file = stdout())

# We need to order the sample IDs in: M matrix (2080), gene expression matrix (938) and samples metadata (2080)
# following the order of the gene expression matrix, as it contains the lowest dimensions, and the rest of samples
# can be in any order but ordered between M matrix and samples metadata. Steps:

cat("\nSplit M matrix in matched with mRNA and non-matched\n", file = stdout())
match.meth.col <- which(colnames(data.matrices.list$meth)%in%colnames(data.matrices.list$mRNA))
nomatch.meth.col <- which(!colnames(data.matrices.list$meth)%in%colnames(data.matrices.list$mRNA))
M.matrix.match <- data.matrices.list$meth[,match.meth.col]
M.matrix.nomatch <- data.matrices.list$meth[,nomatch.meth.col]
dim(M.matrix.match)
dim(M.matrix.nomatch)

cat("\nOrder matched M.matrix with the same order as mRNA matrix\n", file = stdout())
match.ids <- match(colnames(M.matrix.match),colnames(data.matrices.list$mRNA))
M.matrix.match <- M.matrix.match[,match.ids]
cat("\nM matrix\n", file = stdout())
head(colnames(M.matrix.match))
cat("\nmRNA\n", file = stdout())
head(colnames(data.matrices.list$mRNA))

cat("\nJoin again the two M matrices ordered\n", file = stdout())
M.matrix.ordered <- cbind(M.matrix.match,M.matrix.nomatch)
dim(M.matrix.ordered)

cat("\nOrder samples metadata by the M.matrix order\n", file = stdout())
rownames(MOFA.covariates.order) <- MOFA.covariates.order$shareid
sample.ids <- match(colnames(M.matrix.ordered), MOFA.covariates.order$shareid)
MOFA.covariates.order <- MOFA.covariates.order[sample.ids,]
# 2055 23
samples <- "914_2055"

cat("\nCheck all have the same order\n", file = stdout())
cat("\nM matrix\n", file = stdout())
colnames(M.matrix.ordered)[1:20]
cat("\nmRNA\n", file = stdout())
colnames(data.matrices.list$mRNA)[1:20]
cat("\nSamples metadata\n", file = stdout())
MOFA.covariates.order$shareid[1:20]

cat("\nWhere are the outliers?\n", file = stdout())
# mRNA doesn't have them
which(colnames(data.matrices.list$mRNA)%in%tsne.batch.outliers)
# Should be from 938 onwards
which(colnames(M.matrix.ordered)%in%tsne.batch.outliers)

cat("\nCreate again the list data matrices now ordered\n", file = stdout())
data.matrices.list.ordered <- list(mRNA = data.matrices.list$mRNA, meth = M.matrix.ordered)

lapply(data.matrices.list.ordered, function(matrix) {
  print(class(matrix))
  print(paste("Class of rows and columns: ",class(rownames(matrix))," and ",class(colnames(matrix)), sep = ""))
  print(matrix[1:5,1:5])
  print(dim(matrix))
  table(colnames(matrix)%in%rownames(MOFA.covariates.order))
})

write.csv(MOFA.covariates.order, file = paste("MOFA/MOFA_covariates_order_",groups,"_",samples,".csv", sep = ""),
          row.names = TRUE)

#save(data.matrices.list, file = "MOFA/data_matrices.RData")

cat("\n\n... 8) CREATING MultiAssayExperiment OBJECT  ...\n\n", file = stdout())

mae <- MultiAssayExperiment(
  experiments = data.matrices.list 
)
cat("\n\n // MultiAssayExperiment //\n", file = stdout())
mae

cat("\n\n... 9) CREATING MOFA OBJECT  ...\n\n", file = stdout())

# Build the MOFA object
# With groups
if (groups == "groups") {
  MOFAobject <- create_mofa(mae, groups=MOFA.covariates.order$cvd)
}
# Without groups
if (groups == "nogroups") {
  MOFAobject <- create_mofa(mae)
}

cat("\n\n // MOFAobject //\n", file = stdout())
print(MOFAobject)

# If groups
# cat("\n\n Confirm CVD group from MOFA object corresponds to CVD group from MOFA covariates\n\n", file = stdout())
# head(MOFAobject@samples_metadata)
# table(MOFAobject@samples_metadata$sample[1:217]%in%MOFA.covariates.order[MOFA.covariates.order$cvd == 1,"shareid"])

cat("\n\n Plot overview of training Data \n", file = stdout())
png(paste("MOFA/1.create_mofa_object/data_overview_",groups,"_",samples,".png", sep = ""), height=3600, width=6000, res=600, units="px")
plot_data_overview(MOFAobject)
dev.off()

cat("\n\n ... SAVING MOFA OBJECT ... \n\n", file = stdout())

save(MOFAobject, file = paste("MOFA/1.create_mofa_object/MOFAobject_",groups,"_",samples,".RData", sep = ""))

date()

cat("\n\nComputational time\n\n", file = stdout())
end_time <- Sys.time()
end_time - start_time


cat("\n\n ######################## 1. MOFA OBJECT SCRIPT ENDS ######################## \n\n", file = stdout())




