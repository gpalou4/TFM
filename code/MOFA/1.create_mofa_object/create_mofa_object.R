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

cat("\n\n ######################## FINISHED LOADING LIBRARIES ######################## \n\n", file = stdout())

cat("\n\n ######################## 1.MOFA OBJECT SCRIPT BEGINS ######################## \n\n", file = stdout())

cat("\n\n... SETTING WORKING DIRECTORY ...\n\n", file = stdout())
setwd("/projects/regicor/Guille/TFM")
getwd()

groups <- "nogroups"
samples <- "914"

cat("\n\n... 1) LOADING DATA ...\n\n", file = stdout())

cat("\n\n... 1.1) LOADING GENE EXPRESSION DATA: ExpressionSet (its matrix) ...\n\n", file = stdout())
ESet.name <- load(paste("transcriptomics/QC/6.batch_effect/ESet_norm_filt_batch_938.RData", sep = ""))
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

cat("\n\n... 2.1) LOADING METHYLATION SAMPLES METADATA THAT MATCHES WITH MRNA (1522) ...\n\n", file = stdout())
cpgs.samples.metadata <- read.csv(file = "methylation/phenotype/samples_metadata/Analysis/samples_metadata_merged_filt.csv",
                             header=TRUE, stringsAsFactors = F, sep=",")
head(cpgs.samples.metadata)
dim(cpgs.samples.metadata)
# 1522 29

cat("\n\n... 2.2) LOADING METHYLATION SAMPLES METADATA THAT CONTAINS SVA (2080) ...\n\n", file = stdout())
cpgs.samples.metadata.sva <- read.csv(file = "methylation/phenotype/samples_metadata/Analysis/samples_metadata_merged_sva.csv",
                                  header=TRUE, stringsAsFactors = F, sep=",")
head(cpgs.samples.metadata.sva)
dim(cpgs.samples.metadata.sva)
# 2080 30

cat("\n\n... 2.3) LOADING mRNA SAMPLES METADATA (1522)  ...\n\n", file = stdout())
mrna.samples.metadata <- read.csv(file = "transcriptomics/phenotype/Analysis/samples_metadata_merged.csv",
                                      header=TRUE, stringsAsFactors = F, sep=",")
mrna.samples.metadata$cel_files <- substr(mrna.samples.metadata$cel_files,1,nchar(mrna.samples.metadata$cel_files)-4)

head(mrna.samples.metadata)
dim(mrna.samples.metadata)
# 1522 30

cat("\n\n... 3) SUBSETTING SAMPLES METADATA TO MATCH THE SAME INDIVIDUALS OF THE DATA MATRICES ...\n\n", file = stdout())

cat("\n\n... 3.1) SUBSET OF mRNA SAMPLES METADATA (1522) THAT MATCHES TO GENE EXPRESSION MATRIX (1110) ...\n\n", file = stdout())

mrna.samples.metadata.subset <- mrna.samples.metadata[mrna.samples.metadata$cel_files%in%colnames(ESet),]
head(mrna.samples.metadata.subset)
dim(mrna.samples.metadata.subset)
# 938 30

cat("\n\n... 3.2) SUBSET OF METH SAMPLES METADATA (1522) THAT MATCHES TO GENE EXPRESSION MATRIX (1110) ...\n\n", file = stdout())

# Order both samples metadata by the shareid
cpgs.samples.metadata <- cpgs.samples.metadata[order(cpgs.samples.metadata$shareid),]
mrna.samples.metadata <- mrna.samples.metadata[order(mrna.samples.metadata$shareid),]
# Add mRNA sample names as a new column (cel_files)
cpgs.samples.metadata$cel_files <- mrna.samples.metadata$cel_files
head(cpgs.samples.metadata)
dim(cpgs.samples.metadata)
# 1522 30 (29+1)

# Perform the subsetting
cpgs.samples.metadata.subset <- cpgs.samples.metadata[cpgs.samples.metadata$cel_files%in%colnames(ESet),]
head(cpgs.samples.metadata.subset)
dim(cpgs.samples.metadata.subset)
# 938 30

cat("\n\n... 4) MERGE METH AND mRNA SAMPLES METADATA ...\n\n", file = stdout())

cat("\n\n... 4.1) ADD SVA TO METHYLATION SAMPLES METADATA ...\n\n", file = stdout())

cpgs.samples.metadata.sva.subset <- cpgs.samples.metadata.sva[cpgs.samples.metadata.sva$shareid%in%cpgs.samples.metadata.subset$shareid,]
cpgs.samples.metadata.subset$sva1 <- cpgs.samples.metadata.sva.subset$sva1
head(cpgs.samples.metadata.subset)
dim(cpgs.samples.metadata.subset)
# 938 31 (30+1)

cat("\n\n... 4.2) MERGE MRNA AND METHYLATION SAMPLES METADATA ...\n\n", file = stdout())

# Change column names to differentiate mRNA and METH batch columns
common.columns <- c("ID_MRNA","shareid", "sex", "age", "tot_chol", "hdl_chol", "ldl_chol", "trig", "sbp", "dbp", "glucose", "weight", "height", "waist_u", "smoke", "CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran", "cvd", "prev_ex8_cvd_date")
names(cpgs.samples.metadata.subset) <- c("ID_METH","shareid", "sex", "age", "tot_chol", "hdl_chol", "ldl_chol", "trig", "sbp", "dbp", "glucose", "weight", "height", "waist_u", "smoke", "CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran", "cvd", "prev_ex8_cvd_date", "SAMPID", "Sample.Well_METH", "Plate_METH", "Image.date_METH", "Slide_METH", "Array_METH", "ID_MRNA", "sva1_METH")
names(mrna.samples.metadata.subset) <- c("shareid","ID_MRNA", "Batch_MRNA", "Row_MRNA", "Column_MRNA", "Plate_location_MRNA", "dbGaP_Sample_ID", "BODY_SITE", "ANALYTE_TYPE", "sex", "age", "tot_chol", "hdl_chol", "ldl_chol", "trig", "sbp", "dbp", "glucose", "weight", "height", "waist_u", "smoke", "CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran", "cvd", "prev_ex8_cvd_date")
# Merge
MOFA.samples.metadata <- merge(cpgs.samples.metadata.subset, mrna.samples.metadata.subset, by = common.columns, sort = TRUE)
head(MOFA.samples.metadata)
dim(MOFA.samples.metadata)
# 938 38

cat("\n\n... 4.3) REMOVE COLUMNS THAT ARE NOT COVARIATES FOR THE MOFA MODEL (KEEP THE OUTCOME OF INTEREST: 'CVD') ...\n\n", file = stdout())

# We remove ldl_chol as in SVA (colinearity)
# Include cvd because we will save the csv and be used later on during the plots
# but WE WILL REMOVE THE COLUMN BEFORE CREATING THE MOFA MODEL
covariates <- c("cvd","sex","age","tot_chol","hdl_chol","trig","sbp","dbp","glucose","weight","height","waist_u","smoke","CD8T","CD4T","NK","Bcell","Mono","Gran","Array_METH","sva1_METH","Batch_MRNA")
MOFA.covariates <- MOFA.samples.metadata[,covariates]

cat("\n\n... SAVING NEW SAMPLES METADATA MERGED AND COVARIATES THAT WILL BE USED ON MOFA ...\n\n", file = stdout())

write.csv(MOFA.samples.metadata, file = paste("MOFA/MOFA_samples_metadata_",groups,"_",samples,".csv",sep = ""),
          row.names = FALSE)
write.csv(MOFA.covariates, file = paste("MOFA/MOFA_covariates_",groups,"_",samples,".csv",sep = ""),
          row.names = FALSE)

cat("\n\n... 5) SUBSETTING SAMPLES WITH THE NEW SAMPLES METADATA  ...\n\n", file = stdout())

# NOTE: here we use MOFA.samples.metadata (and not MOFA.covariates) because it contains the ID's for both METH and mRNA
# necessary to do the subsetting

cat("\n\n... 5.1) SUBSETTING METHYLATION SAMPLES  ...\n\n", file = stdout())

M.matrix.subset <- M.matrix[,colnames(M.matrix)%in%MOFA.samples.metadata$ID_METH]
cat("\n\n // M values matrix subset //\n", file = stdout())
dim(M.matrix.subset)
# 420764 938
M.matrix.subset[1:5,1:5]

cat("\nChanging sample names by their corresponding shareid\n", file = stdout())
sample.names.index <- match(colnames(M.matrix.subset), MOFA.samples.metadata$ID_METH)
sample.names.shareid <- MOFA.samples.metadata$shareid[sample.names.index]
colnames(M.matrix.subset) <- sample.names.shareid

cat("\n\n // M values matrix subset, names changed to IDs //\n", file = stdout())
dim(M.matrix.subset)
# 420764 938
M.matrix.subset[1:5,1:5]

cat("\nRemove samples from batch 15, because they seem to be outliers in TSNE plots\n", file = stdout())
tsne.batch.outliers <- MOFA.samples.metadata[MOFA.samples.metadata$Batch_MRNA == 15, "shareid"]
tsne.batch.outliers
length(tsne.batch.outliers)
# 24

M.matrix.subset <- M.matrix.subset[,!colnames(M.matrix.subset)%in%tsne.batch.outliers]
cat("\n\n // M values matrix subset batch 15 filtered //\n", file = stdout())
M.matrix.subset[1:5,1:5]
dim(M.matrix.subset)
# 420764 914

cat("\n\n... 5.2) SUBSETTING GENE EXPRESSION SAMPLES  ...\n\n", file = stdout())

ESet.subset <- ESet[,colnames(ESet)%in%MOFA.samples.metadata$ID_MRNA]
cat("\n\n // ExpressionSet subset //\n", file = stdout())
ESet.subset[1:5,1:5]
dim(ESet.subset)
# 20567 938

cat("\nChanging sample names by their corresponding shareid\n", file = stdout())

sample.names.index <- match(colnames(ESet.subset), MOFA.samples.metadata$ID_MRNA)
sample.names.shareid <- MOFA.samples.metadata$shareid[sample.names.index]
colnames(ESet.subset) <- sample.names.shareid

cat("\n\n // ExpressionSet subset, names changed to IDs //\n\n", file = stdout())
ESet.subset[1:5,1:5]
dim(ESet.subset)
# 20567 938

cat("\nConfirm both datasets contain the same samples\n", file = stdout())
table(colnames(ESet.subset)%in%colnames(M.matrix.subset))

# cat("\nRemove samples from batch 15, because they seem to be outliers in TSNE plots\n", file = stdout())
# tsne.batch.outliers <- MOFA.samples.metadata[MOFA.samples.metadata$Batch_MRNA == 15, "shareid"]
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
# 938 22

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

# We need to order the sample IDs in: M matrix (938), gene expression matrix (914) and samples metadata (938)
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
#rownames(MOFA.covariates.order) <- MOFA.covariates.order$shareid
sample.ids <- match(colnames(M.matrix.ordered), rownames(MOFA.covariates.order))
MOFA.covariates.order <- MOFA.covariates.order[sample.ids,]

cat("\nCheck all have the same order\n", file = stdout())
cat("\nM matrix\n", file = stdout())
colnames(M.matrix.ordered)[1:20]
cat("\nmRNA\n", file = stdout())
colnames(data.matrices.list$mRNA)[1:20]
cat("\nSamples metadata\n", file = stdout())
rownames(MOFA.covariates.order)[1:20]

cat("\nAre the outliers correctly ordered?\n", file = stdout())
# mRNA doesn't have them
which(colnames(data.matrices.list$mRNA)%in%colnames(M.matrix.nomatch))
# Should be the last 24
which(colnames(M.matrix.ordered)%in%colnames(M.matrix.nomatch))

# Create again the data list matrices
cat("\nCreate again the list data matrices now ordered\n", file = stdout())
data.matrices.list.ordered <- list(mRNA = data.matrices.list$mRNA, meth = M.matrix.ordered)

lapply(data.matrices.list.ordered, function(matrix) {
  print(class(matrix))
  print(paste("Class of rows and columns: ",class(rownames(matrix))," and ",class(colnames(matrix)), sep = ""))
  print(matrix[1:5,1:5])
  print(dim(matrix))
  table(colnames(matrix)%in%rownames(MOFA.covariates.order))
})

# save(data.matrices.list, file = "MOFA/data_matrices_list.RData")
write.csv(MOFA.covariates.order, file = paste("MOFA/MOFA_covariates_order_",groups,"_",samples,".csv",sep = ""),
          row.names = TRUE)

cat("\n\n... 8) CREATING MultiAssayExperiment OBJECT  ...\n\n", file = stdout())

mae <- MultiAssayExperiment(
  experiments = data.matrices.list, 
)

cat("\n\n // MultiAssayExperiment //\n", file = stdout())
mae

cat("\n\n... 9) CREATING MOFA OBJECT  ...\n\n", file = stdout())

# Build the MOFA object
# With groups
if (groups == "groups") {
  MOFAobject <- create_mofa(mae, groups=MOFA.covariates.order$cvd)
}

if (groups == "nogroups") {
  MOFAobject <- create_mofa(mae)
}

cat("\n\n // MOFAobject //\n", file = stdout())
print(MOFAobject)

cat("\n\n Plot overview of training Data \n", file = stdout())
png(paste("MOFA/1.create_mofa_object/data_overview_",groups,"_",samples,".png",sep = ""), height=3600, width=6000, res=600, units="px")
plot_data_overview(MOFAobject)
dev.off()

cat("\n\n ... SAVING MOFA OBJECT ... \n\n", file = stdout())


save(MOFAobject, file = paste("MOFA/1.create_mofa_object/MOFAobject_",groups,"_",samples,".RData",sep = ""))

date()

cat("\n\nComputational time\n\n", file = stdout())
end_time <- Sys.time()
end_time - start_time


cat("\n\n ######################## 1. MOFA OBJECT SCRIPT ENDS ######################## \n\n", file = stdout())




