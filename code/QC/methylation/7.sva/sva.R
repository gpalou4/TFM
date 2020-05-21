rm(list=ls())
start_time <- Sys.time()
date()

cat("\n\n ######################## START LOADING LIBRARIES ######################## \n\n", file = stdout())

library(sva)
cat("\n\n sva LOADED \n\n", file = stdout())

cat("\n\n ######################## FINISHED LOADING LIBRARIES ######################## \n\n", file = stdout())

sessionInfo()

cat("\n\n ######################## 7. SVA CONTROL SCRIPT BEGINS ######################## \n\n", file = stdout())

cat("\n\n... SETTING WORKING DIRECTORY ...\n\n", file = stdout())
setwd("/projects/regicor/Guille/TFM/methylation")
getwd()

cat("\n\n... LOADING FULL SAMPLES METADATA ...\n\n", file = stdout())

samples.metadata <- read.csv(file="phenotype/samples_metadata/Analysis/samples_metadata_merged.csv",
                             header=TRUE, stringsAsFactors = FALSE, sep=",")
# 2080 individuals
dim(samples.metadata)
head(samples.metadata)
names(samples.metadata)

cat("\n\n... LOADING M VALUES MATRIX FOR THE 2596 INDIVIDUALS ...\n\n", file = stdout())
mvals.name <- load("QC/6.normalization/m_norm_4SD.RData")
mvals <- get(mvals.name)
# 2596 individuals
cat("\n// M values //\n", file = stdout())
class(mvals)
mvals[1:5,1:5]
dim(mvals)
cat("\nna.omit\n", file = stdout())
dim(na.omit(mvals))

# cat("\n\n... LOADING M VALUES MATRIX FOR THE 2596 INDIVIDUALS WITHOUT SD FILTERING ...\n\n", file = stdout())
# mvals.name <- load("QC/6.normalization/m_norm.RData")
# mvals <- get(mvals.name)
# # 2596 individuals
# cat("\n// M values without SD //\n", file = stdout())
# class(mvals)
# mvals[1:5,1:5]
# dim(mvals)
# cat("\nna.omit\n", file = stdout())
# dim(na.omit(mvals))
# cat("\nis.na\n", file = stdout())
# table(rowSums(is.na(mvals)) > 0)

# Match M vals matrix sample names with sample names from samples metadata
cat("\n\n... MATCHING M VALS MATRIX SAMPLE NAMES (2596) WITH SAMPLE NAMES FROM SAMPLES METADATA (2080)  ...\n\n", file = stdout())
mvals.filt <- mvals[,colnames(mvals)%in%samples.metadata$LABID]
cat("\n// M values filtered //\n", file = stdout())
mvals.filt[1:5,1:5]
dim(mvals.filt)
cat("\nna.omit\n", file = stdout())
dim(na.omit(mvals.filt))
cat("\nis.na\n", file = stdout())
table(rowSums(is.na(mvals.filt)) > 0)

cat("\nCheck that Sample names in M values matrix are equal to sample names from the samples metadata\n", file = stdout())

# Three different ways for checking it (Identical does not work, I think because it compares two objects being equal,
# even if they contain the same values? But in theory both objects are equal...)

if(identical(colnames(mvals.filt),samples.metadata$LABID) == TRUE)
{
  print("Identical: Samples ID correctly matched between M values and samples metadata")
}

if(all(colnames(mvals.filt)%in%samples.metadata$LABID))
{
  print("All: Samples ID correctly matched between M values and samples metadata")
}

if(sum(colnames(mvals.filt)%in%samples.metadata$LABID)==length(samples.metadata$LABID))
{
  print("Sum: Samples ID correctly matched between M values and samples metadata")
}

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

cat("\n\n... REMOVING CPGS LOCATED IN BOTH X AND Y CHROMOSOMES ...\n\n", file = stdout())
# mvals.filt.Y <- mvals.filt[-which(Illumina.450k.manifest$CHR=="Y"),]
# mvals.filt.XY <- mvals.filt.Y[-which(Illumina.450k.manifest$CHR=="X"),]
mvals.filt.Y <- mvals.filt[!rownames(mvals.filt)%in%cpgs.chrY,]
mvals.filt.XY <- mvals.filt.Y[!rownames(mvals.filt.Y)%in%cpgs.chrX,]
cat("\n// M values XY filtered //\n", file = stdout())
mvals.filt.XY[1:5,1:5]
dim(mvals.filt.XY)

cat("\n\n... 1) CREATING THE MODELS ...\n\n", file = stdout())

# Outcome of interest
outcome <- "cvd"
# Numerical covariates
# It's possible LDL introduces colinearity so we will remove it
# num.covariates <- c("age","tot_chol","hdl_chol","ldl_chol","sbp","dbp","trig","glucose","CD8T",
#                     "CD4T","NK","Bcell","Mono","Gran")
num.covariates <- c("age","tot_chol","hdl_chol","sbp","dbp","trig","glucose","CD8T",
                    "CD4T","NK","Bcell","Mono","Gran")
# Character/categorical covariates (including batch effect)
chr.covariates <- c("smoke","Plate","sex")

mod <- paste("model.matrix(~as.factor(", outcome,")")
if (!is.null(num.covariates)) {for (covi in num.covariates) mod = paste(mod, "+", covi) }
if (!is.null(chr.covariates)) {for (covi in chr.covariates) mod = paste(mod, "+", "factor(", covi, ")") }
mod <- paste(mod, ", data=samples.metadata)", sep="")
cat("\n // Full model // \n", file = stdout())
mod

mod0 <- NULL
if (!is.null(num.covariates)) {for (covi in num.covariates) mod0 = paste(mod0, "+", covi) }
if (!is.null(chr.covariates)) {for (covi in chr.covariates) mod0 = paste(mod0, "+", "factor(", covi, ")") }
mod0 <- substring(mod0, 3, nchar(mod0))
mod0 <- paste("model.matrix(~", mod0, ", data=samples.metadata)", sep="")
cat("\n // Null model // \n", file = stdout())
mod0

part1 <- mod
part2 <- mod0

mod <- eval(parse(text = mod))
mod0 <- eval(parse(text = mod0))

cat("\nFull model matrix\n", file = stdout())
head(mod)
cat("\nNull model matrix\n", file = stdout())
head(mod0)

cat("\n\n... 2) DOING SVA ...\n\n", file = stdout())

cat("\nM values matrix XY filtered dimensions without NA's\n", file = stdout())
dim(na.omit(mvals.filt.XY))
cat("\nSamples metadata dimensions\n", file = stdout())
dim(samples.metadata)
cat("\nMod dimensions\n", file = stdout())
dim(mod)
cat("\nMod0 dimensions\n", file = stdout())
dim(mod0)

n.sv = num.sv(na.omit(mvals.filt.XY), mod, method="leek")
cat("\nSurrogate variables found\n", file = stdout())
n.sv

cat("\n\n... 3) CONCATENATING SVA RESULTS IN SAMPLES METADATA ...\n\n", file = stdout())

sva.obj <- sva(as.matrix(na.omit(mvals.filt.XY)), mod, mod0, n.sv = n.sv)
cat("\nSVA object\n", file = stdout())
head(sva.obj$sv)

samples.metadata.sva <- cbind(samples.metadata, sva.obj$sv)

if (n.sv == 1) {
  
  names(samples.metadata.sva)[30] <- "sva1"
  
}

if (n.sv == 2) {
  
  names(samples.metadata.sva)[30] <- "sva1"
  names(samples.metadata.sva)[31] <- "sva2"
  
}


write.csv(samples.metadata.sva, file = "phenotype/samples_metadata/Analysis/samples_metadata_merged_sva.csv",
          row.names = FALSE)

date()

cat("\n\nComputational time\n\n", file = stdout())
end_time <- Sys.time()
end_time - start_time


cat("\n\n ######################## 7. SVA CONTROL SCRIPT ENDS ######################## \n\n", file = stdout())