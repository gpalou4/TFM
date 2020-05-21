rm(list=ls())
start_time <- Sys.time()
date()

cat("\n\n ######################## START LOADING LIBRARIES ######################## \n\n", file = stdout())

library(MASS)      
library(lmtest)    
library(parallel)
library(sandwich)
cat("\n\n NONE LOADED \n\n", file = stdout())

cat("\n\n ######################## FINISHED LOADING LIBRARIES ######################## \n\n", file = stdout())

sessionInfo()

cat("\n\n ######################## 10. MOST VARIABLE CPGS SCRIPT BEGINS ######################## \n\n", file = stdout())

cat("\n\n... SETTING WORKING DIRECTORY ...\n\n", file = stdout())
setwd("/projects/regicor/Guille/TFM/methylation")
getwd()

cat("\n\n... LOADING FULL SAMPLES METADATA, SVA INCLUDED...\n\n", file = stdout())

samples.metadata <- read.csv(file="phenotype/samples_metadata/Analysis/samples_metadata_merged_sva.csv",
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

# Match M vals matrix sample names with sample names from samples metadata
cat("\n\n... MATCHING M VALS MATRIX SAMPLE NAMES (2596) WITH SAMPLE NAMES FROM SAMPLES METADATA (2080)  ...\n\n", file = stdout())
mvals.filt <- mvals[,colnames(mvals)%in%samples.metadata$LABID]
cat("\n// M values filtered //\n", file = stdout())
mvals.filt[1:5,1:5]
dim(mvals.filt)
cat("\nna.omit\n", file = stdout())
dim(na.omit(mvals.filt))
#cat("\nis.na\n", file = stdout())
#table(rowSums(is.na(mvals.filt)) > 0)

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

cat("\nSort the samples names (rows) from the samples metadata with the same order as in the M values matrix (columns)\n", file = stdout())

rownames(samples.metadata) <- samples.metadata$LABID
samples.metadata.ordered <- samples.metadata[colnames(mvals.filt),]
# 2080 30
dim(samples.metadata.ordered)
head(samples.metadata.ordered)
names(samples.metadata.ordered)

cat("\n\n... REMOVING CPGS FROM CHR X AND Y ...\n\n", file = stdout())

cat("\n\n... LOADING CPGS METADATA FROM ILLUMINA 450KS MANIFEST ...\n\n", file = stdout())
Illumina.450k.manifest <- read.csv("phenotype/phg000492.v1.FHS_DNAMethylation.marker-info.MULTI/HumanMethylation450_15017482_v.1.1.csv",
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
mvals.filt.Y <- mvals.filt[!rownames(mvals.filt)%in%cpgs.chrY,]
mvals.filt.XY <- mvals.filt.Y[!rownames(mvals.filt.Y)%in%cpgs.chrX,]
cat("\n// M values XY filtered //\n", file = stdout())
mvals.filt.XY[1:5,1:5]
dim(mvals.filt.XY)

cat("\n\n... 1) CREATING THE CHARACTER VECTORS OF THE VARIABLES FOR THE REGRESSION MODEL ...\n\n", file = stdout())

cat("\nOutcome of interest (Y)\n", file = stdout())
outcome <- "cvd"
outcome
cat("\nNumerical covariates (Xcov)\n", file = stdout())
# It's possible LDL introduces colinearity so we will remove it
# num.covariates <- c("age","tot_chol","hdl_chol","ldl_chol","sbp","dbp","trig","glucose","CD8T",
#                     "CD4T","NK","Bcell","Mono","Gran")
num.covariates <- c("age","tot_chol","hdl_chol","sbp","dbp","trig","glucose","CD8T",
                    "CD4T","NK","Bcell","Mono","Gran")
num.covariates
cat("\nCharacter/categorical covariates including batch effect (Xcov)\n", file = stdout())
chr.covariates <- c("smoke","Plate","sex")
chr.covariates

cat("\n\n... 2) OBTAINING THE VALUES OF THE VARIABLES FOR THE REGRESSION MODEL ...\n\n", file = stdout())

cat("\nOutcome of interest (Y)\n", file = stdout())
samples.metadata.ordered.outcome <- samples.metadata.ordered[,outcome]

cat("\n Covariates (Xcov)\n", file = stdout())
Xcov <- "cbind("
if (!is.null(num.covariates)) {for (covi in num.covariates) Xcov=paste(Xcov, paste("samples.metadata$", covi,",", sep=""), sep="") }
if (!is.null(chr.covariates)) {for (covi in chr.covariates) Xcov=paste(Xcov, paste("factor(samples.metadata.ordered$", covi, "),", sep=""), sep="") }
Xcov <- substring(Xcov, 1, nchar(Xcov)-1)
Xcov <- paste(Xcov,")", sep="")
Xcov
# Apply the cbind
Xcov <- eval(parse(text=Xcov))
head(Xcov)

cat("\n\n... 2) ROBUST LINEAR REGRESSION (EWAS) ...\n\n", file = stdout())

#Robust linear regression function
f.RLM.Adjusted.Robust = function(Y, X, Xcov) {
  mod = rlm(Y~X+Xcov, maxit=200)
  cf = coeftest(mod, vcov=vcovHC(mod, type="HC0"))
  coef = cf[2,"Estimate"]
  se = cf[2,"Std. Error"]
  p.value = cf[2,"Pr(>|z|)"]
  out = c(coef, se, p.value)
  names(out) = c("coef","se", "pvalue")
  return(out)
}

#Doing EWAS (Robust linear regression for each CpG)
EWAS.results <- lapply(rownames(mvals.filt.XY), function(i){
  #print(head(rownames(mvals.filt.XY)))
  #print(i)
  #X <- as.matrix(mvals.filt.XY[i,])
  X <- t(mvals.filt.XY)[,i]
  #print(class(X))
  #print(head(X))
  f.RLM.Adjusted.Robust(as.numeric(samples.metadata.ordered.outcome), X, Xcov)[1:3]
})

cat("\n\n... SAVING EWAS RESULTS ...\n\n", file = stdout())

EWAS.results <- matrix(unlist(EWAS.results), ncol=3, byrow=T)
colnames(EWAS.results) <- c("Coefficient","SE", "Pvalue")
rownames(EWAS.results) <- cpg
EWAS.results <- as.data.frame(EWAS.results)
EWAS.results$cpg <- rownames(EWAS.results)
EWAS.results <- EWAS.results[order(EWAS.results[,3]),]
EWAS.results <- EWAS.results[,c(4,1:3)]

write.table(EWAS.results, file="QC/10.reducing_dimensions/EWAS/EWAS_results.csv", row.names=F, col.names=T, sep=";", quote=F)
save(EWAS.results, file="QC/10.reducing_dimensions/EWAS/EWAS_results.RData")

date()

cat("\n\nComputational time\n\n", file = stdout())
end_time <- Sys.time()
end_time - start_time


cat("\n\n ######################## 10. MOST VARIABLE CPGS SCRIPT ENDS ######################## \n\n", file = stdout())