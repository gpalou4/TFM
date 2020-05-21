rm(list=ls())
start_time <- Sys.time()
date()

cat("\n\n ######################## START LOADING LIBRARIES ######################## \n\n", file = stdout())

cat("\n\n NONE LOADED \n\n", file = stdout())

cat("\n\n ######################## FINISHED LOADING LIBRARIES ######################## \n\n", file = stdout())

sessionInfo()

cat("\n\n ######################## 10. CPGS SD RANKING SCRIPT BEGINS ######################## \n\n", file = stdout())

cat("\n\n... SETTING WORKING DIRECTORY ...\n\n", file = stdout())
setwd("/projects/regicor/Guille/TFM/methylation")
getwd()

cat("\n\n... LOADING M VALUES MATRIX FOR THE 2596 INDIVIDUALS ...\n\n", file = stdout())
mvals.name <- load("QC/6.normalization/m_norm_4SD.RData")
mvals <- get(mvals.name)
# 2596 individuals
cat("\n// M values //\n", file = stdout())
class(mvals)
mvals[1:5,1:5]
dim(mvals)

cat("\n\n... REMOVING CpGs FROM CHR X AND Y ...\n\n", file = stdout())

cat("\n\n... LOADING CpGs METADATA FROM ILLUMINA 450KS MANIFEST ...\n\n", file = stdout())
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
mvals.filt.Y <- mvals[!rownames(mvals)%in%cpgs.chrY,]
mvals.filt.XY <- mvals.filt.Y[!rownames(mvals.filt.Y)%in%cpgs.chrX,]
cat("\n// M values XY filtered //\n", file = stdout())
mvals.filt.XY[1:5,1:5]
dim(mvals.filt.XY)

cat("\n\n... APPLYING SD FOR EACH CPG ...\n\n", file = stdout())

cpgs.sd <- apply(mvals.filt.XY, 1, sd, na.rm = T)
cat("\nLength and head\n", file = stdout())
print(length(cpgs.sd))
print(head(cpgs.sd))

cat("\nSort and filter first 20.000 cpgs\n", file = stdout())
cpgs.sd.sorted <- sort(cpgs.sd, decreasing = TRUE)
cpgs.sd.sorted.filt <- cpgs.sd.sorted[1:20000]

cat("\nFirst 30 cpgs\n", file = stdout())
cpgs.sd.sorted.filt[1:30]
cat("\nlast 30 cpgs (of the 20.000 filtered)\n", file = stdout())
cpgs.sd.sorted.filt[19970:20000]

cat("\n\n... SAVING RANKED LIST OF CPGS BY SD ...\n\n", file = stdout())

save(cpgs.sd.sorted.filt, file = "QC/10.reducing_dimensions/SD/cpgs_SD_ranked.RData")

date()

cat("\n\nComputational time\n\n", file = stdout())
end_time <- Sys.time()
end_time - start_time


cat("\n\n ######################## 10. CPGS SD RANKING SCRIPT ENDS ######################## \n\n", file = stdout())