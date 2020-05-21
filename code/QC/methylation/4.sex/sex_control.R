rm(list=ls())

date()

cat("\n\n ######################## START LOADING LIBRARIES ######################## \n\n", file = stdout())

library(minfi)
cat("\n\n minfi LOADED \n\n", file = stdout())
# source("https://bioconductor.org/biocLite.R")
# biocLite("IlluminaHumanMethylation450kanno.ilm10b2.hg19")
# library(IlluminaHumanMethylation450kanno.ilm10b2.hg19)
# cat("\n\n IlluminaHumanMethylation450kanno.ilm10b2.hg19 LOADED \n\n", file = stdout())

cat("\n\n ######################## FINISHED LOADING LIBRARIES ######################## \n\n", file = stdout())

cat("\n\n ######################## 4.sex CONTROL SCRIPT BEGINS ######################## \n\n", file = stdout())


cat("\n\n... SETTING WORKING DIRECTORY ...\n\n", file = stdout())
setwd("/projects/regicor/Guille/TFM/methylation/QC")
getwd()

cat("\n\n... LOADING Methylset.filt OBJECT ...\n\n", file = stdout())
Methylset.filt.name <- load("2.filtering/Methylset_2620_filt2.RData")
Methylset.filt <- get(Methylset.filt.name)
cat("\n// Methylset.filt //\n", file = stdout())
Methylset.filt

cat("\n\n... 1) OBTAINING BETA VALUES FROM CHROMOSOME X ...\n\n", file = stdout())

cat("\nCreate GenomicMethylSet object\n", file = stdout())
GenomicMethylSet <- mapToGenome(Methylset.filt)
cat("\n// GenomicMethylSet //\n", file = stdout())
getClass(GenomicMethylSet)
cat("\nSave GenomicMethylSet object\n", file = stdout())
save(GenomicMethylSet, file="3b.GenomicMethylset/GenomicMethylSet_2620.RData")

### BORRAR
# cat("\n\n... LOADING GenomicMethylSet OBJECT ...\n\n", file = stdout())
# GenomicMethylSet.name <- load("3b.GenomicMethylset/GenomicMethylSet.RData")
# #system("cat /proc/meminfo")
# GenomicMethylSet <- get(GenomicMethylSet.name)
# cat("\n// GenomicMethylSet //\n", file = stdout())
# getClass(GenomicMethylSet)
### BORRAR

cat("\nObtain Betas from ChrX\n", file = stdout())
chrX.meth <- minfi::getMeth(GenomicMethylSet)[getAnnotation(GenomicMethylSet)$chr=='chrX',]
chrX.unmeth <- getUnmeth(GenomicMethylSet)[getAnnotation(GenomicMethylSet)$chr=='chrX',]
cat("\nCalculate Beta values\n", file = stdout())
beta.values <- chrX.meth / (chrX.meth + chrX.unmeth + 100)
cat("\n // Betas // \n", file = stdout())
beta.values[1:5,1:5]
dim(beta.values)

save(beta.values, file = "4.sex/chrX_beta_values.RData")

cat("\n\n... 2) CREATING CLASSICAL MULTIDIMENSIONAL SCALING (MDS) PLOT ...\n\n", file = stdout())

# ### BORRAR
# beta.values.name <- load(file="4.sex/chrX_beta_values.RData")
# beta.values <- get(beta.values.name)
# ### BORRAR

cat("\n\n... LOADING SAMPLES METADATA ...\n\n", file = stdout())

# samples.metadata: 2620 individuals
samples.metadata <- read.csv("../phenotype/samples_metadata/QC/samples_metadata_merged_QC_2620.csv")
cat("\n // Samples metadata // \n", file = stdout())
dim(samples.metadata)
head(samples.metadata)
cat("\nMatch samples from metadata and betas\n", file = stdout())

if (all(colnames(beta.values)%in%samples.metadata$LABID == TRUE)) {
  
  print("Sample names in beta.values matches correctly sample names in samples.metadata")
  
} else {
  
  print("Sample names in beta.values do not match sample names in samples.metadata")
  print(".. Subsetting ...")
  samples.metadata <- samples.metadata[samples.metadata$LABID%in%colnames(beta.values),]
  beta.values <- beta.values[,colnames(beta.values)%in%samples.metadata$LABID]
  # Does it match?
  cat("\nSamples metadata dimensions\n", file = stdout())
  dim(samples.metadata)
  cat("\nBetas dimensions\n", file = stdout())
  dim(beta.values)
  
}

cat("\nCreate MDS\n", file = stdout())
md <- cmdscale(dist(t(beta.values)), 2)

cat("\nMerge MDS dataframe with samples metadata\n", file = stdout())
md.df <- data.frame(md, stringsAsFactors=F)
md.df$LABID <- rownames(md.df)
md.metadata.merged <- merge(samples.metadata, md.df, by="LABID")
rownames(md.metadata.merged) <- md.metadata.merged$LABID
md.metadata.merged <- md.metadata.merged[rownames(md.df),]
identical(rownames(md.metadata.merged), rownames(md.df))

cat("\nCreate and save the MDS plot\n", file = stdout())
png("4.sex/sex_MDS.png", height=3600, width=6000, res=600, units="px")
plot(md, pch=ifelse(md.metadata.merged$SEX==1, "M", "F"), col=ifelse(md.metadata.merged$SEX==1, "red", "blue"))
dev.off()

cat("\n\n... 3) REMOVING FEMALE AND MALE OUTLIERS ...\n\n", file = stdout())

cat("\nNumber of outlier Males\n", file = stdout())
dim(md.metadata.merged[md.metadata.merged$SEX==1 & md.metadata.merged$X1 > 0,])[1]
cat("\nNumber of outlier Females\n", file = stdout())
dim(md.metadata.merged[md.metadata.merged$SEX==2 & md.metadata.merged$X1 < 0,])[1]
cat("\nRemove outliers for the MDS plot\n", file = stdout())
md.metadata.merged.filt <- md.metadata.merged[!(md.metadata.merged$SEX==1 & md.metadata.merged$X1 > 0) & !(md.metadata.merged$SEX==2 & md.metadata.merged$X1 < 0),]
cat("\nRemove outliers from the Methylset.filt object\n", file = stdout())
Methylset.filt2 <- Methylset.filt[,md.metadata.merged.filt$LABID]

cat("\n\n... 4) MDS PLOT WITHOUT OUTLIERS ...\n\n", file = stdout())

md.out <- md[rownames(md)%in%md.metadata.merged.filt$LABID,]

cat("\nCreate and save the MDS plot without outliers\n", file = stdout())
png("4.sex/sex_MDS_no_outliers.png", height=3600, width=6000, res=600, units="px")
plot(md.out, pch=ifelse(md.metadata.merged.filt$SEX==1, "M", "F"), col=ifelse(md.metadata.merged.filt$SEX==1, "red", "blue"))
dev.off()

write.table(md.metadata.merged.filt, file="4.sex/sex_MDS_data.csv", row.names=F, col.names=T, sep=",", quote=F)

cat("\n\n... SAVING Methylset.filt2 OBJECT WITHOUT OUTLIERS ...\n\n", file = stdout())

cat("\n// Methylset.filt2 //\n", file = stdout())
Methylset.filt2

save(Methylset.filt2, file="4.sex/Methylset_2620_filt3.RData")

date()

cat("\n\n ######################## 3.SEX CONTROL SCRIPT ENDS ######################## \n\n", file = stdout())




