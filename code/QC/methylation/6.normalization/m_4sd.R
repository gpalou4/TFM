rm(list=ls())

date()

cat("\n\n ######################## 6.M VALUES OUTLIERS REMOVAL SCRIPT BEGINS ######################## \n\n", file = stdout())

cat("\n\n... SETTING WORKING DIRECTORY ...\n\n", file = stdout())
setwd("/projects/regicor/Guille/TFM/methylation/QC")
getwd()

# betas.norm.name <- load("6.normalization/betas_norm.RData")
# betas.norm <- get(betas.norm.name)
# dim(betas.norm)
# betas.norm <- t(betas.norm)
# dim(betas.norm)
# betas.norm[1:3,1:3]

cat("\n\n... LOADING NORMALIZED M VALUES ...\n\n", file = stdout())
mvals.norm.name <- load("6.normalization/m_norm.RData")
mvals.norm <- get(mvals.norm.name)
dim(mvals.norm)
mvals.norm <- t(mvals.norm)
dim(mvals.norm)
mvals.norm[1:3,1:3]

cat("\n\n... OUTLIERS IDENTIFICATION ...\n\n", file = stdout())

system.time({
  
  cat("\nCalculating missing data for each CpG\n", file = stdout())
  # 0 means no missing data
  cpgs.miss <- colSums(is.na(mvals.norm))
  # Individuals ID
  ind.names <- rownames(mvals.norm)
  # CpGs names
  cpgs.names <- colnames(mvals.norm)
  cat("\nCalculating mean and SD for each CpG\n", file = stdout())
  cpgs.mean <- colMeans(mvals.norm, na.rm = T)
  cpgs.sd <- apply(mvals.norm, 2, sd, na.rm = T)
  cat("\nSetting upper and lower boundaries for remove CpGs outliers\n", file = stdout())
  # that is, CpGs with a M/B value above 4SD or below 4SD of the mean
  lower <- cpgs.mean-4*cpgs.sd
  upper <- cpgs.mean+4*cpgs.sd
  mvals.norm <- t(mvals.norm)
  cat("\nTransforming outliers +-4 SD from the mean to NA's\n", file = stdout())
  mvals.norm <- ifelse(mvals.norm < lower | mvals.norm > upper, NA, mvals.norm)
  mvals.norm <- t(mvals.norm)
  cat("\nCalculating missing data for each CpG again with outliers being now NA\n", file = stdout())
  cpgs.miss.out <- colSums(is.na(mvals.norm))
  cat("\nCalculating SD by CpG without outliers\n", file = stdout())
  cpgs.sd.out <- apply(mvals.norm, 2, sd, na.rm=TRUE)
  
})

cat("\n\n... CHECKS ...\n\n", file = stdout())

cat("\nSD before outliers transformation into NA's\n", file = stdout())
range(cpgs.sd)
cat("\nSD after outliers transformation into NA's\n", file = stdout())
range(cpgs.sd.out)

cat("\nMissing data (NA's) before outliers transformation\n", file = stdout())
sum(cpgs.miss)
cat("\nMissing data (NA's) after outliers transformation\n", file = stdout())
sum(cpgs.miss.out)
cat("\nNumber of outliers transformed to NA\n", file = stdout())
sum(cpgs.miss.out - cpgs.miss)

cat("\n\n... SAVING M VALUES ...\n\n", file = stdout())

mvals.norm <- t(mvals.norm)
mvals.norm <- as.data.frame(mvals.norm, stringsAsFactors = F,)
mvals.norm[1:5,1:5]
dim(mvals.norm)

save(mvals.norm, file="6.normalization/m_norm_4SD.RData")

date()

cat("\n\n ######################## 6.M VALUES OUTLIERS REMOVAL SCRIPT ENDS ######################## \n\n", file = stdout())

