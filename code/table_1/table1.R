rm(list=ls())
start_time <- Sys.time()
date()

cat("\n\n ######################## START LOADING LIBRARIES ######################## \n\n", file = stdout())

cat("\n\n ######################## FINISHED LOADING LIBRARIES ######################## \n\n", file = stdout())

cat("\n\n ######################## RESULTS: TABLE 1 SCRIPT BEGINS ######################## \n\n", file = stdout())

setwd("/projects/regicor/Guille/TFM")

# Descriptive characteristics table

## FOR 2055 INDIVIDUALS 

groups <- "nogroups"
samples <- "914_2055"

MOFA.covariates <- read.csv(file = paste("MOFA/MOFA_covariates_order_",groups,"_",samples,".csv",sep = ""),
                            header=TRUE, stringsAsFactors = F, sep=",")
# Metadata has to contain the columns 'sample' (shareid) and 'group' (cvd) for MOFA to work
colnames(MOFA.covariates)[1] <- "sample"
# Change sample names from integers to characters
MOFA.covariates$sample <- as.character(MOFA.covariates$sample)
MOFA.covariates[1:3,]
dim(MOFA.covariates)
# 2055 24

# Add another covariable: CHD incidence
cov <- read.table("/projects/regicor/data/FHS/phenotype/phs000007.v29.pht003316.v6.p10.c1.vr_survcvd_2014_a_1023s.HMB-IRB-MDS.txt", header=T, stringsAsFactors=F, sep="\t")
head(cov)
dim(cov)
chd.cov <- cov[cov$shareid%in%MOFA.covariates$sample,c("shareid","chd")]
colnames(chd.cov)[1] <- "sample"
MOFA.covariates <- merge(MOFA.covariates, chd.cov, by = "sample")
dim(MOFA.covariates)
#2055 25

# Calculate BMI
bmi <- MOFA.covariates$weight/((MOFA.covariates$height/100)**2)
MOFA.covariates$BMI <- bmi

# Interested covariates
# continous variables for mean + sd
cov.cont.mean <- c("age", "tot_chol", "hdl_chol", "sbp", "dbp", "waist_u", "BMI")
# continous variables for median + IQR
cov.cont.median <- c("trig", "glucose")
cov.cat <- c("sex","smoke", "cvd", "chd")

cov.cont.df.mean <- c()
cov.cont.df.median <- c()
# Mean + sd
for (i in cov.cont.mean) {
  cov.name <- i
  cov.mean <- mean(MOFA.covariates[,i])
  cov.sd <- sd(MOFA.covariates[,i])
  cov.na <- sum(is.na(MOFA.covariates[,i]))
  cov.cont.df.mean <- as.data.frame(rbind(cov.cont.df.mean, cbind(cov.name,cov.mean, cov.sd, cov.na)))
}

# Median + IRQ
for (i in cov.cont.median) {
  cov.name <- i
  cov.median <- median(MOFA.covariates[,i])
  cov.Q3 <- quantile(MOFA.covariates[,i])[4]
  cov.Q1 <- quantile(MOFA.covariates[,i])[2]
  cov.na <- sum(is.na(MOFA.covariates[,i]))
  cov.cont.df.median <- as.data.frame(rbind(cov.cont.df.median, cbind(cov.name, cov.median, cov.Q1, cov.Q3, cov.na)))
}

write.csv(cov.cont.df.mean, file = "manuscript/Results/tables/cov_cont_mean_2055.csv")
write.csv(cov.cont.df.median, file = "manuscript/Results/tables/cov_cont_median_2055.csv")


cov.cat.df <- c()
for (i in cov.cat) {
  
  cov.name <- i
  if (i == "smoke" | i == "sex"){
    cov.number <- as.vector(table(MOFA.covariates[,cov.name])[1])
    cov.percentage <- as.vector(table(MOFA.covariates[,cov.name])[1]/dim(MOFA.covariates)[1]*100)
  }
  else {
    cov.number <- as.vector(table(MOFA.covariates[,cov.name])[2])
    cov.percentage <- as.vector(table(MOFA.covariates[,cov.name])[2]/dim(MOFA.covariates)[1]*100)
  }
  cov.na <- sum(is.na(MOFA.covariates[,cov.name]))
  
  cov.cat.df <- as.data.frame(rbind(cov.cat.df, cbind(cov.name,cov.number, cov.percentage, cov.na)))
  
}
write.csv(cov.cat.df, file = "manuscript/Results/tables/cov_cat_2055.csv")


## FOR 914 INDIVIDUALS 

groups <- "nogroups"
samples <- "914"

MOFA.covariates.mRNA <- read.csv(file = paste("MOFA/MOFA_covariates_order_",groups,"_",samples,".csv",sep = ""),
                            header=TRUE, stringsAsFactors = F, sep=",")
# Metadata has to contain the columns 'sample' (shareid) and 'group' (cvd) for MOFA to work
colnames(MOFA.covariates.mRNA)[1] <- "sample"
# Change sample names from integers to characters
MOFA.covariates.mRNA$sample <- as.character(MOFA.covariates.mRNA$sample)
MOFA.covariates.mRNA[1:3,]
dim(MOFA.covariates.mRNA)
# 914 24

# Add another covariable: CHD incidence
cov <- read.table("/projects/regicor/data/FHS/phenotype/phs000007.v29.pht003316.v6.p10.c1.vr_survcvd_2014_a_1023s.HMB-IRB-MDS.txt", header=T, stringsAsFactors=F, sep="\t")
head(cov)
dim(cov)
chd.cov <- cov[cov$shareid%in%MOFA.covariates.mRNA$sample,c("shareid","chd")]
colnames(chd.cov)[1] <- "sample"
MOFA.covariates.mRNA <- merge(MOFA.covariates.mRNA, chd.cov, by = "sample")
dim(MOFA.covariates.mRNA)
# 914 25

# Calculate BMI
bmi <- MOFA.covariates.mRNA$weight/((MOFA.covariates.mRNA$height/100)**2)
MOFA.covariates.mRNA$BMI <- bmi

# Interested covariates
# continous variables for mean + sd
cov.cont.mean <- c("age", "tot_chol", "hdl_chol", "sbp", "dbp", "waist_u", "BMI")
# continous variables for median + IQR
cov.cont.median <- c("trig", "glucose")
cov.cat <- c("sex","smoke", "cvd", "chd")

cov.cont.df.mean <- c()
cov.cont.df.median <- c()
# Mean + sd
for (i in cov.cont.mean) {
  cov.name <- i
  cov.mean <- mean(MOFA.covariates.mRNA[,i])
  cov.sd <- sd(MOFA.covariates.mRNA[,i])
  cov.na <- sum(is.na(MOFA.covariates.mRNA[,i]))
  cov.cont.df.mean <- as.data.frame(rbind(cov.cont.df.mean, cbind(cov.name,cov.mean, cov.sd, cov.na)))
}

# Median + IRQ
for (i in cov.cont.median) {
  cov.name <- i
  cov.median <- median(MOFA.covariates.mRNA[,i])
  cov.Q3 <- quantile(MOFA.covariates.mRNA[,i])[4]
  cov.Q1 <- quantile(MOFA.covariates.mRNA[,i])[2]
  cov.na <- sum(is.na(MOFA.covariates.mRNA[,i]))
  cov.cont.df.median <- as.data.frame(rbind(cov.cont.df.median, cbind(cov.name, cov.median, cov.Q1, cov.Q3, cov.na)))
}

write.csv(cov.cont.df.mean, file = "manuscript/Results/tables/cov_cont_mean_914.csv")
write.csv(cov.cont.df.median, file = "manuscript/Results/tables/cov_cont_median_914.csv")


cov.cat.df <- c()
for (i in cov.cat) {
  
  cov.name <- i
  if (i == "smoke" | i == "sex"){
    cov.number <- as.vector(table(MOFA.covariates.mRNA[,cov.name])[1])
    cov.percentage <- as.vector(table(MOFA.covariates.mRNA[,cov.name])[1]/dim(MOFA.covariates.mRNA)[1]*100)
  }
  else {
    cov.number <- as.vector(table(MOFA.covariates.mRNA[,cov.name])[2])
    cov.percentage <- as.vector(table(MOFA.covariates.mRNA[,cov.name])[2]/dim(MOFA.covariates.mRNA)[1]*100)
  }
  cov.na <- sum(is.na(MOFA.covariates.mRNA[,cov.name]))
  
  cov.cat.df <- as.data.frame(rbind(cov.cat.df, cbind(cov.name,cov.number, cov.percentage, cov.na)))
  
}
write.csv(cov.cat.df, file = "manuscript/Results/tables/cov_cat_914.csv")


## Significant differences between the two groups (meth vs mrna) for each covariate

mRNA.samples <- MOFA.covariates.mRNA$sample
MOFA.covariates$omic <- "Meth"
MOFA.covariates.copy <- MOFA.covariates
MOFA.covariates.copy$omic <- "mRNA"
MOFA.covariates <- rbind(MOFA.covariates,MOFA.covariates.copy[MOFA.covariates.copy$sample%in%mRNA.samples,])
table(MOFA.covariates$omic)

# T-test for means
cov.cont.df.ttest <- c()
for (i in cov.cont.mean) {
  cov.name <- i
  #cov.ttest <- t.test(MOFA.covariates[MOFA.covariates$sample%in%mRNA.samples,i], MOFA.covariates[,i])
  cov.ttest <- t.test(eval(parse(text=cov.name)) ~ omic, data = MOFA.covariates, paired = FALSE)
  cov.cont.df.ttest <- as.data.frame(rbind(cov.cont.df.ttest, cbind(cov.name,round(cov.ttest$p.value,4), "T-test")))
}

# Kruskal-Wallis test for medians
for (i in cov.cont.median) {
  cov.name <- i
  cov.ttest <- kruskal.test(eval(parse(text=cov.name)) ~ omic, data = MOFA.covariates)
  cov.cont.df.ttest <- as.data.frame(rbind(cov.cont.df.ttest, cbind(cov.name,round(cov.ttest$p.value,4), "Kruskal-Wallis test")))
}

colnames(cov.cont.df.ttest) <- c("Covariate","p-value","Test")
write.csv(cov.cont.df.ttest, file = "manuscript/Results/tables/cov_cont_ttest.csv")

cov.cat.df.ttest <- c()
for (i in cov.cat) {
  cov.name <- i
  table.cat <- cbind(table(MOFA.covariates[MOFA.covariates$omic == "Meth",i]),table(MOFA.covariates[MOFA.covariates$omic == "mRNA",i]))
  cov.fishert <- fisher.test(table.cat)
  cov.cat.df.ttest <- as.data.frame(rbind(cov.cat.df.ttest, cbind(cov.name,round(cov.fishert$p.value,4))))
}

colnames(cov.cat.df.ttest) <- c("Covariate","Fisher Test")
write.csv(cov.cat.df.ttest, file = "manuscript/Results/tables/cov_cat_Fisher_test.csv")













