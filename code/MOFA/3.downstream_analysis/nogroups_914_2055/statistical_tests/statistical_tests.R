rm(list=ls())
start_time <- Sys.time()
date()

cat("\n\n ######################## START LOADING LIBRARIES ######################## \n\n", file = stdout())

library(ggplot2)
cat("\n\n ggplot2 LOADED \n\n", file = stdout())
library(MOFA2)
cat("\n\n MOFA LOADED \n\n", file = stdout())


cat("\n\n ######################## FINISHED LOADING LIBRARIES ######################## \n\n", file = stdout())

cat("\n\n ######################## 3. MOFA DOWNSTREAM ANALYSIS: STATISTICAL TESTS SCRIPT BEGINS ######################## \n\n", file = stdout())

cat("\n\n... SETTING WORKING DIRECTORY ...\n\n", file = stdout())
setwd("/projects/regicor/Guille/TFM")
getwd()

cat("\n\n... 1) LOADING DATA ...\n\n", file = stdout())


cat("\n\n... 1.1) LOADING SAMPLES METADATA ...\n\n", file = stdout())

groups <- "nogroups"
samples <- "914_2055"

MOFA.covariates <- read.csv(file = paste("MOFA/MOFA_covariates_order_",groups,"_",samples,".csv", sep=""),
                                  header=TRUE, stringsAsFactors = F, sep=",")
# Metadata has to contain the columns 'sample' (shareid) and 'group' (cvd) for MOFA to work
colnames(MOFA.covariates)[1] <- "sample"
# Change sample names from integers to characters
MOFA.covariates$sample <- as.character(MOFA.covariates$sample)
MOFA.covariates[1:3,]
dim(MOFA.covariates)
# 2055 24

cat("\n\n... 1.2) LOADING MOFA TRAINED OBJECT ...\n\n", file = stdout())

MOFA.model <- load_model(paste("MOFA/2.training_model/MOFA_trained_",groups,"_",samples,".hdf5", sep = ""))
MOFA.model

cat("\n\n... 2) ADDING SAMPLES METADATA (COVARIABLES) TO THE TRAINED MODEL OBJECT ...\n\n", file = stdout())

samples_metadata(MOFA.model) <- MOFA.covariates[-2]

cat("\n\n... 3) INSPECTING THE MOFA OBJECT ...\n\n", file = stdout())

cat("\nSlot names of the MOFA model\n", file = stdout())
slotNames(MOFA.model)
cat("\nMetadata stored in the MOFA model\n", file = stdout())
head(MOFA.model@samples_metadata)

cat("\n\n... 4) STATISTICAL TESTS BETWEEN INDIVIDUALS WITH CVD AND NO CVD ... \n\n", file = stdout())

cat("\n\n... 4.1) T-TEST AND MANN-WHITNEY TEST FOR INTERESTED FACTORS  ...\n\n", file = stdout())

testFactor <- function(factors, factor, covariable, group, MOFA_samples_metadata) {

  factors.df <- as.data.frame(factors)

  # Create a dataframe with the factor values and the covariable
  factor.df <- data.frame(interested_factor = eval(parse(text=paste("factors.df$Factor",factor,sep=""))),
                          covariable = as.factor(eval(parse(text=paste("MOFA_samples_metadata$",covariable,sep="")))),
                          row.names = rownames(factors.df))
  colnames(factor.df)[1] <- paste("Factor",factor,sep="")
  colnames(factor.df)[2] <- covariable
  print(head(factor.df))

  # Check histogram
  png(paste("MOFA/3.downstream_analysis/",groups,"_",samples,"/statistical_tests/factors_histograms/factor",factor,"_",covariable,"_",groups,"_groups_hist.png", sep = ""), height=3600, width=6000, res=600, units="px")
  histogram <- ggplot(factor.df, aes(eval(parse(text=paste("Factor",factor,sep=""))))) +
    geom_histogram(fill = "white", color = "grey30", bins = 30) +
    ggtitle(paste("Histogram for factor ",factor, sep = "")) +
    labs(y="Count", x = "Values") +
    theme (plot.title = element_text(hjust = 0.5))
  print(histogram)
  dev.off()

  # Check boxplot
  png(paste("MOFA/3.downstream_analysis/",groups,"_",samples,"/statistical_tests/factors_boxplots/factor",factor,"_",covariable,"_",groups,"_groups_boxplot.png", sep = ""), height=3600, width=6000, res=600, units="px")
  boxplot <- ggplot(factor.df, aes(as.factor(eval(parse(text=covariable))), eval(parse(text=paste("Factor",factor,sep=""))))) +
    geom_boxplot() +
    ggtitle(paste("Boxplot for factor ",factor, sep = "")) +
    labs(y="Factor values", x = covariable) +
    theme (plot.title = element_text(hjust = 0.5))
  print(boxplot)
  dev.off()

  # T-test for mean (parametric test)
  t.test.res <- t.test(eval(parse(text=paste("Factor",factor,sep=""))) ~ eval(parse(text=covariable)), data = factor.df, paired = FALSE)

  # Mann–Whitney–Wilcoxon for means (non-parametric test)
  wilcox.test.res <- wilcox.test(eval(parse(text=paste("Factor",factor,sep=""))) ~ eval(parse(text=covariable)), data = factor.df, paired = FALSE, conf.int = 0.95)

  results.df <- data.frame(factor = paste("Factor",factor,sep=""), covariable = covariable, groups = group,
                           t_test = t.test.res$p.value,
                           m_whitney_wilcoxon = wilcox.test.res$p.value)

  return(results.df)

}

factors <- get_factors(MOFA.model,
                       groups = "all",
                       factors = "all"
)

lapply(factors,dim)
lapply(factors,head)

all.factors.res <- c()
for (i in 1:30) {
  test.factor.res <- testFactor(factors = factors$group1, factor = i, covariable = "cvd", group = "all", MOFA_samples_metadata = MOFA.model@samples_metadata)
  all.factors.res <- rbind(all.factors.res, test.factor.res)
}

write.csv(all.factors.res, file = paste("MOFA/3.downstream_analysis/",groups,"_",samples,"/statistical_tests/factors_cvd_tests.csv", sep = ""),
          row.names = FALSE)

cat("\n\n... 4.2) T-TEST AND MANN-WHITNEY TEST FOR TOP 20 FEATURE WEIGTHS FROM HEATMAPS + CORRELATION BETWEEN FEATURES ...\n\n", file = stdout())

corrFeatures <- function(matrix, factor, omic) {

  # Spearman Correlation
  features.corr.R = cor(na.omit(matrix), method = c("spearman"))
  features.corr.R

  # Visualize
  library(corrplot)
  library(RColorBrewer)
  col<- colorRampPalette(c("red", "white", "blue"))(20)

  # With circles
  png(paste("MOFA/3.downstream_analysis/",groups,"_",samples,"/statistical_tests/",omic,"_correlation_features_1_factor_",factor,".png", sep = ""), height=3600, width=6000, res=600, units="px")
  corrplot(features.corr.R,
           order="hclust",
           col = col,
           type = "upper",
           tl.srt=45)
  dev.off()

  # Adding significant p-values
  cor.mtest <- function(mat, ...) {
    mat <- as.matrix(mat)
    n <- ncol(mat)
    p.mat<- matrix(NA, n, n)
    diag(p.mat) <- 0
    for (i in 1:(n - 1)) {
      for (j in (i + 1):n) {
        tmp <- cor.test(mat[, i], mat[, j], method = c("spearman"))
        p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
      }
    }
    colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
    p.mat
  }

  p.mat <- cor.mtest(matrix)

  # Alternative plot coloring significant correlations (in blank unsignificant ones)
  col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))

  png(paste("MOFA/3.downstream_analysis/",groups,"_",samples,"/statistical_tests/",omic,"_correlation_features_2_factor_",factor,".png", sep = ""), height=3600, width=6000, res=600, units="px")
  corrplot(features.corr.R, method="color", col=col(200),
           type="upper", order="hclust",
           addCoef.col = "black", # Add coefficient of correlation
           tl.col="black", tl.srt=45, #Text label color and rotation
           # Combine with significance
           p.mat = p.mat, sig.level = 0.01, insig = "blank",
           # hide correlation coefficient on the principal diagonal
           diag=FALSE, number.cex = 0.45, tl.cex = 1
  )
  dev.off()

}

testFeaturesFactor <- function(features.matrix, factor, omic) {

  heatmap.features.tests <- data.frame()

  for (i in 1:length(colnames(features.matrix))) {

    features.factor.df <- data.frame(feature = features.matrix[,i], cvd = as.factor(MOFA.model@samples_metadata$cvd),
                                     row.names = rownames(features.matrix))
    if (omic == "mRNA") {
      colnames(features.factor.df)[1] <- "mRNA"
      label <- "Transcript expression"
    }
    else if (omic == "meth") {
      colnames(features.factor.df)[1] <- "cpg"
      label <- "CpG M value"
    }

    feature <- colnames(features.factor.df)[1]

    # Histogram
    png(paste("MOFA/3.downstream_analysis/",groups,"_",samples,"/statistical_tests/features_histograms/hist_",colnames(features.matrix)[i],".png", sep = ""), height=3600, width=6000, res=600, units="px")
    histogram <- ggplot(features.factor.df, aes(eval(parse(text=feature)))) +
      geom_histogram(fill = "white", color = "grey30", bins = 30) +
      ggtitle(paste("Histogram for ",colnames(features.matrix)[i]," (Factor ",factor,")", sep = "")) +
      labs(y = "Count", x = label) +
      theme (plot.title = element_text(hjust = 0.5))
    print(histogram)
    dev.off()

    # Boxplot
    png(paste("MOFA/3.downstream_analysis/",groups,"_",samples,"/statistical_tests/features_boxplots/boxplot_",colnames(features.matrix)[i],".png", sep = ""), height=3600, width=6000, res=600, units="px")
    boxplot <- ggplot(features.factor.df, aes(as.factor(cvd), eval(parse(text=feature)))) +
      geom_boxplot() +
      ggtitle(paste("Boxplot for ",colnames(features.matrix)[i]," (Factor ",factor,")", sep = "")) +
      labs(y = label, x = "CVD") +
      theme (plot.title = element_text(hjust = 0.5))
    print(boxplot)
    dev.off()

    # T-test for mean (parametric test)
    t.test.res <- t.test(eval(parse(text=feature)) ~ cvd, data = features.factor.df, paired = FALSE)

    # Mann–Whitney–Wilcoxon for means (non-parametric test)
    wilcox.test.res <- wilcox.test(eval(parse(text=feature)) ~ cvd, data = features.factor.df, paired = FALSE, conf.int = 0.95)

    feature.name <- colnames(features.matrix)[i]
    results.df <- data.frame(feature_name = feature.name, factor = factor, t_test = t.test.res$p.value,
                             m_whitney_wilcoxon = wilcox.test.res$p.value)
    heatmap.features.tests <- rbind(heatmap.features.tests, results.df)

  }
  return(heatmap.features.tests)
}

cat("\n\n... 4.2.1) METHYLATION ...\n\n", file = stdout())

# Factor 11

factor <- 11
meth.30features.heatmap.name <- load(file = paste("MOFA/3.downstream_analysis/",groups,"_",samples,"/variance_decomposition/meth_30features_heatmap.RData", sep = ""))
meth.30features.heatmap <- get(meth.30features.heatmap.name)
cat("\n// Top 30 features from Factor of interest heatmap for Methylation //\n", file = stdout())
meth.30features.heatmap[[factor]]

cat("\nM matrix subset for that specific features\n", file = stdout())
MOFA.data <- get_data(MOFA.model)
M.matrix <- t(MOFA.data$meth$group1[rownames(MOFA.data$meth$group1)%in%meth.30features.heatmap[[factor]],])
dim(M.matrix)
M.matrix[1:5,1:5]

cat("\n\nPerforming histogram, boxplot, t-test and Mann-Whitney-Wilcoxon test for each feature\n\n", file = stdout())
meth.heatmap.features.tests.1 <- testFeaturesFactor(features.matrix = M.matrix, factor = factor, omic = "meth")

cat("\nResults for methylation features tests\n", file = stdout())
meth.heatmap.features.tests.1

cat("\nCorrelations between methylation features\n", file = stdout())
corrFeatures(M.matrix, factor = factor, omic = "meth")

# Factor 23

factor <- 23
meth.30features.heatmap.name <- load(file = paste("MOFA/3.downstream_analysis/",groups,"_",samples,"/variance_decomposition/meth_30features_heatmap.RData", sep = ""))
meth.30features.heatmap <- get(meth.30features.heatmap.name)
cat("\n// Top 20 features from Factor of interest heatmap for Methylation //\n", file = stdout())
meth.30features.heatmap[[factor]]

cat("\nM matrix subset for that specific features\n", file = stdout())
MOFA.data <- get_data(MOFA.model)
M.matrix <- t(MOFA.data$meth$group1[rownames(MOFA.data$meth$group1)%in%meth.30features.heatmap[[factor]],])
dim(M.matrix)
M.matrix[1:5,1:5]

cat("\n\nPerforming histogram, boxplot, t-test and Mann-Whitney-Wilcoxon test for each feature\n\n", file = stdout())
meth.heatmap.features.tests.2 <- testFeaturesFactor(features.matrix = M.matrix, factor = factor, omic = "meth")

cat("\nResults for methylation features tests\n", file = stdout())
meth.heatmap.features.tests.2

cat("\nCorrelations between methylation features\n", file = stdout())
corrFeatures(M.matrix, factor = factor, omic = "meth")


cat("\n\n... 4.2.2) MRNA ...\n\n", file = stdout())

# Factor 11

factor <- 11
mRNA.20features.heatmap.name <- load(file = paste("MOFA/3.downstream_analysis/",groups,"_",samples,"/variance_decomposition/mRNA_30features_heatmap.RData", sep = ""))
mRNA.20features.heatmap <- get(mRNA.20features.heatmap.name)
cat("\n// Top 30 features from Factor 11 heatmap for mRNA //\n", file = stdout())
mRNA.20features.heatmap[[factor]]

cat("\nGene expression matrix subset for that specific features\n", file = stdout())
MOFA.data <- get_data(MOFA.model)
mRNA.matrix <- t(MOFA.data$mRNA$group1[rownames(MOFA.data$mRNA$group1)%in%mRNA.20features.heatmap[[factor]],])
dim(mRNA.matrix)
mRNA.matrix[1:5,1:5]

cat("\n\nPerforming histogram, boxplot, t-test and Mann-Whitney-Wilcoxon test for each feature\n\n", file = stdout())
mRNA.heatmap.features.tests.1 <- testFeaturesFactor(features.matrix = mRNA.matrix, factor = factor, omic = "mRNA")

cat("\nResults for mRNA features\n", file = stdout())
mRNA.heatmap.features.tests.1

cat("\nCorrelations between mRNA features\n", file = stdout())
corrFeatures(M.matrix, factor = factor, omic = "mRNA")

# Factor 23

factor <- 23
mRNA.20features.heatmap.name <- load(file = paste("MOFA/3.downstream_analysis/",groups,"_",samples,"/variance_decomposition/mRNA_30features_heatmap.RData", sep = ""))
mRNA.20features.heatmap <- get(mRNA.20features.heatmap.name)
cat("\n// Top 30 features from Factor 11 heatmap for mRNA //\n", file = stdout())
mRNA.20features.heatmap[[factor]]

cat("\nGene expression matrix subset for that specific features\n", file = stdout())
MOFA.data <- get_data(MOFA.model)
mRNA.matrix <- t(MOFA.data$mRNA$group1[rownames(MOFA.data$mRNA$group1)%in%mRNA.20features.heatmap[[factor]],])
dim(mRNA.matrix)
mRNA.matrix[1:5,1:5]

cat("\n\nPerforming histogram, boxplot, t-test and Mann-Whitney-Wilcoxon test for each feature\n\n", file = stdout())
mRNA.heatmap.features.tests.2 <- testFeaturesFactor(features.matrix = mRNA.matrix, factor = factor, omic = "mRNA")

cat("\nResults for mRNA features\n", file = stdout())
mRNA.heatmap.features.tests.2

cat("\nCorrelations between mRNA features\n", file = stdout())
corrFeatures(M.matrix, factor = factor, omic = "mRNA")

# # Factor 23
#
# cat("\n\n... 4.2.2) MRNA ...\n\n", file = stdout())
# factor <- 23
# mRNA.20features.heatmap.name <- load(file = paste("MOFA/3.downstream_analysis/",groups,"_",samples,"/variance_decomposition/mRNA_30features_heatmap.RData", sep = ""))
# mRNA.20features.heatmap <- get(mRNA.20features.heatmap.name)
# cat("\n// Top 30 features from Factor 23 heatmap for mRNA //\n", file = stdout())
# mRNA.20features.heatmap[[factor]]
#
# cat("\nGene expression matrix subset for that specific features\n", file = stdout())
# MOFA.data <- get_data(MOFA.model)
# mRNA.matrix <- t(MOFA.data$mRNA$group1[rownames(MOFA.data$mRNA$group1)%in%mRNA.20features.heatmap[[factor]],])
# dim(mRNA.matrix)
# mRNA.matrix[1:5,1:5]
#
# cat("\n\nPerforming histogram, boxplot, t-test and Mann-Whitney-Wilcoxon test for each feature\n\n", file = stdout())
# mRNA.heatmap.features.tests.2 <- data.frame()
#
# for (i in 1:length(colnames(mRNA.matrix))) {
#
#   factor.transcript.cvd <- data.frame(transcript = mRNA.matrix[,i], cvd = as.factor(MOFA.model@samples_metadata$cvd),
#                                       row.names = rownames(mRNA.matrix))
#
#   png(paste("MOFA/3.downstream_analysis/",groups,"_",samples,"/statistical_tests/features_histograms/hist_",colnames(mRNA.matrix)[i],".png", sep = ""), height=3600, width=6000, res=600, units="px")
#   histogram <- ggplot(factor.transcript.cvd, aes(transcript)) +
#     geom_histogram(fill = "white", color = "grey30", bins = 30) +
#     ggtitle(paste("Histogram for ",colnames(mRNA.matrix)[i], sep = "")) +
#     labs(y="Count", x = "Transcript expression") +
#     theme (plot.title = element_text(hjust = 0.5))
#   print(histogram)
#   dev.off()
#
#   png(paste("MOFA/3.downstream_analysis/",groups,"_",samples,"/statistical_tests/features_boxplots/boxplot_",colnames(mRNA.matrix)[i],".png", sep = ""), height=3600, width=6000, res=600, units="px")
#   boxplot <- ggplot(factor.transcript.cvd, aes(as.factor(cvd), transcript)) +
#     geom_boxplot() +
#     ggtitle(paste("Boxplot for ",colnames(mRNA.matrix)[i], sep = "")) +
#     labs(y="Transcript expression", x = "CVD") +
#     theme (plot.title = element_text(hjust = 0.5))
#   print(boxplot)
#   dev.off()
#
#   # T-test for mean (parametric test)
#   t.test.res <- t.test(transcript ~ cvd, data = factor.transcript.cvd, paired = FALSE)
#
#   # Mann–Whitney–Wilcoxon for means (non-parametric test)
#   wilcox.test.res <- wilcox.test(transcript ~ cvd, data = factor.transcript.cvd, paired = FALSE, conf.int = 0.95)
#
#   transcript.name <- colnames(mRNA.matrix)[i]
#   results.df <- data.frame(feature_name = transcript.name, factor = 12, t_test = t.test.res$p.value,
#                            m_whitney_wilcoxon = wilcox.test.res$p.value)
#   print(results.df)
#   mRNA.heatmap.features.tests.2 <- rbind(mRNA.heatmap.features.tests.2, results.df)
# }
#
# cat("\nResults for mRNA features\n", file = stdout())
# mRNA.heatmap.features.tests.2


a <- read.csv("MOFA/3.downstream_analysis/nogroups_914_2055/statistical_tests/heatmap_features_tests.csv")

head(a)

# p-values FDR adjust
meth.heatmap.features.tests.1$t_test_adj <- p.adjust(meth.heatmap.features.tests.1$t_test, method = "fdr")
meth.heatmap.features.tests.1$m_whitney_wilcoxon_adj <- p.adjust(meth.heatmap.features.tests.1$m_whitney_wilcoxon, method = "fdr")
meth.heatmap.features.tests.2$t_test_adj <- p.adjust(meth.heatmap.features.tests.2$t_test, method = "fdr")
meth.heatmap.features.tests.2$m_whitney_wilcoxon_adj <- p.adjust(meth.heatmap.features.tests.2$m_whitney_wilcoxon, method = "fdr")
mRNA.heatmap.features.tests.1$t_test_adj <- p.adjust(mRNA.heatmap.features.tests.1$t_test, method = "fdr")
mRNA.heatmap.features.tests.1$m_whitney_wilcoxon_adj <- p.adjust(mRNA.heatmap.features.tests.1$m_whitney_wilcoxon, method = "fdr")
mRNA.heatmap.features.tests.2$t_test_adj <- p.adjust(mRNA.heatmap.features.tests.2$t_test, method = "fdr")
mRNA.heatmap.features.tests.2$m_whitney_wilcoxon_adj <- p.adjust(mRNA.heatmap.features.tests.2$m_whitney_wilcoxon, method = "fdr")

total.heatmap.features.tests <- rbind(meth.heatmap.features.tests.1, meth.heatmap.features.tests.2,
                                      mRNA.heatmap.features.tests.1, mRNA.heatmap.features.tests.2)

write.csv(total.heatmap.features.tests, file = paste("MOFA/3.downstream_analysis/",groups,"_",samples,"/statistical_tests/heatmap_features_tests.csv", sep = ""),
          col.names = TRUE, row.names = FALSE, sep = ",")


cat("\n\n... 5) STATISTICAL TESTS BETWEEN INDIVIDUALS WITHIN CVD SUBGROUPS ... \n\n", file = stdout())

# We make the CVD subgroup variable based on the factor 11 only
factor <- 11

# Subset of samples_metadata to match only individuals with CVD
MOFA.cov.cvd <- MOFA.model@samples_metadata[MOFA.model@samples_metadata$cvd == 1,]

cat("\nAll factor values of the interested factor\n", file = stdout())
factors <- get_factors(MOFA.model,
                       groups = "all",
                       factors = "all"
)
lapply(factors,dim)
lapply(factors,head)

# Subset of factor values from the individuals with CVD
factors.cvd <- factors$group1[rownames(factors$group1)%in%MOFA.cov.cvd$sample,]
# Same order of sample ID from factors and samples_metadata
sample.ids.index <- match(MOFA.cov.cvd$sample, rownames(factors.cvd))
factors.cvd <- factors.cvd[sample.ids.index,]
head(factors.cvd)
dim(factors.cvd)

cat("\nDefine the two groups that split the CVD group into subgroups to a new covariable using the factor 11 values\n", file = stdout())
cat("\nFactor of interest values for CVD group \n", file = stdout())
summary(factors.cvd[,factor])
# min -0.9, max 7.8, let's use a threshold of median
table(factors.cvd[,factor] > summary(factors.cvd[,factor])[3])
# Convert TRUE to 1 and FALSE to 0
cvd.subgroups <- ifelse(factors.cvd[,factor] > summary(factors.cvd[,factor])[3], 1,0)
cat("\nNew 2 subgroups from within CVD group\n", file = stdout())
table(cvd.subgroups)

# Add the new subgrouping variable
MOFA.cov.cvd$cvd_subgroups <- cvd.subgroups
# Add another covariable: CHD incidence
cov <- read.table("/projects/regicor/data/FHS/phenotype/phs000007.v29.pht003316.v6.p10.c1.vr_survcvd_2014_a_1023s.HMB-IRB-MDS.txt", header=T, stringsAsFactors=F, sep="\t")
head(cov)
dim(cov)
chd.cov <- cov[cov$shareid%in%MOFA.cov.cvd$sample,c("shareid","chd")]
colnames(chd.cov)[1] <- "sample"
MOFA.cov.cvd <- merge(MOFA.cov.cvd, chd.cov, by = "sample")

# cat("\n\n... 5.1) T-TEST FOR INTERESTED FACTORS FOR CHD ... \n\n", file = stdout())
# 
# all.factors.res <- c()
# for (i in 1:30) {
#   test.factor.res <- testFactor(factors = factors.cvd, factor = i, covariable = "chd", group = "cvd", MOFA_samples_metadata = MOFA.cov.cvd)
#   all.factors.res <- rbind(all.factors.res, test.factor.res)
# }
# 
# write.csv(all.factors.res, file = paste("MOFA/3.downstream_analysis/",groups,"_",samples,"/statistical_tests/factors_chd_cvd_group_tests.csv", sep = ""),
#           row.names = FALSE)

cat("\n\n... 5.2) T-TEST AND CHI-SQUARE FOR ALL COVARIABLES ... \n\n", file = stdout())

testCov <- function(covariable, type, MOFA_samples_metadata, factor){
  
  if (type == "categorical") {
    
    if (covariable == "smoke") {
      MOFA_samples_metadata <- MOFA_samples_metadata[!MOFA_samples_metadata$smoke == "f1-5",]
    }
    
    test <- "chi_square"
    
    # Change class to factor
    MOFA_samples_metadata[,which(colnames(MOFA_samples_metadata) == covariable)] <- as.factor(eval(parse(text=paste("MOFA_samples_metadata$",covariable,sep=""))))
    print(class(eval(parse(text=paste("MOFA_samples_metadata$",covariable,sep="")))))
    
    # Check barplot
    png(paste("MOFA/3.downstream_analysis/",groups,"_",samples,"/statistical_tests/cov_plots/",covariable,"_cvd_subgroups_barplot.png", sep = ""), height=3600, width=6000, res=600, units="px")
    table <- table(as.factor(eval(parse(text=paste("MOFA_samples_metadata$",covariable,sep="")))), as.factor(eval(parse(text="MOFA_samples_metadata$cvd_subgroups"))))
    barplot(table,
            main = paste("Number of individuals by ",covariable, " and CVD subgroups", sep = ""),
            xlab = "CVD subgroups",
            beside = TRUE,
            col = 1:length(rownames(table)))  
    legend("topright", 
           cex = 0.75,
           fill = 1:length(rownames(table)), ncol = 2, legend = rownames(table))
    dev.off()
    
    chisq.test.res <- chisq.test(table)
    mean_group0 <- "NULL"
    mean_group1 <- "NULL"
    p_val <- chisq.test.res$p.value
  }
  else if (type == "continous") {
    
    test <- "t_test"
    # Change class to numeric
    MOFA_samples_metadata[,which(colnames(MOFA_samples_metadata) == covariable)] <- as.numeric(eval(parse(text=paste("MOFA_samples_metadata$",covariable,sep=""))))
    print(class(eval(parse(text=paste("MOFA_samples_metadata$",covariable,sep="")))))
    
    # Check histogram
    png(paste("MOFA/3.downstream_analysis/",groups,"_",samples,"/statistical_tests/cov_plots/",covariable,"_cvd_subgroups_histogram.png", sep = ""), height=3600, width=6000, res=600, units="px")
    histogram <- ggplot(MOFA_samples_metadata, aes(eval(parse(text=covariable)))) +
      geom_histogram(fill = "white", color = "grey30", bins = 30) +
      ggtitle(paste("Histogram for factor ",factor, sep = "")) +
      labs(y="Count", x = covariable) +
      theme (plot.title = element_text(hjust = 0.5))
    print(histogram)
    dev.off()
    
    # Check boxplot
    png(paste("MOFA/3.downstream_analysis/",groups,"_",samples,"/statistical_tests/cov_plots/",covariable,"_cvd_subgroups_boxplot.png", sep = ""), height=3600, width=6000, res=600, units="px")
    boxplot <- ggplot(MOFA_samples_metadata, aes(as.factor(cvd_subgroups) , eval(parse(text=covariable)))) +
      geom_boxplot() +
      ggtitle(paste("Boxplot for factor ",factor, sep = "")) +
      labs(y=covariable, x = MOFA_samples_metadata$cvd_subgroups) +
      theme (plot.title = element_text(hjust = 0.5))
    print(boxplot)
    dev.off()
    
    # T-test for mean (parametric test)
    t.test.res <- t.test(eval(parse(text=covariable)) ~ as.factor(cvd_subgroups), data = MOFA_samples_metadata, paired = FALSE)
    mean_group0 <- round(as.numeric(t.test.res$estimate)[1],6)
    mean_group1 <- round(as.numeric(t.test.res$estimate)[2],6)
    p_val <- t.test.res$p.value
    
  }

  results.df <- data.frame(covariable = covariable, type = type,
                           mean_group0 = mean_group0,
                           mean_group1 = mean_group1,
                           test = test, p_val = p_val,
                           stringsAsFactors = FALSE)

 return(results.df)
  
}

### NO HAY DOS GRUPOS SEPARADOS EN CVD ASÍ QUE NO PODEMOS HACER ESTO, PERO LO HACEMOS VA

cov.classes <- c("categorical","categorical","categorical","continous","continous","continous","continous",
                 "continous","continous","continous","continous","continous","continous","categorical",
                 "continous","continous","continous","continous","continous","continous","categorical",
                 "continous","categorical","categorical","categorical", "categorical")
all.results.df <- c()

for (iteration in 1:(length(names(MOFA.cov.cvd)))) {

  print(iteration)
  covariable <- names(MOFA.cov.cvd)[iteration]

  if (covariable%in%c("sample","group","cvd","cvd_subgroups")) {
    next
  }
  print(covariable)
  results.df <- testCov(covariable = covariable, type = cov.classes[iteration], MOFA_samples_metadata = MOFA.cov.cvd , factor = 11)
  all.results.df <- rbind(all.results.df,results.df)

}

cat("\n Tests results\n", file = stdout())
all.results.df
write.csv(all.results.df, file = paste("MOFA/3.downstream_analysis/",groups,"_",samples,"/statistical_tests/cov_within_cvd_subgroups_tests.csv", sep = ""),
          row.names = FALSE)
date()

cat("\n\nComputational time\n\n", file = stdout())
end_time <- Sys.time()
end_time - start_time

cat("\n\n ######################## 3. MOFA DOWNSTREAM ANALYSIS: STATISTICAL TESTS SCRIPT ENDS ######################## \n\n", file = stdout())



# pvals <- c(7.16E-05, 4.95E-05, 4.02E-06, 5.40E-06, 0.000796716, 0.000108133, 0.000270022, 0.000115567, 0.000291015, 2.16E-05, 4.83E-06, 1.55E-05, 1.08E-05, 8.25E-05, 1.63E-05, 3.68E-05, 5.67E-06, 4.64E-05, 5.37E-05, 2.64E-05)
# p.adjust(pvals, method = "bonferroni")



