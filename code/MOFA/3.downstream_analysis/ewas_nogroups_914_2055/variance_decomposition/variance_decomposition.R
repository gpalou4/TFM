rm(list=ls())
start_time <- Sys.time()
date()

cat("\n\n ######################## START LOADING LIBRARIES ######################## \n\n", file = stdout())

library(ggplot2)
cat("\n\n ggplot2 LOADED \n\n", file = stdout())
library(MOFA2)
cat("\n\n MOFA LOADED \n\n", file = stdout())


cat("\n\n ######################## FINISHED LOADING LIBRARIES ######################## \n\n", file = stdout())

cat("\n\n ######################## 3. MOFA DOWNSTREAM ANALYSIS: VARIANCE DECOMPOSITION SCRIPT BEGINS ######################## \n\n", file = stdout())

cat("\n\n... SETTING WORKING DIRECTORY ...\n\n", file = stdout())
setwd("/projects/regicor/Guille/TFM")
getwd()

cat("\n\n... LOADING MOFA TRAINED OBJECT ...\n\n", file = stdout())

groups <- "nogroups"
samples <- "914_2055"
# interested factors
factors <- c(9,19,21,27)

MOFA.model <- load_model(paste("MOFA/2.training_model/MOFA_trained_ewas_",groups,"_",samples,".hdf5",sep = ""))
MOFA.model

cat("\n\n... LOADING SAMPLES METADATA ...\n\n", file = stdout())

MOFA.covariates <- read.csv(file = paste("MOFA/MOFA_covariates_order_",groups,"_",samples,".csv",sep = ""),
                                  header=TRUE, stringsAsFactors = F, sep=",")
# Metadata has to contain the columns 'sample' (shareid) and 'group' (cvd) for MOFA to work
colnames(MOFA.covariates)[1] <- "sample"
# Change sample names from integers to characters
MOFA.covariates$sample <- as.character(MOFA.covariates$sample)
MOFA.covariates[1:3,]
dim(MOFA.covariates)
# 2055 24


groups <- "ewas_nogroups"
samples <- "914_2055"

cat("\n\n... 2) ADDING SAMPLES METADATA (COVARIABLES) TO THE TRAINED MODEL OBJECT ...\n\n", file = stdout())

samples_metadata(MOFA.model) <- MOFA.covariates[,-2]

cat("\n\n... 3) INSPECTING THE MOFA OBJECT ...\n\n", file = stdout())

cat("\nSlot names of the MOFA model\n", file = stdout())
slotNames(MOFA.model)
cat("\nMetadata stored in the MOFA model\n", file = stdout())
head(MOFA.model@samples_metadata)

cat("\n\n... 4) VARIANCE DECOMPOSITION ...\n\n", file = stdout())

cat("\n\n... 4.1) CORRELATIONS BETWEEN FACTORS AND COVARIABLES ...\n\n", file = stdout())

png(paste("MOFA/3.downstream_analysis/",groups,"_",samples,"/variance_decomposition/factors_correlation.png", sep = ""), height=3600, width=6000, res=600, units="px")
corr <- plot_factor_cor(MOFA.model)
dev.off()

#MOFA.model@samples_metadata$sex[1:10] <- NA

names(MOFA.model@samples_metadata) <- c("sample","CVD","sex","age","total-C","HDL-C","triglycerides","SBP","DBP","glucose","weight",
                                        "height","waist","smoke","CD8T","CD4T","NK","Bcell","monocytes","granulocytes","batch meth","SVA meth","batch mRNA","group")

png(paste("MOFA/3.downstream_analysis/",groups,"_",samples,"/variance_decomposition/cov_factors_correlation_R.png", sep = ""), height=3600, width=6000, res=600, units="px")
corr_factors_cov <- correlate_factors_with_covariates(MOFA.model,
                                  covariates = names(MOFA.model@samples_metadata)[-c(1,24)],
                                  factors = "all",
                                  groups = "all",
                                  abs = FALSE,
                                  plot = c("r"),
                                  alpha = 0.05,
                                  return_data = FALSE,
                                  transpose = FALSE,
                                  tl.srt=45,
                                  number.cex = 0.25, 
                                  tl.cex = 0.75,
                                  cex.main = 0.01
)

# Change names again
samples_metadata(MOFA.model) <- MOFA.covariates[,-2]

corr_factors_cov[1:5,1:5]
dev.off()

cat(paste("\nFactor ",factors[1], "\n", sep = ""), file = stdout())
sort(corr_factors_cov[factors[1],])
cat(paste("\nFactor ",factors[2], "\n", sep = ""), file = stdout())
sort(corr_factors_cov[factors[2],])

### P-values for correlation coefficients

## Get factors
all_factors <- get_factors(MOFA.model,
                           groups = "all",
                           factors = "all"
)

# Z <- get_factors(object, factors = factors, groups = groups, as.data.frame=FALSE)
# Z <- do.call(rbind, Z)
all_factors$group1[1:5,1:5]

## Get covariates

covariates <- MOFA.model@samples_metadata[,-24]

# convert character columns to factors
cols <- which(sapply(covariates,class)=="character")
if (length(cols>=1)) {
  covariates[cols] <- lapply(covariates[cols], as.factor)
}

# convert all columns to numeric
cols <- which(sapply(covariates,class)!="numeric")
if (length(cols>=1)) {
  cols.factor <- which(sapply(covariates,class)=="factor")
  covariates[cols] <- lapply(covariates[cols], as.numeric)
  warning("There are non-numeric values in the covariates data.frame, converting to numeric...")
  covariates[cols] <- lapply(covariates[cols], as.numeric)
}
stopifnot(all(sapply(covariates,class)=="numeric"))

cor <- psych::corr.test(all_factors$group1, covariates, method = "pearson", adjust = "BH")
cat("\n\n Correlation coefficients \n\n", file = stdout())
cor$r
cat("\n\n P-values for the correlation coefficients \n\n", file = stdout())
cor$p

# DA ERROR

# png(paste("MOFA/3.downstream_analysis/",groups,"_",samples,"/variance_decomposition/cov_factors_correlation_heatmap.png", sep = ""), height=3600, width=6000, res=600, units="px")
# correlate_factors_with_covariates(MOFA.model,
#                                   covariates = names(MOFA.model@samples_metadata)[-24],
#                                   factors = "all",
#                                   groups = "all",
#                                   abs = FALSE,
#                                   plot = c("log_pval"),
#                                   alpha = 0.05,
#                                   return_data = FALSE,
#                                   transpose = FALSE
# )
# dev.off()

cat("\n\n... 4.2) VARIANCE EXPLAINED ...\n\n", file = stdout())

cat("\n\n Plot variance decomposition \n\n", file = stdout())
png(paste("MOFA/3.downstream_analysis/",groups,"_",samples,"/variance_decomposition/variance_decomposition_1.png", sep = ""), height=3600, width=6000, res=600, units="px")
plot.variance.1 <- plot_variance_explained(MOFA.model, x="group", y="factor", xlab = "")
plot.variance.1
dev.off()
plot.variance.1$data
cat(paste("\nVariance explained for each omic by Factor ",factors[1], "\n\n", sep = ""), file = stdout())
factor1.var <- plot.variance.1$data[plot.variance.1$data$factor == paste("Factor",factors[1],sep=""),]
factor1.var
sum(factor1.var[,"value"])

cat(paste("\nVariance explained for each omic by Factor ",factors[2], "\n\n", sep = ""), file = stdout())
factor2.var <- plot.variance.1$data[plot.variance.1$data$factor == paste("Factor",factors[2],sep=""),]
factor2.var
sum(factor2.var[,"value"])

png(paste("MOFA/3.downstream_analysis/",groups,"_",samples,"/variance_decomposition/variance_decomposition_2.png", sep = ""), height=3600, width=6000, res=600, units="px")
p <- plot_variance_explained(MOFA.model, x="view", y="group")
p + theme(axis.text.x = element_text(color="black", angle=40, vjust=1, hjust=1))
dev.off()

png(paste("MOFA/3.downstream_analysis/",groups,"_",samples,"/variance_decomposition/variance_decomposition_3.png", sep = ""), height=3600, width=6000, res=600, units="px")
plot.variance.3 <- plot_variance_explained(MOFA.model, x="group", y="factor", plot_total = T)
plot.variance.3[[2]]
# cat(paste("\nTotal variance explained for each omic by Factor \n",factors[2], sep = ""), file = stdout())
plot.variance.3[[2]]$data
dev.off()

cat("\n\n... 4.3) CLUSTERING FOR EACH FACTOR ...\n\n", file = stdout())

cat("\n\n... 4.3.1) ALL FACTORS ...\n\n", file = stdout())

cat("\nChecking an outlier from factor 1...\n", file = stdout())

# factor.values <- all_factors$group1[,1]
# summary(factor.values)
# factor.values[which(factor.values > 15 | factor.values < -15)]
# outliers1 <- which(factor.values > 15 | factor.values < -15 )
# 
# cat("\nChecking outliers from factor 6...\n", file = stdout())
# 
# factor.values <- all_factors$group1[,6]
# summary(factor.values)
# factor.values[which(factor.values > 15 | factor.values < -15)]
# outliers2 <- which(factor.values > 15 | factor.values < -15 )
# 
# cat("\nChecking outliers from factor 15...\n", file = stdout())
# 
# factor.values <- all_factors$group1[,15]
# summary(factor.values)
# factor.values[which(factor.values == max(factor.values))]
# outliers3 <- which(factor.values == max(factor.values))
# 
# outliers <- as.vector(c(outliers1,outliers2,outliers3))
# cat("\n...Outliers...\n", file = stdout())
# outliers
# 
# #samples.filter <- samples.filter[which(factor.values != max(factor.values))]
# 
# samples.filter <- rownames(all_factors$group1)[-outliers]
# length(samples.filter)
# MOFA.model.filt <- subset_samples(MOFA.model, samples.filter)
# MOFA.model.filt
#
#MOFA.model.filt@samples_metadata$cvd[1:10] <- NA
# 
# 
# png(paste("MOFA/3.downstream_analysis/",groups,"_",samples,"/variance_decomposition/factors1-15_sex_cvd.png", sep = ""), height=3600, width=6000, res=600, units="px")
# p <- plot_factor(MOFA.model.filt,
#             factor = 1:15,
#             color_by = "cvd",
#             shape_by = "sex"
# )
# p + theme(axis.title.x=element_blank(),
#           axis.text.x=element_blank(),
#           axis.ticks.x=element_blank())
# dev.off()
# 
# png(paste("MOFA/3.downstream_analysis/",groups,"_",samples,"/variance_decomposition/factors16-30_sex_cvd.png", sep = ""), height=3600, width=6000, res=600, units="px")
# p <- plot_factor(MOFA.model.filt,
#             factor = 16:30,
#             color_by = "cvd",
#             shape_by = "sex"
# )
# p + theme(axis.title.x=element_blank(),
#           axis.text.x=element_blank(),
#           axis.ticks.x=element_blank())
# dev.off()
# # For the violin plot you cannot use quantitative variables as it makes the plot for each group
# 
# png(paste("MOFA/3.downstream_analysis/",groups,"_",samples,"/variance_decomposition/factors1-15_cvd_violin.png", sep = ""), height=3600, width=6000, res=600, units="px")
# p <- plot_factor(MOFA.model.filt,
#             factor = 1:15,
#             color_by = "cvd",
#             dot_size = 1.3,      # change dot size
#             dodge = TRUE,           # dodge points with different colors
#             legend = TRUE,          # remove legend
#             add_violin = TRUE,      # add violin plots,
#             color_violin = TRUE,
#             scale = TRUE
#             #violin_alpha = 0.01  # transparency of violin plots
# ) # + scale_color_manual(values=colors) + scale_fill_manual(values=colors)
# p + theme(axis.title.x=element_blank(),
#           axis.text.x=element_blank(),
#           axis.ticks.x=element_blank())
# dev.off()
# 
# png(paste("MOFA/3.downstream_analysis/",groups,"_",samples,"/variance_decomposition/factors16-30_cvd_violin.png", sep = ""), height=3600, width=6000, res=600, units="px")
# p <- plot_factor(MOFA.model.filt,
#             factor = 16:30,
#             color_by = "cvd",
#             dot_size = 1.3,      # change dot size
#             dodge = TRUE,           # dodge points with different colors
#             legend = TRUE,          # remove legend
#             add_violin = TRUE,      # add violin plots,
#             color_violin = TRUE,
#             scale = TRUE
#             #violin_alpha = 0.01  # transparency of violin plots
# ) # + scale_color_manual(values=colors) + scale_fill_manual(values=colors)
# p + theme(axis.title.x=element_blank(),
#           axis.text.x=element_blank(),
#           axis.ticks.x=element_blank())
# dev.off()
# 
# # png("MOFA/3.downstream_analysis/factor1_factor2.png", sep = ""), height=3600, width=6000, res=600, units="px")
# # plot_factors(MOFA.model,
# #              factors = c(1,2),
# #              color_by = "cvd"
# # )
# # dev.off()
# 
# cat("\n\n... 4.3.2) SPECIFIC FACTORS COLORED BY COVARIABLES ...\n\n", file = stdout())
# 
# library(ggpubr)
# library(ggsignif)
# 
# 
# ## Violin for factor 9
# 
# png(paste("MOFA/3.downstream_analysis/",groups,"_",samples,"/variance_decomposition/factor9_cvd_violin.png", sep = ""), height=3600, width=6000, res=600, units="px")
# p <- plot_factor(MOFA.model,
#             factor = 9,
#             color_by = "cvd",
#             dot_size = 1.3,      # change dot size
#             dodge = TRUE,           # dodge points with different colors
#             legend = TRUE,          # remove legend
#             add_violin = TRUE,      # add violin plots,
#             color_violin = TRUE,
#             scale = TRUE
#             #violin_alpha = 0.01  # transparency of violin plots
# ) # + scale_color_manual(values=colors) + scale_fill_manual(values=colors)
# p + theme(axis.title.x=element_blank(),
#           axis.text.x=element_blank(),
#           axis.ticks.x=element_blank())+
#   stat_compare_means(method = "t.test") +
#   geom_boxplot(alpha=0.5, position=position_dodge(width=1), show.legend = FALSE, width = 0.3, outlier.shape = NA)
# dev.off()
# 
# ## Violin for factor 19
# 
# png(paste("MOFA/3.downstream_analysis/",groups,"_",samples,"/variance_decomposition/factor19_cvd_violin.png", sep = ""), height=3600, width=6000, res=600, units="px")
# p <- plot_factor(MOFA.model,
#             factor = 19,
#             color_by = "cvd",
#             dot_size = 1.3,      # change dot size
#             dodge = TRUE,           # dodge points with different colors
#             legend = TRUE,          # remove legend
#             add_violin = TRUE,      # add violin plots,
#             color_violin = TRUE,
#             scale = TRUE
#             #violin_alpha = 0.01  # transparency of violin plots
# ) # + scale_color_manual(values=colors) + scale_fill_manual(values=colors)
# #str(p)
# p + theme(axis.title.x=element_blank(),
#           axis.text.x=element_blank(),
#           axis.ticks.x=element_blank()) +
#   stat_compare_means(method = "t.test") + 
#   geom_boxplot(alpha=0.5, position=position_dodge(width=1), show.legend = FALSE, width = 0.3, outlier.shape = NA)
# dev.off()
# 
# ## Violin for factor 21
# 
# png(paste("MOFA/3.downstream_analysis/",groups,"_",samples,"/variance_decomposition/factor21_cvd_violin.png", sep = ""), height=3600, width=6000, res=600, units="px")
# p <- plot_factor(MOFA.model,
#             factor = 21,
#             color_by = "cvd",
#             dot_size = 1.3,      # change dot size
#             dodge = TRUE,           # dodge points with different colors
#             legend = TRUE,          # remove legend
#             add_violin = TRUE,      # add violin plots,
#             color_violin = TRUE,
#             scale = TRUE
#             #violin_alpha = 0.01  # transparency of violin plots
# ) # + scale_color_manual(values=colors) + scale_fill_manual(values=colors)
# p + theme(axis.title.x=element_blank(),
#           axis.text.x=element_blank(),
#           axis.ticks.x=element_blank()) +
#   stat_compare_means(method = "t.test") +
#   geom_boxplot(alpha=0.5, position=position_dodge(width=1), show.legend = FALSE, width = 0.3, outlier.shape = NA)
# dev.off()
# 
# ## Violin for factor 27
# 
# png(paste("MOFA/3.downstream_analysis/",groups,"_",samples,"/variance_decomposition/factor27_cvd_violin.png", sep = ""), height=3600, width=6000, res=600, units="px")
# p <- plot_factor(MOFA.model,
#             factor = 27,
#             color_by = "cvd",
#             dot_size = 1.3,      # change dot size
#             dodge = TRUE,           # dodge points with different colors
#             legend = TRUE,          # remove legend
#             add_violin = TRUE,      # add violin plots,
#             color_violin = TRUE,
#             scale = TRUE
#             #violin_alpha = 0.01  # transparency of violin plots
# ) # + scale_color_manual(values=colors) + scale_fill_manual(values=colors)
# p + theme(axis.title.x=element_blank(),
#           axis.text.x=element_blank(),
#           axis.ticks.x=element_blank())+
#   stat_compare_means(method = "t.test") +
#   geom_boxplot(alpha=0.5, position=position_dodge(width=1), show.legend = FALSE, width = 0.3, outlier.shape = NA)
# dev.off()
# 
# ##Scatterplots
# 
# #library(ggpmisc)
# library(ggpubr)
# 
# scatterplot_factor_cov <- function(factor,covariate){
#   
#   if (factor == 27) {
#     factor.values <- all_factors$group1[,27]
#     # print("head of factor.values")
#     # print(head(factor.values))
#     outliers <- which(factor.values > 6 | factor.values < -6 )
#     # print("outliers")
#     # print(outliers)
#     # print("numeric outliers")
#     # print(as.numeric(outliers))
#     # print("printing factor value of 187")
#     # print(factor.values[187])
#     # print("printing factor value of 14549")
#     # print(factor.values[14549])
#     factor.values <- factor.values[-as.numeric(outliers)]
#     samples.filter <- rownames(all_factors$group1)[outliers]
#     # print("samples.filter")
#     # print(as.numeric(samples.filter))
#     # print(length(factor.values))
#     # print("covariates$sample")
#     # print(head(covariates$sample))
#     # print("covariates$sample%in%as.numeric(samples.filter)")
#     # print(table(covariates$sample%in%as.numeric(samples.filter)))
#     covariates <- covariates[-as.numeric(outliers),]
#     # print(dim(covariates))
#   }
#   else {
#   factor.values <- all_factors$group1[,factors[factor]]
#   }
#   
#   df.scatterp <- data.frame(factor_values = factor.values, covariate = covariates[,covariate])
#   png(paste("MOFA/3.downstream_analysis/",groups,"_",samples,"/variance_decomposition/factor_cov/factor",factors[factor],"_",covariate,"_scatterplot.png", sep = ""), height=3600, width=6000, res=600, units="px")
#   sp <- ggscatter(df.scatterp, x = "factor_values", y = "covariate", add = "reg.line",
#                   add.params = list(color = "blue", fill = "lightgray"),
#                   conf.int = TRUE,
#                   shape = 1,
#                   ellipse.alpha = 0.3)
#   print(sp + stat_cor( aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
#                  label.x = 0 ) + 
#     labs(y=paste(covariate, sep = ""), x = paste("Factor", factors[factor], sep ="")))
#   
#   dev.off()
# }

factors <- 1:30

cat("\nFACTOR 9\n", file = stdout())

factor <- 9
covariates.f9 <- c("age", "tot_chol","hdl_chol","glucose","sbp","dbp","CD4T",
                "CD8T","NK","Bcell")
covariates.f9.qual <- c("cvd","sex")

# for (covariate in covariates.f9) {
#   scatterplot_factor_cov(factor = factor, covariate = covariate)
# }

cat("\nFACTOR 19\n", file = stdout())

factor <- 19
covariates.f19 <- c("age","hdl_chol","weight","height","sbp","dbp","Bcell")
covariates.f19.qual <- c("cvd","sex","smoke")

# for (covariate in covariates.f19) {
#   scatterplot_factor_cov(factor = factor, covariate = covariate)
# }

cat("\nFACTOR 21\n", file = stdout())

factor <- 21
covariates.f21 <- c("hdl_chol","tot_chol","weight","height","CD8T")
covariates.f21.qual <- c("cvd","sex")

# for (covariate in covariates.f21) {
#   scatterplot_factor_cov(factor = factor, covariate = covariate)
# }

cat("\nFACTOR 27\n", file = stdout())

cat("\nChecking outliers from factor 27...\n", file = stdout())

# factor.values <- all_factors$group1[,27]
# summary(factor.values)
# factor.values[which(factor.values > 6 | factor.values < -6)]
# outliers <- which(factor.values > 6 | factor.values < -6 )
# cat("\n...Outliers...\n", file = stdout())
# outliers

# samples.filter <- rownames(factors$group1)[-outliers]
# length(samples.filter)
# MOFA.model.filt <- subset_samples(MOFA.model, samples.filter)
# MOFA.model.filt

factor <- 27
covariates.f27 <- c("age","Gran","Bcell","Mono")
covariates.f27.qual <- c("cvd","sex")

# for (covariate in covariates.f27) {
#   scatterplot_factor_cov(factor = factor, covariate = covariate)
# }

# # cat("\n\n... 4.4) VISUALIZATION OF COMBINATION OF FACTORS ...\n\n", file = stdout())
# #
# # png(paste("MOFA/3.downstream_analysis/",groups,"_",samples,"/variance_decomposition/factors1-30_pairs.png", sep = ""), height=3600, width=6000, res=600, units="px")
# # plot_factors(MOFA.model,
# #              factors = 1:30,
# #              color_by = "cvd"
# # )
# # dev.off()

plotFactorCov <- function(MOFA_model,factor,color_by,shape_by,corr_factor_cov, corr_factor_cov_pval){

  p <- plot_factor(MOFA_model,
                   factor = factor,
                   color_by = color_by,
                   shape_by = shape_by,
                   dodge = TRUE
  ) + scale_shape_manual(values=c(21,25))
  print(p)
  png(paste("MOFA/3.downstream_analysis/",groups,"_",samples,"/variance_decomposition/factor_cov/factor",factor,"_",color_by,".png", sep = ""), height=3600, width=6000, res=600, units="px")
  print(p + ggtitle(paste("R² = ",corr_factor_cov*corr_factor_cov, ", p = ",corr_factor_cov_pval,sep="")) +
          theme(plot.title = element_text(hjust = 0.5),
                axis.title.x=element_blank(),
                axis.text.x=element_blank(),
                axis.ticks.x=element_blank()) +
          labs(fill = color_by))
  dev.off()
}
# 
# cat("\nFACTOR 9\n", file = stdout())
# 
# shape_by <- NULL
# 
# colnames(corr_factors_cov) <- names(MOFA.model@samples_metadata)[-c(1,24)]
# 
# factor <- 9
# color_by <- "cvd"
# corr_factor_cov <- round(corr_factors_cov[factor,color_by],3)
# corr_factor_cov_pval <- signif(cor$p[factor,color_by], digits = 2)
# sort(corr_factors_cov[factor,])
# print("Correlation p-value")
# print(corr_factor_cov_pval)
# plotFactorCov(MOFA_model = MOFA.model, factor = factor, color_by = color_by, shape_by = shape_by, corr_factor_cov = corr_factor_cov, corr_factor_cov_pval = corr_factor_cov_pval)
# 
# color_by <- "sex"
# corr_factor_cov <- round(corr_factors_cov[factor,color_by],3)
# corr_factor_cov_pval <- signif(cor$p[factor,color_by], digits = 2)
# plotFactorCov(MOFA_model = MOFA.model, factor = factor, color_by = color_by, shape_by = shape_by, corr_factor_cov = corr_factor_cov, corr_factor_cov_pval = corr_factor_cov_pval)
# 
# # shape_by <- "cvd"
# color_by <- "age"
# corr_factor_cov <- round(corr_factors_cov[factor,color_by],3)
# plotFactorCov(MOFA_model = MOFA.model, factor = factor, color_by = color_by, shape_by = shape_by, corr_factor_cov = corr_factor_cov, corr_factor_cov_pval = corr_factor_cov_pval)
# 
# color_by <- "tot_chol"
# corr_factor_cov <- round(corr_factors_cov[factor,color_by],3)
# plotFactorCov(MOFA_model = MOFA.model, factor = factor, color_by = color_by, shape_by = shape_by, corr_factor_cov = corr_factor_cov, corr_factor_cov_pval = corr_factor_cov_pval)
# 
# color_by <- "hdl_chol"
# corr_factor_cov <- round(corr_factors_cov[factor,color_by],3)
# plotFactorCov(MOFA_model = MOFA.model, factor = factor, color_by = color_by, shape_by = shape_by, corr_factor_cov = corr_factor_cov, corr_factor_cov_pval = corr_factor_cov_pval)
# 
# color_by <- "glucose"
# corr_factor_cov <- round(corr_factors_cov[factor,color_by],3)
# plotFactorCov(MOFA_model = MOFA.model, factor = factor, color_by = color_by, shape_by = shape_by, corr_factor_cov = corr_factor_cov, corr_factor_cov_pval = corr_factor_cov_pval)
# 
# color_by <- "sbp"
# corr_factor_cov <- round(corr_factors_cov[factor,color_by],3)
# plotFactorCov(MOFA_model = MOFA.model, factor = factor, color_by = color_by, shape_by = shape_by, corr_factor_cov = corr_factor_cov, corr_factor_cov_pval = corr_factor_cov_pval)
# 
# color_by <- "dbp"
# corr_factor_cov <- round(corr_factors_cov[factor,color_by],3)
# plotFactorCov(MOFA_model = MOFA.model, factor = factor, color_by = color_by, shape_by = shape_by, corr_factor_cov = corr_factor_cov, corr_factor_cov_pval = corr_factor_cov_pval)
# 
# color_by <- "CD4T"
# corr_factor_cov <- round(corr_factors_cov[factor,color_by],3)
# plotFactorCov(MOFA_model = MOFA.model, factor = factor, color_by = color_by, shape_by = shape_by, corr_factor_cov = corr_factor_cov, corr_factor_cov_pval = corr_factor_cov_pval)
# 
# color_by <- "CD8T"
# corr_factor_cov <- round(corr_factors_cov[factor,color_by],3)
# plotFactorCov(MOFA_model = MOFA.model, factor = factor, color_by = color_by, shape_by = shape_by, corr_factor_cov = corr_factor_cov, corr_factor_cov_pval = corr_factor_cov_pval)
# 
# color_by <- "NK"
# corr_factor_cov <- round(corr_factors_cov[factor,color_by],3)
# plotFactorCov(MOFA_model = MOFA.model, factor = factor, color_by = color_by, shape_by = shape_by, corr_factor_cov = corr_factor_cov, corr_factor_cov_pval = corr_factor_cov_pval)
# 
# color_by <- "Bcell"
# corr_factor_cov <- round(corr_factors_cov[factor,color_by],3)
# plotFactorCov(MOFA_model = MOFA.model, factor = factor, color_by = color_by, shape_by = shape_by, corr_factor_cov = corr_factor_cov, corr_factor_cov_pval = corr_factor_cov_pval)
# 
# 
# cat("\nFACTOR 19\n", file = stdout())
# 
# shape_by <- NULL
# 
# factor <- 19
# color_by <- "cvd"
# corr_factor_cov <- round(corr_factors_cov[factor,color_by],3)
# corr_factor_cov_pval <- signif(cor$p[factor,color_by], digits = 2)
# sort(corr_factors_cov[factor,])
# plotFactorCov(MOFA_model = MOFA.model, factor = factor, color_by = color_by, shape_by = shape_by, corr_factor_cov = corr_factor_cov, corr_factor_cov_pval = corr_factor_cov_pval)
# 
# color_by <- "sex"
# corr_factor_cov <- round(corr_factors_cov[factor,color_by],3)
# corr_factor_cov_pval <- signif(cor$p[factor,color_by], digits = 2)
# plotFactorCov(MOFA_model = MOFA.model, factor = factor, color_by = color_by, shape_by = shape_by, corr_factor_cov = corr_factor_cov, corr_factor_cov_pval = corr_factor_cov_pval)
# 
# color_by <- "smoke"
# corr_factor_cov <- round(corr_factors_cov[factor,color_by],3)
# corr_factor_cov_pval <- signif(cor$p[factor,color_by], digits = 2)
# plotFactorCov(MOFA_model = MOFA.model, factor = factor, color_by = color_by, shape_by = shape_by, corr_factor_cov = corr_factor_cov, corr_factor_cov_pval = corr_factor_cov_pval)

# shape_by <- "cvd"
# color_by <- "age"
# corr_factor_cov <- round(corr_factors_cov[factor,color_by],3)
# plotFactorCov(MOFA_model = MOFA.model, factor = factor, color_by = color_by, shape_by = shape_by, corr_factor_cov = corr_factor_cov, corr_factor_cov_pval = corr_factor_cov_pval)
# 
# color_by <- "hdl_chol"
# corr_factor_cov <- round(corr_factors_cov[factor,color_by],3)
# plotFactorCov(MOFA_model = MOFA.model, factor = factor, color_by = color_by, shape_by = shape_by, corr_factor_cov = corr_factor_cov, corr_factor_cov_pval = corr_factor_cov_pval)
# 
# color_by <- "weight"
# corr_factor_cov <- round(corr_factors_cov[factor,color_by],3)
# plotFactorCov(MOFA_model = MOFA.model, factor = factor, color_by = color_by, shape_by = shape_by, corr_factor_cov = corr_factor_cov, corr_factor_cov_pval = corr_factor_cov_pval)
# 
# color_by <- "height"
# corr_factor_cov <- round(corr_factors_cov[factor,color_by],3)
# plotFactorCov(MOFA_model = MOFA.model, factor = factor, color_by = color_by, shape_by = shape_by, corr_factor_cov = corr_factor_cov, corr_factor_cov_pval = corr_factor_cov_pval)
# 
# color_by <- "smoke"
# corr_factor_cov <- round(corr_factors_cov[factor,color_by],3)
# plotFactorCov(MOFA_model = MOFA.model, factor = factor, color_by = color_by, shape_by = shape_by, corr_factor_cov = corr_factor_cov, corr_factor_cov_pval = corr_factor_cov_pval)
# 
# color_by <- "sbp"
# corr_factor_cov <- round(corr_factors_cov[factor,color_by],3)
# plotFactorCov(MOFA_model = MOFA.model, factor = factor, color_by = color_by, shape_by = shape_by, corr_factor_cov = corr_factor_cov, corr_factor_cov_pval = corr_factor_cov_pval)
# 
# color_by <- "dbp"
# corr_factor_cov <- round(corr_factors_cov[factor,color_by],3)
# plotFactorCov(MOFA_model = MOFA.model, factor = factor, color_by = color_by, shape_by = shape_by, corr_factor_cov = corr_factor_cov, corr_factor_cov_pval = corr_factor_cov_pval)
# 
# color_by <- "Bcell"
# corr_factor_cov <- round(corr_factors_cov[factor,color_by],3)
# plotFactorCov(MOFA_model = MOFA.model, factor = factor, color_by = color_by, shape_by = shape_by, corr_factor_cov = corr_factor_cov, corr_factor_cov_pval = corr_factor_cov_pval)
# 
# cat("\nFACTOR 21\n", file = stdout())
# 
# factor <- 21
# color_by <- "cvd"
# shape_by <- "sex"
# corr_factor_cov <- round(corr_factors_cov[factor,color_by],3)
# corr_factor_cov_pval <- signif(cor$p[factor,color_by], digits = 2)
# sort(corr_factors_cov[factor,])
# plotFactorCov(MOFA_model = MOFA.model, factor = factor, color_by = color_by, shape_by = shape_by, corr_factor_cov = corr_factor_cov, corr_factor_cov_pval = corr_factor_cov_pval)
# 
# shape_by <- NULL
# 
# color_by <- "sex"
# corr_factor_cov <- round(corr_factors_cov[factor,color_by],3)
# corr_factor_cov_pval <- signif(cor$p[factor,color_by], digits = 2)
# plotFactorCov(MOFA_model = MOFA.model, factor = factor, color_by = color_by, shape_by = shape_by, corr_factor_cov = corr_factor_cov, corr_factor_cov_pval = corr_factor_cov_pval)
# 
# shape_by <- "cvd"
# color_by <- "tot_chol"
# corr_factor_cov <- round(corr_factors_cov[factor,color_by],3)
# plotFactorCov(MOFA_model = MOFA.model, factor = factor, color_by = color_by, shape_by = shape_by, corr_factor_cov = corr_factor_cov, corr_factor_cov_pval = corr_factor_cov_pval)
# 
# color_by <- "hdl_chol"
# corr_factor_cov <- round(corr_factors_cov[factor,color_by],3)
# plotFactorCov(MOFA_model = MOFA.model, factor = factor, color_by = color_by, shape_by = shape_by, corr_factor_cov = corr_factor_cov, corr_factor_cov_pval = corr_factor_cov_pval)
# 
# color_by <- "weight"
# corr_factor_cov <- round(corr_factors_cov[factor,color_by],3)
# plotFactorCov(MOFA_model = MOFA.model, factor = factor, color_by = color_by, shape_by = shape_by, corr_factor_cov = corr_factor_cov, corr_factor_cov_pval = corr_factor_cov_pval)
# 
# color_by <- "height"
# corr_factor_cov <- round(corr_factors_cov[factor,color_by],3)
# plotFactorCov(MOFA_model = MOFA.model, factor = factor, color_by = color_by, shape_by = shape_by, corr_factor_cov = corr_factor_cov, corr_factor_cov_pval = corr_factor_cov_pval)
# 
# color_by <- "CD8T"
# corr_factor_cov <- round(corr_factors_cov[factor,color_by],3)
# plotFactorCov(MOFA_model = MOFA.model, factor = factor, color_by = color_by, shape_by = shape_by, corr_factor_cov = corr_factor_cov, corr_factor_cov_pval = corr_factor_cov_pval)
# 
# cat("\nFACTOR 27\n", file = stdout())
# 
# cat("\nChecking outliers from factor 27...\n", file = stdout())
# 
# factor.values <- all_factors$group1[,27]
# summary(factor.values)
# factor.values[which(factor.values > 6 | factor.values < -6)]
# outliers <- which(factor.values > 6 | factor.values < -6 )
# cat("\n...Outliers...\n", file = stdout())
# outliers
# 
# samples.filter <- rownames(all_factors$group1)[-outliers]
# length(samples.filter)
# MOFA.model.filt <- subset_samples(MOFA.model, samples.filter)
# MOFA.model.filt
# 
# shape_by <- NULL
# 
# factor <- 27
# color_by <- "cvd"
# corr_factor_cov <- round(corr_factors_cov[factor,color_by],3)
# corr_factor_cov_pval <- signif(cor$p[factor,color_by], digits = 2)
# sort(corr_factors_cov[factor,])
# plotFactorCov(MOFA_model = MOFA.model.filt, factor = factor, color_by = color_by, shape_by = shape_by, corr_factor_cov = corr_factor_cov, corr_factor_cov_pval = corr_factor_cov_pval)
# 
# color_by <- "sex"
# corr_factor_cov <- round(corr_factors_cov[factor,color_by],3)
# corr_factor_cov_pval <- signif(cor$p[factor,color_by], digits = 2)
# plotFactorCov(MOFA_model = MOFA.model.filt, factor = factor, color_by = color_by, shape_by = shape_by, corr_factor_cov = corr_factor_cov, corr_factor_cov_pval = corr_factor_cov_pval)
#
# shape_by <- "cvd"
# 
# color_by <- "age"
# corr_factor_cov <- round(corr_factors_cov[factor,color_by],3)
# plotFactorCov(MOFA_model = MOFA.model.filt, factor = factor, color_by = color_by, shape_by = shape_by, corr_factor_cov = corr_factor_cov)
# 
# color_by <- "Gran"
# corr_factor_cov <- round(corr_factors_cov[factor,color_by],3)
# plotFactorCov(MOFA_model = MOFA.model.filt, factor = factor, color_by = color_by, shape_by = shape_by, corr_factor_cov = corr_factor_cov)
# 
# color_by <- "Bcell"
# corr_factor_cov <- round(corr_factors_cov[factor,color_by],3)
# plotFactorCov(MOFA_model = MOFA.model.filt, factor = factor, color_by = color_by, shape_by = shape_by, corr_factor_cov = corr_factor_cov)
# 
# color_by <- "Mono"
# corr_factor_cov <- round(corr_factors_cov[factor,color_by],3)
# plotFactorCov(MOFA_model = MOFA.model.filt, factor = factor, color_by = color_by, shape_by = shape_by, corr_factor_cov = corr_factor_cov)
# 
# # cat("\n\n... 4.4) VISUALIZATION OF COMBINATION OF FACTORS ...\n\n", file = stdout())
# #
# # png(paste("MOFA/3.downstream_analysis/",groups,"_",samples,"/variance_decomposition/factors1-30_pairs.png", sep = ""), height=3600, width=6000, res=600, units="px")
# # plot_factors(MOFA.model,
# #              factors = 1:30,
# #              color_by = "cvd"
# # )
# # dev.off()

cat("\n\n... 4.5) VISUALIZATION OF FEATURES WEIGTHS FOR ALL FACTORS ...\n\n", file = stdout())

# Top 10 features plots

# pdf(paste("MOFA/3.downstream_analysis/",groups,"_",samples,"/variance_decomposition/mRNA_10features_factors1-30.pdf", sep = ""))
# for (factor in 1:30) {
#   #png(paste(paste("MOFA/3.downstream_analysis/",groups,"_",samples,"/variance_decomposition/mRNA_features_factor",factor,".png", sep = ""),sep = ""), height=3600, width=6000, res=600, units="px")
#   p <- plot_weights(MOFA.model,
#                view = "mRNA",
#                factor = factor,
#                nfeatures = 10,     # Top number of features to highlight
#                scale = T           # Scale weights from -1 to 1
#   )
#   print(p)
#   #dev.off()
# }
# dev.off()
# 
### Only for factors 9,19,21,27
# factors <- c(9,19,21,27)
# 
# for (factor in factors) {
#   png(paste(paste("MOFA/3.downstream_analysis/",groups,"_",samples,"/variance_decomposition/meth_8features_factor",factor,".png", sep = ""),sep = ""), height=3600, width=6000, res=600, units="px")
#   p <- plot_weights(MOFA.model,
#                       view = "meth",
#                       factor = factor,
#                       nfeatures = 8,     # Top number of features to highlight
#                       scale = T           # Scale weights from -1 to 1
#   )
#   print(p)
# }
# dev.off()
# 
# factor <- 9
# png(paste(paste("MOFA/3.downstream_analysis/",groups,"_",samples,"/variance_decomposition/meth_8features_factor",factor,".png", sep = ""),sep = ""), height=3600, width=6000, res=600, units="px")
# p <- plot_weights(MOFA.model,
#                   view = "meth",
#                   factor = factor,
#                   nfeatures = 8,     # Top number of features to highlight
#                   scale = T           # Scale weights from -1 to 1
# )
# print(p)
# 
# factor <- 19
# png(paste(paste("MOFA/3.downstream_analysis/",groups,"_",samples,"/variance_decomposition/meth_8features_factor",factor,".png", sep = ""),sep = ""), height=3600, width=6000, res=600, units="px")
# p <- plot_weights(MOFA.model,
#                   view = "meth",
#                   factor = factor,
#                   nfeatures = 8,     # Top number of features to highlight
#                   scale = T           # Scale weights from -1 to 1
# )
# print(p)
# 
# factor <- 21
# png(paste(paste("MOFA/3.downstream_analysis/",groups,"_",samples,"/variance_decomposition/meth_8features_factor",factor,".png", sep = ""),sep = ""), height=3600, width=6000, res=600, units="px")
# p <- plot_weights(MOFA.model,
#                   view = "meth",
#                   factor = factor,
#                   nfeatures = 8,     # Top number of features to highlight
#                   scale = T           # Scale weights from -1 to 1
# )
# print(p)
# 
# factor <- 27
# png(paste(paste("MOFA/3.downstream_analysis/",groups,"_",samples,"/variance_decomposition/meth_8features_factor",factor,".png", sep = ""),sep = ""), height=3600, width=6000, res=600, units="px")
# p <- plot_weights(MOFA.model,
#                   view = "meth",
#                   factor = factor,
#                   nfeatures = 8,     # Top number of features to highlight
#                   scale = T           # Scale weights from -1 to 1
# )
# print(p)
# 
## For all factors
# 
# pdf(paste("MOFA/3.downstream_analysis/",groups,"_",samples,"/variance_decomposition/meth_10features_factors1-30.pdf", sep = ""))
# for (factor in 1:30) {
#   #png(paste(paste("MOFA/3.downstream_analysis/",groups,"_",samples,"/variance_decomposition/meth_features_factor",factor,".png", sep = ""),sep = ""), height=3600, width=6000, res=600, units="px")
#   p <- plot_weights(MOFA.model,
#                view = "meth",
#                factor = factor,
#                nfeatures = 10,     # Top number of features to highlight
#                scale = T           # Scale weights from -1 to 1
#   )
#   print(p)
#   #dev.off()
# }
# dev.off()

# # Top 30 features weigths plots
# 
# ## Only for factors 9,19,21,27
# for (factor in factors) {
#   png(paste(paste("MOFA/3.downstream_analysis/",groups,"_",samples,"/variance_decomposition/meth_30features_weights_factor",factor,".png", sep = ""),sep = ""), height=3600, width=6000, res=600, units="px")
#   p <- plot_top_weights(MOFA.model,
#                         view = "meth",
#                         factor = factor,
#                         abs = T,
#                         nfeatures = 30,     # Top number of features to highlight
#                         scale = T           # Scale weights from -1 to 1
#   )
#   print(p)
#   #dev.off()
# }
# dev.off()

## For all factors

# pdf(paste("MOFA/3.downstream_analysis/",groups,"_",samples,"/variance_decomposition/meth_30features_weigths_factors1-30.pdf", sep = ""))
# top30.meth.features.weigths <- lapply(1:30, function(factor) {
#   #png(paste("MOFA/3.downstream_analysis/",groups,"_",samples,"/variance_decomposition/meth_features_factor14_weigths.png", sep = ""), height=3600, width=6000, res=600, units="px")
#   p <- plot_top_weights(MOFA.model,
#                    view = "meth",
#                    factor = factor,
#                    abs = T,
#                    nfeatures = 30,     # Top number of features to highlight
#                    scale = T           # Scale weights from -1 to 1
#   )
#   print(p)
#   as.character(p$data$feature)
#   #dev.off()
#
# })
# dev.off()
# top30.meth.features.weigths

# pdf(paste("MOFA/3.downstream_analysis/",groups,"_",samples,"/variance_decomposition/mRNA_30features_weigths_factors1-30.pdf", sep = ""))
# top30.mRNA.features.weigths <- lapply(1:30, function(factor) {
#   #png(paste("MOFA/3.downstream_analysis/",groups,"_",samples,"/variance_decomposition/mRNA_features_factor14_weigths.png", sep = ""), height=3600, width=6000, res=600, units="px")
#   p <- plot_top_weights(MOFA.model,
#                    view = "mRNA",
#                    factor = factor,
#                    abs = T,
#                    nfeatures = 30,     # Top number of features to highlight
#                    scale = T           # Scale weights from -1 to 1
#   )
#   print(p)
#   as.character(p$data$feature)
#   #dev.off()
# })
# dev.off()
# top30.mRNA.features.weigths

cat("\n\n... 5) VISUALISATION OF PATTERNS IN THE DATA ...\n\n", file = stdout())

cat("\n\n... 5.1) HEATMAP FOR EACH FACTOR...\n\n", file = stdout())

# # Change covariate values by quantiles
# library(gtools)
# MOFA.covariates <- read.csv(file = paste("MOFA/MOFA_covariates_order_",groups,"_",samples,".csv", sep = ""),
#                             header=TRUE, stringsAsFactors = F, sep=",")
# # Metadata has to contain the columns 'sample' (shareid) and 'group' (cvd) for MOFA to work
# colnames(MOFA.covariates)[1] <- "sample"
# # Change sample names from integers to characters
# MOFA.covariates$sample <- as.character(MOFA.covariates$sample)
# MOFA.covariates[1:3,]
# dim(MOFA.covariates)
# # 2055 24
# for (cov in names(MOFA.covariates)[-c(1:4,15,22)]) {
#   MOFA.covariates[,cov] <- quantcut(MOFA.covariates[,cov])
# }
# write.csv(MOFA.covariates, file = paste("MOFA/MOFA_covariates_order_",groups,"_",samples,"_quantiles.csv", sep = ""))

groups <- "nogroups"
samples <- "914_2055"

MOFA.covariates <- read.csv(file = paste("MOFA/MOFA_covariates_order_",groups,"_",samples,"_quantiles.csv", sep = ""),
                            header=TRUE, stringsAsFactors = F, sep=",")
MOFA.covariates$sample <- as.character(MOFA.covariates$sample)

samples_metadata(MOFA.model) <- as.data.frame(apply(MOFA.covariates[,-c(1,3)],2,factor))

## Heatmap for factors 9,19,21,27
factor <- 9
cov.interest <- c("sex","age","CD4T","CD8T","NK")
for (cov in cov.interest) {
  png(paste("MOFA/3.downstream_analysis/",groups,"_",samples,"/variance_decomposition/meth_features_factor",factor,"_",cov,"_heatmap.png", sep = ""), height=3600, width=6000, res=600, units="px")
  p <- plot_data_heatmap(MOFA.model,
                         view = "meth",         # view of interest
                         factor = factor,             # factor of interest
                         features = 30,          # number of features to plot (they are selected by loading)
                         # extra arguments that are passed to the `pheatmap` function
                         cluster_rows = TRUE, cluster_cols = TRUE,
                         show_rownames = TRUE, show_colnames = FALSE,
                         main = paste("Heatmap for methylation, factor: ",factor, sep = ""),
                         annotation_samples = c(cov,"cvd") #annotation_colors = list("CVD"=colors)
  )
  meth.features.heatmap.f9 <- p$tree_row$labels
  print(p)
  dev.off()  
}

factor <- 19
cov.interest <- c("sex","age","Bcell")
for (cov in cov.interest) {
  png(paste("MOFA/3.downstream_analysis/",groups,"_",samples,"/variance_decomposition/meth_features_factor",factor,"_",cov,"_heatmap.png", sep = ""), height=3600, width=6000, res=600, units="px")
  p <- plot_data_heatmap(MOFA.model,
                         view = "meth",         # view of interest
                         factor = factor,             # factor of interest
                         features = 30,          # number of features to plot (they are selected by loading)
                         # extra arguments that are passed to the `pheatmap` function
                         cluster_rows = TRUE, cluster_cols = TRUE,
                         show_rownames = TRUE, show_colnames = FALSE,
                         main = paste("Heatmap for methylation, factor: ",factor, sep = ""),
                         annotation_samples = c(cov,"cvd") #annotation_colors = list("CVD"=colors)
  )
  meth.features.heatmap.f19 <- p$tree_row$labels
  dev.off()
}

factor <- 21
cov.interest <- c("sex","CD8T","height","tot_chol")
for (cov in cov.interest) {
  png(paste("MOFA/3.downstream_analysis/",groups,"_",samples,"/variance_decomposition/meth_features_factor",factor,"_",cov,"_heatmap.png", sep = ""), height=3600, width=6000, res=600, units="px")
  p <- plot_data_heatmap(MOFA.model,
                         view = "meth",         # view of interest
                         factor = factor,             # factor of interest
                         features = 30,          # number of features to plot (they are selected by loading)
                         # extra arguments that are passed to the `pheatmap` function
                         cluster_rows = TRUE, cluster_cols = TRUE,
                         show_rownames = TRUE, show_colnames = FALSE,
                         main = paste("Heatmap for methylation, factor: ",factor, sep = ""),
                         annotation_samples = c(cov,"cvd") #annotation_colors = list("CVD"=colors)
  )
  meth.features.heatmap.f21 <- p$tree_row$labels
  dev.off()
}

factor <- 27
cov.interest <- c("Gran","Bcell","Mono")
for (cov in cov.interest) {
  png(paste("MOFA/3.downstream_analysis/",groups,"_",samples,"/variance_decomposition/meth_features_factor",factor,"_",cov,"_heatmap.png", sep = ""), height=3600, width=6000, res=600, units="px")
  p <- plot_data_heatmap(MOFA.model,
                         view = "meth",         # view of interest
                         factor = factor,             # factor of interest
                         features = 30,          # number of features to plot (they are selected by loading)
                         # extra arguments that are passed to the `pheatmap` function
                         cluster_rows = TRUE, cluster_cols = TRUE,
                         show_rownames = TRUE, show_colnames = FALSE,
                         main = paste("Heatmap for methylation, factor: ",factor, sep = ""),
                         annotation_samples = c(cov,"cvd") #annotation_colors = list("CVD"=colors)
  )
  meth.features.heatmap.f27 <- p$tree_row$labels
  dev.off()
}

stop("dew")

meth.features.heatmap <- list(meth.features.heatmap.f9,meth.features.heatmap.f19,meth.features.heatmap.f21,meth.features.heatmap.f27)

save(meth.features.heatmap, file = paste("MOFA/3.downstream_analysis/",groups,"_",samples,"/variance_decomposition/meth_30features_heatmap.RData", sep = ""))

## Heatmap for all factors

# pdf(paste("MOFA/3.downstream_analysis/",groups,"_",samples,"/variance_decomposition/meth_features_factors1-30_heatmap.pdf", sep = ""))
# meth.features.heatmap <- lapply(c(1:14, 17:30), function(factor) {
# 
#   p <- plot_data_heatmap(MOFA.model,
#                          view = "meth",         # view of interest
#                          factor = factor,             # factor of interest
#                          features = 30,          # number of features to plot (they are selected by loading)
#                          # extra arguments that are passed to the `pheatmap` function
#                          cluster_rows = TRUE, cluster_cols = TRUE,
#                          show_rownames = TRUE, show_colnames = FALSE,
#                          main = paste("Heatmap for methylation, factor: ",factor, sep = ""),
#                          annotation_samples = c("sex","cvd") #annotation_colors = list("CVD"=colors)
#   )
#   p$tree_row$labels
# 
# })
# dev.off()
# meth.features.heatmap
# save(meth.features.heatmap, file = paste("MOFA/3.downstream_analysis/",groups,"_",samples,"/variance_decomposition/meth_30features_heatmap.RData", sep = ""))
# 
# pdf(paste("MOFA/3.downstream_analysis/",groups,"_",samples,"/variance_decomposition/mRNA_features_factors1-30_heatmap.pdf", sep = ""))
# mRNA.features.heatmap <- lapply(1:30, function(factor) {
# 
#   p <- plot_data_heatmap(MOFA.model,
#                          view = "mRNA",         # view of interest
#                          factor = factor,             # factor of interest
#                          features = 30,          # number of features to plot (they are selected by loading)
#                          # extra arguments that are passed to the `pheatmap` function
#                          cluster_rows = TRUE, cluster_cols = TRUE,
#                          show_rownames = TRUE, show_colnames = FALSE,
#                          main = paste("Heatmap for mRNA, factor: ",factor, sep = ""),
#                          annotation_samples = c("sex","cvd") #annotation_colors = list("CVD"=colors)
#   )
#   p$tree_row$labels
# 
# })
# dev.off()
# mRNA.features.heatmap
# save(mRNA.features.heatmap, file = paste("MOFA/3.downstream_analysis/",groups,"_",samples,"/variance_decomposition/mRNA_30features_heatmap.RData", sep = ""))

cat("\n\n... 5.2) SCATTERPLOTS ...\n\n", file = stdout())

# png(paste("MOFA/3.downstream_analysis/",groups,"_",samples,"/variance_decomposition/meth_features_factor14_scatterplot.png", sep = ""), height=3600, width=6000, res=600, units="px")
# plot_data_scatter(MOFA.model,
#                   view = "meth",         # view of interest
#                   factor = 14,             # factor of interest
#                   #groups = "all",
#                   features = 5,           # number of features to plot (they are selected by loading)
#                   sign = "all",
#                   add_lm = TRUE,          # add linear regression
#                   color_by = "cvd",
#                   shape_by = "sex",
#                   legend = TRUE,
#                   alpha = 1,
#                   dot_size = 1,
#                   text_size = 5,
#                   imputed = FALSE
# )
# dev.off()
#
# png(paste("MOFA/3.downstream_analysis/",groups,"_",samples,"/variance_decomposition/mRNA_features_factor14_scatterplot.png", sep = ""), height=3600, width=6000, res=600, units="px")
# plot_data_scatter(MOFA.model,
#                   view = "mRNA",         # view of interest
#                   factor = 14,             # factor of interest
#                   #groups = "all",
#                   features = 5,           # number of features to plot (they are selected by loading)
#                   sign = "all",
#                   add_lm = TRUE,          # add linear regression
#                   color_by = "cvd",
#                   shape_by = "sex",
#                   legend = TRUE,
#                   alpha = 1,
#                   dot_size = 1,
#                   text_size = 5,
#                   imputed = FALSE
# )
# dev.off()

cat("\n\n... 5.3) T-SNE ...\n\n", file = stdout())

model <- run_tsne(MOFA.model)

for (i in 1:length(colnames(MOFA.model@samples_metadata))) {

  if ( colnames(MOFA.model@samples_metadata)[i] == "sample" | colnames(MOFA.model@samples_metadata)[i] == "shareid") {
    next
  }
  print(i)
  print(colnames(MOFA.model@samples_metadata)[i])
  png(paste("MOFA/3.downstream_analysis/",groups,"_",samples,"/variance_decomposition/TSNE/TSNE_",colnames(MOFA.model@samples_metadata)[i],".png", sep = ""), height=3600, width=6000, res=600, units="px")
  p <- plot_dimred(model,
                   method = "TSNE",
                   color_by = colnames(MOFA.model@samples_metadata)[i],
                   legend = TRUE,
                   show_missing = TRUE
  )
  print(p)
  dev.off()

}

png(paste("MOFA/3.downstream_analysis/",groups,"_",samples,"/variance_decomposition/TSNE/TSNE_factor2.png", sep = ""), height=3600, width=6000, res=600, units="px")
p <- plot_dimred(model,
                 method = "TSNE",
                 color_by = "Factor2",
                 legend = TRUE
)
print(p)
dev.off()

png(paste("MOFA/3.downstream_analysis/",groups,"_",samples,"/variance_decomposition/TSNE/TSNE_factor2.png", sep = ""), height=3600, width=6000, res=600, units="px")
p <- plot_dimred(model,
                 method = "TSNE",
                 color_by = "Factor11",
                 legend = TRUE
)
print(p)
dev.off()


date()

cat("\n\nComputational time\n\n", file = stdout())
end_time <- Sys.time()
end_time - start_time

cat("\n\n ######################## 3. MOFA DOWNSTREAM ANALYSIS: VARIANCE DECOMPOSITION SCRIPT ENDS ######################## \n\n", file = stdout())





