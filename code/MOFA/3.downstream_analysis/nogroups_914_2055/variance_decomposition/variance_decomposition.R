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
# Factor 11
factors <- c(11,23)

MOFA.model <- load_model(paste("MOFA/2.training_model/MOFA_trained_",groups,"_",samples,".hdf5",sep = ""))
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
# 2080 24

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

png(paste("MOFA/3.downstream_analysis/",groups,"_",samples,"/variance_decomposition/cov_factors_correlation_R.png", sep = ""), height=3600, width=6000, res=600, units="px")
corr_factors_cov <- correlate_factors_with_covariates(MOFA.model,
                                  covariates = names(MOFA.model@samples_metadata)[-24],
                                  factors = "all",
                                  groups = "all",
                                  abs = FALSE,
                                  plot = c("r"),
                                  alpha = 0.05,
                                  return_data = FALSE,
                                  transpose = FALSE
)
corr_factors_cov[1:5,1:5]
dev.off()

cat(paste("\nFactor ",factors[1], "\n", sep = ""), file = stdout())
sort(corr_factors_cov[factors[1],])
cat(paste("\nFactor ",factors[2], "\n", sep = ""), file = stdout())
sort(corr_factors_cov[factors[2],])

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
plot.variance.1 <- plot_variance_explained(MOFA.model, x="group", y="factor")
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

factors <- get_factors(MOFA.model,
                       groups = "all",
                       factors = "all"
)

factor.values <- factors$group1[,1]
summary(factor.values)
factor.values[which(factor.values == max(factor.values))]
outlier1 <- which(factor.values == max(factor.values))
outlier1


cat("\nChecking an outlier from factor 19...\n", file = stdout())

factor.values <- factors$group1[,19]
summary(factor.values)
factor.values[which(factor.values == max(factor.values))]
outlier2 <- which(factor.values == max(factor.values))
outliers <- as.vector(c(outlier1,outlier2))
cat("\n...Outliers...\n", file = stdout())
outliers

#samples.filter <- samples.filter[which(factor.values != max(factor.values))]

samples.filter <- rownames(factors$group1)[-outliers]
length(samples.filter)
MOFA.model.filt <- subset_samples(MOFA.model, samples.filter)
MOFA.model.filt

png(paste("MOFA/3.downstream_analysis/",groups,"_",samples,"/variance_decomposition/factors1-15_sex_cvd.png", sep = ""), height=3600, width=6000, res=600, units="px")
plot_factor(MOFA.model.filt,
            factor = 1:15,
            color_by = "cvd",
            shape_by = "sex"
)
dev.off()

png(paste("MOFA/3.downstream_analysis/",groups,"_",samples,"/variance_decomposition/factors16-30_sex_cvd.png", sep = ""), height=3600, width=6000, res=600, units="px")
plot_factor(MOFA.model.filt,
            factor = 16:30,
            color_by = "cvd",
            shape_by = "sex"
)
dev.off()
# For the violin plot you cannot use quantitative variables as it makes the plot for each group

png(paste("MOFA/3.downstream_analysis/",groups,"_",samples,"/variance_decomposition/factors1-15_cvd_violin.png", sep = ""), height=3600, width=6000, res=600, units="px")
plot_factor(MOFA.model.filt,
            factor = 1:15,
            color_by = "cvd",
            dot_size = 1.3,      # change dot size
            dodge = TRUE,           # dodge points with different colors
            legend = TRUE,          # remove legend
            add_violin = TRUE,      # add violin plots,
            color_violin = TRUE,
            scale = TRUE
            #violin_alpha = 0.01  # transparency of violin plots
) # + scale_color_manual(values=colors) + scale_fill_manual(values=colors)
dev.off()

png(paste("MOFA/3.downstream_analysis/",groups,"_",samples,"/variance_decomposition/factors16-30_cvd_violin.png", sep = ""), height=3600, width=6000, res=600, units="px")
plot_factor(MOFA.model.filt,
            factor = 16:30,
            color_by = "cvd",
            dot_size = 1.3,      # change dot size
            dodge = TRUE,           # dodge points with different colors
            legend = TRUE,          # remove legend
            add_violin = TRUE,      # add violin plots,
            color_violin = TRUE,
            scale = TRUE
            #violin_alpha = 0.01  # transparency of violin plots
) # + scale_color_manual(values=colors) + scale_fill_manual(values=colors)
dev.off()

# png("MOFA/3.downstream_analysis/factor1_factor2.png", sep = ""), height=3600, width=6000, res=600, units="px")
# plot_factors(MOFA.model,
#              factors = c(1,2),
#              color_by = "cvd"
# )
# dev.off()

cat("\n\n... 4.3.2) SPECIFIC FACTORS COLORED BY COVARIABLES ...\n\n", file = stdout())

##Violin for factor 11

png(paste("MOFA/3.downstream_analysis/",groups,"_",samples,"/variance_decomposition/factor11_cvd_violin.png", sep = ""), height=3600, width=6000, res=600, units="px")
plot_factor(MOFA.model.filt,
            factor = 11,
            color_by = "cvd",
            dot_size = 1.3,      # change dot size
            dodge = TRUE,           # dodge points with different colors
            legend = TRUE,          # remove legend
            add_violin = TRUE,      # add violin plots,
            color_violin = TRUE,
            scale = TRUE
            #violin_alpha = 0.01  # transparency of violin plots
) # + scale_color_manual(values=colors) + scale_fill_manual(values=colors)
dev.off()

## Violin for factor 23

png(paste("MOFA/3.downstream_analysis/",groups,"_",samples,"/variance_decomposition/factor23_cvd_violin.png", sep = ""), height=3600, width=6000, res=600, units="px")
plot_factor(MOFA.model.filt,
            factor = 23,
            color_by = "cvd",
            dot_size = 1.3,      # change dot size
            dodge = TRUE,           # dodge points with different colors
            legend = TRUE,          # remove legend
            add_violin = TRUE,      # add violin plots,
            color_violin = TRUE,
            scale = TRUE
            #violin_alpha = 0.01  # transparency of violin plots
) # + scale_color_manual(values=colors) + scale_fill_manual(values=colors)
dev.off()

plotFactorCov <- function(factor,color_by,shape_by,corr_factor_cov){

  p <- plot_factor(MOFA.model,
                   factor = factor,
                   color_by = color_by,
                   shape_by = shape_by
  )
  png(paste("MOFA/3.downstream_analysis/",groups,"_",samples,"/variance_decomposition/factor_cov/factor",factor,"_",color_by,".png", sep = ""), height=3600, width=6000, res=600, units="px")
  print(p + ggtitle(paste("Factor",factor," - ",color_by," clustering. R²: ",corr_factor_cov,sep="")) +
    theme(plot.title = element_text(hjust = 0.5)) +
    labs(fill = color_by))
  dev.off()
}

cat("\nFACTOR 2\n", file = stdout())
factor <- 2
color_by <- "cvd"
shape_by <- "sex"
corr_factor_cov <- round(corr_factors_cov[factor,color_by],3)
sort(corr_factors_cov[factor,])
plotFactorCov(factor = factor, color_by = color_by, shape_by = shape_by, corr_factor_cov = corr_factor_cov)

color_by <- "CD8T"
corr_factor_cov <- round(corr_factors_cov[factor,color_by],3)
plotFactorCov(factor = factor, color_by = color_by, shape_by = shape_by, corr_factor_cov = corr_factor_cov)

cat("\nFACTOR 4\n", file = stdout())
factor <- 4
color_by <- "CD4T"
shape_by <- "cvd"
corr_factor_cov <- round(corr_factors_cov[factor,color_by],3)
sort(corr_factors_cov[factor,])
plotFactorCov(factor = factor, color_by = color_by, shape_by = shape_by, corr_factor_cov = corr_factor_cov)

color_by <- "Bcell"
corr_factor_cov <- round(corr_factors_cov[factor,color_by],3)
plotFactorCov(factor = factor, color_by = color_by, shape_by = shape_by, corr_factor_cov = corr_factor_cov)

color_by <- "CD8T"
corr_factor_cov <- round(corr_factors_cov[factor,color_by],3)
plotFactorCov(factor = factor, color_by = color_by, shape_by = shape_by, corr_factor_cov = corr_factor_cov)

color_by <- "Gran"
corr_factor_cov <- round(corr_factors_cov[factor,color_by],3)
plotFactorCov(factor = factor, color_by = color_by, shape_by = shape_by, corr_factor_cov = corr_factor_cov)

color_by <- "Mono"
corr_factor_cov <- round(corr_factors_cov[factor,color_by],3)
plotFactorCov(factor = factor, color_by = color_by, shape_by = shape_by, corr_factor_cov = corr_factor_cov)

color_by <- "NK"
corr_factor_cov <- round(corr_factors_cov[factor,color_by],3)
plotFactorCov(factor = factor, color_by = color_by, shape_by = shape_by, corr_factor_cov = corr_factor_cov)

color_by <- "age"
corr_factor_cov <- round(corr_factors_cov[factor,color_by],3)
plotFactorCov(factor = factor, color_by = color_by, shape_by = shape_by, corr_factor_cov = corr_factor_cov)

cat("\nFACTOR 5\n", file = stdout())
factor <- 5
color_by <- "CD8T"
shape_by <- "cvd"
corr_factor_cov <- round(corr_factors_cov[factor,color_by],3)
sort(corr_factors_cov[factor,])
plotFactorCov(factor = factor, color_by = color_by, shape_by = shape_by, corr_factor_cov = corr_factor_cov)

cat("\nFACTOR 6\n", file = stdout())
factor <- 6
color_by <- "sex"
shape_by <- "cvd"
corr_factor_cov <- round(corr_factors_cov[factor,color_by],3)
sort(corr_factors_cov[factor,])
plotFactorCov(factor = factor, color_by = color_by, shape_by = shape_by, corr_factor_cov = corr_factor_cov)

color_by <- "weight"
corr_factor_cov <- round(corr_factors_cov[factor,color_by],3)
plotFactorCov(factor = factor, color_by = color_by, shape_by = shape_by, corr_factor_cov = corr_factor_cov)

color_by <- "height"
corr_factor_cov <- round(corr_factors_cov[factor,color_by],3)
plotFactorCov(factor = factor, color_by = color_by, shape_by = shape_by, corr_factor_cov = corr_factor_cov)

color_by <- "Bcell"
corr_factor_cov <- round(corr_factors_cov[factor,color_by],3)
plotFactorCov(factor = factor, color_by = color_by, shape_by = shape_by, corr_factor_cov = corr_factor_cov)

cat("\nFACTOR 8\n", file = stdout())
factor <- 8
color_by <- "sex"
shape_by <- "cvd"
corr_factor_cov <- round(corr_factors_cov[factor,color_by],3)
sort(corr_factors_cov[factor,])
plotFactorCov(factor = factor, color_by = color_by, shape_by = shape_by, corr_factor_cov = corr_factor_cov)

color_by <- "height"
corr_factor_cov <- round(corr_factors_cov[factor,color_by],3)
plotFactorCov(factor = factor, color_by = color_by, shape_by = shape_by, corr_factor_cov = corr_factor_cov)

color_by <- "weight"
corr_factor_cov <- round(corr_factors_cov[factor,color_by],3)
plotFactorCov(factor = factor, color_by = color_by, shape_by = shape_by, corr_factor_cov = corr_factor_cov)

color_by <- "hdl_chol"
corr_factor_cov <- round(corr_factors_cov[factor,color_by],3)
plotFactorCov(factor = factor, color_by = color_by, shape_by = shape_by, corr_factor_cov = corr_factor_cov)

cat("\nFACTOR 9\n", file = stdout())
factor <- 9
color_by <- "age"
shape_by <- "cvd"
corr_factor_cov <- round(corr_factors_cov[factor,color_by],3)
sort(corr_factors_cov[factor,])
plotFactorCov(factor = factor, color_by = color_by, shape_by = shape_by, corr_factor_cov = corr_factor_cov)

color_by <- "CD8T"
corr_factor_cov <- round(corr_factors_cov[factor,color_by],3)
plotFactorCov(factor = factor, color_by = color_by, shape_by = shape_by, corr_factor_cov = corr_factor_cov)

color_by <- "NK"
corr_factor_cov <- round(corr_factors_cov[factor,color_by],3)
plotFactorCov(factor = factor, color_by = color_by, shape_by = shape_by, corr_factor_cov = corr_factor_cov)

color_by <- "Bcell"
corr_factor_cov <- round(corr_factors_cov[factor,color_by],3)
plotFactorCov(factor = factor, color_by = color_by, shape_by = shape_by, corr_factor_cov = corr_factor_cov)

color_by <- "Gran"
corr_factor_cov <- round(corr_factors_cov[factor,color_by],3)
plotFactorCov(factor = factor, color_by = color_by, shape_by = shape_by, corr_factor_cov = corr_factor_cov)

cat("\nFACTOR 10\n", file = stdout())
factor <- 10
color_by <- "sex"
shape_by <- "cvd"
corr_factor_cov <- round(corr_factors_cov[factor,color_by],3)
sort(corr_factors_cov[factor,])
plotFactorCov(factor = factor, color_by = color_by, shape_by = shape_by, corr_factor_cov = corr_factor_cov)

color_by <- "weight"
corr_factor_cov <- round(corr_factors_cov[factor,color_by],3)
plotFactorCov(factor = factor, color_by = color_by, shape_by = shape_by, corr_factor_cov = corr_factor_cov)

color_by <- "height"
corr_factor_cov <- round(corr_factors_cov[factor,color_by],3)
plotFactorCov(factor = factor, color_by = color_by, shape_by = shape_by, corr_factor_cov = corr_factor_cov)

color_by <- "Bcell"
corr_factor_cov <- round(corr_factors_cov[factor,color_by],3)
plotFactorCov(factor = factor, color_by = color_by, shape_by = shape_by, corr_factor_cov = corr_factor_cov)

cat("\nFACTOR 11\n", file = stdout())
factor <- 11
color_by <- "cvd"
shape_by <- "sex"
corr_factor_cov <- round(corr_factors_cov[factor,color_by],3)
sort(corr_factors_cov[factor,])
plotFactorCov(factor = factor, color_by = color_by, shape_by = shape_by, corr_factor_cov = corr_factor_cov)

shape_by <- "cvd"

color_by <- "tot_chol"
corr_factor_cov <- round(corr_factors_cov[factor,color_by],3)
plotFactorCov(factor = factor, color_by = color_by, shape_by = shape_by, corr_factor_cov = corr_factor_cov)

color_by <- "hdl_chol"
corr_factor_cov <- round(corr_factors_cov[factor,color_by],3)
plotFactorCov(factor = factor, color_by = color_by, shape_by = shape_by, corr_factor_cov = corr_factor_cov)

color_by <- "weight"
corr_factor_cov <- round(corr_factors_cov[factor,color_by],3)
plotFactorCov(factor = factor, color_by = color_by, shape_by = shape_by, corr_factor_cov = corr_factor_cov)

color_by <- "height"
corr_factor_cov <- round(corr_factors_cov[factor,color_by],3)
plotFactorCov(factor = factor, color_by = color_by, shape_by = shape_by, corr_factor_cov = corr_factor_cov)

color_by <- "CD8T"
corr_factor_cov <- round(corr_factors_cov[factor,color_by],3)
plotFactorCov(factor = factor, color_by = color_by, shape_by = shape_by, corr_factor_cov = corr_factor_cov)

color_by <- "CD4T"
corr_factor_cov <- round(corr_factors_cov[factor,color_by],3)
plotFactorCov(factor = factor, color_by = color_by, shape_by = shape_by, corr_factor_cov = corr_factor_cov)

color_by <- "age"
corr_factor_cov <- round(corr_factors_cov[factor,color_by],3)
plotFactorCov(factor = factor, color_by = color_by, shape_by = shape_by, corr_factor_cov = corr_factor_cov)

color_by <- "Mono"
corr_factor_cov <- round(corr_factors_cov[factor,color_by],3)
plotFactorCov(factor = factor, color_by = color_by, shape_by = shape_by, corr_factor_cov = corr_factor_cov)

color_by <- "Gran"
corr_factor_cov <- round(corr_factors_cov[factor,color_by],3)
plotFactorCov(factor = factor, color_by = color_by, shape_by = shape_by, corr_factor_cov = corr_factor_cov)

color_by <- "glucose"
corr_factor_cov <- round(corr_factors_cov[factor,color_by],3)
plotFactorCov(factor = factor, color_by = color_by, shape_by = shape_by, corr_factor_cov = corr_factor_cov)

color_by <- "sex"
corr_factor_cov <- round(corr_factors_cov[factor,color_by],3)
plotFactorCov(factor = factor, color_by = color_by, shape_by = shape_by, corr_factor_cov = corr_factor_cov)

cat("\nFACTOR 12\n", file = stdout())
factor <- 12
color_by <- "Array_METH"
shape_by <- "cvd"
corr_factor_cov <- round(corr_factors_cov[factor,color_by],3)
sort(corr_factors_cov[factor,])
plotFactorCov(factor = factor, color_by = color_by, shape_by = shape_by, corr_factor_cov = corr_factor_cov)

cat("\nFACTOR 14\n", file = stdout())
factor <- 14
color_by <- "sex"
shape_by <- "cvd"
corr_factor_cov <- round(corr_factors_cov[factor,color_by],3)
sort(corr_factors_cov[factor,])
plotFactorCov(factor = factor, color_by = color_by, shape_by = shape_by, corr_factor_cov = corr_factor_cov)

color_by <- "height"
corr_factor_cov <- round(corr_factors_cov[factor,color_by],3)
plotFactorCov(factor = factor, color_by = color_by, shape_by = shape_by, corr_factor_cov = corr_factor_cov)

color_by <- "CD8T"
corr_factor_cov <- round(corr_factors_cov[factor,color_by],3)
plotFactorCov(factor = factor, color_by = color_by, shape_by = shape_by, corr_factor_cov = corr_factor_cov)

cat("\nFACTOR 23\n", file = stdout())
factor <- 23
color_by <- "cvd"
shape_by <- "sex"
corr_factor_cov <- round(corr_factors_cov[factor,color_by],3)
sort(corr_factors_cov[factor,])
plotFactorCov(factor = factor, color_by = color_by, shape_by = shape_by, corr_factor_cov = corr_factor_cov)

shape_by <- "cvd"

color_by <- "sex"
corr_factor_cov <- round(corr_factors_cov[factor,color_by],3)
plotFactorCov(factor = factor, color_by = color_by, shape_by = shape_by, corr_factor_cov = corr_factor_cov)

color_by <- "age"
corr_factor_cov <- round(corr_factors_cov[factor,color_by],3)
plotFactorCov(factor = factor, color_by = color_by, shape_by = shape_by, corr_factor_cov = corr_factor_cov)

color_by <- "tot_chol"
corr_factor_cov <- round(corr_factors_cov[factor,color_by],3)
plotFactorCov(factor = factor, color_by = color_by, shape_by = shape_by, corr_factor_cov = corr_factor_cov)

color_by <- "hdl_chol"
corr_factor_cov <- round(corr_factors_cov[factor,color_by],3)
plotFactorCov(factor = factor, color_by = color_by, shape_by = shape_by, corr_factor_cov = corr_factor_cov)

color_by <- "glucose"
corr_factor_cov <- round(corr_factors_cov[factor,color_by],3)
plotFactorCov(factor = factor, color_by = color_by, shape_by = shape_by, corr_factor_cov = corr_factor_cov)

color_by <- "Gran"
corr_factor_cov <- round(corr_factors_cov[factor,color_by],3)
plotFactorCov(factor = factor, color_by = color_by, shape_by = shape_by, corr_factor_cov = corr_factor_cov)

color_by <- "sbp"
corr_factor_cov <- round(corr_factors_cov[factor,color_by],3)
plotFactorCov(factor = factor, color_by = color_by, shape_by = shape_by, corr_factor_cov = corr_factor_cov)

color_by <- "dbp"
corr_factor_cov <- round(corr_factors_cov[factor,color_by],3)
plotFactorCov(factor = factor, color_by = color_by, shape_by = shape_by, corr_factor_cov = corr_factor_cov)

color_by <- "CD4T"
corr_factor_cov <- round(corr_factors_cov[factor,color_by],3)
plotFactorCov(factor = factor, color_by = color_by, shape_by = shape_by, corr_factor_cov = corr_factor_cov)

color_by <- "CD8T"
corr_factor_cov <- round(corr_factors_cov[factor,color_by],3)
plotFactorCov(factor = factor, color_by = color_by, shape_by = shape_by, corr_factor_cov = corr_factor_cov)

color_by <- "Bcell"
corr_factor_cov <- round(corr_factors_cov[factor,color_by],3)
plotFactorCov(factor = factor, color_by = color_by, shape_by = shape_by, corr_factor_cov = corr_factor_cov)

color_by <- "Mono"
corr_factor_cov <- round(corr_factors_cov[factor,color_by],3)
plotFactorCov(factor = factor, color_by = color_by, shape_by = shape_by, corr_factor_cov = corr_factor_cov)

# cat("\n\n... 4.4) VISUALIZATION OF COMBINATION OF FACTORS ...\n\n", file = stdout())
#
# png(paste("MOFA/3.downstream_analysis/",groups,"_",samples,"/variance_decomposition/factors1-30_pairs.png", sep = ""), height=3600, width=6000, res=600, units="px")
# plot_factors(MOFA.model,
#              factors = 1:30,
#              color_by = "cvd"
# )
# dev.off()

cat("\n\n... 4.5) VISUALIZATION OF FEATURES WEIGTHS FOR ALL FACTORS ...\n\n", file = stdout())

# Top 10 features plots

pdf(paste("MOFA/3.downstream_analysis/",groups,"_",samples,"/variance_decomposition/mRNA_10features_factors1-30.pdf", sep = ""))
for (factor in 1:30) {
  #png(paste(paste("MOFA/3.downstream_analysis/",groups,"_",samples,"/variance_decomposition/mRNA_features_factor",factor,".png", sep = ""),sep = ""), height=3600, width=6000, res=600, units="px")
  p <- plot_weights(MOFA.model,
               view = "mRNA",
               factor = factor,
               nfeatures = 10,     # Top number of features to highlight
               scale = T           # Scale weights from -1 to 1
  )
  print(p)
  #dev.off()
}
dev.off()

## Only for factor 11
factor <- 11
png(paste(paste("MOFA/3.downstream_analysis/",groups,"_",samples,"/variance_decomposition/meth_8features_factor",factor,".png", sep = ""),sep = ""), height=3600, width=6000, res=600, units="px")
p <- plot_weights(MOFA.model,
                    view = "meth",
                    factor = factor,
                    nfeatures = 8,     # Top number of features to highlight
                    scale = T           # Scale weights from -1 to 1
)
print(p)
dev.off()


## For all factors

pdf(paste("MOFA/3.downstream_analysis/",groups,"_",samples,"/variance_decomposition/meth_10features_factors1-30.pdf", sep = ""))
for (factor in 1:30) {
  #png(paste(paste("MOFA/3.downstream_analysis/",groups,"_",samples,"/variance_decomposition/meth_features_factor",factor,".png", sep = ""),sep = ""), height=3600, width=6000, res=600, units="px")
  p <- plot_weights(MOFA.model,
               view = "meth",
               factor = factor,
               nfeatures = 10,     # Top number of features to highlight
               scale = T           # Scale weights from -1 to 1
  )
  print(p)
  #dev.off()
}
dev.off()

# Top 30 features weigths plots


## For factor 11 only
factor <- 11
png(paste(paste("MOFA/3.downstream_analysis/",groups,"_",samples,"/variance_decomposition/meth_30features_weights_factor",factor,".png", sep = ""),sep = ""), height=3600, width=6000, res=600, units="px")
p <- plot_top_weights(MOFA.model,
                        view = "meth",
                        factor = factor,
                        abs = T,
                        nfeatures = 30,     # Top number of features to highlight
                        scale = T           # Scale weights from -1 to 1
)
print(p)
as.character(p$data$feature)
dev.off()


## For all factors

pdf(paste("MOFA/3.downstream_analysis/",groups,"_",samples,"/variance_decomposition/meth_30features_weigths_factors1-30.pdf", sep = ""))
top30.meth.features.weigths <- lapply(1:30, function(factor) {
  #png(paste("MOFA/3.downstream_analysis/",groups,"_",samples,"/variance_decomposition/meth_features_factor14_weigths.png", sep = ""), height=3600, width=6000, res=600, units="px")
  p <- plot_top_weights(MOFA.model,
                   view = "meth",
                   factor = factor,
                   abs = T,
                   nfeatures = 30,     # Top number of features to highlight
                   scale = T           # Scale weights from -1 to 1
  )
  print(p)
  as.character(p$data$feature)
  #dev.off()

})
dev.off()
top30.meth.features.weigths

pdf(paste("MOFA/3.downstream_analysis/",groups,"_",samples,"/variance_decomposition/mRNA_30features_weigths_factors1-30.pdf", sep = ""))
top30.mRNA.features.weigths <- lapply(1:30, function(factor) {
  #png(paste("MOFA/3.downstream_analysis/",groups,"_",samples,"/variance_decomposition/mRNA_features_factor14_weigths.png", sep = ""), height=3600, width=6000, res=600, units="px")
  p <- plot_top_weights(MOFA.model,
                   view = "mRNA",
                   factor = factor,
                   abs = T,
                   nfeatures = 30,     # Top number of features to highlight
                   scale = T           # Scale weights from -1 to 1
  )
  print(p)
  as.character(p$data$feature)
  #dev.off()
})
dev.off()
top30.mRNA.features.weigths

cat("\n\n... 5) VISUALISATION OF PATTERNS IN THE DATA ...\n\n", file = stdout())

cat("\n\n... 5.1) HEATMAP FOR EACH FACTOR...\n\n", file = stdout())

## Heatmap for factor 11

factor <- 11

png(paste("MOFA/3.downstream_analysis/",groups,"_",samples,"/variance_decomposition/meth_30_features_heatmap_factor11.png", sep = ""), height=3600, width=6000, res=600, units="px")
p <- plot_data_heatmap(MOFA.model,
                       view = "meth",         # view of interest
                       factor = 11,             # factor of interest
                       features = 30,          # number of features to plot (they are selected by loading)
                       # extra arguments that are passed to the `pheatmap` function
                       cluster_rows = TRUE, cluster_cols = TRUE,
                       show_rownames = TRUE, show_colnames = FALSE,
                       main = paste("Heatmap for methylation, factor: ",factor, sep = ""),
                       annotation_samples = c("sex","cvd") #annotation_colors = list("CVD"=colors)
)
dev.off()


## Heatmap for all factors

pdf(paste("MOFA/3.downstream_analysis/",groups,"_",samples,"/variance_decomposition/meth_features_factors1-30_heatmap.pdf", sep = ""))
meth.features.heatmap <- lapply(1:30, function(factor) {

  p <- plot_data_heatmap(MOFA.model,
                         view = "meth",         # view of interest
                         factor = factor,             # factor of interest
                         features = 30,          # number of features to plot (they are selected by loading)
                         # extra arguments that are passed to the `pheatmap` function
                         cluster_rows = TRUE, cluster_cols = TRUE,
                         show_rownames = TRUE, show_colnames = FALSE,
                         main = paste("Heatmap for methylation, factor: ",factor, sep = ""),
                         annotation_samples = c("sex","cvd") #annotation_colors = list("CVD"=colors)
  )
  p$tree_row$labels

})
dev.off()
meth.features.heatmap
save(meth.features.heatmap, file = paste("MOFA/3.downstream_analysis/",groups,"_",samples,"/variance_decomposition/meth_30features_heatmap.RData", sep = ""))

pdf(paste("MOFA/3.downstream_analysis/",groups,"_",samples,"/variance_decomposition/mRNA_features_factors1-30_heatmap.pdf", sep = ""))
mRNA.features.heatmap <- lapply(1:30, function(factor) {

  p <- plot_data_heatmap(MOFA.model,
                         view = "mRNA",         # view of interest
                         factor = factor,             # factor of interest
                         features = 30,          # number of features to plot (they are selected by loading)
                         # extra arguments that are passed to the `pheatmap` function
                         cluster_rows = TRUE, cluster_cols = TRUE,
                         show_rownames = TRUE, show_colnames = FALSE,
                         main = paste("Heatmap for mRNA, factor: ",factor, sep = ""),
                         annotation_samples = c("sex","cvd") #annotation_colors = list("CVD"=colors)
  )
  p$tree_row$labels

})
dev.off()
mRNA.features.heatmap
save(mRNA.features.heatmap, file = paste("MOFA/3.downstream_analysis/",groups,"_",samples,"/variance_decomposition/mRNA_30features_heatmap.RData", sep = ""))

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





