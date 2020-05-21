rm(list=ls())
start_time <- Sys.time()
date()

cat("\n\n ######################## START LOADING LIBRARIES ######################## \n\n", file = stdout())

library(oligo)
cat("\n\n oligo LOADED \n\n", file = stdout())
library(genefilter)
cat("\n\n genefilter LOADED \n\n", file = stdout())

cat("\n\n ######################## FINISHED LOADING LIBRARIES ######################## \n\n", file = stdout())

cat("\n\n ######################## 5. FILTERING SCRIPT BEGINS ######################## \n\n", file = stdout())

cat("\n\n... SETTING WORKING DIRECTORY ...\n\n", file = stdout())
setwd("/projects/regicor/Guille/TFM/transcriptomics")
getwd()

cat("\n\n... LOADING NORMALIZED ExpressionSet OBJECT ...\n\n", file = stdout())
ESet.name <- load("QC/3.normalization/ESet_norm.RData")
ESet <- get(ESet.name)
cat("\n\n // ExpressionSet //\n", file = stdout())
getClass(ESet)

cat("\n\n 1) FILTERING OUTLIER SAMPLES \n\n", file = stdout())

cat("\n\n 1.1) LOADING OUTLIER SAMPLES \n\n", file = stdout())

cat("\n\n... LOADING BOXPLOT OUTLIER SAMPLES ...\n\n", file = stdout())

boxplot.outlier.samples.name <- load("QC/2.plots/boxplots/boxplot_outlier_samples.RData")
boxplot.outlier.samples <- get(boxplot.outlier.samples.name)
cat("\n\n // BOXPLOT outliers samples //\n", file = stdout())
names(boxplot.outlier.samples) <- paste(names(boxplot.outlier.samples),"-",sep="")
unlist(boxplot.outlier.samples)

cat("\n\n... LOADING DENSITY PLOTS OUTLIER SAMPLES ...\n\n", file = stdout())

# density.outlier.samples.name <- load("QC/4.norm_plots/density_plots/density_outlier_samples.RData")
# density.outlier.samples <- get(density.outlier.samples.name)
# cat("\n\n // DENSITY PLOTS outliers samples //\n", file = stdout())
# density.outlier.samples

cat("\n\n... LOADING NUSE PLOTS OUTLIER SAMPLES ...\n\n", file = stdout())

NUSE.outlier.samples.name <- load("QC/2.plots/NUSE_plots/NUSE_outlier_samples.RData")
NUSE.outlier.samples <- get(NUSE.outlier.samples.name)
cat("\n\n // NUSE outliers samples //\n", file = stdout())
names(NUSE.outlier.samples) <- paste(names(NUSE.outlier.samples),"-",sep="")
unlist(NUSE.outlier.samples)

cat("\n\n... LOADING RLE PLOTS SAMPLES ...\n\n", file = stdout())

RLE.outlier.samples.name <- load("QC/2.plots/RLE_plots/RLE_outlier_samples.RData")
RLE.outlier.samples <- get(RLE.outlier.samples.name)
cat("\n\n // RLE outliers samples //\n", file = stdout())
names(RLE.outlier.samples) <- paste(names(RLE.outlier.samples),"-",sep="")
unlist(RLE.outlier.samples)

cat("\n\n... MA plots manually ...\n\n", file = stdout())

MA.outlier.samples <- c("KAW101122_FHS_OFF_13C10","KAW101122_FHS_OFF_13D11","KAW101122_FHS_OFF_13E09",
  "KAW101122_FHS_OFF_13F07","KAW101122_FHS_OFF_13G12", "KAW101129_FHS_OFF_14C07",
  "KAW101129_FHS_OFF_14E07", "KAW101129_FHS_OFF_14H08", "KAW101215_FHS_OFF_16B09",
  "KAW101215_FHS_OFF_16C10", "KAW101215_FHS_OFF_16D11", "KAW101215_FHS_OFF_16E08",
  "KAW101215_FHS_OFF_16E09", "KAW101215_FHS_OFF_16E10", "KAW101215_FHS_OFF_16E11",
  "KAW101227_FHS_OFF_18A09", "KAW101227_FHS_OFF_18B08", "KAW101227_FHS_OFF_18E07",
  "KAW110104_FHS_OFF_19A04", "KAW110104_FHS_OFF_19G02", "KAW110125_FHS_OFF_22A07",
  "KAW110125_FHS_OFF_22A10", "KAW110125_FHS_OFF_22B11", "KAW110125_FHS_OFF_22E08",
  "KAW110125_FHS_OFF_22E09", "KAW110207_FHS_OFF_24A07", "KAW110207_FHS_OFF_24D10",
  "KAW110209_FHS_OFF_25C11", "KAW110209_FHS_OFF_25D09", "KAW110209_FHS_OFF_25E08",
  "KAW110209_FHS_OFF_25E12", "KAW110209_FHS_OFF_25H08", "KAW110222_FHS_OFF_26A11",
  "KAW110222_FHS_OFF_26H09", "KAW110316_FHS_OFF_28A11", "KAW110321_FHS_OFF_27H12",
  "NRC101220_FHS_OFF_17C07", "NRC101220_FHS_OFF_17H08", "PLX110324_FHS_OFF_30D09",
  "PLX110324_FHS_OFF_30D12", "RWX100909_FHS_OFF_07B05", "RWX100909_FHS_OFF_07C01",
  "RWX100909_FHS_OFF_07E03", "RWX100912_FHS_OFF_07A09", "RWX100912_FHS_OFF_07F12",
  "RWX100912_FHS_OFF_07H12", "RWX100917_FHS_OFF_08B08", "RWX100917_FHS_OFF_08C09",
  "RWX100917_FHS_OFF_08G11", "RWX101101_FHS_OFF_10H06", "RWX101108_FHS_OFF_11A04",
  "RWX101108_FHS_OFF_11H01", "RWX101109_FHS_OFF_11D09", "RWX101220_FHS_OFF_17A06",
  "RWX101220_FHS_OFF_17B04", "RWX101220_FHS_OFF_17D02", "RWX110104_FHS_OFF_19B12")

MA.outlier.samples.medians <- c(-0.714,-0.549,-0.427,
           -0.567,-0.559,-0.0651,
            0.695,0.591,0.576,
            0.51,0.255,0.69,
            0.756, 0.562, 0.468,
           -0.0931, 0.138, 0.636,
           0.498, 0.34, 0.55,
           0.0935, -0.143, 0.41,
           0.535, -0.023, -0.115,
           0.363, 0.864, 0.921,
           0.625, 0.993, -0.537,
           -0.205, 0.126, 0.116,
           1.27, -0.0157, -0.446,
           -0.61, -1.06, -1.12, 
           -0.865, -1.39, -1.56,
           -1.66, -1.55, -1.65,
           -1.67, -0.893, -1.4,
           -1.22, -0.593, 0.787, 
           -1.15, -1.14, 0.711)

# outlier <- "RWX110104_FHS_OFF_19B12"
# outlier%in%unlist(RLE.outlier.samples)
# outlier%in%unlist(NUSE.outlier.samples)
# outlier%in%unlist(boxplot.outlier.samples)



cat("\n\n 1.2) CONCATENATING IN THE SAME VECTOR \n\n", file = stdout())
all.outlier.samples <- as.vector(c(unlist(RLE.outlier.samples), unlist(NUSE.outlier.samples), unlist(boxplot.outlier.samples), MA.outlier.samples))
cat("\nNumber of unique outliers (bad quality samples in at least 1 plot) \n", file = stdout())
length(unique(all.outlier.samples))
cat("\nNumber of outliers that are at least bad quality samples in 2 plots\n", file = stdout())
outlier.samples.duplicated <- unique(all.outlier.samples[which(duplicated(all.outlier.samples))])
length(outlier.samples.duplicated)

outlier.samples.4plots <- c()
outlier.samples.3plots <- c()
outlier.samples.2plots <- c()
#a <- c()

for (i in 1:length(outlier.samples.duplicated)) {
  
    # Outlier in 4 plots
    if (outlier.samples.duplicated[i]%in%unlist(RLE.outlier.samples) && outlier.samples.duplicated[i]%in%unlist(NUSE.outlier.samples) 
        && outlier.samples.duplicated[i]%in%unlist(boxplot.outlier.samples) && outlier.samples.duplicated[i]%in%MA.outlier.samples) {
      outlier.samples.4plots <- c(outlier.samples.4plots, outlier.samples.duplicated[i])
    }
    # Outlier in 3 plots
    else if (outlier.samples.duplicated[i]%in%unlist(RLE.outlier.samples) && outlier.samples.duplicated[i]%in%unlist(NUSE.outlier.samples) 
        && outlier.samples.duplicated[i]%in%unlist(boxplot.outlier.samples) ||
        outlier.samples.duplicated[i]%in%unlist(RLE.outlier.samples) && outlier.samples.duplicated[i]%in%unlist(NUSE.outlier.samples) 
        && outlier.samples.duplicated[i]%in%MA.outlier.samples ||
        outlier.samples.duplicated[i]%in%unlist(boxplot.outlier.samples) && outlier.samples.duplicated[i]%in%unlist(NUSE.outlier.samples) 
        && outlier.samples.duplicated[i]%in%MA.outlier.samples ||
        outlier.samples.duplicated[i]%in%unlist(RLE.outlier.samples) && outlier.samples.duplicated[i]%in%MA.outlier.samples
        && outlier.samples.duplicated[i]%in%unlist(boxplot.outlier.samples)
        ) {
      outlier.samples.3plots <- c(outlier.samples.3plots, outlier.samples.duplicated[i])
    }
  
    # Outlier in 2 plots that are: RLE + NUSE or NUSE + boxplot or RLE + boxplot (i.e. no MA + boxplot/RLE/NUSE)
    else if ((outlier.samples.duplicated[i]%in%unlist(RLE.outlier.samples) && outlier.samples.duplicated[i]%in%unlist(NUSE.outlier.samples)) 
           || (outlier.samples.duplicated[i]%in%unlist(NUSE.outlier.samples) && outlier.samples.duplicated[i]%in%unlist(boxplot.outlier.samples))
           || (outlier.samples.duplicated[i]%in%unlist(RLE.outlier.samples) && outlier.samples.duplicated[i]%in%unlist(boxplot.outlier.samples))) {
    outlier.samples.2plots <- c(outlier.samples.2plots, outlier.samples.duplicated[i])
    }
  # else if (outlier.samples.duplicated[i]%in%unlist(boxplot.outlier.samples) && outlier.samples.duplicated[i]%in%MA.outlier.samples) {
  #   a <- c(a, outlier.samples.duplicated[i])
  #   }
}
cat("\nNumber of outliers that are bad quality samples in 2 plots (MA plot not taken into account)\n", file = stdout())
length(outlier.samples.2plots)
cat("\nNumber of outliers that are bad quality samples in 3 plots\n", file = stdout())
length(outlier.samples.3plots)
cat("\nNumber of outliers that are bad quality samples in 4 plots\n", file = stdout())
length(outlier.samples.4plots)

cat("\nTotal number of outliers that are bad quality samples\n", file = stdout())
all.outlier.samples <- c(outlier.samples.2plots, outlier.samples.3plots, outlier.samples.4plots)
length(all.outlier.samples)

# sort(substr(outlier.samples.4plots, nchar(outlier.samples.4plots)-4,nchar(outlier.samples.4plots)))
# outlier <- "27H12"
# outlier%in%substr(MA.outlier.samples, nchar(MA.outlier.samples)-4,nchar(MA.outlier.samples))
# outlier%in%substr(unlist(RLE.outlier.samples), nchar(unlist(RLE.outlier.samples))-4,nchar(unlist(RLE.outlier.samples)))
# outlier%in%substr(unlist(NUSE.outlier.samples), nchar(unlist(NUSE.outlier.samples))-4,nchar(unlist(NUSE.outlier.samples)))
# outlier%in%substr(unlist(boxplot.outlier.samples), nchar(unlist(boxplot.outlier.samples))-4,nchar(unlist(boxplot.outlier.samples)))

# cat("\n % of RLE outliers that are also outliers in NUSE\n", file = stdout())
# table(unlist(RLE.outlier.samples)%in%unlist(NUSE.outlier.samples))
# mean(unlist(RLE.outlier.samples)%in%unlist(NUSE.outlier.samples))
# cat("\n RLE outliers that are NOT NUSE outliers\n", file = stdout())
# unlist(RLE.outlier.samples)[which(!unlist(RLE.outlier.samples)%in%unlist(NUSE.outlier.samples))]

cat("\n\n 1.3) REMOVING OUTLIER SAMPLES FROM THE EXPRESSION SET \n\n", file = stdout())

ESet.filt.samples <- ESet[,!rownames(pData(ESet))%in%all.outlier.samples]
cat("\n\n // ExpressionSet filtered //\n", file = stdout())
getClass(ESet.filt.samples)
assayDataElement(ESet.filt.samples,'exprs')[1:5,1:5]

cat("\n\n 2) FILTERING OUTLIER LOWLY EXPRESSED GENES \n\n", file = stdout())

## CRITERIA ##
# Transcripts that do not have expression values lower than a given threshold in at least 
# as many samples as the smallest experimental group (140 CVD) are excluded.

cat("\n\n 2.1) GENE EXPRESSION THRESHOLD \n\n", file = stdout())

# Can be the minimum value of all genes, or one chosen by us (let's say 4).
#exp.cutoff <- min(exprs(ESet.filt.samples))
exp.cutoff <- 4
exp.cutoff

cat("\n\n 2.2) SAMPLES NUMBER THRESHOLD \n\n", file = stdout())

samples.metadata <- read.csv(file = "phenotype/Analysis/samples_metadata_merged.csv",
                                 header=TRUE, stringsAsFactors = F, sep=",")
samples.cutoff <- min(table(samples.metadata$cvd))
samples.cutoff

cat("\n\n 2.3) USING GENEFILTER \n\n", file = stdout())
f1 <- kOverA(samples.cutoff, exp.cutoff)
ffun <- filterfun(f1)
wh1 <- genefilter(exprs(ESet.filt.samples), ffun)
cat("\nNumber of transcripts that pass the filtering process\n", file = stdout())
sum(wh1)

cat("\n\n 2.4) USING OUR FUNCTIONS \n\n", file = stdout())

ESet.filt.samples.medians <- rowMedians(Biobase::exprs(ESet.filt.samples))

pdf("QC/5.filtering/lowly_exp_genes.pdf")
hist_res <- hist(ESet.filt.samples.medians, 100, col = "cornsilk", freq = FALSE, 
                 main = "Histogram of the median intensities",
                 border = "antiquewhite4",
                 xlab = "Median intensities")
abline(v = exp.cutoff, col = "coral4", lwd = 2)
dev.off()

genes.filter <- apply(Biobase::exprs(ESet.filt.samples), 1, 
                           function(x) { sum(x > exp.cutoff) >= samples.cutoff } 
                           )
cat("\nNumber of transcripts that pass the filtering process\n", file = stdout())
table(genes.filter)

cat("\n\n 2.5) REMOVING OUTLIER GENES FROM THE EXPRESSION SET \n\n", file = stdout())
ESet.filt.samples.genes <- ESet.filt.samples[genes.filter,]
#ESet.filt.samples.genes <- subset(ESet.filt.samples, genes.filter)

cat("\n\n // ExpressionSet filtered //\n", file = stdout())
getClass(ESet.filt.samples.genes)
assayDataElement(ESet.filt.samples.genes,'exprs')[1:5,1:5]

cat("\n\n ... SAVING FILTERED EXPRESSION SET ...\n\n", file = stdout())

save(ESet.filt.samples.genes, file = "QC/5.filtering/ESet_norm_filt.RData")

date()

cat("\nComputational time\n", file = stdout())
end_time <- Sys.time()
print(end_time - start_time)

cat("\n\n ######################## 5. FILTERING SCRIPT ENDS ######################## \n\n", file = stdout())

