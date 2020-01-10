### Installation

# en terminal
#sudo apt-get install libcurl4-openssl-dev
#sudo apt-get install python-virtualenv
# install.packages("RCurl")
# BiocManager::install("MultiAssayExperiment")
# install.packages("reticulate")
# py_install("mofapy", envname = "r-reticulate", method="auto")
# BiocManager::install("MOFA")
# BiocManager::install("MOFAdata")

### Loading libraries

library(reticulate)
library(ggplot2)
library(MultiAssayExperiment)
library(MOFA)
library(MOFAdata)

### 1) Load data and create MOFA object

# Load data
# import list with mRNA, Methylation, Drug Response and Mutation data. 
data("CLL_data") 

# check dimensionalities, samples are columns, features are rows
lapply(CLL_data, dim)
lapply(CLL_data, head)

# Load sample metadata: Sex and Diagnosis
data("CLL_covariates")
head(CLL_covariates)

# Create MultiAssayExperiment object 
mae_CLL <- MultiAssayExperiment(
  experiments = CLL_data, 
  colData = CLL_covariates
)

# Build the MOFA object
MOFAobject <- createMOFAobject(mae_CLL)
MOFAobject

# Overview of training Data
plotDataOverview(MOFAobject)

### 2) Fit the MOFA model

### 2.1) Data options

# scaleViews: logical indicating whether to scale views to have unit variance. 
# As long as the scale of the different data sets is not too high, this is not required. Default is FALSE.

# removeIncompleteSamples: logical indicating whether to remove samples that are not profiled in all omics. 
# The model can cope with missing assays, so this option is not required. Default is FALSE.

DataOptions <- getDefaultDataOptions()
DataOptions

### 2.2) Model options

# numFactors: number of factors (default is 0.5 times the number of samples). By default, the model will only 
# remove a factor if it explains exactly zero variance in the data. You can increase this threshold on minimum 
# variance explained by setting TrainOptions$dropFactorThreshold to a value higher than zero.

# likelihoods: likelihood for each view. Usually we recommend gaussian for continuous data, bernoulli for binary 
# data and poisson for count data. By default, the model tries to guess it from the data.

# sparsity: do you want to use sparsity? This makes the interpretation easier so it is recommended (Default is TRUE).

ModelOptions <- getDefaultModelOptions(MOFAobject)
ModelOptions$numFactors <- 25
ModelOptions

### 2.3) Training options

# maxiter: maximum number of iterations. Ideally set it large enough and use the convergence criterion TrainOptions$tolerance.

# tolerance: convergence threshold based on change in the evidence lower bound. For an exploratory run you can use a 
# value between 1.0 and 0.1, but for a “final” model we recommend a value of 0.01.

# DropFactorThreshold: hyperparameter to automatically learn the number of factors based on a minimum variance explained 
# criteria. Factors explaining less than DropFactorThreshold fraction of variation in all views will be removed. 
# For example, a value of 0.01 means that factors that explain less than 1% of variance in all views will be discarded. 
# By default this it zero, meaning that all factors are kept unless they explain no variance at all.

TrainOptions <- getDefaultTrainOptions()

# Automatically drop factors that explain less than 2% of variance in all omics
TrainOptions$DropFactorThreshold <- 0.02
TrainOptions$seed <- 2017
TrainOptions$verbose <- TRUE
TrainOptions

### 2.4) Prepare MOFA

MOFAobject <- prepareMOFA(
  MOFAobject, 
  DataOptions = DataOptions,
  ModelOptions = ModelOptions,
  TrainOptions = TrainOptions
)

# Optionally, we can choose to regress out some (technical) covariates before training, using a simple linear model. 
# For example, here we can choose to remove the effect of sex. Ideally, all undesired sources of variation should be 
# removed a priori from the model. The reason ebing that, if strong technical factors exist, the model will “focus” 
# on capturing the variability driven by the technical factors, and small sources of biological variability could be missed.
# (Note: uncomment and running the function below will lead to a slight modification of the results)

# MOFAobject <- regressCovariates(
#   object = MOFAobject,
#   views = c("Drugs","Methylation","mRNA"),
#   covariates = MOFAobject@InputData$Gender
# )

### 2.5) Run MOFA

# This step can take some time (around 15 min with default parameters)

MOFAobject <- runMOFA(MOFAobject)
MOFAobject

### 3) Analyse a trained MOFA model

### 3.1) Part 1: Disentangling the heterogeneity: calculation of variance explained by each factor in each view

# Calculation of variance explained by each factor in each view. This is probably the most important plot that MOFA 
# generates, as it summarises the entire heterogeneity of the dataset in a single figure. Here we can see in which 
# view a factor explains variation which can guide further characterisation of the factors by investigating the weights 
# in those views.

# This is done by calculateVarianceExplained (to get the numerical values) and plotVarianceExplained (to get the plot). 
# The resulting figure gives an overview of which factors are active in which view(s). If a factor is active in more than 
# one view, this means that is capturing shared signal (co-variation) between features of different data modalities.
# Here, for example Factor 1 is active in all data modalities, while Factor 4 is specific to mRNA.

# Calculate the variance explained (R2) per factor in each view 
r2 <- calculateVarianceExplained(MOFAobject)
r2$R2Total
## POR QUE NO SUMA 1 TODO?
# Variance explained by each factor in each view
head(r2$R2PerFactor)

# Plot it
plotVarianceExplained(MOFAobject)

### 3.2) Part 2: Characterisation of individual factors

# Inspection of top features with highest loadings: the loading is a measure of feature importance, so features with 
# high loading are the ones driving the heterogeneity captured by the factor.

# Feature set enrichment analysis (where set annotations are present, e.g. gene sets for mRNA views).

# Ordination of samples by factors to reveal clusters and/or gradients: this is similar to what is traditionally done 
# with Principal Component Analysis or t-SNE.

# Other analyses, including imputation of missing values and clustering of samples are also available. See below for 
# a short illustration of these functionalities. In addition, the factors can be used in further analyses, for example 
# as predictors, e.g. for predicting clinical outcome or classifying patients, or to control for unwanted sources of 
# variation. Vignettes illustrating this are coming soon.

### 3.2.1) Inspection of top weighted features in the active views

# To get an overview of the weights across all factors in a given view you can use the plotWeightsHeatmap function. 
# For example, here we plot all weights from all factors in the Mutations data:

plotWeightsHeatmap(
  MOFAobject, 
  view = "Mutations", 
  factors = 1:5,
  show_colnames = FALSE
)

# We observe that Factors 1 and 2 have large non-zero weights. To explore a given factor in more detail we can plot 
# all weights for a single factor using the plotWeights function. For example, here we plot all weights from Factor 1 
# in the Mutation data. With nfeatures we can set how many features should be labelled

plotWeights(
  MOFAobject, 
  view = "Mutations", 
  factor = 1, 
  nfeatures = 5
)

# `manual` let’s you specify features manually to be labelled in the plot

plotWeights(
  MOFAobject, 
  view = "Mutations", 
  factor = 1, 
  nfeatures = 5,
  manual = list(c("BRAF"),c("MED12")),
  color_manual = c("red","blue")
  
)

plotTopWeights(
  MOFAobject, 
  view="Mutations", 
  factor=1
)

# Again, features with large weight in a given factor means that they follow the pattern of covariation 
# associated with the factor. For example, here the factor is associated with the B-cell of tumour’s origin, 
# consistent with a large weight on the IGHV status (see our manuscript for more details).

# From the previous plots, we can clearly see that Factor 1 is associated to IGHV status. As the variance decomposition 
# above told us that this factor is also relevant on all the other data modalities we can investigate its weights on 
# other modalities, e.g. mRNA, to make connections of the IGHV-linked axes of variation to other molecular layers.

plotTopWeights(
  MOFAobject, 
  view = "mRNA", 
  factor = 1
)

# Finally, instead of looking at an “abstract” weight, it is useful to observe the coordinated heterogeneity of the top 
# features in the original data. This can be done using the plotDataHeatmap function. In this plot samples (in rows) are 
# ordered according to their value on the factor (here Factor 1). Here, this shows clear patterns of the samples’ gene 
# expression in the 20 top weighted genes along the factor.

plotDataHeatmap(
  MOFAobject, 
  view = "mRNA", 
  factor = 1, 
  features = 20, 
  show_rownames = FALSE
)


### 3.2.2) Feature set enrichment analysis in the active views

# Sometimes looking at the loadings of single features can be challenging, and often the combination of signal from 
# functionally related sets of features (i.e. gene ontologies) is required.

# Here we implemented a function for feature set enrichment analysis method (runEnrichmentAnalysis) derived from the PCGSE package.

# The input of this function is a MOFA trained model (MOFAmodel), the factors for which to perform feature set enrichment 
# (a character vector), the feature sets (a binary matrix) and a set of options regarding how the analysis should be performed, 
# see also documentation of runEnrichmentAnalysis

# We illustrate the use of this function using the reactome annotations, which are contained in the package. Depending on your 
# data other gene or feature sets might be useful and you can prepare your customized feature sets and specify it using the 
# feature.sets argument of the function.

# Load reactome annotations
data("reactomeGS") # binary matrix with feature sets in rows and features in columns
reactomeGS[1:5,1:5]
str(reactomeGS)

# perform enrichment analysis
gsea <- runEnrichmentAnalysis(
  MOFAobject,
  view = "mRNA",
  feature.sets = reactomeGS,
  alpha = 0.01
)

# The next step is to visualise the results of the Gene Set Enrichment Analysis. There are several ways:
  
# Plot the number of enriched gene sets per factor

plotEnrichmentBars(gsea, alpha=0.01)

# From this we find enriched at a FDR of 1% gene sets on Factors 3-6 and 8. To look into which gene sets these 
# are we can choose a factor of interest and visualize the most enriched gene sets as follows:
  
# Plot the top enriched pathways for every factor

interestingFactors <- 4:5

fseaplots <- lapply(interestingFactors, function(factor) {
  plotEnrichment(
    MOFAobject,
    gsea,
    factor = factor,
    alpha = 0.01,
    max.pathways = 10 # The top number of pathways to display
  )
})

cowplot::plot_grid(fseaplots[[1]], fseaplots[[2]],
                   ncol = 1, labels = paste("Factor", interestingFactors))

# This shows us that Factor 4 is capturing variation related to immune response (possibly due to T-cell contamination of the samples) 
# and Factor 5 is related to differences in stress response, as discussed in our paper.

# 3.2.3) Ordination of samples by factors to reveal clusters and gradients in the sample space

# Samples can be visualized along factors of interest using the plotFactorScatter function. We can use features included 
# in the model (such as IGHV or trisomy12) to color or shape the samples by. Alternatively, external covariates can also be used 
# for this purpose.

plotFactorScatter(
  MOFAobject,
  factors = 1:2,
  color_by = "IGHV",      # color by the IGHV values that are part of the training data
  shape_by = "trisomy12"  # shape by the trisomy12 values that are part of the training data
)

# Here we find again a clear separation of samples based on their IGHV status (color) along Factor 1 and by the absence 
# or prescence of trisomy 12 (shape) along Factor 2 as indicated by the corresponding factor weights in the Mutations view.

# An overview of pair-wise sctterplots for all or a subset of factors is produced by the plotFactorScatters function

plotFactorScatters(
  MOFAobject,
  factors = 1:3,
  color_by = "IGHV"
)

# A single factor can be visualised using the plotFactorBeeswarm function

plotFactorBeeswarm(
  MOFAobject,
  factors = 1,
  color_by = "IGHV"
)

### 3.2.4) Customized analysis

# For customized exploration of weights and factors, you can directly fetch the variables from the model using ‘get’ 
# functions: getWeights, getFactors and getTrainData:

MOFAweights <- getWeights(
  MOFAobject, 
  views = "all", 
  factors = "all", 
  as.data.frame = TRUE    # if TRUE, it outputs a long dataframe format. If FALSE, it outputs a wide matrix format
)
head(MOFAweights)
dim(MOFAweights)

MOFAfactors <- getFactors(
  MOFAobject, 
  factors = c(1,2),
  as.data.frame = FALSE   # if TRUE, it outputs a long dataframe format. If FALSE, it outputs a wide matrix format
)
head(MOFAfactors)
dim(MOFAfactors)

MOFAtrainData <- getTrainData(
  MOFAobject,
  as.data.frame = TRUE, 
  views = "Mutations"
)
head(MOFAtrainData)

### 4) Further functionalities
### 4.1) Prediction of views

# With the predict function, full views can be predicted based on the MOFA model with all or a subset of factors.

predictedDrugs <- predict(
  MOFAobject,
  view = "Drugs",
  factors = "all"
)[[1]]


# training data (incl. missing values)
drugData4Training <- getTrainData(MOFAobject, view="Drugs")[[1]]
pheatmap::pheatmap(drugData4Training[1:40,1:20],
                   cluster_rows = FALSE, cluster_cols = FALSE,
                   show_rownames = FALSE, show_colnames = FALSE)
# predicted data
pheatmap::pheatmap(predictedDrugs[1:40,1:20],
                   cluster_rows = FALSE, cluster_cols = FALSE,
                   show_rownames = FALSE, show_colnames = FALSE)

### 4.2) Imputation...

### 4.3) Clustering of samples based on latent factors

# Samples can be clustered according to their values on some or all latent factors using the clusterSamples function. 
# Clusters can for example be visualised using the plotFactorScatters function

set.seed(1234)
clusters <- clusterSamples(
  MOFAobject, 
  k = 2,        # Number of clusters for the k-means function
  factors = 1   # factors to use for the clustering
)

plotFactorScatter(
  MOFAobject, 
  factors = 1:2, 
  color_by = clusters
)

