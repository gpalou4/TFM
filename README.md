# Introduction
This is the repository for the manuscript *"DNA methylation and gene expression integration in cardiovascular disease"*

Authors: Guillermo Palou Márquez, Isaac Subirana, Lara Nonell, Alba Fernández-Sanlés, Roberto Elosua

**Abstract**:

*Background:* The integration of different layers of omics information is an opportunity to tackle the complexity of cardiovascular diseases (CVD) and to identify new predictive biomarkers and potential therapeutic targets. Our aim was to integrate DNA methylation and gene expression data in an effort to identify biomarkers related to cardiovascular disease risk in a community-based population. We accessed data from the Framingham Offspring Study, a cohort study with data on DNA methylation (Infinium HumanMethylation450 BeadChip; Illumina) and gene expression (Human Exon 1.0 ST Array; Affymetrix). Using the MOFA2 R package, we integrated
these data to identify biomarkers related to the risk of presenting a cardiovascular event.

*Results:* Four independent latent factors (9, 19, 21 –only in women– and 27), driven by DNA methylation, were associated with cardiovascular disease independently of classical risk factors and cell type counts. In a sensitivity analysis, we also identified factor 21 as associated with CVD in women. Factors 9, 21 and 27 were also associated with coronary heart disease risk. Moreover, in a replication effort in an independent study three of the genes included in factor 27 were also present in a factor identified to be associated with myocardial infarction (CDC42BPB, MAN2A2 and RPTOR). Factor 9 was related to age and cell type proportions; factor 19 was related to age and B cells count; factor 21 pointed to human immunodeficiency virus infection-related pathways and inflammation; and factor 27, related to lifestyle factors such as alcohol consumption, smoking, and body mass index. Inclusion of factor 21 (only in women) improved the discriminative and reclassification capacity of the Framingham classical risk function and factor 27 improved its discrimination.

*Conclusions*: Unsupervised multi-omics data integration methods have the potential to provide insights into the pathogenesis of cardiovascular diseases. We identified four independent factors (one only in women) pointing to inflammation, endothelium homeostasis, visceral fat, cardiac remodeling, and lifestyles as key players in the determination of cardiovascular risk. Moreover, two of these factors improved the predictive capacity of a classical risk function.

# Repository organization

This repository contains three folders:

-**Code** --> All the R code used for the quality control of the data, MOFA integration and statistical / prediction analysis.

-**Manuscript** --> Supplementary material of the manuscript

-**Notes** --> Weekly Notebook (not updated) for progress reviews. Started on 11/12/2019

# Code 

The _Code_ folder is divided in the following structure:

## 1. QC:

Quality control for both omics data. For a more detailed and visual overview check the supplementary figure S1 (for methylation) and S2 (for gene expression), provided in the *Manuscript* folder.

-   _gene_expression_: Quality control for the gene expression data. Steps performed: 

    1. Reading the input files.
    2. Visualizing several plots (NUSE, RLE, MA, boxplots and density plots) for sample filtering.
    3. Normalization with `oligo` package.
    4. Redoing the plots with the normalized data.
    5. Filtering low-expressed genes.
    6. Batch effect check with MDS plot and removal with `comBat`

-   _methylation_: Quality control for the methylation data. Steps performed: 

    1. Reading the input files.
    2. Filtering of bad quality probes/samples and cross-reactive probes using `Minfi` package.
    3. Sex control , removing samples with a different matched-reported sex. M-values determination
    4. Normalization with `Dasen` package, removal of outliers CpGs and Blood Cell count with `FlowSorted.Blood.450k` package.
    5. Batch effect check with MDS plot.
    6. Identification of one surrogate variable via `SVA` package.
    7. Reducing dimensions of the methylation matrix with two strategies:
 
        -EWAS --> Selecting the top 20.000 CpGs with lowest p-value in association with CVD (EWAS).
        
        -SD --> Selecting the top 20.000 CpGs with highest variability using standard deviation.
    
## 2. MOFA

MOFA integration analysis of the two omic layers (gene expression and DNA methylation) using `MOFA2` R package.

-   _create_mofa_object_: Integrating both Omic datasets and the samples covariates into a MOFA object

-   _training_model_: Training of the model with previously build MOFA object.

    1. training_model.R --> Training model of the main analysis
    2. REGICOR_validation/training_model.R --> Replication analysis in the independent cohort REGICOR.

-   _downstream_analysis_: Main analysis, including variance decomposition of the identified latent factors, statistical analysis and prediction improvement. 

    1. ewas_nogroups_914_2055 --> Main analysis (explained in the next section in more detail)
    2. MOFA_samples_metadata_nogroups_914_2055.csv --> Samples metadata for the main analysis
    3. REGICOR_validation --> Same analysis but for the replication in REGICOR
    4. cpgs_shared.R --> Calculation of the number of CpGs shared between the EWAS (main) and the SD (sensitivity) strategies.

-   _samples_covariates_: Different samples metadata files for different steps of the analysis

### 2.1 Main analysis

Found in MOFA/downstream_analysis/ewas_nogroups_914_2055. We have three folders:

- 1) variance_decomposition --> It includes the following: correlation between factors and covariates; correlation between factors; variance explained by each omic and factor; individual or grouped violin plots for all factors; feature weights (from both omics) for all factors; heatmaps only for the interesting factors; scatterplots (for quantitative variables) or boxplots (for qualitative variables) of factors vs covariates; TSNE of all individuals colored by different covariates.
The different resulting plots are either on the same folder or on their respective folders.

-

-


