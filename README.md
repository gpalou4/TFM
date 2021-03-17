# Introduction
This is the repository for the manuscript *"DNA methylation and gene expression integration in cardiovascular disease"*

Authors: Guillermo Palou Márquez, Isaac Subirana, Lara Nonell, Alba Fernández-Sanlés, Roberto Elosua

Guillermo Palou-Márquez 1,2,3 , Isaac Subirana 1,4 , Lara Nonell 5 , Alba Fernández-Sanlés 1,6 † , Roberto Elosua 1,7,8 †

1 Cardiovascular Epidemiology and Genetics Research Group, REGICOR Study group, Hospital
del Mar Medical Research Institute (IMIM), Barcelona, Spain.
2 Pompeu Fabra University (UPF), Barcelona, Spain.
3 Institute for Research in Biomedicine (IRB Barcelona), The Barcelona Institute of Science and
Technology, Barcelona, Spain
4 CIBER Epidemiology and Public Health (CIBERESP), Spain.
5 MARGenomics, Hospital del Mar Medical Research Institute (IMIM), Barcelona, Spain,
6 MRC Integrative Epidemiology Unit at the University of Bristol, Bristol, UK
7 CIBER Cardiovascular Diseases (CIBERCV), Spain.
8 Medicine Department, Medical School, University of Vic-Central University of Catalonia (UVic-
UCC), Vic, Spain.
† These two authors equally contributed to this work.

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
    2. _REGICOR_validation_/training_model.R --> Replication analysis in the independent cohort REGICOR.

-   _downstream_analysis_: Main analysis, including variance decomposition of the identified latent factors, statistical analysis and prediction improvement. 

    1. _ewas_nogroups_914_2055_ --> Main analysis (explained in the next section in more detail)
    2. MOFA_samples_metadata_nogroups_914_2055.csv --> Samples metadata for the main analysis
    3. _REGICOR_validation_ --> Same analysis but for the replication in REGICOR
    4. cpgs_shared.R --> Calculation of the number of CpGs shared between the EWAS (main) and the SD (sensitivity) strategies.

-   _samples_covariates_: Different samples metadata files for different steps of the analysis

### 2.1 Main analysis

Found in _MOFA/downstream_analysis/ewas_nogroups_914_2055_. We have three folders:

- 1) _variance_decomposition_: It includes the following: correlation between factors and covariates; correlation between factors; variance explained by each omic and factor; individual or grouped violin plots for all factors; feature weights (from both omics) for all factors; heatmaps only for the interesting factors; scatterplots (for quantitative variables) or boxplots (for qualitative variables) of factors vs covariates; TSNE of all individuals colored by different covariates.
The different resulting plots are either on the same folder or on their respective folders.

- 2) _statistical_tests_: It includes the following: t-tests and Mann-Whitney tests for all factors vs CVD (factors_cvd_tests.csv) and CHD incidence (factors_chd_tests.csv); T-test (heatmap_features_T_tests.csv) or cox regression (heatmap_features_cox_regression.csv) for the top 30 CpGs leading the heatmaps using CVD as outcome and adjusted by cell type proportions and one surrogate variable; simple histograms (factors_histograms) and boxplots (factors_boxplots) for all factors vs CVD; simple histograms (features_histograms) and boxplots (features_boxplots) for the top 30 features of the interested factors plus their correlations (features_correlations).

- 3) _prediction_: It includes the following:
    1. Cox regression (discrimination analysis) --> Done for all the factors with the `survival` R package for either CVD (cox_results_cvd.csv) or CHD (cox_results_chd.csv).
    2. NRI and clinical NRI (reclassification analysis) --> Done for the significantly associated factors with the `Hmisc` R package. Results can be found in the _results_reclassication_ folder.
    3. Kaplan meier plots for the interesting factors.
    4. _factor_cov_ --> It includes the samples metadata including the factor values, separately for each factor (temp files just to speed the analyses).


