# TFM project notebook

This notebook is simple a file where I write everything I do during my master thesis, on a weekly basis, since I obtained the data to analyze (11/12/2019).

### October-November

Reading articles, defining the project, writing the project proposal (you can find it on Reports folder), making other collaborations analysis within the group, learning the quality control steps and EWAS analysis of the methylion data, and beginning to understand/propose the quality control for the gene expression array data.
As for the collaboration within the group, my work is to analyze differentially methylated regions (DMRs) from the REGICOR population (+ validation in FHS and EPIC populations). The traits analyzed and associated to CVDs are: diabetes, glucose and metabolic syndrome.

### Week 9-13 December 2019

#### 11/12/2019: dbGAP approval and QC gene exp papers

We have already received approval by `dbGAP` to download the data. The next step is to know how exactly download it and where to find the gene expression array files, plus the metadata (phenotype). There are many many files, so that is not an easy task to do.
I have requested access to some ports in order to download it. While I wait the approval by the head of IT department I have beginning to search for some papers where they have already used FHS gene expression data, to look for the quality control steps.

> Paper 1: Whole blood gene expression and interleukin-6 levels. \
*The raw data were quantile-normalized and log2 transformed, followed by summarization using Robust Multi-array Average. The gene annotations were obtained from Affymetrix NetAffx Analysis Center (version 31). We excluded transcript clusters that were not mapped to RefSeq transcripts, resulting in 17,873 distinct transcripts (17,324 distinct genes) for downstream analysis.*

> Paper 2: Integromic Analysis of Genetic Variation and Gene Expression Identifies Networks for Cardiovascular Disease Phenotypes \
*mRNA expression profiling was assessed with the use of the Affymetrix Human Exon 1.0 ST GeneChip platform (Affymetrix Inc, Santa Clara, CA), which contains >5.5 million probes targeting the expression of 17 873 genes. The Robust Multi-Array Average package20 was used to normalize the gene expression values and remove any technical or spurious background variation. Linear regression models were used to adjust for technical covariates (batch, first principal component, and resid-ual of all probeset mean values).*

#### 12/12/2019: New computer + data download

Spent several hours installing the new computer with the informaticians, linux Virtual Machine, Rstudio, etc.
The IT department will be in charged of download the data for security reasons. Meanwhile I will focus on finish my collaboration within the group.

#### 13/12/2019: Collaborations --> DMRs Ven Diagrams + Juan R. G. meeting

IMIM collaboration: Made two Ven Diagrams (BMI and No BMI) for the DMRs results, to know the overlap between the traits (diabetes, glucose and metabolic syndrome)
ISGlobal collaboration: Meeting with Juan R. González to update the R scripts.

### Week 16-20 December 2019

#### 16/12/2019: Collaborations --> Lab meeting + Supplementary tables

IMIM collaboration: Lab meeting to discuss the paper and whether it is necessary or not to make the functional annotation of the CpGs/DMRs. I also began to work with the supplementary tables for the DMRs results. I run some scripts for FHS validation of EWAS results.

#### 17/12/2019: Supplementary tables

IMIM collaboration: Finished the supplementary tables.

#### 18/12/2019: Manuscript + IPA functional analysis

IMIM collaboration: Wrote a little bit in the manuscript and added the DMRs ven Diagrams into supplementary material. I also started to use the IPA software to perform pathway enrichment analysis.

#### 19/12/2019 IPA + finished manuscript

IMIM collaboration: I finished to write in the manuscript and performed the enrichmen pathway analysis for all the cases. The results weren't coherents so Alba and I decided to discard this analysis. 

#### 20/12/2019

The head of the IT department gave us the authorization to download the dbGAP data. We were able to do it but we couldn't find the gene expression nor the methylation data, only the samples' metadata. Therefore, Roberto wrote them an email asking where is the data, as it is missing. 

### CHRISTMAS HOLIDAYS 20/12/2019 - 07/01/2020

### Week 7-10 January 2020

#### 7/1/2020: MOFA R libraries

dbGAP have still not answered, thus we have decided to search for other data to work meanwhile.
I have begun to set up the R environment in my new Linux computer (download libraries, etc). I will firstly use the Rstudio from the computer with the MOFA example to see how it works in the following days.

#### 8/1/2020: MOFA

Roberto has contacted with the coordinator for the dbGAP access for our specific data and has forwarded the email to the help technicians. We are still waiting for their answer.
I have finished to download all R libraries and begun to use the MOFA example.

#### 9/1/2020: MOFA + download dbGAP data

We have finally solved the issues and we are currently downloading the data. I have found the annotation data for the mRNA expression data, to be used during the analysis.

#### 10/1/2020: Decrypting the data and Methylation QC

I have started decrypting the data but it takes so long.
Thus, meanwhile, I have started to do the pipeline for the quality control for the methylation data with Alba's EPIC dataset. The filtering step is already done and working.

### Week 13-17 January 2020

#### 13/1/2020: Looking at the downloaded files

Finished decrypting the data.
Started to look at all the files and searched what files are needed.

#### 14/1/2020: Looking at the downloaded files

I finally have had an idea of what kind of files there are, where is the metadata, etc. However, I don't find the file containing the batches, so we have written to dbGAP Helper.
I compared the samplesIDs from mRNA's metadata with the samplesIDs from methylation. Around 2300 ID's are matched, from the OFF_CSS files (1900 OFF, and 400 CSS aprox).
We need to upload the IDATs (methylation raw data) and CEL (gene expression raw data) to the cluster, but we need IT's help.

#### 15/1/2020: Methylation QC pipeline

continue QC meth pipeline --> sex filtering, searching sex metadata, obtained a dataframe with sex and samples IDs. Everything works so far with the EPIC data until step 4) Sex. B and M values obtained (step 5).

#### 16/1/2020:

No work

#### 17/1/2020: Methylation QC pipeline

Finished the pipeline. But there are still a few things to do. I have created a a list of cross-reactive probes to remove them ("Discovery of cross-reactive probes and polymorphic CpGs in the Illumina Infinium HumanMethylation450 microarray", + Illumina Manifest), but I can't find the CpGs to remove from the Illumina manifest, instead I have searched for a most recent paper ("Comprehensive characterization, annotation and innovative use of Infinium DNA methylation BeadChip probes. Zhou et al. 2016), that contains some CpGs. Chen has 30k CpGs, and Zhou 60ks, but 60ks uniquely between them.

I also can begin now to try to upload just a few samples from my data and use them to run all the pipeline and change whatever I need to change.
I have started doing the input reading part.

### Week 20-24 January 2020

#### 20/1/2020: gene expression QC pipeline

I have started to make the Quality Control for the gene expression data with 5 samples, and doing all the plots by batches (created batches because I don't have the real metadata yet).

#### 21/1/2020: Methylation QC pipeline

I finally was able to integrate metadata with the IDATs files, there was an error happening: "Slide" and "Array" fields in the metadata file must be filled in all samples for `minfi` to work properly when reading the data. I had to filled these fields because most of the samples didn't have information. I tried this with an example of 5 individuals and it works.

#### 22/1/2020: Methylation QC pipeline

I have adapted the pipeline to the 450k subset real data (5 samples). All the steps works so far. Now I would like to try it out with all the samples but I am still waiting for the data to be uploaded on the cluster...

#### 23/1/2020: gene expression batch effect

Checked papers about how to deal with the batch effect in gene expression. `ComBat` is the best method.
> Removing Batch Effects in Analysis of Expression Microarray Data: An Evaluation of Six Batch Adjustment Methods
*Of the batch adjustment methods we evaluated, we found that the empirical Bayes algorithm implemented in ComBat was best able to reduce and remove batch effects while increasing precision and accuracy*
However, it turns out `MOFA` can deal with batch effects by inserting covariates into the model. It is not ideal, authors say, but it is okey. We shall discuss which method to use.
I found the batch metadata from within a file.

#### 24/1/2020: TFM proposal comments meeting + Methylation normalization

During the meeting we discussed how to deal with the comments presented by the reviewers on the thesis proposal. 
1) REVIEWER C: Batch effect on Methylation data --> We have decided to perform the normalization recommended, with `Dasen`, and not to do the standardization. We will also use SVA to help deal with the batch effect.
2) REVIEWER E: About the relative degree, we will probably try to download the genomic data, perform a PCA, and use the first 10 PCs as covariates in the MOFA model. To avoid overfitting we can't have an independent cohort, but we might use a bootstrap strategy and perform MOFA multiple times on diferent subsamples and create an average of the results somehow.
3) REVIEWER D: Missing information of the individuals. MOFA can deal with missing information, it just ignores it. And there are no missing individuals in our dataset.
I have already done the normalization with `Dasen` for the methylation data.

### Week 27-31 January 2020

#### 27/1/2020: Singularity

Alba helped me to create a singularity image through the web application. But there were a lot of errors.

#### 28/1/2020: Singularity + real metadata gene expression

I finally solved the errors and created a singularity image.
I also swapped the invented metadata of the gene expression for the real one.

#### 29/1/2020: Upload data to the cluster

As I have been waiting so long for the data to be uploaded I decided to manually do it myself. I spent the whole day uploading the gene expression raw files.

#### 30/1/2020: Upload data to the cluster

I spent all day uploading the methylation raw files.

#### 31/1/2020: Upload data to the cluster + Alba PhD defense

I finished to upload all raw files.
Alba PhD defense.

### Week 3-7 February 2020

#### 3/2/2020: Methylation QC pipeline

I have begun running the methylation QC pipeline. There is an error uploading the raw data, I think it must be because of RAM memory problems.

#### 4/2/2020: Methylation QC pipeline

I tried to solve the error but it persists.
The cluster did not work correctly.

#### 5/2/2020: Methylation QC pipeline

I reduced the number of individuals to analyze by matching gene expression raw files with methylation raw files. Now only 1786 individuals can be analyzed (from 2600 initially) --> RAM error solved

#### 6/2/2020: Methylation QC pipeline

Alba sent me the samples metadata with the covariates and all the stuff.
I continued with the QC pipeline, several errors encountered.

#### 7/2/2020: Methylation QC pipeline

I corrected the sex script.
The pipeline works until normalization with Dasen.
Meged samples metadata from Alba and my files.

### Week 10-14 February 2020

#### 10/2/2020: Methylation QC pipeline

Rerunning all the pipeline because some individuals have missing CVD data (from Alba's metadata file). Now 1786 --> 1399 individuals to analyze.
Meanwhile, I continue the pipeline with the previous files, now doing SVA analysis.
Cluster has stopped to work, I cannot work anymore today.

#### 11/2/2020: Matching gene expression and methylation samples metadata

I finally matched correctly the samples metadata for both omics, gene expression and methylation --> 1549 individuals shared
I have begun the methylation QC pipeline with the 1549 individuals.

#### 12/2/2020: Both QC pipelines

Continuing with methylation QC pipeline
Began the gene expression QC pipeline (reading input)
Alba said it is better to perform the QC with ALL samples and after, do the matching. So, I created two new samples metadata files, one for each omics (2620 for methylation and 1818 for gene expression) and I'll begin to run methylation QC pipeline again, hopefully for the last time.

#### 13/2/2020: Methylation QC pipeline

Started the methylation QC pipeline for the 2620 individuals, this time without getting a RAM memory error.

#### 14/2/2020: Methylation QC pipeline

Followed the methylation QC pipeline. I had an error in filterting step so I divided the step in two substeps and it worked.

### Week 17-21 February 2020

#### 17/2/2020: Both QC pipelines

There was an error trying to upload the mRNA data to a ExpressionFeatureSet object.
Methylation QC worked fine until normalization steps.

#### 18/2/2020: Both QC pipelines

The mRNA error persisted, it seems is a RAM problem.
Cluster errors, couldn't load M values from methylation data. I kept with the QC methylation pipeline, doing SVA and batch effect plots.

#### 19/2/2020: Both QC pipelines

The mRNA error persisted, it seems is a RAM problem.
Removal of X/Y chrs CpGs for SVA, and batch effect plots, changing the scripts.

#### 20/2/2020: Both QC pipelines

Didn't go

#### 21/2/2020: Both QC pipelines
 
Talked to the group about the mRNA computational problem.
Problem with NA's (after SD outliers converted to NA), implies only ~30.000 CpGs are used for SVA. We can't do much more!
Finishing batch effect plots

### Week 24-28 February 2020

#### 24/2/2020: mRNA computational

I estimated the computational cost in RAM needed to upload the 1820 samples from mRNA. We have 255GB, we need around 400 GB. Solution?

#### 25/2/2020: QC mRNA

Lab meeting
Doing scripts for the plots before normalization

#### 26/2/2020: parallel computing + gene expression QC

Alternative to upgrade the cluster memory: parallel computing --> I talked with Joan Marc, he said it's difficult but feasible and he has to check it. I will persist reminding him.
A merge of 2 ExpressionFeatureSets doesn't work, the function accepts ExpressionSets (so, after the QC..., not very convenient)
Correcting plots scripts

#### 27/2/2020: gene expression QC

A lot of problems with obtaining a threshold for the plots, in order to remove bad quality samples.
NUSE prolems solved, plot density thresholds still unkown.

#### 28/2/2020: gene expression QC

Normalized plots of boxplots and MA. (RLE, NUSE not necessary,as it is already normalized)
Filtering samples using the thresholds
Filtering genes (with some functions + genefilter).

### Week 2-6 March 2020: parallel + batch effect QC mrna + begin mofa

#### 2/3/2020:

No work

#### 3/3/2020: batch effect QC mrna + begin MOFA

Finished batch effect plots.
Samples metadata subset and match between mRNA and methylation data. It is already done, getting ready to start MOFA tomorrow.
I rerun SVA for the QC methylation.

#### 4/3/2020: Begin MOFA

QC meth: SVA found 1 sva variable
MOFA: singularity MOFA image, lleva muchos tries hasta que funciona. De momento no workea.

#### 5/3/2020: MOFA

MOFA: started writing the MOFA scripts using an in-house image. The training function doesn't work because of the image.

#### 6/3/2020:  Meeting

Meeting + couldn't work that day

### Week 9-13 March 2020: MOFA

#### 9/3/2020: MOFA image + training

Tried several images, none of them working properly. Then I used the Dockerfile the MOFA authors provide but still didn't work because of the reticulate package.
I wrote the question on the SLACK MOFA group.
Finally solved it, I created my own Dockerfile using theirs, + added some ubuntu dependencies that were needed for the Reticulate package.
The training function is already running.

#### 10/3/2020: Reduction Methylation dimensions

Several strategies are to be done to reduce methylation matrix dimensions:
1) Obtain the 20.000 most variable CpGs, ranked by SD.
2) Perform an EWAS to CVD and obtain the 20.000 CpGs with the lowest p-value
3) Obtain DMRs
4) Merge CpGs from the same gene
5) Select CpGs that we know are already related to CVD

I did 1) and 2)

#### 11/3/2020: MOFA analysis

COVID-19 --> WORKING FROM HOME
MOFA plots + looking the results

#### 12/3/2020: MOFA analysis

MOFA plots + looking the results


#### 13/3/2020: MOFA analysis

More plots, more results...

### Week 16-20 March 2020: COVID-19 ALARM STATE + MOFA analysis

Due to the COVID-19 situation I will be working from home (since 11th march).
This week I performed two MOFA models (with groups and without groups) and corrected (chatting with the MOFA author Ricard) some plots that were not working, among other things.

Meeting (19/3/2020):

    1) Check CpGs in interested factors (10/13 groups, 14 nogroups) and merge into genes. Also check CpGs/mRNA and look for differences in expression between CVD and no CVD (Alba's code).
    2) Obain factor values per individual and perform an average and median between CVD and no CVD people. T-test + U Man-Wittney.
    3) GSEA: pick most expressed transcripts + perform GSEA nogroups
    4) Perform models without covariables.
    5) Correlation between factors and covariables and perform clustering of the factors using that covariables.
    6) Upload the image
    7) EWAS: Alba will check

### Week 23-27 March 2020: MOFA analysis

#### 23/3/2020: MOFA analysis

I did the following: 

1) MOFA models with no covariables 
2) Checked top20  CpGs/mRNA from heatmap in interested factors (10/13 in groups, 14 in nogroups) and looked for differences between CVD and no CVD, using
T-test and Mann-Whitney test. 
3) Identified clusters on heatmap in these factors of interest, and check EWAS catalog + affymetrix annotation (not finished)
4) T-test and Mann-Whitney test for factors 10/13 (done between CVD-no CVD, but should it be between individuals of the same group?, i.e. 2 tests?)
results weird?

#### 24/3/2020: MOFA analysis

1) MOFA models with no covariables: I asked Ricard about that. It seems MOFA does not use covariables. Answer: "covariates cannot be used for model training (see FAQ point 5.3). You can either
    (1) regress them out before training if they are undesired variation"
    (2) add them to metadata of the MOFA object after training if they are useful to characterise the factors"
2) Finished the annotation of features from the clusters heatmap.
3) Performed a correlation between features from factors 10,13 and 14. Matrix of Spearman correlations, diagonal of 1.

***Ricard: just as a reminder: when using the multi-group framework the features are centered to zero-mean per group. That means that you are assuming that the two groups are like “batches”. The aim is not to find a factor that essentially separates the two groups, but rather to find out which  sources of variation are common or unique to the different groups***


#### 25/3/2020: MOFA analysis

1) Features that overlap between factors
2) Clustering plots of factors-covariable written in the cluster script (not in local)

#### 26/3/2020: MOFA analysis + Lara meeting

Lara meeting, about mRNA QC:

    1) ArrayQualityMetrics: 3+/6 plots mal == outlier.
    2) Por encima Percentil 90 y por debajo percentil 10. Estudiar distribución mediana de invidiuos, individuos con ediana por encima percentil 95 o por debajo de percentil 5 son outliers, Hacer eso para cada plot (boxplot, NUSE...). Outliers en 3+/6 plots son muestras a eliminar.
    3) Estructura cluster entre diferentes batch. Hierarchichal clustering
    4) limpiar batch con Combat
    5) Hacer mi MDS plot pero con batch degradado

1) For the groups model: On factor 10, obtained a categorical variable that splits CVD group in two. Performed t-test (for continous variables) and
chi-square tests (categorical variables) between the obtained variable and covariables.

#### 27/3/2020: MOFA analysis + QC mRNA

Perform the same as yesterday but with factor 14.

### Week 30-03 March-April 2020: MOFA analysis + mRNA QC

#### 30/3/2020: MOFA analysis

I performed a t-test between factor 10,14 values from CVD group and CHD, and finished the tests between CVD subgroups and all covariables.

#### 31/3/2020: QC mRNA

I did again the QC of the mRNA as Lara suggested: distribution of medians for every plot (NUSE, RLE, boxplot and MA); individuals with a median above 95% quantile or below 5% quantile are considered suspicious. If the same sample is an outlier in minimum 2/4 plots, this sample is excluded. MA-plots were visualized manually as median data cannot be extracted from the function value. I checked, one by one, that all outlier samples (90) are real outliers in their respective plots
I also used ComBat to remove batch effect and plotted the MDS plot before the removal with a gradient legend.

#### 1/4/2020: MOFA analysis

I created two new MOFA models:
-With the corrected individuals from the QC from yesterday, 934 x 934.
-With all samples (using missing individuals), 934 individuals with gene expression data and 2080 individuals with methylation data. 

#### 2/4/2020: MOFA analysis: prediction and survival analysis

Performed the prediction analysis for the Factor 14 from the old model (nogroups_984): Cox regression, kaplan meier plot and discrimination analysis. Reclassification gave me error.

#### 3/4/2020: MOFA analysis

Meeting + correction of scripts from prediction analysis, and new models.
I created two submodels (grups and nogroups) for each model, so in total there are 4 models: groups_934, nogroups_934, groups_934_2080, nogroups_934_2080.

### Week 06-10 April 2020: MOFA analysis with new models

#### 6/4/2020: MOFA analysis

Changed and adapted the scripts for the new models. I finished the first two: nogroups_934, nogroups_934_2080

#### 7/4/2020: MOFA analysis

Same as yesterday but for the groups models. The groups_934_2080 is incorrect, I have rerunned it again, and I will finish it tomorrw
I performed the survival analysis, etc, for the models.

#### 8/4/2020: MOFA analysis

There has been another problem with groups_934_2080 model, related with sample names ordering. I have spent all day fixing it, and now it's running again.
Doing the singularity image for the survival analysis.
Related with the previous error, I found out that the batch effect was not performed in the 938 models due to the samples ordering.  

#### 9/4/2020: MOFA analysis

I had to update the models: 938 models corrected by batch effect; removed mRNA data from mRNA batch 15 (were clustered weirdly in the MOFA plots); and removed
the genes from ChrX/Y from gene expression data.

#### 10/4/2020: MOFA analysis

I created the models again and began to run and correct the scripts accordingly.

### Week 13-17 April 2020: MOFA analysis with updated models + manuscript writing

#### 13/4/2020: MOFA analysis

Analysis of the updated models. The mRNA batch 15 surprisingly was still producing weird results, eventhough we only had methylation data from that individuals. We decided to remove the complete batch from both methylation and gene expression (already done) data. Thus, I began to update, create and run the new models.

#### 14/4/2020: MOFA analysis + manuscript writing (methods)

New models: nogroups_914, groups_914, nogroups_914_2055, groups_914_2055.
Still running the models.
Began to write the materials and methods section of the manuscript

#### 15/4/2020: MOFA analysis + manuscript writing (methods)

Analysis of the updated models
Meeting

#### 16/4/2020: MOFA analysis + manuscript writing (methods)

#### 17/4/2020: MOFA analysis + manuscript writing (methods)

### Week 20-24 April 2020: MOFA analysis with updated models + manuscript writing

#### 20/4/2020: MOFA analysis + manuscript writing (methods)
#### 21/4/2020: MOFA analysis + manuscript writing (methods)

Ricard meeting