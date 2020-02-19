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

Finished the pipeline. But there are still a few things to do. I have created a a list of cross-reactive probes to remove them (Chen 2013 et al. + Illumina Manifest), but I can't find the CpGs to remove from the Illumina manifest, instead I have searched for a most recent paper ("Comprehensive characterization, annotation and innovative use of Infinium DNA methylation BeadChip probes. Zhou et al. 2016), that contains some CpGs. Chen has 30k CpGs, and Zhou 60ks, but 60ks uniquely between them.

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

mRNA input --> Error
meth --> all ok until norm + 4sd

#### 18/2/2020: Both QC pipelines

mRNA input --> still error bc of RAM... what do we do?
meth --> mrna_match error, not loading Mvalues!
     --> once completed I will start sva
     --> Started batch effect plots

#### 19/2/2020: Both QC pipelines

mRNA input --> still error bc of RAM... what do we do?

meth --> sva started
    --> removed X/Y chrs cpgs to do sva and batch effect
    --> chanching script of batch effect plots
 