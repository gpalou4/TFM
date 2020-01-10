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
Thus, meanwhile, I have started to do the pipeline for the quality control for the methylation data with Alba's EPIC dataset.

### Week 13-17 January 2020

#### 13/1/2020:

