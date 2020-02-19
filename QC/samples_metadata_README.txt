METHYLATION samples metadata:

Initial metadata files:
samples_metadata_covariates: file obtained directly by Alba Fernandez, containing covariates information for 2123 samples. 
The samples missing (2620-2123 = 497) are the ones having NA's in CVD incidence.
phs000724.v2.pht004247.v1.p9.c1.Framingham_DNA_Methylation_Sample_Attributes_II.HMB-IRB-MDS_clean.txt: contains batch information for 2620 individuals.

Initial raw data files:
2620 Raw data sample files in raw_data folder

Steps followed to process the samples metadata:

1) Matched raw data file names (2620) with the batch metadata file: 2620 in 2620 --> all matched
2) Matched raw data file names (2620) with the covariates metadata file: 2123 in 2620 --> All matched, but 497 raw data files do not have its corresponding metadata, thus,
will have to be removed from the analysis.
3) Merge both metadata files based on the `LABID`.

Results: merged_samples_metadata.csv file with 2123 individuals and 29 columns of metadata information.


GENE EXPRESSION samples metadata:

Initial metadata files:
samples_metadata_covariates (same as before): file obtained directly by Alba Fernandez, containing covariates information for 2123 samples.
phs000363.v17.pht002944.v4.p11.c1.Framingham_SABRe_Sample_Attributes.HMB-IRB-MDS_clean.txt: contains ID information for 5373 individuals (OFFs + CCS + GENIII). This file is necessary
for mergin batch and covariates information using the ID, because batch information does not have this ID.
MasterFile_Gene_OFF_2446.txt: contains batch information for 2446 individuals (OFFs)

Initial raw data files:
2311 Raw data sample files in raw_data folder
- 419 are Case-Control samples (we discard them)
- 1892 are from the offspring generation

Steps followed to process the samples metadata:

1) Matched raw data file names with the batch metadata file: 2446 in 1892 --> 1818 TRUE, 628 FALSE. The 74 missing individuals (1892-1818) are not OFFs (e.g. FFS_EX12...)).
The rest will probably be from CCS.
2) Matched raw data file names with the IDs metadata file: 5373 in 1892 --> 1818 TRUE, 9008 FALSE. Same as before, good!
3) Merge both metadata files based on the `cel_file` (name of the .CEL files) --> 1818 individuals
4) Matched raw data file names with covariates metadata file: 2123 in 1818: 1549 TRUE, 574 FALSE. Ok, the 269 missing individuals (1818-1549), will be individuals with CVD incidence missing
because they do not match with this metadata, the original having 2620 individuals.
5) Merge previous merged dataframe with the covariates metadata file --> Results in 1549 individuals

Results: merged_samples_metadata.csv file with 1549 individuals and 30 columns of metadata information.


MATCH BETWEEN BOTH OMICS samples metadata:

Now we must match both samples metadata files from both omics: 2123 in 1549 --> 1549 TRUE, 574 FALSE. Thus, we remove 574 individuals from the methylation samples metadata and
we finally have a samples metadata for each omic with the same samples.

Methylation: 2123 individuals --> 1549 individuals (same ID)
Gene expression: 1549 individuals (same ID)


