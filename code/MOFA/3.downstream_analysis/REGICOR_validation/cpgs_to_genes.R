rm(list=ls())
start_time <- Sys.time()
date()

cat("\n\n... SETTING WORKING DIRECTORY ...\n\n", file = stdout())
setwd("/projects/regicor/Guille/TFM")
getwd()

# Get top 30 CpGs from REGICOR factors
meth.features.heatmap.name <- load("MOFA/3.downstream_analysis/regicor/variance_decomposition/meth_30features_heatmap.RData")
meth.features.heatmap.regicor <- get(meth.features.heatmap.name)

# Get top 30 CpGs from FOS factors
meth.features.heatmap.name <- load("MOFA/3.downstream_analysis/ewas_nogroups_914_2055/variance_decomposition/meth_30features_heatmap.RData")
meth.features.heatmap.fos <- get(meth.features.heatmap.name)

# Get EPIC Illumina Manifest
Illumina.EPIC.manifest <- read.csv("/projects/regicor/data/REGICOR/methylation/epic/manifest/MethylationEPIC_v-1-0_B4.csv",
                                   skip = 7, stringsAsFactors=F, header=T, sep=",", quote="")
rownames(Illumina.EPIC.manifest) <- Illumina.EPIC.manifest$IlmnID
head(Illumina.EPIC.manifest)
dim(Illumina.EPIC.manifest)

# Get 450k Illumina Manifest
Illumina.450k.manifest <- read.csv("methylation/phenotype/phg000492.v1.FHS_DNAMethylation.marker-info.MULTI/HumanMethylation450_15017482_v.1.1.csv",
                                   skip = 7, stringsAsFactors=F, header=T, sep=",", quote="")
rownames(Illumina.450k.manifest) <- Illumina.450k.manifest$IlmnID
head(Illumina.450k.manifest)
dim(Illumina.450k.manifest)

# Check directly the CpGs from FOS in REGICOR

for (i in 1:30) {
  print(paste("REGICOR factor: ",i, sep = ""))
  regicor.cpgs <- meth.features.heatmap.regicor[[i]]
  for (f in 1:4) {
    print(paste("###### FOS factor: ",f, sep = ""))
    fos.cpgs <- meth.features.heatmap.fos[[f]]
    print(fos.cpgs[fos.cpgs%in%regicor.cpgs])
    
  }
  print("###################")
}


# Transform CpGs into Genes for FOS
for (i in 1:4) {
meth.features.heatmap.fos.genes <- Illumina.450k.manifest[rownames(Illumina.450k.manifest)%in%meth.features.heatmap.fos[[i]],c("IlmnID","UCSC_RefGene_Name")]
meth.features.heatmap.fos.genes <- meth.features.heatmap.fos.genes[order(rownames(meth.features.heatmap.fos.genes)),]
meth.features.heatmap.fos[[i]] <- meth.features.heatmap.fos.genes$UCSC_RefGene_Name
meth.features.heatmap.fos[[i]] <- meth.features.heatmap.fos[[i]][!meth.features.heatmap.fos[[i]] == ""]
}

# Transform CpGs into Genes for REGICOR
for (i in 1:30) {
  meth.features.heatmap.regicor.genes <- Illumina.EPIC.manifest[rownames(Illumina.EPIC.manifest)%in%meth.features.heatmap.regicor[[i]],c("IlmnID","UCSC_RefGene_Name")]
  meth.features.heatmap.regicor.genes <- meth.features.heatmap.regicor.genes[order(rownames(meth.features.heatmap.regicor.genes)),]
  meth.features.heatmap.regicor[[i]] <- meth.features.heatmap.regicor.genes$UCSC_RefGene_Name
  meth.features.heatmap.regicor[[i]] <- meth.features.heatmap.regicor[[i]][!meth.features.heatmap.regicor[[i]] == ""]
}

write.csv(meth.features.heatmap.regicor, file = paste("MOFA/3.downstream_analysis/regicor/cpgs_to_genes.RData", sep = ""),
          row.names = FALSE)

lapply(meth.features.heatmap.regicor, function(x) write.table( data.frame(x), "MOFA/3.downstream_analysis/regicor/cpgs_to_genes.csv"  , append= T, sep=',' ))


# Now check genes from FOS in REGICOR

for (i in 1:30) {
  print(paste("REGICOR factor: ",i, sep = ""))
  regicor.genes <- strsplit(as.character(meth.features.heatmap.regicor[[i]]), ";")
  regicor.genes.unique <- unique(unlist(regicor.genes))
  for (f in 1:4) {
    print(paste("###### FOS factor: ",f, sep = ""))
    fos.genes <- strsplit(as.character(meth.features.heatmap.fos[[f]]), ";")
    fos.genes.unique <- unique(unlist(fos.genes))
    print(fos.genes.unique[fos.genes.unique%in%regicor.genes.unique])
    
  }
  print("###################")
}


# if (length(a[[26]]) > 1) {
#   
#   for (i in 1:length(a[[26]])) {
#     
#     print(a[[26]][i]%in%meth.features.heatmap.regicor[[13]])
#     
#   }
#   
# } else {
#   
#   print(a[[26]]%in%meth.features.heatmap.regicor[[13]])
#   
#   
# }





