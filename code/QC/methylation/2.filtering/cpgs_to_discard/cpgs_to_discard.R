rm(list=ls())

date()

library(minfi)
cat("\n\n minfi LOADED \n\n", file = stdout())

cat("\n\n... SETTING WORKING DIRECTORY ...\n\n", file = stdout())
setwd("/projects/regicor/Guille/TFM/methylation/QC")
getwd()

cat("\n\n... LOADING RGChannelSetExtended OBJECT ...\n\n", file = stdout())
RGCSext <- load("1.read_input/RGChannelSetExtended.RData")
RGCSext <- get(RGCSext)
Methylset <- preprocessRaw(RGCSext)
cat("\n // Methylset before removing \\ \n", file = stdout())
Methylset

# Chen et al. 2013: Around 30ks CpGs
chen2013.cpgs <- read.csv("2.filtering/cpgs_to_discard/chen2013_cpgs_to_discard.csv", 
                          header=TRUE, stringsAsFactors = F, sep=";")
cat("\n Chen 2013 CpGs \n", file = stdout())
length(chen2013.cpgs$TargetID)

# Zhou et al. 2016: Around 60ks CpGs
zhou2016.cpgs <- readRDS("2.filtering/cpgs_to_discard/zhou2016_cpgs_to_discard.rds")
cat("\n Zhou 2016 CpGs \n", file = stdout())
table(zhou2016.cpgs@elementMetadata$MASK_general)[2]

zhou.cpgs.names <- names(zhou2016.cpgs[which(zhou2016.cpgs@elementMetadata$MASK_general),])
chen.cpgs.names <- chen2013.cpgs$TargetID

# 60ks Unique CpGs to remove
cpgs.discard <- unique(c(zhou.cpgs.names, chen.cpgs.names))
cat("\n Total CpGs to discard \n", file = stdout())
length(cpgs.discard)

# cpgs.discard.name <- load("2.filtering/cpgs_to_discard/cpgs_to_discard.txt")
# cpgs.discard <- get(cpgs.discard.name)
# head(cpgs.discard)

write(cpgs.discard, file = "2.filtering/cpgs_to_discard/cpgs_to_discard.txt")

Methylset <- Methylset[-which(Methylset@NAMES%in%cpgs.discard),]
cat("\n // Methylset after removing \\ \n", file = stdout())
Methylset
