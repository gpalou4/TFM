
plotbyBatch <- function(ExonFeatureSet, batches, directory, plot.type){
  
  # Check that plot.type is provided
  plots.avail <- c("boxplot", "density", "RLE", "NUSE")
  plots.match <- charmatch(plot.type, plots.avail, nomatch = 0)
  if (plots.match == 0)
    stop("plot.type must be 'boxplot' or 'density' \n")

    lapply(batches, function(batch) {
    
    # Obtain samples names for a given batch
    samples.names <- rownames(ExonFeatureSet@phenoData@data[ExonFeatureSet@phenoData@data$PLATE==batch,]) 
    # Filter a subset of the EFSet with the given samples names
    EFSet.subset <- ExonFeatureSet[,samples.names]
    # BOXPLOT
    if (plots.match == 1) {
      pdf(paste(directory,"/boxplot_batch_",batch,".pdf", sep = ""))
      oligo::boxplot(EFSet.subset, target = "core", transfo=log2, nsample=10000,
                   main = paste("Boxplot for batch ",batch, sep = ""))
      dev.off()
    }
    # DENSITY PLOT
    else if (plots.match == 2) {
      pdf(paste(directory,"/density_batch_",batch,".pdf", sep = ""))
      oligo::hist(EFSet.subset, target = "core",
                  main = paste("Density plot for batch ",batch, sep = ""))
      dev.off()
    }
    # RLE PLOT
    else if (plots.match == 3) {
      fit.PLM <- fitProbeLevelModel(EFSet.subset)
      pdf(paste(directory,"/RLE_batch_",batch,".pdf", sep = ""))
      oligo::RLE(fit.PLM, target = "core",
                  main = paste("Relative Log Expression plot for batch ",batch, sep = ""))
      dev.off()
    } # NUSE PLOT
    else if (plots.match == 4) {
      fit.PLM <- fitProbeLevelModel(EFSet.subset)
      pdf(paste(directory,"/NUSE_batch_",batch,".pdf", sep = ""))
      oligo::NUSE(fit.PLM, target = "core",
                 main = paste("Normalized Unscaled Standard Errors plot for batch ",batch, sep = ""))
      dev.off()
    }
      
  })
  cat("PLOTS DRAWN", file = stdout())
}


