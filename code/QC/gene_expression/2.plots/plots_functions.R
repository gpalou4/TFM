
plotbyBatch <- function(Set, batches, directory, plot.type){
  
  # Check that plot.type is provided
  plots.avail <- c("boxplot", "MA", "RLE", "NUSE")
  plots.match <- charmatch(plot.type, plots.avail, nomatch = 0)
  if (plots.match == 0)
    stop("plot.type must be 'boxplot', 'MA', 'RLE', or 'NUSE' \n")
  # Check the type of Set provided
  sets.avail <- c("ExonFeatureSet", "ExpressionSet")
  set.class <- class(Set)[1]
  sets.match <- charmatch(set.class, sets.avail, nomatch = 0)
  if (sets.match == 0)
    stop("Set must be 'ExonFeatureSet' or 'ExpressionSet' \n")
  
  outlier.samples <- lapply(batches, function(batch) {
    
    # Obtain samples names for a given batch
    samples.names <- rownames(Set@phenoData@data[Set@phenoData@data$Batch==batch,]) 
    # Filter a subset of the EFSet with the given samples names
    Set.subset <- Set[,samples.names]
    # ORDENAR COLUMNAS ALPHABETICALLY
    # BOXPLOT
    if (plots.match == 1) {
      # Making the boxplot
      pdf(paste(directory,"/boxplot_batch_",batch,".pdf", sep = ""))
      oligo.boxplot <- oligo::boxplot(Set.subset, target = "core", transfo=log2, nsample=10000,
                           main = paste("Boxplot for batch ",batch, sep = ""), las = 2,
                           #names = substr(samples.names, nchar(samples.names)-4, nchar(samples.names)),
                           ylab = "log2-Intensity", cex.axis = 0.5)
      
      print(paste("BATCH NUMBER ",batch, sep = ""))
      print("Medians")
      
      # ExonFeatureSet
      if (sets.match == 1) {  
        sample.medians <- oligo.boxplot[[1]]$stats[3,]
        sample.names <- oligo.boxplot[[1]]$names
      } # ExpressionSet
        else if (sets.match == 2) {
          sample.medians <- oligo.boxplot$stats[3,]
          sample.names <- oligo.boxplot$names
        }
      
      # Medians
      print(sort(sample.medians))
      print("Mean for all the samples medians")
      mean <- mean(sample.medians)
      print(mean)
      
      # Cut-offs
      Q <- as.vector(quantile(sample.medians, probs=c(.05, .95), na.rm = FALSE))
      higher.cutoff <- Q[2]
      lower.cutoff <- Q[1]
      print("Thresholds to filter")
      print(paste("95% quantile (Higher cut-off): ",higher.cutoff, sep = ""))
      print(paste("5% quantile (Lower cut-off): ",lower.cutoff, sep = ""))
      
      # Outliers
      outlier.samples.index <- which(sample.medians > higher.cutoff | sample.medians < lower.cutoff)
      print("Outlier samples")
      print(sample.medians[outlier.samples.index])
      outlier.samples.names <- sample.names[outlier.samples.index]
      print("Outlier samples names")
      print(outlier.samples.names)
      
      # Add cut-off values lines on the plot
      abline(h = higher.cutoff, col = "red", lwd = 1.5)
      abline(h = lower.cutoff, col = "red", lwd = 1.5)
      
      # save plot
      dev.off()
    
      # Return outliers to a list
      outlier.samples.names

      
    }
    # 
    else if (plots.match == 2) {

      MAplot <- MAplot(Set.subset)
      # print(paste("BATCH NUMBER ",batch, sep = ""))
      # print("sample.means")
      # sample.means <- apply(oligo.density[[1]][[1]], 2, FUN = mean)
      # print(sample.means)
      # print("mean for all the samples")
      # mean <- mean(apply(oligo.density[[1]][[1]], 2, FUN = mean))
      # print(mean)
      # outlier.samples.index <- which(as.matrix(sample.means > higher.cutoff + mean | sample.means < mean + lower.cutoff))
      # print("outlier.samples.index")
      # print(outlier.samples.index)      
      # outlier.samples.names <- names(sample.means)[outlier.samples.index]
      # print("outlier.samples.names")
      # print(outlier.samples.names)
      print(MAplot)
      
      # Add mean and outliers line
      abline(v = mean, col = "blue", lwd = 2)
      abline(v = higher.cutoff + mean, col = "red", lwd = 2)
      abline(v = lower.cutoff + mean, col = "red", lwd = 2)
      
      # save plot
      dev.off()
      
      # Return outliers to a list
      outlier.samples.names
      
    }
    # RLE PLOT
    else if (plots.match == 3) {
      
      if (sets.match != 1) {
        stop("Set object must be an 'ExonFeatureSet' \n")
      }
      
      fit.PLM <- fitProbeLevelModel(Set.subset)
      pdf(paste(directory,"/RLE_batch_",batch,".pdf", sep = ""))
      oligo::RLE(fit.PLM, target = "core",
                     main = paste("Relative Log Expression plot for batch ",batch, sep = ""), las = 2,
                     names = substr(samples.names, nchar(samples.names)-4, nchar(samples.names)),
                     cex.axis = 0.5)
      # Medians
      print(paste("BATCH NUMBER ",batch, sep = ""))
      print("Medians")
      sample.medians <- oligo::RLE(fit.PLM, type='stats')[3,]
      print(sort(sample.medians))
      
      # Cut-offs
      Q <- as.vector(quantile(sample.medians, probs=c(.05, .95), na.rm = FALSE))
      higher.cutoff <- Q[2]
      lower.cutoff <- Q[1]
      print("Thresholds to filter")
      print(paste("95% quantile (Higher cut-off): ",higher.cutoff, sep = ""))
      print(paste("5% quantile (Lower cut-off): ",lower.cutoff, sep = ""))
      
      # Outliers
      outlier.samples.index <- which(sample.medians > higher.cutoff | sample.medians < lower.cutoff)
      print("Outlier samples medians")
      print(sample.medians[outlier.samples.index])
      outlier.samples.names <- names(sample.medians)[outlier.samples.index]
      print("Outlier samples names")
      print(outlier.samples.names)
  
      # Add cut-off values lines on the plot
      abline(h = higher.cutoff, col = "red", lwd = 1.5)
      abline(h = lower.cutoff, col = "red", lwd = 1.5)
      
      # save plot
      dev.off()
      
      # Return outliers to a list
      outlier.samples.names
      
    } # NUSE PLOT
    else if (plots.match == 4) {
      
      if (sets.match != 1) {
        stop("Set object must be an 'ExonFeatureSet' \n")
      }
      
      fit.PLM <- fitProbeLevelModel(Set.subset)
      
      pdf(paste(directory,"/NUSE_batch_",batch,".pdf", sep = ""))
      oligo::NUSE(fit.PLM, target = "core",
                  main = paste("Normalized Unscaled Standard Errors plot for batch ",batch, sep = ""), las = 2,
                  names = substr(samples.names, nchar(samples.names)-4, nchar(samples.names)),
                  cex.axis = 0.5)
      
      print(paste("BATCH NUMBER ",batch, sep = ""))
      
      ######## IMPORTANT ######## 
      # The function NUSE(fit.PLM, type='stats') does not work because of this error:
      # Error in quantile.default(newX[, i], ...) : 
      #   missing values and NaN's not allowed if 'na.rm' is FALSE
      # Basically, when calculating the normalized unscaled standard errors for each sample, some NA's can appear
      # and then, the quantile CANNOT BE performed --> thus error. We are going to update the function, and omit the NA's
      # in order to obtain the quantiles and extract the median.
      # This cannot be done: oligo::NUSE(fit.PLM, type='stats'))
      # Let's calculate the NUSE values first
      NUSE.values <- sweep(se(fit.PLM), 1, rowMedians(se(fit.PLM)), '/')
      # Now the quantiles, only changing: NUSE --> na.omit(NUSE)
      NUSE.quantiles <- apply(na.omit(NUSE.values), 2, quantile, c(0, .25, .50, .75, 1))
      print("Medians")
      sample.medians <- NUSE.quantiles[3,]
      print(sort(sample.medians))
      ######## IMPORTANT ######## 
      
      # Cut-offs
      Q <- as.vector(quantile(sample.medians, probs=c(.05, .95), na.rm = FALSE))
      higher.cutoff <- Q[2]
      lower.cutoff <- Q[1]
      print("Thresholds to filter")
      print(paste("95% quantile (Higher cut-off): ",higher.cutoff, sep = ""))
      print(paste("5% quantile (Lower cut-off): ",lower.cutoff, sep = ""))
      
      # Outliers
      outlier.samples.index <- which(sample.medians > higher.cutoff | sample.medians < lower.cutoff)
      print("Outlier samples medians")
      print(sample.medians[outlier.samples.index])
      outlier.samples.names <- names(sample.medians)[outlier.samples.index]
      print("Outlies samples names")
      print(outlier.samples.names)
      
      # Add cut-off values lines on the plot
      abline(h = higher.cutoff, col = "red", lwd = 1.5)
      abline(h = lower.cutoff, col = "red", lwd = 1.5)
      
      # save plot
      dev.off()
      
      # Return outliers to a list
      outlier.samples.names
    }
    
  })
  cat("\n\n ... ALL PLOTS DRAWN ...\n\n", file = stdout())
  return(outlier.samples)
}


