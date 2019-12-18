#' @title identifies empty droplets using DropletUtils::emptyDrops
#' @description distinguishes between droplets containing cells and ambient RNA in a droplet-based single-cell RNA sequencing experiment.
#' @param sce SingleCellExperiment object. Must contain a raw counts matrix.
#' @param lower Numeric. The lower bound on the total UMI count, at or below which all barcodes are assumed to correspond to empty droplets.
#' @return SingleCellExperiment object with the emptyDrops output stored in coldata: "Total", "LogProb", "PValue", "Limited", "FDR"; and a column "sample".
#' @examples
#' sce <- runDropletUtilsEmptyDrops(sce=sce, lower=1000)
#' @export
#' @import DropletUtils
runDropletUtilsEmptyDrops <-
  function(sce, lower = 1000) {
    # check if required packages aren't installed
    list.of.packages <- c("DropletUtils", "SingleCellExperiment")
    new.packages <-
      list.of.packages[!(list.of.packages %in% installed.packages()[, "Package"])]
    
    # install packages if needed
    if (length(new.packages))
      install.packages(new.packages)
    
    # load the packages
    library(DropletUtils)
    library(SingleCellExperiment)
    library(tidyverse)
    
    sce.matrix <-
      tryCatch(
        expr = {
          as.matrix(assay(sce))
        },
        warning = function(w) {
          print("Warning: ", w)
        },
        error = function(e) {
          stop("Error: sce assay cannot be converted to matrix. \n\r", e)
        }
      )
    if (!is.numeric(sce.matrix)) {
      stop("Error: sce assay is not numeric")
    }
    rm(sce.matrix)
    if (is.null(sce@colData$sample)) {
      sce@colData$sample <- 1
    }
    
    # get all sample ids in the sce
    allNamesSamples <- sort(unique(sce@colData$sample))
    
    for (i in allNamesSamples) {
      sceSIds <- (which(sce@colData$sample == i))
      sceSamp <- sce[, sceSIds]
      
      # run DropletsUtils::emptyDrops on the count matrix
      output <- emptyDrops(m = assay(sceSamp), lower = lower)
      sceSamp@colData <- output
      sceSamp@colData$sample <- rep(i, nrow(sceSamp@colData))
      
      # if iterating for the first time create an empty df with colnames of the emptyDrops output + "sample" column
      if (which(allNamesSamples == i) == 1) {
        names <- colnames(output)
        names[length(names) + 1] <- "sample"
        colData <-
          setNames(data.frame(matrix(
            ncol = ncol(output) + 1, nrow = 0
          )), names)
      }
      # append the emptyDrops output to colData df
      colData <- rbind(colData, sceSamp@colData)
      
    }
    sce@colData <- colData
    
    return (sce)
  }

# test function. Calling runDropletUtilsEmptyDrops on unfiltered pbmc3k data.
testRunDropletUtilsEmptyDrops <- function() {
  library(SingleCellExperiment)
  
  sce <- readRDS("data/testUnfilteredSce.rds")
  
  #calling runDropletUtilsEmptyDrops
  sce_test <- runDropletUtilsEmptyDrops(sce = sce, lower = 1000)
  
}

