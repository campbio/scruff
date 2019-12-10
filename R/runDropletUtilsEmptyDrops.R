#' @title identifies empty droplets using DropletUtils::emptyDrops
#' @description distinguishes between droplets containing cells and ambient RNA in a droplet-based single-cell RNA sequencing experiment.
#' @param sce SingleCellExperiment object. Must contain a raw counts matrix.
#' @param lower Numeric. The lower bound on the total UMI count, at or below which all barcodes are assumed to correspond to empty droplets.
#' @return SingleCellExperiment object with the emptyDrops output stored in coldata: emptyDroplet_Total, emptyDroplet_LogProb, emptyDroplet_PValue, emptyDroplet_Limited, emptyDroplet_FDR.
#' @examples
#' sce <- runDropletUtilsEmptyDrops(sce=sce, lower=100)
#' @export
#' @import DropletUtils
runDropletUtilsEmptyDrops <-
  function(sce, lower = 10) {
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
    
    sce.matrix <- 
      tryCatch(expr = {
        as.matrix(assay(sce))}, 
        warning = function(w) {
          print("Warning: ", w)
        },
        error = function(e) {
          stop("Error: sce assay cannot be converted to matrix. \n\r", e)
        })
    if (!is.numeric(sce.matrix)) {
      stop("Error: sce assay is not numeric")
    }
    rm(sce.matrix)
    
    # run DropletsUtils::emptyDrops on the count matrix, store the output to the sce@CalData, and return the updated sce.
    output <-
      emptyDrops(m = assay(sce),
                 lower = lower)
    colnames(output) <-
      paste0("emptyDroplet_", colnames(output))
    sce@colData <-
      cbind(sce@colData, output)
    
    return (sce)
  }
    
# test function. Calling runDropletUtilsEmptyDrops on unfiltered pbmc3k data.
testRunDropletUtilsEmptyDrops <- function() {
  library(SingleCellExperiment)
  
  sce <- readRDS("data/testUnfilteredSce.rds")
  
  #calling runDropletUtilsEmptyDrops
  sce_test <- runDropletUtilsEmptyDrops(sce=sce, lower = 10)
  
}

#allNamesSamples <- sort(unique(sce@colData$sample))

#for (num in allNamesSamples){
#  sceIdx <- (which(sce@colDatasample == num))
  
  
#}
