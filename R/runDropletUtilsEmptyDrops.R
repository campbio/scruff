#' @title identifies empty droplets using DropletUtils::emptyDrops
#' @description distinguishes between droplets containing cells and ambient RNA in a droplet-based single-cell RNA sequencing experiment.
#' @param sce SingleCellExperiment object. Must contain a raw counts matrix.
#' @param lower Numeric. The lower bound on the total UMI count, at or below which all barcodes are assumed to correspond to empty droplets.
#' @param assayType Character. The name of the counts matrix within the input
#'  singleCellExperiment object. It could be "count", "raw", "raw counts", "counts" or NULL.
#' @return SingleCellExperiment object with the emptyDrops output stored in coldata: emptyDroplet_Total, emptyDroplet_LogProb, emptyDroplet_PValue, emptyDroplet_Limited, emptyDroplet_FDR.
#' @examples
#' sce <- runDropletUtilsEmptyDrops(sce=sce, lower=100, assayType = "counts")
#' @export
#' @import DropletUtils
runDropletUtilsEmptyDrops <-
  function(sce,
           lower = 100,
           assayType = NULL) {
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
    
    # check if assayType exists in the sce
    if (!is.null(assayType) &&
        !assayType %in% c("count", "raw", "raw counts", "counts")) {
      stop("Wrong assay type. EmptyDrops takes raw counts only.")
    }
    # A helper function. It runs DropletsUtils::emptyDrops on the count matrix, stores the output to the sce@CalData, and returns the updated sce.
    runSingleAssay <- function(sce, lower) {
      output <-
        emptyDrops(m = assay(sce),
                   lower = lower)
      colnames(output) <-
        paste0("emptyDroplet_", colnames(output))
      sce@colData <-
        cbind(sce@colData, output)
      return (sce)
      
    }
    # if the user does not specify assayType, checking if sce assay is a matrix, then call the helper function.
    if (is.null(assayType)) {
      if (is.matrix(assay(sce))) {
        sce_updated <- runSingleAssay(sce, lower)
      } else {
        stop("Assay format is wrong. Assay should be type of matrix.")
      }
    }
    # if the user specifies the assayType, checking if it's in the sce and is a matrix. Calling the helper function on the specified assay.
    else if (assayType %in% names(assays(sce))) {
      if (is.matrix(assays(sce)[[assayType]])) {
        sce_updated <-
          runSingleAssay(assays(sce)[[assayType]], lower)
      } else {
        stop("Assay of assayType should be type of matrix.")
      }
    }  else {
      stop("Assay type not found is assays.")
    }
    return (sce_updated)
  }

# test function. Calling runDropletUtilsEmptyDrops on unfiltered pbmc3k data.
testRunDropletUtilsEmptyDrops() {
  library(SingleCellExperiment)
  
  pbmc3k <- readRDS("data/pbmc3k_unfiltered.rds")
  # testing
  #creating SCE
  sce <-
    SingleCellExperiment(assays = list(counts = pbmc3k, logcounts = pbmc3k))
  #calling runDropletUtilsEmptyDrops
  sce_test <- runDropletUtilsEmptyDrops(sce, lower = 1000)
  
}
