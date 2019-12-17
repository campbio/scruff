#' @title Identify empty droplets using \link[DropletUtils]{emptyDrops}.
#' @description Run \link[DropletUtils]{emptyDrops} on the count matrix in the
#'  provided \link[SingleCellExperiment]{SingleCellExperiment} object.
#'  Distinguish between droplets containing cells and ambient RNA in a
#'  droplet-based single-cell RNA sequencing experiment. This function assumes
#'  the provided
#'  \link[SingleCellExperiment]{SingleCellExperiment} object contains cells
#'  from only one sample and any neccesary subsetting of cells happened
#'  beforehand.
#' @param sce A \link[SingleCellExperiment]{SingleCellExperiment} object.
#'  Must contain an unfiltered raw counts matrix.
#' @param lower Numeric. The lower bound on the total UMI count, at or below
#'  which all barcodes are assumed to correspond to empty droplets.
#' @param ... Additional arguments to be passed to
#'  \link[DropletUtils]{emptyDrops}.
#' @return A \link[SingleCellExperiment]{SingleCellExperiment} object with the
#'  \link[DropletUtils]{emptyDrops} output table appended to the
#'  \link[SingleCellExperiment]{colData} slot. The columns include
#'  \emph{EmptyDrops_Total}, \emph{EmptyDrops_LogProb},
#'  \emph{EmptyDrops_PValue}, \emph{EmptyDrops_Limited}, \emph{EmptyDrops_FDR}.
#'  Please refer to the documentation of \link[DropletUtils]{emptyDrops} for
#'  details.
#' @examples
#' # The following unfiltered PBMC_1k_v3 data were downloaded from
#' # https://support.10xgenomics.com/single-cell-gene-expression/datasets/3.0.0
#' # /pbmc_1k_v3
#' # Only the top 10 cells with most counts and the last 10 cells with non-zero
#' # counts are included in this example.
#' # This example only serves as an proof of concept and a tutoriol on how to
#' # run the function. The results should not be
#' # used for drawing scientific conclusions.
#' data(emptyDropsSceExample, package = "scruff")
#' sce <- runDropletUtilsEmptyDrops(sce = emptyDropsSceExample, lower = 100)
#' @import DropletUtils
#' @export
runDropletUtilsEmptyDrops <- function(sce, lower = 10, ...) {
    cts <- tryCatch(expr = {
        as(SingleCellExperiment::counts(sce), "dgTMatrix")},
        warning = function(w) {
            warning("Warning: ", w)
        },
        error = function(e) {
            stop("Error: The count matrix in the provided SCE object cannot be",
                " converted to 'dgTMatrix'.\n", e)
        })

    # run DropletUtils::emptyDrops on the count matrix, store the output to
    # colData(sce), and return the updated sce.
    output <- DropletUtils::emptyDrops(m = cts, lower = lower, ...)
    colnames(output) <- paste0("EmptyDrops_", colnames(output))
    SingleCellExperiment::colData(sce) <-
        cbind(SingleCellExperiment::colData(sce), output)

    return(sce)
}
