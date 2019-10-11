.readBarcodesCellRanger <- function(path,
    header = FALSE,
    colname = "cell_barcode") {

    res <- data.table::fread(path, header = header)
    if (ncol(res) == 1) {
        colnames(res) <- colname
    } else {
        warning("'barcodes' file contains >1 columns!",
            " The column names are kept as is.")
    }
    return(res)
}


.readFeaturesCellRanger <- function(path,
    header = FALSE,
    colnames = c("feature_ID", "feature_name", "feature_type")) {

    res <- data.table::fread(path, header = header)
    if (ncol(res) == 2) {
        colnames(res) <- colnames[1:2]
    } else if (ncol(res) == 3) {
        colnames(res) <- colnames
    } else {
        warning("'features' file contains <2 or >3 columns!",
            " The column names are kept as is.")
    }
    return(res)
}


.readMatrixCellRanger <- function(path, gzipped = TRUE) {
    if (isTRUE(gzipped)) {
        path <- gzfile(path)
    }
    res <- Matrix::readMM(path)
    return(res)
}


# dir <- "outs/filtered_feature_bc_matrix/"
.constructSCEFromCellRangerOutputs <- function(dir,
    sample,
    matrixFileName = "matrix.mtx.gz",
    featuresFileName = "features.tsv.gz",
    barcodesFileName = "barcodes.tsv.gz",
    gzipped = TRUE) {

    cb <- .readBarcodesCellRanger(file.path(dir, barcodesFileName))
    fe <- .readFeaturesCellRanger(file.path(dir, featuresFileName))
    ma <- .readMatrixCellRanger(file.path(dir, matrixFileName),
        gzipped = gzipped)
    expr <- as.matrix(ma)

    coln <- paste(sample, cb[[1]], sep = "_")
    # colnames(expr) <- coln
    rownames(expr) <- fe[[1]]

    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = expr))
    SummarizedExperiment::rowData(sce) <- fe
    SummarizedExperiment::colData(sce) <- S4Vectors::DataFrame(cb,
        column_name = coln,
        sample = sample,
        row.names = coln)

    return(sce)
}


.getCellRangerDir <- function(cellRangerDir, sample, cellRangerOuts) {
    path <- file.path(cellRangerDir, sample, cellRangerOuts)
    return (path)
}


#' @title Construct SCE object from Cell Ranger output
#' @description Read the filtered barcodes, features, and matrices for all
#'  samples from a single run of Cell Ranger output. Combine them into one big
#'  \link[SingleCellExperiment]{SingleCellExperiment} object.
#' @param cellRangerDir The root directory where Cell Ranger was run. This
#'  folder should contain sample specific folders.
#' @param samples A vector of sample names. Must be the same as the folder
#'  names of the samples. The cells in the final SCE object will be ordered by
#'  samples.
#' @param cellRangerOuts The intermidiate path to filtered feature count files
#'  saved in sparce matrix format. Reference genome names might need to be
#'  appended for CellRanger version below 3.0.0 if reads were mapped to
#'  multiple genomes when running Cell Ranger pipeline. Default
#'  \code{"outs/filtered_feature_bc_matrix/"}.
#' @param matrixFileName Filename for the Market Exchange Format (MEX) sparse
#'  matrix file. Default \emph{matrix.mtx.gz}.
#' @param featuresFileName Filename for the feature annotation file. It can be
#'  \emph{features.tsv.gz} or \emph{genes.tsv}. Default
#'  \emph{features.tsv.gz}.
#' @param barcodesFileName Filename for the cell barcode list file. Default
#'  \emph{barcodes.tsv.gz}.
#' @param gzipped Boolean. \code{TRUE} if the Cell Ranger output files
#'  (barcodes.tsv, features.tsv, and matrix.mtx) were
#'  gzip compressed. \code{FALSE} otherwise. This is true after Cell Ranger
#'  3.0.0 update. Default \code{TRUE}.
#' @return A \code{SingleCellExperiment} object containing the combined count
#'  matrix, the feature annotations, and the cell annotation.
#' @examples
#' # Example #1
#' # The following filtered feature, cell, and matrix files were downloaded from
#' # https://support.10xgenomics.com/single-cell-gene-expression/datasets/
#' # 3.0.0/hgmm_1k_v3
#' # The top 10 hg19 & mm10 genes are included in this example.
#' # Only the first 20 cells are included.
#' sce <- importCellRanger(
#'     cellRangerDir = system.file("extdata", package = "scruff"),
#'     samples = "hgmm_1k_v3_20x20")
#' @export
importCellRanger <- function(
    cellRangerDir,
    samples,
    cellRangerOuts = "outs/filtered_feature_bc_matrix/",
    matrixFileName = "matrix.mtx.gz",
    featuresFileName = "features.tsv.gz",
    barcodesFileName = "barcodes.tsv.gz",
    gzipped = TRUE) {

    res <- vector("list", length = length(samples))
    coldata <- vector("list", length = length(samples))

    for (i in seq_along(samples)) {
        dir <- .getCellRangerDir(cellRangerDir, samples[i], cellRangerOuts)
        scei <- .constructSCEFromCellRangerOutputs(dir,
            sample = samples[i],
            matrixFileName = matrixFileName,
            featuresFileName = featuresFileName,
            barcodesFileName = barcodesFileName,
            gzipped = gzipped)
        res[[i]] <- scei
    }

    sce <- do.call(BiocGenerics::cbind, res)
    return(sce)
}
