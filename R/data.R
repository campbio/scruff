#' A vector of example cell barcodes.
#'
#' A vector containing 96 predefined cell barcodes which will be used for demultiplexing the example FASTQ files.
#'
#' @format A vector of cell barcode sequences. Cell barcodes for this study (van den Brink, et al.) are of length 8.
"barcodeExample"

#' Example SingleCellExperiment Object
#'
#' An example SingleCellExperiment object containing count matrix, cell and gene annotations, and all QC metrics for mouse mitochonrial genes generated from example FASTQ reads.
#'
#' @format A SingleCellExperiment object.
"sceExample"

#' Example GAlignments Object
#'
#' An example GAlignments object containing read alignment information for cell "vandenBrink_b1_cell_0095" of example FASTQ files. Used as an example for \code{rview} function.
#'
#' @format A GAlignments object.
"bamExample"
