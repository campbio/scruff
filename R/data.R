#' A vector of example cell barcodes.
#'
#' A vector containing 96 predefined cell barcodes which will be used for
#'  demultiplexing the example FASTQ files.
#'
#' @format A vector of cell barcode sequences. Cell barcodes for this study
#'  (van den Brink, et al.) are of length 8.
"barcodeExample"

#' Example SingleCellExperiment Object
#'
#' An example SingleCellExperiment object containing count matrix, cell and
#'  gene annotations, and all QC metrics for mouse mitochonrial genes generated
#'  from example FASTQ reads.
#'
#' @format A \code{SingleCellExperiment} object.
"sceExample"

#' Example GAlignments Object
#'
#' An example GAlignments object containing read alignment information for cell
#'  "vandenBrink_b1_cell_0095" of example FASTQ files. Used as an example for
#'  \code{rview} function.
#'
#' @format A \code{GAlignments} object.
"bamExample"

#' Cell barcode whitelist (737K-august-2016.txt)
#'
#' A barcode whitelist is the list of all known barcode sequences that have
#'  been included in the assay kit and are available during library preparation.
#'  There are roughly 737,000 cell barcodes in the whitelist
#'  (737K-august-2016.txt) for Cell Ranger's
#'  Single Cell 3' and V(D)J applications.
#'
#' @format A \code{data.table} object.
"validCb"
