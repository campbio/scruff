#' A vector of example cell barcodes.
#'
#' A vector containing 96 predefined cell barcodes which will be used for demultiplexing the example FASTQ files.
#'
#' @format A vector of cell barcode sequences. Cell barcodes for this plate are of length 8 in the example dataset (van den Brink, et al.).
"barcodeExample"

#' Example \code{ShortReadQ} objects of FASTQ files.
#'
#' A list of \code{ShortReadQ} objects containing 4 example FASTQ files from 4 cells o a publicly available dataset (van den Brink, et al.). Each object contains 20,000 reads which are a subset of original dataset. Contains paired end read 1 (R1, UMI + cell barcodes) and read 2 (R2, genomic sequence) reads. Read library was generated using CEL-Seq protocol.
#'
#' @format A list of 4 \code{ShortReadQ} \code{S4} objects each containing 20,000 reads.
"fastqExample"

#' Example annotation table for demultiplexing cell barcodes
#' 
#' A \code{data.table} object of the annotation table for example FASTQ files.
#' 
#' @format A \code{data.table} with 5 columns:
#' \describe{
#' \item{\code{project}}{Project name. Project names must be identical across samples.}
#' \item{\code{sample}}{Sample name}
#' \item{\code{lane}}{Flow cell lane number. If FASTQ files from multiple lanes are concatenated, any placeholder would be sufficient, e.g. "L001".}
#' \item{\code{read1_path}}{Path to read1 FASTQ file. This is the read file with UMI and cell barcode information.}
#' \item{\code{read2_path}}{Path to read2 FASTQ file. This read file contains genomic sequences.}
#' }
"annotationExample"

#' GRCm38 mitochondrial DNA fasta sequence
#'
#' An \code{data.table} containing GRCm38 mitochondrial (MT) DNA sequence.
"GRCm38MitochondrialFasta"

#' GRCm38 mitochondrial gene sets (GTF) file
#' 
#' A \code{data.table} object of GRCm38 mitochondrial (MT) GTF file.
"GRCm38MitochondrialGTF"