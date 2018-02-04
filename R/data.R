#' Example barcodes table.
#'
#' A \code{data.table} object containing 96 barcodes for demultiplexing example fastq files.
#'
#' @format A vector of barcode sequences. In the example dataset (van den Brink, et al.), cell barcodes are of length 8.
"examplebc"

#' Example fastq \code{ShortReadQ} objects
#'
#' A list of \code{ShortReadQ} objects of example fastq reads (van den Brink, et al.). Contains read 1 (R1, UMI + cell barcodes) and read 2 (R2, genomic sequence) reads. Read library was generated using CEL-SEQ protocol.
#'
#' @format A list of 4 \code{ShortReadQ} \code{S4} objects each containing 20,000 reads.
"examplefastq"

#' Example annotation table
#' 
#' A \code{data.table} object of the annotation table for example data.
#' 
#' @format A \code{data.table} with 5 columns:
#' \describe{
#' \item{\code{project}}{Project name. It is advised to keep project names identical across samples.}
#' \item{\code{sample}}{Sample name}
#' \item{\code{lane}}{Flow cell lane number. If fastq files from multiple lanes are concatenated after sequencing, any placeholder would be sufficient, e.g. "L001".}
#' \item{\code{read1_path}}{Path to read1 fastq file}
#' \item{\code{read2_path}}{Path to read2 fastq file}
#' }
"exampleannot"

#' GRCm38 mitochondrial dna fasta sequence
#'
#' An \code{data.table} containing GRCm38 mitochondrial (MT) DNA sequence.
"GRCm38_MT_fa"

#' GRCm38 mitochondrial gtf file
#' 
#' A \code{data.table} object of GRCm38 mitochondrial gtf file.
"GRCm38_MT_gtf"