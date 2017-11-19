#' Example barcodes.
#'
#' 96 barcodes for demultiplexing example fastq files.
#'
#' @format A vector of barcode sequences.
"examplebc"

#' Example fastq \code{ShortReadQ} objects
#'
#' A list of \code{ShortReadQ} objects of example fastq reads. Contains R1 and R2 reads
#' from Example_01 and R1 and R2 reads from Example_02. Read library is generated
#' using CEL-SEQ2 protocol.
#'
#' @format A list of 4 \code{ShortReadQ} \code{S4} objects.
"examplefastq"

#' Example annotation table
#' 
#' A \code{data.table} object of the annotation table for example data. 
#' 
#' @format A \code{data.table} with 5 columns:
#' \describe{
#' \item{\code{project}}{Project name}
#' \item{\code{sample}}{Sample name}
#' \item{\code{lane}}{Flow cell lane number}
#' \item{\code{read1_path}}{Path to read1 fastq file}
#' \item{\code{read2_path}}{Path to read2 fastq file}
#' }
"exampleannot"

#' GRCh38 mitochondrial dna fasta sequence
#'
#' An vector containing GRCh38 chrMT fasta sequence.
"GRCh38_MT"

#' GRCh38 mitochondrial gtf file
#' 
#' A \code{data.table} object of GRCh38 mitochondrial gtf file.
"GRCh38_MT_gtf"