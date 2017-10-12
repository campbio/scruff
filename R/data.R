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
#' A \code{data.table} object of example annotation table. 
#' 
#' @format A \code{data.table} with 6 columns:
#' \describe{
#' \item{\code{project}}{Project name}
#' \item{\code{id}}{Sample ID}
#' \item{\code{num}}{Sample number based on the order that samples are listed during fastq generation}
#' \item{\code{lane}}{Lane number}
#' \item{\code{read}}{The read. R1 means Read 1 and R2 means Read 2}
#' \item{\code{dir}}{Directory to the fastq file}
#' }
"exampleannot"

#' GRCh38 mitochondial dna fasta sequence
#'
#' An \strong{R} object (vector) containing GRCh38 chrMT fasta sequence.
"GRCh38_MT"