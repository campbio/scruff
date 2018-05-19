#' Run scruff pipeline
#' 
#' Run the \code{scruff} pipeline. This function performs all \code{demultiplex}, \code{alignRsubread}, and \code{countUMI} functions. Write demultiplex information, alignment information, and expression table in output directories. QC table is also generated.
#' 
#' @param fastqAnnot An annotation data table or data frame that contains information about input fastq files. For example, see \code{?annotationExample}.
#' @param bc A vector of pre-determined cell barcodes. For example, see \code{?barcodeExample}.
#' @param index Path to the \code{Rsubread} index of the reference genome. For generation of Rsubread indices, please refer to \code{buildindex} function in \code{Rsubread} package.
#' @param reference Can be in one of the following 2 formats:
#' \itemize{
#'  \item{"Path to the reference GTF file."}{The TxDb obeject of the GTF file will be generated and saved in the same directory as the GTF file with ".sqlite" suffix.}
#'  \item{"A TxDb object."}{A TxDb object contains feature information about the reference genome. For more detail, please refer to \code{GenomicFeatures} package.}}
#' @param bcStart Integer or vector of integers containing the cell barcode start positions (inclusive, one-based numbering).
#' @param bcStop Integer or vector of integers containing the cell barcode stop positions (inclusive, one-based numbering).
#' @param bcEdit Maximally allowed edit distance for barcode correction. Barcodes with mismatches equal or fewer than this will be assigned a corrected barcode if the inferred barcode matches uniquely in the provided predetermined barcode list.
#' @param umiStart Integer or vector of integers containing the start positions (inclusive, one-based numbering) of UMI sequences.
#' @param umiStop Integer or vector of integers containing the stop positions (inclusive, one-based numbering) of UMI sequences.
#' @param keep Read trimming. Read length or number of nucleotides to keep for read 2 (the read that contains transcript sequence information). Longer reads will be clipped at 3' end. Shorter reads will not be affected.
#' @param minQual Minimally acceptable Phred quality score for barcode and UMI sequences. Phread quality scores are calculated for each nucleotide in the sequence. Sequences with at least one nucleotide with score lower than this will be filtered out. Default is \strong{10}.
#' @param yieldReads The number of reads to yield when drawing successive subsets from a fastq file, providing the number of successive records to be returned on each yield. This parameter is passed to the \code{n} argument of the \code{FastqStreamer} function in \emph{ShortRead} package. Default is \strong{1e06}.
#' @param alignmentFileFormat Format of sequence alignment results. \strong{"BAM"} or \strong{"SAM"}. Default is \strong{"BAM"}.
#' @param demultiplexOutDir Output folder path for demultiplex results. Demultiplexed cell specifc FASTQ files will be stored in folders in this path, respectively. \strong{Make sure the folder is empty.} Default is \code{"./Demultiplex"}.
#' @param alignmentOutDir Output directory for alignment results. Sequence alignment maps will be stored in folders in this directory, respectively. \strong{Make sure the folder is empty.} Default is \code{"./Alignment"}.
#' @param countUmiOutDir Output directory for UMI counting results. UMI corrected count matrix will be stored in this directory. Default is \code{"./Count"}.
#' @param demultiplexSummaryPrefix Prefix for demultiplex summary file. Default is \code{"demultiplex"}.
#' @param alignmentSummaryPrefix Prefix for alignment summary file. Default is \code{"alignment"}.
#' @param exprPrefix Prefix for UMI count matrix filename. Default is \code{"countUMI"}.
#' @param logfilePrefix Prefix for log file. Default is current date and time in the format of \code{format(Sys.time(), "\%Y\%m\%d_\%H\%M\%S")}.
#' @param overwrite Boolean indicating whether to overwrite the output directory. Default is \strong{FALSE}.
#' @param verbose Boolean indicating whether to print log messages. Useful for debugging. Default to \strong{FALSE}.
#' @param cores Number of cores to use for parallelization. Default is \code{max(1, parallel::detectCores() / 2)}.
#' @param threads \strong{Do not change}. Number of threads/CPUs used for mapping for each core. Refer to \code{align} function in \code{Rsubread} for details. Default is \strong{1}. It should not be changed in most cases.
#' @return A \code{data.table} object of UMI filtered count matrix.
#' @export
scruff <- function(fastqAnnot,
                   bc,
                   index,
                   reference,
                   bcStart,
                   bcStop,
                   bcEdit = 0,
                   umiStart,
                   umiStop,
                   keep,
                   minQual = 10,
                   yieldReads = 1e6,
                   alignmentFileFormat = "BAM",
                   demultiplexOutDir = "./Demultiplex",
                   alignmentOutDir = "./Alignment",
                   countUmiOutDir = "./Count",
                   demultiplexSummaryPrefix = "demultiplex",
                   alignmentSummaryPrefix = "alignment",
                   exprPrefix = "countUMI",
                   logfilePrefix = format(Sys.time(), "%Y%m%d_%H%M%S"),
                   overwrite = FALSE,
                   verbose = FALSE,
                   cores = max(1, parallel::detectCores() / 2),
                   threads = 1) {
  
  # run pipeline
  message(paste(Sys.time(), "Start running scruff ..."))
  print(match.call(expand.dots = TRUE))
  
  de <- demultiplex(
    fastqAnnot = fastqAnnot,
    bc = bc,
    bcStart = bcStart,
    bcStop = bcStop,
    bcEdit = bcEdit,
    umiStart = umiStart,
    umiStop = umiStop,
    keep = keep,
    minQual = minQual,
    yieldReads = yieldReads,
    outDir = demultiplexOutDir,
    summaryPrefix = demultiplexSummaryPrefix,
    overwrite = overwrite,
    verbose = verbose,
    cores = cores,
    logfilePrefix = logfilePrefix
  )
  
  al <- alignRsubread(
    fastqPaths = de[!(is.na(cell_num)), fastq_path],
    index = index,
    format = alignmentFileFormat,
    outDir = alignmentOutDir,
    cores = cores,
    threads = threads,
    summaryPrefix = alignmentSummaryPrefix,
    overwrite = overwrite,
    verbose = verbose,
    logfilePrefix = logfilePrefix
  )
  
  co <- countUMI(
    alignment = al$Samples,
    reference = reference,
    format = alignmentFileFormat,
    outDir = countUmiOutDir,
    cores = cores,
    outputPrefix = exprPrefix,
    verbose = verbose,
    logfilePrefix = logfilePrefix
  )
  
  message(paste(Sys.time(), "Finished running scruff ..."))
  
  return (list(demultiplex = de,
               alignment = al,
               expression = co))
}



