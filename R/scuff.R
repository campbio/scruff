#' Run scuff pipeline
#' 
#' Run the \code{scuff} pipeline including \code{demultiplex}, \code{align.rsubread}, and \code{count.umi}. Write expression table in output directory.
#' 
#' @param fastq Can be in one of the following formats: \enumerate{
#'   \item An annotation data table or data frame that contains information about input fastq files. For example, see \code{?exampleannot}.
#'   \item The directory to fastq files. }
#' @param bc A vector of cell barcodes determined from experimental design. For example, see \code{?examplebc}.
#' @param index Directory to the \code{Rsubread} index of reference sequences. For generation of Rsubread indices, please refer to \code{buildindex} function in \code{Rsubread} package.
#' @param features Directory to the gtf reference file. For generation of TxDb objects from gtf files, please refer to \code{makeTxDbFromGFF} function in \code{GenomicFeatures} package.
#' @param umi.pos An integer vector of length 2 consisting of the start and end index of umi sequences (one-based numbering). Default is \code{c(1, 5)}.
#' @param bc.pos An integer vector of length 2 consisting of the start and end index of barcodes (one-based numbering). Default is \code{c(6, 11)}.
#' @param keep Read length or number of nucleotides to keep for read that contains transcript sequence information. Longer reads will be clipped at 3' end. Default is \strong{50}.
#' @param min.qual Minimal acceptable Phred quality score for barcode and umi sequences. Phread quality scores are calculated for each nucleotide in the sequence. Sequences with at least one nucleotide with score lower than this will be filtered out. Default is \strong{10}.
#' @param alignment.file.format Format of sequence alignment results. \strong{"BAM"} or \strong{"SAM"}. Default is \strong{"BAM"}.
#' @param demultiplex.out.dir Output directory for demultiplexing results. Demultiplexed fastq files will be stored in folders in this directory, respectively. \strong{Make sure the folder is empty.} Default is \code{"../Demultiplex"}.
#' @param alignment.out.dir Output directory for alignment results. Sequence alignment maps will be stored in folders in this directory, respectively. \strong{Make sure the folder is empty.} Default is \code{"../Alignment"}.
#' @param umi.out.dir Output directory for UMI counting results. Expression table will be stored in this directory. Default is \code{"../Count"}.
#' @param demultiplex.summary.prefix Prefix for demultiplex summary file. Default is \code{"demultiplex"}.
#' @param alignment.summary.prefix Prefix for alignment summary file. Default is \code{"alignment"}.
#' @param expr.prefix Prefix for expression table filename. Default is \code{"countUMI"}.
#' @param logfile.prefix Prefix for log file. Default is current date and time in the format of \code{format(Sys.time(), "\%Y\%m\%d_\%H\%M\%S")}.
#' @param overwrite Overwrite the output directory. Default is \strong{FALSE}.
#' @param verbose Print log messages. Useful for debugging. Default to \strong{FALSE}.
#' @param cores Number of cores used for parallelization. Default is \code{max(1, parallel::detectCores() - 1)}.
#' @param threads Number of threads/CPUs used for mapping for each core. Refer to \code{align} function in \code{Rsubread} for details. Default is \strong{1}.
#' @export
scuff <- function(fastq,
                  bc,
                  index,
                  features,
                  bc.pos = c(6, 11),
                  umi.pos = c(1, 5),
                  keep = 50,
                  min.qual = 10,
                  alignment.file.format = "BAM",
                  demultiplex.out.dir = "../Demultiplex",
                  alignment.out.dir = "../Alignment",
                  umi.out.dir = "../Count",
                  demultiplex.summary.prefix = "demultiplex",
                  alignment.summary.prefix = "alignment",
                  expr.prefix = "count",
                  logfile.prefix = format(Sys.time(), "%Y%m%d_%H%M%S"),
                  overwrite = FALSE,
                  verbose = FALSE,
                  cores = max(1, parallel::detectCores() - 1),
                  threads = 1) {
  # run pipeline
  message(paste(Sys.time(), "Start running scuff ..."))
  message(match.call(expand.dots = TRUE))
  
  dem <- demultiplex(
    fastq = fastq,
    bc = bc,
    bc.pos = bc.pos,
    umi.pos = umi.pos,
    keep = keep,
    min.qual = min.qual,
    out.dir = demultiplex.out.dir,
    summary.prefix = demultiplex.summary.prefix,
    overwrite = overwrite,
    cores = cores,
    verbose = verbose,
    logfile.prefix = logfile.prefix
  )
  
  fastq.dir <- dem[!(is.na(cell_num)), dir]
  
  ali <- align.rsubread(
    fastq.dir = fastq.dir,
    index = index,
    format = alignment.file.format,
    out.dir = alignment.out.dir,
    cores = cores,
    threads = threads,
    summary.prefix = alignment.summary.prefix,
    overwrite = overwrite,
    verbose = verbose,
    logfile.prefix = logfile.prefix
  )
  
  expr <- count.umi(
    alignment = ali,
    features = features,
    format = alignment.file.format,
    out.dir = umi.out.dir,
    cores = cores,
    output.prefix = expr.prefix,
    verbose = verbose,
    logfile.prefix = logfile.prefix
  )
  
  message(paste(Sys.time(), "Finished running scuff ..."))
  
  return (expr)
}



