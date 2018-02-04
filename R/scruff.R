#' Run scruff pipeline
#' 
#' Run the \code{scruff} pipeline including \code{demultiplex}, \code{align.rsubread}, and \code{count.umi} functions. Write demultiplex information, alignment information, and expression table in output directories.
#' 
#' @param fastq.annot An annotation data table or data frame that contains information about input fastq files. For example, see \code{?exampleannot}.
#' @param bc A vector of pre-determined cell barcodes. For example, see \code{?examplebc}.
#' @param index Directory to the \code{Rsubread} index of reference sequences. For generation of Rsubread indices, please refer to \code{buildindex} function in \code{Rsubread} package.
#' @param features Directory to the gtf reference file. For generation of TxDb objects from gtf files, please refer to \code{makeTxDbFromGFF} function in \code{GenomicFeatures} package.
#' @param bc.start Integer or vector of integers containing the cell barcode start positions (inclusive, one-based numbering).
#' @param bc.stop Integer or vector of integers containing the cell barcode stop positions (inclusive, one-based numbering).
#' @param bc.edit Maximally allowed edit distance for barcode correction. Barcodes with mismatches equal or fewer than this will be assigned a corrected barcode if the inferred barcode matches uniquely in the provided predetermined barcode list.
#' @param umi.start Integer or vector of integers containing the start positions (inclusive, one-based numbering) of UMI sequences.
#' @param umi.stop Integer or vector of integers containing the stop positions (inclusive, one-based numbering) of UMI sequences.
#' @param keep Read trimming. Read length or number of nucleotides to keep for the read that contains transcript sequence information. Longer reads will be clipped at 3' end.
#' @param min.qual Minimally acceptable Phred quality score for barcode and umi sequences. Phread quality scores are calculated for each nucleotide in the sequence. Sequences with at least one nucleotide with score lower than this will be filtered out. Default is \strong{10}.
#' @param yield.reads The number of reads to yield when drawing successive subsets from a fastq file, providing the number of successive records to be returned on each yield. Default is \strong{1e6}.
#' @param alignment.file.format Format of sequence alignment results. \strong{"BAM"} or \strong{"SAM"}. Default is \strong{"BAM"}.
#' @param demultiplex.out.dir Output directory for demultiplexing results. Demultiplexed fastq files will be stored in folders in this directory, respectively. \strong{Make sure the folder is empty.} Default is \code{"./Demultiplex"}.
#' @param alignment.out.dir Output directory for alignment results. Sequence alignment maps will be stored in this directory, respectively. \strong{Make sure the folder is empty.} Default is \code{"./Alignment"}.
#' @param umi.out.dir Output directory for UMI counting results. UMI corrected expression table will be stored in this directory. Default is \code{"./Count"}.
#' @param demultiplex.summary.prefix Prefix for demultiplex summary file. Default is \code{"demultiplex"}.
#' @param alignment.summary.prefix Prefix for alignment summary file. Default is \code{"alignment"}.
#' @param expr.prefix Prefix for expression table filename. Default is \code{"countUMI"}.
#' @param logfile.prefix Prefix for log file. Default is current date and time in the format of \code{format(Sys.time(), "\%Y\%m\%d_\%H\%M\%S")}.
#' @param overwrite Overwrite the output directory. Default is \strong{FALSE}.
#' @param verbose Print log messages. Useful for debugging. Default to \strong{FALSE}.
#' @param cores Number of cores to use for parallelization. Default is \code{max(1, parallel::detectCores() / 2)}.
#' @param threads \strong{Do not change}. Number of threads/CPUs used for mapping for each core. Refer to \code{align} function in \code{Rsubread} for details. Default is \strong{1}.
#' @export
scruff <- function(fastq.annot,
                   bc,
                   index,
                   features,
                   bc.start,
                   bc.stop,
                   bc.edit = 1,
                   umi.start,
                   umi.stop,
                   keep,
                   min.qual = 10,
                   yield.reads = 1e6,
                   alignment.file.format = "BAM",
                   demultiplex.out.dir = "./Demultiplex",
                   alignment.out.dir = "./Alignment",
                   umi.out.dir = "./Count",
                   demultiplex.summary.prefix = "demultiplex",
                   alignment.summary.prefix = "alignment",
                   expr.prefix = "countUMI",
                   logfile.prefix = format(Sys.time(), "%Y%m%d_%H%M%S"),
                   overwrite = FALSE,
                   verbose = FALSE,
                   cores = max(1, parallel::detectCores() / 2),
                   threads = 1) {
  
  # run pipeline
  message(paste(Sys.time(), "Start running scruff ..."))
  message(match.call(expand.dots = TRUE))
  
  de <- demultiplex(
    fastq.annot = fastq.annot,
    bc = bc,
    bc.start = bc.start,
    bc.stop = bc.stop,
    bc.edit = bc.edit,
    umi.start = umi.start,
    umi.stop = umi.stop,
    keep = keep,
    min.qual = min.qual,
    yield.reads = yield.reads,
    out.dir = demultiplex.out.dir,
    summary.prefix = demultiplex.summary.prefix,
    overwrite = overwrite,
    verbose = verbose,
    cores = cores,
    logfile.prefix = logfile.prefix
  )
  
  al <- align.rsubread(
    fastq.paths = de[!(is.na(cell_num)), fastq_path],
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
  
  co <- count.umi(
    alignment = al$Samples,
    features = features,
    format = alignment.file.format,
    out.dir = umi.out.dir,
    cores = cores,
    output.prefix = expr.prefix,
    verbose = verbose,
    logfile.prefix = logfile.prefix
  )
  
  message(paste(Sys.time(), "Finished running scruff ..."))
  
  return (list(demultiplex = de,
               alignment = al,
               expression = co))
}



