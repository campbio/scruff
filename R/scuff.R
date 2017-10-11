#' Run scuff pipeline
#' 
#' Run the \code{scuff} pipeline including \code{demultiplex}, \code{align.rsubread}, and \code{count.umi}. Write expression table in output directory.
#' 
#' @param fastq An annotation data table or data frame of input fastq files.
#' @param bc A vector of cell barcodes known from experimental design.
#' @param index Rsubread index for reference sequences.
#' @param features 
#' @param umi.pos Start and end index number of umi sequences.
#' @param bc.pos Start and end index number of barcodes.
#' @param keep Read length or number of nucleotides to keep. Longer reads will be clipped at 3' end.
#' @param bc.qual Minimal acceptable quality score for barcode sequence.
#' @param alignmentFileFormat
#' @param demultiplexOutDir 
#' @param alignmentOutDir
#' @param umiCountOutDir
#' @param demultiplexSummaryPrefix
#' @param alignmentSummaryPrefix
#' @param exprPrefix
#' @param logfile.prefix
#' @param overwrite Overwrite the output directory. Default is TRUE.
#' @param verbose
#' @param cores Number of cores to use.
#' @param threads Number of threads to run for each core. Default is \strong{16}.
#' @export
scuff = function(fastq,
                 bc,
                 index,
                 features,
                 bc.pos = c(6, 11),
                 umi.pos = c(1, 5),
                 keep = 50,
                 bc.qual = 10,
                 alignmentFileFormat = "BAM",
                 demultiplexOutDir = "../Demultiplex",
                 alignmentOutDir = "../Alignment",
                 umiCountOutDir = "../Count",
                 demultiplexSummaryPrefix = "demultiplex",
                 alignmentSummaryPrefix = "alignment",
                 exprPrefix = "count",
                 logfile.prefix = format(Sys.time(), "%Y%m%d_%H%M%S"),
                 overwrite = FALSE,
                 verbose = FALSE,
                 cores = max(1, parallel::detectCores() - 1),
                 threads = 1)
{
  # run pipeline
  message(paste(Sys.time(), "Start running scuff ..."))
  message()
  
  dem = demultiplex(
    fastq = fastq,
    bc = bc,
    bc.pos = bc.pos,
    umi.pos = umi.pos,
    keep = keep,
    bc.qual = bc.qual,
    out.dir = demultiplexOutDir,
    summary.prefix = demultiplexSummaryPrefix,
    overwrite = overwrite,
    cores = cores,
    verbose = verbose,
    logfile.prefix = logfile.prefix
  )
  
  fastq.dir = demultiplex.res[!(is.na(cell_num)), dir]
  
  ali = align.rsubread(
    fastq.dir = fastq.dir,
    index = index,
    format = format,
    out.dir = alignmentOutDir,
    cores = cores,
    threads = threads,
    summary.prefix = alignmentSummaryPrefix,
    overwrite = overwrite,
    verbose = verbose,
    logfile.prefix = logfile.prefix
  )
  
  expr = count.umi(
    alignment = ali,
    features = features,
    format = format,
    out.dir = umiCountOutDir,
    cores = cores,
    output.prefix = exprPrefix,
    verbose = verbose,
    logfile.prefix = logfile.prefix
  )
  
  message(paste(Sys.time(), "Finished running scuff ..."))
  
  return (expr)
}



