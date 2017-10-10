#' Run scuff pipeline
#' 
#' Run the \code{scuff} pipeline including \code{demultiplex}, \code{align.rsubread}, and \code{count.umi}. Write expression table in output directory.
#' 
#' @param fastq An annotation data table or data frame of input fastq files.
#' @param index Rsubread index for reference sequences.
#' @param bc A vector of cell barcodes known from experimental design.
#' @param annot Directory of the gene annotation file.
#' @param clip.at Read length or number of nucleotides to keep. Longer reads will be clipped at 3' end.
#' @param umi.pos Start and end index number of umi sequences.
#' @param bc.pos Start and end index number of barcodes.
#' @param bc.qual Minimal acceptable quality score for barcode sequence.
#' @param out Output directory.
#' @param overwrite Overwrite the output directory. Default is TRUE.
#' @param ncore Number of cores to use.
#' @param threads Number of threads to run for each core. Default is \strong{16}.
#' @export
scuff = function(fastq, bc, index, features, bc.pos = c(6, 11), umi.pos = c(1, 5), keep = 50,
                 bc.qual = 10, alignmentFileFormat = "BAM",
                 exprPrefix = "count",
                 demultiplexOutDir = "../Demultiplex", 
                 alignmentOutDir = "../Alignment",
                 umiCountOutDir = "../Count",
                 demultiplexSummaryPrefix = "demultiplex", 
                 alignmentSummaryPrefix = "alignment", 
                 logfile.prefix = format(Sys.time(), "%Y%m%d_%H%M%S"),
                 overwrite = FALSE, 
                  verbose = FALSE,
                 cores = max(1, parallel::detectCores() - 1), threads = 1) 
{
  # run pipeline
  message("Start running scuff ...")
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
    logfile.prefix = logfile.prefix)
  
}



