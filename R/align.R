#' A wrapper to \code{Rsubread} read alignment 
#' 
#' Align cell specific reads to reference genome and write sequence alignment results to output directory. A wrapper to the \code{align} function in \code{Rsubread} package. For details please refer to \code{Rsubread} manual.
#' 
#' @param fastq Can be in one of the following formats: \enumerate{
#'   \item A vcector of directories to input files.
#'   \item The main directory to \code{demultiplex} output (sample specific folders containing fastq files). }
#' @param index Directory to the \code{Rsubread} index of reference sequences. For generation of Rsubread indices, please refer to \code{buildindex} function in \code{Rsubread} package.
#' @param format Format of sequence alignment results. \strong{"BAM"} or \strong{"SAM"}. Default is \strong{"BAM"}.
#' @param out.dir Output directory for alignment results. Sequence alignment maps will be stored in folders in this directory, respectively. Default is \code{"../Alignment"}.
#' @param cores Number of cores used for parallelization. Default is \code{max(1, parallel::detectCores() - 1)}.
#' @param threads Number of threads/CPUs used for mapping. Refer to \code{align} function in \code{Rsubread}. Default is \strong{1}.
#' @param summary.prefix Prefix for summary files. Default is \code{"alignment"}.
#' @param overwrite Overwrite the output directory. Default is \strong{FALSE}.
#' @param verbose Print log messages. Useful for debugging. Default to \strong{FALSE}.
#' @param logfile.prefix Prefix for log file. Default is current date and time in the format of \code{format(Sys.time(), "\%Y\%m\%d_\%H\%M\%S")}.
#' @import data.table foreach
#' @export
align.rsubread <- function(fastq, index, format = "BAM", out.dir = "../Alignment",
                           cores = max(1, parallel::detectCores() - 1), threads = 1,
                           summary.prefix = "alignment", overwrite = FALSE, verbose = FALSE, 
                           logfile.prefix = format(Sys.time(), "%Y%m%d_%H%M%S"))

  
  
align.rsubread <- function(demultiplex.dir, index, format = "BAM", out.dir = "../Alignment",
                           threads = 1, mc.cores = 16) {
  i <- list.files(demultiplex.dir, full.names = F)
  mclapply(i, align.sample, demultiplex.dir, index, out.format, out.dir, nthreads,
           mc.cores = mc.cores)
}


align.sample <- function(i, demultiplex.dir, index, out.format, out.dir, nthreads) {
  print(paste(Sys.time(), "Align sample", i))
  
  # delete results from previous run
  unlink(file.path(out.dir, i), recursive = T)
  dir.create(file.path(out.dir, i), showWarnings = FALSE, recursive = T)
  
  samples <- parse.input.files(file.path(demultiplex.dir, i))
  sink("/dev/null")
  align(index = index,
        readfile1 = samples,
        nthreads = nthreads,
        output_format = out.format,
        output_file = file.path(out.dir, i,
                                paste0(sub(pattern = "(.*?)\\..*$",
                                           replacement = "\\1", basename(samples)),
                                       ".", out.format)))
  sink.reset()
}


sink.reset <- function(){
  for(i in seq_len(sink.number())){
    sink(NULL)
  }
}


report <- function(out.dir, mc.cores = 16) {
  i <- list.files(out.dir)
  mclapply(i, report.sample, out.dir, mc.cores = mc.cores)
}


report.sample <- function(i, out.dir) {
  bams <- list.files(path = file.path(out.dir, i), pattern="*.BAM$", full.names = T)
  map_prob <- propmapped(bams)
  write.table(map_prob, file.path(out.dir, i, paste(Sys.Date(),
                                                    i, "alignmentstats.tab", sep="_")), sep="\t")
}


build.index <- function(fa.dir, prefix) {
  out.dir <- dirname(fa.dir)
  bowtie_build(references = fa.dir, outdir = out.dir, force=T, prefix = prefix)
}


align.wrapper <- function(fastq.dir, GRCh38.index, 
                          out.format = "BAM", out.dir = "../Alignment", nthreads = 16, mc.cores = 16) {
  suppressPackageStartupMessages(library(Rsubread))
  align.rsubread(fastq.dir, GRCh38.index, out.format, out.dir, nthreads, mc.cores)
  report(out.dir)
}




