
#' Count the number of UMIs for each transcript and output expression matrix
#' 
#' Count unique \emph{UMI:transcript} pairs for single cell RNA-sequencing alignment files. Write resulting table to output directory. Columns are samples (cells) and rows are transcript IDs.
#' 
#' @param alignment A character vector of the directories to input alignment files.
#' @param features Directory to the gtf reference file. For generation of TxDb objects from gtf files, please refer to \code{makeTxDbFromGFF} function in \code{GenomicFeatures} package.
#' @param format Format of input sequence alignment files. \strong{"BAM"} or \strong{"SAM"}. Default is \strong{"BAM"}.
#' @param out.dir Output directory for UMI counting results. Expression table will be stored in this directory. \strong{Make sure the folder is empty.} Default is \code{"../Count"}.
#' @param cores Number of cores used for parallelization. Default is \code{max(1, parallel::detectCores() - 1)}.
#' @param output.prefix Prefix for expression table filename. Default is \code{"countUMI"}.
#' @param verbose Print log messages. Useful for debugging. Default to \strong{FALSE}.
#' @param logfile.prefix Prefix for log file. Default is current date and time in the format of \code{format(Sys.time(), "\%Y\%m\%d_\%H\%M\%S")}.
#' @return A expression matrix \code{data.table} containing the raw counts of unique \emph{UMI:transcript} pairs.
#' @import data.table foreach
#' @export
count.umi <- function(alignment, features, format = "BAM", out.dir = "../Count", 
                      cores = max(1, parallel::detectCores() - 1),
                      output.prefix = "countUMI", verbose = FALSE,
                      logfile.prefix = format(Sys.time(), "%Y%m%d_%H%M%S")) {
  message(paste(Sys.time(), "Start UMI counting ..."))
  
  logfile <- paste0(logfile.prefix, "_countUMI_log.txt")
  
  log.messages(Sys.time(), "... Start UMI counting", logfile=logfile, append=FALSE)
  log.messages(Sys.time(), alignment, logfile=logfile, append=TRUE)
  
  if (verbose) {
    print("... Input alignment files:")
    print(alignment)
  }
  
  log.messages(Sys.time(), "... Creating output directory", out.dir,
               logfile=logfile, append=TRUE)
  dir.create(file.path(out.dir), showWarnings = FALSE, recursive = T)
  
  log.messages(Sys.time(), paste("... Loading TxDb file", gtf.db.file),
               logfile=logfile, append=TRUE)
  features = gtf.db.read(features)
  
  # parallelization
  cl <- if (verbose) parallel::makeCluster(cores, outfile = logfile) else parallel::makeCluster(cores)
  doParallel::registerDoParallel(cl)
  
  if (format == "SAM") {
    alignment = foreach::foreach(i = alignment, .verbose = verbose, .combine = c,
                     .multicombine=TRUE, .packages = c("Rsubread")) %dopar% {
                       log.messages(Sys.time(), "... Converting", i, "to BAM format (if not exist)",
                                    logfile=logfile, append=TRUE)
    convert.to.bam(i)
                     }
  }
  
  expr = foreach::foreach(i = alignment, .verbose = verbose, .combine = list,
                          .multicombine=TRUE, .packages = c("Rsubread")) %dopar% {
                            log.messages(Sys.time(), "... UMI counting sample", i,
                                         logfile=logfile, append=TRUE)
    count.umi.unit(i, features, format, out.dir, logfile)
  }
  
  parallel::stopCluster(cl)
  
  expr = base::Reduce(bass:merge, expr)
  
  print(paste(Sys.time(), paste("... Write expression table to", 
                                file.path(out.dir, paste0(format(Sys.time(),
                                                                 "%Y%m%d_%H%M%S"), "_",
                                                          output.prefix, ".tab")))))
  
  fwrite(expr, file.path(out.dir, paste0(format(Sys.time(), "%Y%m%d_%H%M%S"), "_",
                                           output.prefix, ".tab")), sep="\t")
  
  message(paste(Sys.time(), "... UMI counting done!"))
  return(expr)
}


count.umi.unit <- function(i, features, format, out.dir, logfile) {
  log.messages(Sys.time(), "... UMI counting sample", i, logfile=logfile, append=TRUE)

  bfl <- Rsamtools::BamFile(i)
  bamGA <- GenomicAlignments::readGAlignments(bfl, use.names=T)
  names(bamGA) <- data.table::last(data.table::tstrsplit(names(bamGA), ":"))
  ol = GenomicAlignments::findOverlaps(features, bamGA)
  ol.dt <- data.table(gene.id=names(features)[queryHits(ol)],
                      umi=names(bamGA)[subjectHits(ol)],
                      pos=start(bamGA)[subjectHits(ol)],
                      hits=subjectHits(ol))
  
  # remove ambiguous gene alignments
  ol.dt <- ol.dt[!(data.table::duplicated(ol.dt, by="hits") |
                     data.table::duplicated(ol.dt, by="hits", fromLast = TRUE)), ]
  count.umi <- base::table(unique(ol.dt[,.(gene.id, umi)])[,gene.id])
  
  count.umi.dt <- data.table::data.table(gene.id=names(features))
  
  count.umi.dt[[basename(i)]] <- 0
  count.umi.dt[gene.id %in% names(count.umi),
               eval(basename(i)) := as.numeric(count.umi[gene.id])]

  #fwrite(count.umi.dt, file.path(output.dir, paste0(sid, ".tab")), sep="\t")
  #print(paste(Sys.time(), sid, "umi counting finished!"))
  return (count.umi.dt)
}


count.wrapper <- function(alignment.dir = out.dir, gtf.file,
                          if.bam = T, output.dir = "../Count", mc.cores = 16) {
  suppressPackageStartupMessages(library(data.table))
  suppressPackageStartupMessages(library(GenomicAlignments))
  suppressPackageStartupMessages(library(GenomicFeatures))
  suppressPackageStartupMessages(library(Rsamtools))
  suppressPackageStartupMessages(library(gtools))
  
  features <- gtf.db.read(gtf.file)
  count.umi(alignment.dir, features, if.bam, output.dir, mc.cores)
}
