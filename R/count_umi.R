
#' Count the number of UMIs for each transcript and output expression matrix
#' 
#' Count unique \emph{UMI:transcript} pairs for single cell RNA-sequencing alignment files. Write resulting table to output directory. Columns are samples (cells) and rows are transcript IDs.
#' 
#' @param alignment A character vector of the directories to input alignment files.
#' @param features Directory to the gtf reference file. For generation of TxDb objects from gtf files, please refer to \code{makeTxDbFromGFF} function in \code{GenomicFeatures} package.
#' @param format Format of input sequence alignment files. \strong{"BAM"} or \strong{"SAM"}. Default is \strong{"BAM"}.
#' @param out.dir Output directory for UMI counting results. Expression table will be stored in this directory. Default is \code{"../Count"}.
#' @param cores Number of cores used for parallelization. Default is \code{max(1, parallel::detectCores() - 1)}.
#' @param output.prefix Prefix for expression table filename. Default is \code{"countUMI"}.
#' @param verbose Print log messages. Useful for debugging. Default to \strong{FALSE}.
#' @param logfile.prefix Prefix for log file. Default is current date and time in the format of \code{format(Sys.time(), "\%Y\%m\%d_\%H\%M\%S")}.
#' @return A expression matrix \code{data.table} containing the raw counts of unique \emph{UMI:transcript} pairs.
#' @import data.table foreach
#' @export
count.umi <- function(alignment,
                      features,
                      format = "BAM",
                      out.dir = "../Count",
                      cores = max(1, parallel::detectCores() - 1),
                      output.prefix = "countUMI",
                      verbose = FALSE,
                      logfile.prefix = format(Sys.time(), "%Y%m%d_%H%M%S")) {
  
  message(paste(Sys.time(), "Start UMI counting ..."))
  
  logfile <- paste0(logfile.prefix, "_countUMI_log.txt")
  
  if (verbose) {
    log.messages(Sys.time(),
                 "... Start UMI counting",
                 logfile = logfile,
                 append = FALSE)
    log.messages(Sys.time(), alignment, logfile = logfile, append = TRUE)
    print("... Input alignment files:")
    print(alignment)
  } else {
    log.messages(Sys.time(),
                 "... Start UMI counting",
                 logfile = NULL,
                 append = FALSE)
  }
  
  message(paste(Sys.time(),
                "... Creating output directory",
                out.dir))
  dir.create(file.path(out.dir),
             showWarnings = FALSE,
             recursive = TRUE)
  
  print(paste(Sys.time(),
              paste("... Loading TxDb file")))
  features <- suppressPackageStartupMessages(gtf.db.read(features, logfile))
  
  # parallelization
  cl <- if (verbose)
    parallel::makeCluster(cores, outfile = logfile)
  else
    parallel::makeCluster(cores)
  doParallel::registerDoParallel(cl)
  
  if (format == "SAM") {
    alignment <- foreach::foreach(
      i = alignment,
      .verbose = verbose,
      .combine = c,
      .multicombine = TRUE
    ) %dopar% {
      to.bam(i, logfile, overwrite = FALSE, index = TRUE)
    }
  }
  
  expr <- foreach::foreach(
    i = alignment,
    .verbose = verbose,
    .combine = cbind,
    .multicombine = TRUE,
    .packages = c("BiocGenerics", "S4Vectors",
                  "GenomicFeatures", "GenomicAlignments")
  ) %dopar% {
    count.umi.unit(i, features, format, logfile, verbose)
  }
  
  parallel::stopCluster(cl)
  
  expr <- data.table::data.table(expr, keep.rownames = TRUE)
  colnames(expr)[1] <- "gene.id"
  
  print(paste(Sys.time(), paste(
    "... Write expression table to",
    file.path(out.dir, paste0(
      format(Sys.time(),
             "%Y%m%d_%H%M%S"), "_",
      output.prefix, ".tab"
    ))
  )))
  
  data.table::fwrite(expr, file.path(out.dir, paste0(
    format(Sys.time(), "%Y%m%d_%H%M%S"), "_",
    output.prefix, ".tab"
  )), sep = "\t")
  
  message(paste(Sys.time(), "... UMI counting done!"))
  return(expr)
}


count.umi.unit <- function(i, features, format, logfile, verbose) {
  if (verbose) {
    log.messages(Sys.time(),
                 "... UMI counting sample",
                 i,
                 logfile = logfile,
                 append = TRUE)
  }

  bfl <- Rsamtools::BamFile(i)
  bamGA <- GenomicAlignments::readGAlignments(bfl, use.names = T)
  
  # get reads mapped to genome from bamGA
  reads.mapped.to.genome <- length(bamGA[!grepl("ERCC",
                                                GenomeInfoDb::seqnames(bamGA))])
  
  names(bamGA) <- data.table::last(data.table::tstrsplit(names(bamGA), ":"))
  ol <- GenomicAlignments::findOverlaps(features, bamGA)
  ol.dt <- data.table::data.table(
    gene.id = base::names(features)[S4Vectors::queryHits(ol)],
    umi = base::names(bamGA)[S4Vectors::subjectHits(ol)],
    pos = BiocGenerics::start(bamGA)[S4Vectors::subjectHits(ol)],
    hits = S4Vectors::subjectHits(ol)
  )
  
  # remove ambiguous gene alignments
  # reads mapped to genes
  ol.dt <- ol.dt[!(
    base::duplicated(ol.dt, by = "hits") |
      base::duplicated(ol.dt, by = "hits", fromLast = TRUE)
  ), ]
  
  # umi filtering
  count.umi <- base::table(unique(ol.dt[, .(gene.id, umi)])[, gene.id])
  
  # clean up
  count.umi.dt <- data.table::data.table(gene.id = c(names(features),
                                                     "reads_mapped_to_genome",
                                                     "reads_mapped_to_genes"))
  cell <- remove.last.extension(i)
  count.umi.dt[[cell]] <- 0
  count.umi.dt[gene.id == "reads_mapped_to_genome",
               cell] <- reads.mapped.to.genome
  count.umi.dt[gene.id == "reads_mapped_to_genes",
               cell] <- nrow(ol.dt[!grepl("ERCC",
                                          ol.dt[,gene.id]),
                                   ])
  count.umi.dt[gene.id %in% names(count.umi),
               eval(cell) := as.numeric(count.umi[gene.id])]
  count.umi.dt <- data.frame(count.umi.dt,
                             row.names = 1,
                             check.names = FALSE,
                             fix.empty.names = FALSE)
  return (count.umi.dt)
}

