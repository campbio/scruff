#' A wrapper to \code{Rsubread} read alignment 
#' 
#' Align cell specific reads to reference genome and write sequence alignment results to output directory. A wrapper to the \code{align} function in \code{Rsubread} package. For details please refer to \code{Rsubread} manual.
#' 
#' @param fastq.dir A vcector containing directories to input fastq files.
#' @param index Directory to the \code{Rsubread} index of reference sequences. For generation of Rsubread indices, please refer to \code{buildindex} function in \code{Rsubread} package.
#' @param format Format of sequence alignment results. \strong{"BAM"} or \strong{"SAM"}. Default is \strong{"BAM"}.
#' @param out.dir Output directory for alignment results. Sequence alignment maps will be stored in folders in this directory, respectively. \strong{Make sure the folder is empty.} Default is \code{"../Alignment"}.
#' @param cores Number of cores used for parallelization. Default is \code{max(1, parallel::detectCores() - 1)}.
#' @param threads Number of threads/CPUs used for mapping for each core. Refer to \code{align} function in \code{Rsubread} for details. Default is \strong{1}.
#' @param summary.prefix Prefix for alignment summary file. Default is \code{"alignment"}.
#' @param overwrite Overwrite the output directory. Default is \strong{FALSE}.
#' @param verbose Print log messages. Useful for debugging. Default to \strong{FALSE}.
#' @param logfile.prefix Prefix for log file. Default is current date and time in the format of \code{format(Sys.time(), "\%Y\%m\%d_\%H\%M\%S")}.
#' @return A character vector of the directories to output alignment files.
#' @import data.table foreach
#' @export
align.rsubread <- function(fastq.dir, index, format = "BAM", out.dir = "../Alignment",
                           cores = max(1, parallel::detectCores() - 1), threads = 1,
                           summary.prefix = "alignment", overwrite = FALSE, verbose = FALSE,
                           logfile.prefix = format(Sys.time(), "%Y%m%d_%H%M%S")) {
  
  message(paste(Sys.time(), "Start alignment ..."))
  
  logfile <- paste0(logfile.prefix, "_alignment_log.txt")
  
  if (verbose) {
    log.messages(Sys.time(), "... Start alignment", logfile=logfile, append=FALSE)
    log.messages(Sys.time(), fastq.dir, logfile=logfile, append=TRUE)
    print("... Input fastq.dir:")
    print(fastq.dir)
  } else {
    log.messages(Sys.time(), "... Start alignment", logfile=NULL, append=FALSE)
  }
  
  if (overwrite) {
    # delete results from previous run
    log.messages(Sys.time(), "... Delete (if any) existing alignment results",
                 logfile = logfile, append = TRUE)
    unlink(file.path(out.dir), recursive = TRUE)
  } else {
    alignment.dir = getalignmentfiledir(fastq.dir, format, out.dir)
    if (any(file.exists(alignment.dir))) {
      log.messages(paste("Abort.", alignment.dir[which(file.exists(alignment.dir) == TRUE)],
                         "already exists in output directory", file.path(out.dir), "\n"),
                   logfile = logfile, append = TRUE)
      stop("Abort. Try re-running the function by setting overwrite to TRUE\n")
    }
  }
  
  log.messages(Sys.time(), "... Creating output directory", out.dir, 
               logfile=logfile, append=TRUE)
  dir.create(file.path(out.dir), showWarnings = FALSE, recursive = T)
  
  sink(logfile, append = TRUE)
  
  # parallelization
  cl <- if (verbose) parallel::makeCluster(cores, outfile = logfile) else parallel::makeCluster(cores)
  doParallel::registerDoParallel(cl)

  alignmentfiledir = foreach::foreach(i = fastq.dir, .verbose = verbose, .combine = c,
                                    .multicombine=TRUE, .packages = c("Rsubread")) %dopar%    {
    if (verbose) {
      align.rsubread.unit(i, index, format, out.dir, threads, logfile)
    } else {
      suppressMessages(align.rsubread.unit(i, index, format,
                                         out.dir, threads, logfile = NULL))
    }
  }
  
  res.dt = foreach::foreach(i = alignmentfiledir, .verbose = verbose, .combine = rbind, 
                            .multicombine=TRUE, .packages = c("Rsubread")) %dopar% {
    if (verbose) {
      Rsubread::propmapped(i)
    } else {
      suppressMessages(Rsubread::propmapped(i))
    }
  }
  
  parallel::stopCluster(cl)
  
  sink.reset()
  
  print(paste(Sys.time(), paste("... Write alignment summary to", 
              file.path(out.dir, paste0(format(Sys.time(), "%Y%m%d_%H%M%S"), "_",
                                                          summary.prefix, ".tab")))))
  
  fwrite(res.dt, file.path(out.dir, paste0(format(Sys.time(), "%Y%m%d_%H%M%S"), "_",
                                           summary.prefix, ".tab")), sep="\t")
  
  message(paste(Sys.time(), "... Alignment done!"))
  return(alignmentfiledir)
}


align.rsubread.unit <- function(i, index, format, out.dir, threads, logfile) {
  log.messages(Sys.time(), "... mapping sample", i, logfile=logfile, append=TRUE)
  filedir = getalignmentfiledir(i, format, out.dir)
  align(index = index, readfile1 = i, nthreads = threads,
        output_format = format, output_file = filedir)
  return (filedir)
}

