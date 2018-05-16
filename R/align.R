#' A wrapper to \code{Rsubread} read alignment
#' 
#' Align cell specific reads to reference genome and write sequence alignment results to output directory. A wrapper to the \code{align} function in \code{Rsubread} package. For details please refer to \code{Rsubread} manual.
#' 
#' @param fastq.paths A vcector containing the paths to input fastq files.
#' @param index Path to the \code{Rsubread} index of reference sequences. For generation of Rsubread indices, please refer to \code{buildindex} function in \code{Rsubread} package.
#' @param unique Argument passed to \code{align} function in \code{Rsubread} package. Logical indicating if only uniquely mapped reads should be reported. A uniquely mapped read has one single mapping location that has less mis-matched bases than any other candidate locations. If set to \strong{FALSE}, multi-mapping reads will be reported in addition to uniquely mapped reads. Number of alignments reported for each multi-mapping read is determined by the nBestLocations parameter. Default is \strong{TRUE}.
#' @param nBestLocations Argument passed to \code{align} function in \code{Rsubread} package. Numeric value specifying the maximal number of equally-best mapping locations that will be reported for a multi-mapping read. 1 by default. The allowed value is between 1 to 16 (inclusive). In the mapping output, ‘NH’ tag is used to indicate how many alignments are reported for the read and ‘HI’ tag is used for numbering the alignments reported for the same read. This argument is only applicable when unique option is set to \strong{FALSE}.
#' @param format Format of sequence alignment results. \strong{"BAM"} or \strong{"SAM"}. Default is \strong{"BAM"}.
#' @param out.dir Output directory for alignment results. Sequence alignment maps will be stored in folders in this directory, respectively. \strong{Make sure the folder is empty.} Default is \code{"../Alignment"}.
#' @param cores Number of cores used for parallelization. Default is \code{max(1, parallel::detectCores() / 2)}.
#' @param threads Number of threads/CPUs used for mapping for each core. Refer to \code{align} function in \code{Rsubread} for details. Default is \strong{1}. It should not be changed in most cases.
#' @param summary.prefix Prefix for alignment summary file. Default is \code{"alignment"}.
#' @param overwrite Overwrite the output directory. Default is \strong{FALSE}.
#' @param verbose Print log messages. Useful for debugging. Default to \strong{FALSE}.
#' @param logfile.prefix Prefix for log file. Default is current date and time in the format of \code{format(Sys.time(), "\%Y\%m\%d_\%H\%M\%S")}.
#' @return A character vector of the paths to output alignment files.
#' @import data.table foreach
#' @export
align.rsubread <- function(fastq.paths,
                           index,
                           unique = TRUE,
                           nBestLocations = 1,
                           format = "BAM",
                           out.dir = "./Alignment",
                           cores = max(1, parallel::detectCores() / 2),
                           threads = 1,
                           summary.prefix = "alignment",
                           overwrite = FALSE,
                           verbose = FALSE,
                           logfile.prefix = format(Sys.time(),
                                                   "%Y%m%d_%H%M%S")) {
  
  message(paste(Sys.time(), "Start alignment ..."))
  
  logfile <- paste0(logfile.prefix, "_alignment_log.txt")
  
  if (verbose) {
    log.messages(Sys.time(),
                 "... Start alignment",
                 logfile = logfile,
                 append = FALSE)
    log.messages(Sys.time(), fastq.paths, logfile = logfile, append = TRUE)
    print("... Input fastq paths:")
    print(fastq.paths)
  }
  
  if (overwrite) {
    # delete results from previous run
    message(paste(Sys.time(), "... Delete (if any) existing alignment results"))
    unlink(file.path(out.dir), recursive = TRUE)
  } else {
    alignment.paths <- get.alignment.file.paths(fastq.paths, format, out.dir)
    if (any(file.exists(alignment.paths))) {
      message(
        paste(
          "Abort.",
          alignment.paths[which(file.exists(alignment.paths) == TRUE)],
          "already exists in output directory",
          file.path(out.dir),
          "\n"
        )
      )
      stop("Abort. Try re-running the function by setting overwrite to TRUE\n")
    }
  }
  
  print(paste(Sys.time(), "... Creating output directory", out.dir))
  dir.create(file.path(out.dir), showWarnings = FALSE, recursive = TRUE)
  
  #sink(logfile, append = TRUE)
  
  # parallelization
  cl <- if (verbose)
    parallel::makeCluster(cores, outfile = logfile)
  else
    parallel::makeCluster(cores)
  doParallel::registerDoParallel(cl)

  alignment.file.paths <- foreach::foreach(
    i = fastq.paths,
    .verbose = verbose,
    .combine = c,
    .multicombine = TRUE,
    .packages = c("Rsubread")
  ) %dopar% {
    if (verbose) {
      align.rsubread.unit(i,
                          index,
                          unique,
                          nBestLocations,
                          format,
                          out.dir,
                          threads,
                          logfile)
    } else {
      suppressMessages(align.rsubread.unit(i,
                                           index,
                                           unique,
                                           nBestLocations,
                                           format,
                                           out.dir,
                                           threads,
                                           logfile = NULL))
    }
  }
  
  res.dt <- foreach::foreach(
    i = alignment.file.paths,
    .verbose = verbose,
    .combine = rbind,
    .multicombine = TRUE,
    .packages = c("Rsubread")
  ) %dopar% {
    if (verbose) {
      propmapped_wrapper(i)
    } else {
      suppressMessages(propmapped_wrapper(i))
    }
  }
  
  parallel::stopCluster(cl)
  
  res.dt <- data.table::data.table(res.dt)
  
  #sink.reset()
  
  print(paste(Sys.time(), paste(
    "... Write alignment summary to",
    file.path(out.dir, paste0(
      format(Sys.time(), "%Y%m%d_%H%M%S"),
      "_",
      summary.prefix,
      ".tab"
    ))
  )))
  
  fwrite(res.dt, file.path(out.dir, paste0(
    format(Sys.time(), "%Y%m%d_%H%M%S"),
    "_",
    summary.prefix,
    ".tab"
  )), sep = "\t")
  
  message(paste(Sys.time(), "... Alignment done!"))
  return(res.dt)
}


align.rsubread.unit <- function(i,
                                index,
                                unique,
                                nBestLocations,
                                format,
                                out.dir,
                                threads,
                                logfile) {
  
  file.path <- get.alignment.file.paths(i, format, out.dir)
  
  if (file.size(i) == 0) {
    file.create(file.path, showWarnings = FALSE)
    return (file.path)
  } else {
    log.messages(Sys.time(),
                 "... mapping sample",
                 i,
                 logfile = logfile,
                 append = TRUE)
    
    Rsubread::align(
      index = index,
      readfile1 = i,
      unique = unique,
      nBestLocations = nBestLocations,
      nthreads = threads,
      output_format = format,
      output_file = file.path
    )
    return (file.path)
  }
}


propmapped_wrapper <- function(i) {
  if (file.size(i) == 0)
    return (data.frame(Samples = i,
                       NumTotal = 0,
                       NumMapped = 0,
                       PropMapped = NA))
  else
    return (Rsubread::propmapped(i))
}
