#' A wrapper to \code{Rsubread} read alignment function \code{align}
#' 
#' This function is \strong{not} available in Windows environment. Align cell specific reads to reference genome and write sequence alignment results to output directory. A wrapper to the \code{align} function in \code{Rsubread} package. For details please refer to \code{Rsubread} manual.
#' 
#' @param sce A \code{SingleCellExperiment} object of which the \code{colData} slot contains the \strong{fastq_path} column with paths to input cell-specific FASTQ files.
#' @param index Path to the \code{Rsubread} index of the reference genome. For generation of Rsubread indices, please refer to \code{buildindex} function in \code{Rsubread} package.
#' @param unique Argument passed to \code{align} function in \code{Rsubread} package. Boolean indicating if only uniquely mapped reads should be reported. A uniquely mapped read has one single mapping location that has less mis-matched bases than any other candidate locations. If set to \strong{FALSE}, multi-mapping reads will be reported in addition to uniquely mapped reads. Number of alignments reported for each multi-mapping read is determined by the nBestLocations parameter. Default is \strong{TRUE}.
#' @param nBestLocations Argument passed to \code{align} function in \code{Rsubread} package. Numeric value specifying the maximal number of equally-best mapping locations that will be reported for a multi-mapping read. 1 by default. The allowed value is between 1 to 16 (inclusive). In the mapping output, "NH" tag is used to indicate how many alignments are reported for the read and "HI" tag is used for numbering the alignments reported for the same read. This argument is only applicable when unique option is set to \strong{FALSE}.
#' @param format Format of sequence alignment results. \strong{"BAM"} or \strong{"SAM"}. Default is \strong{"BAM"}.
#' @param outDir Output directory for alignment results. Sequence alignment maps will be stored in folders in this directory, respectively. \strong{Make sure the folder is empty.} Default is \code{"./Alignment"}.
#' @param cores Number of cores used for parallelization. Default is \code{max(1, parallel::detectCores() / 2)}.
#' @param threads \strong{Do not change}. Number of threads/CPUs used for mapping for each core. Refer to \code{align} function in \code{Rsubread} for details. Default is \strong{1}. It should not be changed in most cases.
#' @param summaryPrefix Prefix for alignment summary file. Default is \code{"alignment"}.
#' @param overwrite Boolean indicating whether to overwrite the output directory. Default is \strong{FALSE}.
#' @param verbose Boolean indicating whether to print log messages. Useful for debugging. Default to \strong{FALSE}.
#' @param logfilePrefix Prefix for log file. Default is current date and time in the format of \code{format(Sys.time(), "\%Y\%m\%d_\%H\%M\%S")}.
#' @param ... Additional arguments passed to the \code{align} function in \code{Rsubread} package.
#' @return \strong{SingleCellExperiment} object containing the alignment summary information in the \code{colData} slot. Contains a character vector of the paths to output alignment files.
#' @import data.table foreach
#' @export
alignRsubread <- function(sce,
                          index,
                          unique = TRUE,
                          nBestLocations = 1,
                          format = "BAM",
                          outDir = "./Alignment",
                          cores = max(1, parallel::detectCores() / 2),
                          threads = 1,
                          summaryPrefix = "alignment",
                          overwrite = FALSE,
                          verbose = FALSE,
                          logfilePrefix = format(Sys.time(),
                                                 "%Y%m%d_%H%M%S"),
                          ...) {
  
  if (!requireNamespace("Rsubread", quietly = TRUE)) {
    stop(paste("Package \"Rsubread\" needed for \"alignRsubread\"",
               "function to work.",
               "Please install it if you are using Linux or macOS.",
               "The function is not available in Windows environment."),
         call. = FALSE)
  }
  
  message(paste(Sys.time(), "Start alignment ..."))
  print(match.call(expand.dots = TRUE))
  
  fastqPaths <- SummarizedExperiment::colData(sce)$fastq_path
  
  logfile <- paste0(logfilePrefix, "_alignment_log.txt")
  
  if (verbose) {
    .logMessages(Sys.time(),
                 "... Start alignment",
                 logfile = logfile,
                 append = FALSE)
    .logMessages(Sys.time(), fastqPaths, logfile = logfile, append = TRUE)
    message("... Input fastq paths:")
    print(fastqPaths)
  }
  
  if (overwrite) {
    # delete results from previous run
    message(paste(Sys.time(), "... Delete (if any) existing alignment results"))
    unlink(file.path(outDir), recursive = TRUE)
  } else {
    alignmentPaths <- .getAlignmentFilePaths(fastqPaths, format, outDir)
    if (any(file.exists(alignmentPaths))) {
      message(
        paste(
          "Abort.",
          alignmentPaths[which(file.exists(alignmentPaths) == TRUE)],
          "already exists in output directory",
          file.path(outDir),
          "\n"
        )
      )
      stop("Abort. Try re-running the function by setting overwrite to TRUE\n")
    }
  }
  
  message(paste(Sys.time(), "... Creating output directory", outDir))
  dir.create(file.path(outDir), showWarnings = FALSE, recursive = TRUE)
  
  # parallelization
  cl <- if (verbose)
    parallel::makeCluster(cores, outfile = logfile)
  else
    parallel::makeCluster(cores)
  doParallel::registerDoParallel(cl)
  
  alignmentFilePaths <- foreach::foreach(
    i = fastqPaths,
    .verbose = verbose,
    .combine = c,
    .multicombine = TRUE,
    .packages = c("Rsubread")
  ) %dopar% {
    if (verbose) {
      .alignRsubreadUnit(i,
                         index,
                         unique,
                         nBestLocations,
                         format,
                         outDir,
                         threads,
                         logfile,
                         ...)
    } else {
      suppressPackageStartupMessages(.alignRsubreadUnit(i,
                                                        index,
                                                        unique,
                                                        nBestLocations,
                                                        format,
                                                        outDir,
                                                        threads,
                                                        logfile = NULL,
                                                        ...))
    }
  }
  
  resDt <- foreach::foreach(
    i = alignmentFilePaths,
    .verbose = verbose,
    .combine = rbind,
    .multicombine = TRUE,
    .packages = c("Rsubread")
  ) %dopar% {
    if (verbose) {
      propmappedWrapper(i)
    } else {
      suppressPackageStartupMessages(propmappedWrapper(i))
    }
  }
  
  parallel::stopCluster(cl)
  
  resDt <- data.table::data.table(resDt)
  
  message(paste(Sys.time(), paste(
    "... Write alignment summary to",
    file.path(outDir, paste0(
      format(Sys.time(), "%Y%m%d_%H%M%S"),
      "_",
      summaryPrefix,
      ".tab"
    ))
  )))
  
  fwrite(resDt, file.path(outDir, paste0(
    format(Sys.time(), "%Y%m%d_%H%M%S"),
    "_",
    summaryPrefix,
    ".tab"
  )), sep = "\t")
  
  colnames(resDt) <- c("alignment_path",
                       "reads",
                       "aligned_reads_incl_ercc",
                       "fraction_aligned")
  
  message(paste(Sys.time(), "... Add alignment information to SCE colData."))
  SummarizedExperiment::colData(sce) <-
    cbind(SummarizedExperiment::colData(sce), resDt[, -"reads"])
  
  message(paste(Sys.time(), "... Alignment done!"))
  return(sce)
}


.alignRsubreadUnit <- function(i,
                               index,
                               unique,
                               nBestLocations,
                               format,
                               outDir,
                               threads,
                               logfile,
                               ...) {
  
  file.path <- .getAlignmentFilePaths(i, format, outDir)
  
  if (file.size(i) == 0) {
    file.create(file.path, showWarnings = FALSE)
    return (file.path)
  } else {
    .logMessages(Sys.time(),
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
      output_file = file.path,
      ...
    )
    return (file.path)
  }
}


propmappedWrapper <- function(i) {
  if (file.size(i) == 0)
    return (data.frame(Samples = i,
                       NumTotal = 0,
                       NumMapped = 0,
                       PropMapped = NA))
  else
    return (Rsubread::propmapped(i))
}
