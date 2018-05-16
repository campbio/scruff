#' Demultiplex cell barcodes and assign cell specific reads
#' 
#' Demultiplex fastq files and write cell specific reads in compressed fastq format to output directory
#' 
#' @param fastqAnnot An annotation data table or data frame that contains information about input fastq files. For example, see \code{?annotationExample}.
#' @param bc A vector of pre-determined cell barcodes. For example, see \code{?barcodeExample}.
#' @param bcStart Integer or vector of integers containing the cell barcode start positions (inclusive, one-based numbering).
#' @param bcStop Integer or vector of integers containing the cell barcode stop positions (inclusive, one-based numbering).
#' @param bcEdit Maximally allowed edit distance for barcode correction. Barcodes with mismatches equal or fewer than this will be assigned a corrected barcode if the inferred barcode matches uniquely in the provided predetermined barcode list.
#' @param umiStart Integer or vector of integers containing the start positions (inclusive, one-based numbering) of UMI sequences.
#' @param umiStop Integer or vector of integers containing the stop positions (inclusive, one-based numbering) of UMI sequences.
#' @param keep Read trimming. Read length or number of nucleotides to keep for read 2 (the read that contains transcript sequence information). Longer reads will be clipped at 3' end. Shorter reads will not be affected.
#' @param minQual Minimally acceptable Phred quality score for barcode and UMI sequences. Phread quality scores are calculated for each nucleotide in the sequence. Sequences with at least one nucleotide with score lower than this will be filtered out. Default is \strong{10}.
#' @param yieldReads The number of reads to yield when drawing successive subsets from a fastq file, providing the number of successive records to be returned on each yield. This parameter is passed to the \code{n} argument of the \code{FastqStreamer} function in \emph{ShortRead} package. Default is \strong{1e06}.
#' @param outDir Output folder path for demultiplex results. Demultiplexed cell specifc FASTQ files will be stored in folders in this path, respectively. \strong{Make sure the folder is empty.} Default is \code{"./Demultiplex"}.
#' @param summaryPrefix Prefix for demultiplex summary file. Default is \code{"demultiplex"}.
#' @param overwrite Boolean indicating whether to overwrite the output directory. Default is \strong{FALSE}.
#' @param cores Number of cores used for parallelization. Default is \code{max(1, parallel::detectCores() / 2)}.
#' @param verbose Poolean indicating whether to print log messages. Useful for debugging. Default to \strong{FALSE}.
#' @param logfilePrefix Prefix for log file. Default is current date and time in the format of \code{format(Sys.time(), "\%Y\%m\%d_\%H\%M\%S")}.
#' @return Demultiplexed annotation \code{data.table}.
#' @import data.table foreach
#' @export
demultiplex <- function(fastqAnnot,
                        bc,
                        bcStart,
                        bcStop,
                        bcEdit = 0,
                        umiStart,
                        umiStop,
                        keep,
                        minQual = 10,
                        yieldReads = 1e6,
                        outDir = "./Demultiplex",
                        summaryPrefix = "demultiplex",
                        overwrite = FALSE,
                        cores = max(1, parallel::detectCores() / 2),
                        verbose = FALSE,
                        logfilePrefix = format(Sys.time(), "%Y%m%d_%H%M%S")) {
  
  message(paste(Sys.time(), "Start demultiplexing ..."))
  if (overwrite) {
    message(paste(Sys.time(), "All files in", outDir,  "will be deleted ..."))
  }
  if (verbose) {
    print("Input annotation table for FASTQ files:")
    print(fastqAnnot)
  }
  
  logfile <- paste0(logfilePrefix, "_demultiplex_log.txt")
  
  fastqAnnotDt <- data.table::data.table(fastqAnnot)
  barcodeDt <- data.table::data.table("cell_num" = seq_len(length(bc)),
                                      "barcode" = bc)
  sampleId <- fastqAnnotDt[, unique(sample)]
  
  # disable threading in ShortRead package
  nthreads <- .Call(ShortRead:::.set_omp_threads, 1L)
  on.exit(.Call(ShortRead:::.set_omp_threads, nthreads))
  
  # parallelization
  cl <- if (verbose)
    parallel::makeCluster(cores, outfile = logfile)
  else
    parallel::makeCluster(cores)
  doParallel::registerDoParallel(cl)
  
  resDt <- foreach::foreach(
    i = sampleId,
    .verbose = verbose,
    .combine = rbind,
    .multicombine = TRUE,
    .packages = c("data.table", "ShortRead")
  ) %dopar% {
    if (verbose) {
      ## Generate a unique log file name based on given prefix and parameters
      logfile <- paste0(logfilePrefix, "_sample_", i, "_log.txt")
      demultiplexUnit(
        i,
        fastqAnnotDt,
        barcodeDt,
        bcStart,
        bcStop,
        bcEdit,
        umiStart,
        umiStop,
        keep,
        minQual,
        yieldReads,
        outDir,
        summaryPrefix,
        overwrite,
        logfile
      )
    } else {
      suppressPackageStartupMessages(
        demultiplexUnit(
          i,
          fastqAnnotDt,
          barcodeDt,
          bcStart,
          bcStop,
          bcEdit,
          umiStart,
          umiStop,
          keep,
          minQual,
          yieldReads,
          outDir,
          summaryPrefix,
          overwrite,
          logfile = NULL
        )
      )
    }
  }
  parallel::stopCluster(cl)
  
  print(paste(
    Sys.time(),
    paste(
      "... Write demultiplex summary for all samples to",
      file.path(outDir, paste0(
        format(Sys.time(), "%Y%m%d_%H%M%S"),
        "_",
        summaryPrefix,
        ".tab"
      ))
    )
  ))
  
  fwrite(resDt, file = file.path(outDir, paste0(
    format(Sys.time(), "%Y%m%d_%H%M%S"),
    "_",
    summaryPrefix,
    ".tab"
  )), sep = "\t")
  
  message(paste(Sys.time(), "... Demultiplex done!"))
  return(resDt)
}


# demultiplex unit function for one sample (unique sample id)
demultiplexUnit <- function(i,
                            fastq,
                            barcodeDt,
                            bcStart,
                            bcStop,
                            bcEdit,
                            umiStart,
                            umiStop,
                            keep,
                            minQual,
                            yieldReads,
                            outDir,
                            summaryPrefix,
                            overwrite,
                            logfile) {
  
  log.messages(Sys.time(),
               "... demultiplexing sample",
               i,
               logfile = logfile,
               append = FALSE)
  
  sampleMetaDt <- fastq[sample == i, ]
  lanes <- unique(sampleMetaDt[, lane])
  summaryDt <- data.table::copy(barcodeDt)
  summaryDt[, filename := paste0(sampleMetaDt[,
                                              paste(unique(project),
                                                    i, sep = "_")],
                                 "_cell_",
                                 sprintf("%04d", cell_num),
                                 ".fastq.gz")]
  summaryDt[, reads := 0]
  summaryDt[, percent_assigned := 0]
  
  # initialize summary data table
  summaryDt <- data.table::rbindlist(
    list(
      summaryDt,
      list(
        cell_num = c(NA, NA, NA),
        barcode = c(NA, NA, NA),
        filename = c("low_quality", "undetermined", "total"),
        reads = c(0, 0, 0),
        percent_assigned = c(0, 0, 1)
      )
    ),
    use.names = TRUE,
    fill = TRUE,
    idcol = FALSE
  )
  summaryDt[, sample := i]
  
  if (overwrite) {
    # delete existing results
    log.messages(
      Sys.time(),
      "... Delete (if any) existing demultiplex results for sample ",
      i,
      logfile = logfile,
      append = TRUE
    )
    unlink(file.path(outDir, i), recursive = TRUE)
  } else {
    if (any(file.exists(file.path(outDir, i,
                                  summaryDt[!(is.na(cell_num)), filename])))) {
      log.messages(
        paste(
          "Stop.",
          summaryDt[!(is.na(cell_num)), ]
          [which(file.exists(file.path(outDir, i,
                                       summaryDt[!(is.na(cell_num)),
                                                 filename])) == TRUE), filename],
          "already exists in output directory",
          file.path(outDir, i),
          "\n"
        ),
        logfile = logfile,
        append = TRUE
      )
      stop("Stop. Try re-running the function by setting overwrite to TRUE\n")
    }
  }
  
  for (j in lanes) {
    log.messages(Sys.time(), "... Processing Lane", j,
                 logfile = logfile, append=TRUE)
    
    f1 <- sampleMetaDt[lane == j, read1_path]
    f2 <- sampleMetaDt[lane == j, read2_path]
    fq1 <- ShortRead::FastqStreamer(f1, n = yieldReads)
    fq2 <- ShortRead::FastqStreamer(f2, n = yieldReads)
    repeat {
      fqy1 <- ShortRead::yield(fq1)
      fqy2 <- ShortRead::yield(fq2)
      if (length(fqy1) != length(fqy2))
      {
        log.messages(
          Sys.time(),
          "Stop. Unequal number of reads between read1 and read2 fastq files:",
          f1,
          f2,
          logfile = logfile,
          append = TRUE
        )
        stop(paste(
          "Unequal number of reads between read1 and read2 fastq files:",
          f1, f2))
      } else if (length(fqy1) == 0 & length(fqy2) == 0)
        break
      
      summaryDt[filename == "total", reads := reads + length(fqy1)]
      
      umiBcQual <- ""
      umiSeq <- ""
      bcSeq <- ""
      for (k in seq_len(length(umiStart))) {
        umiBcQual <- paste0(umiBcQual,
                            substr(fqy1@quality@quality,
                                   umiStart[k],
                                   umiStop[k]))
        umiSeq <- paste(umiSeq, substr(fqy1@sread,
                                       umiStart[k],
                                       umiStop[k]), sep = "_")
        umiSeq <- .stripLeadingUnderscore(umiSeq)
      }
      
      for (k in seq_len(length(bcStart))) {
        umiBcQual <- paste0(umiBcQual, 
                            substr(fqy1@quality@quality,
                                   bcStart[k],
                                   bcStop[k]))
        bcSeq <- paste(bcSeq, substr(fqy1@sread,
                                     bcStart[k],
                                     bcStop[k]), sep = "_")
        bcSeq <- .stripLeadingUnderscore(bcSeq)
      }
      
      minBasePhred1 <- min(methods::as(Biostrings::PhredQuality(umiBcQual),
                                       "IntegerList"))
      
      fqyDt <- data.table::data.table(
        rname1 = data.table::tstrsplit(fqy1@id, " ")[[1]],
        rname2 = data.table::tstrsplit(fqy2@id, " ")[[1]],
        read1 = as.character(fqy1@sread),
        read2 = substr(fqy2@sread, 1, keep),
        qtring1 = as.character(fqy1@quality@quality),
        qtring2 = substr(fqy2@quality@quality, 1, keep),
        min.phred1 = minBasePhred1,
        length1 = S4Vectors::width(fqy1),
        #umi = substr(fqy1@sread, umi.pos[1], umi.pos[2]),
        umi = umiSeq,
        #barcode = substr(fqy1@sread, bc.pos[1], bc.pos[2])
        # barcodes are separated by "_"
        barcode = bcSeq
      )
      
      # remove low quality and short R1 reads
      
      fqyDt <- fqyDt[min.phred1 >= minQual &
                       length1 >= sum(bcStop - bcStart) + length(bcStart) +
                       sum(umiStop - umiStart) + length(umiStart), ]
      
      if (bcEdit > 0) {
        # cell barcode correction
        fqyDt[, bc_correct := vapply(barcode,
                                     .bcCorrectMem,
                                     chracter(1),
                                     barcodeDt[, barcode],
                                     bcEdit)]
      } else {
        fqyDt[, bc_correct := barcode]
      }
      
      summaryDt[filename == "low_quality",
                reads := reads + length(fqy1) - nrow(fqyDt)]
      
      for (k in barcodeDt[, cell_num]) {
        cellBarcode <- barcodeDt[cell_num == k, barcode]
        cfqDt <- fqyDt[bc_correct == cellBarcode, ]
        
        dir.create(file.path(outDir, i),
                   recursive = TRUE,
                   showWarnings = FALSE)
        # project_sample_"cell"_cellnum.fastq.gz
        outFname <- summaryDt[cell_num == k, filename]
        outFull <- file.path(outDir, i, outFname)
        if (!file.exists(outFull)) {
          file.create(outFull, showWarnings = FALSE)
        }
        
        # if barcode exists in fastq reads
        if (nrow(cfqDt) != 0) {
          fqOut <- ShortRead::ShortReadQ(
            sread = Biostrings::DNAStringSet(cfqDt[, read2]),
            quality = Biostrings::BStringSet(cfqDt[, qtring2]),
            id = Biostrings::BStringSet(cfqDt[, paste0(rname2, ":UMI:",
                                                       umi, ":")])
          )
          # write reads to output file
          ShortRead::writeFastq(fqOut, outFull, mode = "a")
        }
        summaryDt[barcode == cellBarcode, reads := reads + nrow(cfqDt)]
      }
      
      summaryDt[!(is.na(cell_num)), 
                fastq_path := file.path(outDir, sample, filename)]
      
      undeterminedDt <- fqyDt[!(bc_correct %in% barcodeDt[, barcode]), ]
      undeterminedFqOutR1 <- ShortRead::ShortReadQ(
        sread = Biostrings::DNAStringSet(undeterminedDt[, read1]),
        quality = Biostrings::BStringSet(undeterminedDt[, qtring1]),
        id = Biostrings::BStringSet(undeterminedDt[, rname1])
      )
      undeterminedFqOutR2 <- ShortRead::ShortReadQ(
        sread = Biostrings::DNAStringSet(undeterminedDt[, read2]),
        quality = Biostrings::BStringSet(undeterminedDt[, qtring2]),
        id = Biostrings::BStringSet(undeterminedDt[, rname2])
      )
      outFullUndeterminedR1 <- file.path(outDir,
                                         i,
                                         "Undetermined_R1.fastq.gz")
      outFullUndeterminedR2 <- file.path(outDir,
                                         i,
                                         "Undetermined_R2.fastq.gz")
      if (!file.exists(outFullUndeterminedR1)) {
        file.create(outFullUndeterminedR1, showWarnings = FALSE)
      }
      if (!file.exists(outFullUndeterminedR2)) {
        file.create(outFullUndeterminedR2, showWarnings = FALSE)
      }
      ShortRead::writeFastq(undeterminedFqOutR1,
                            outFullUndeterminedR1, mode = "a")
      
      ShortRead::writeFastq(undeterminedFqOutR2,
                            outFullUndeterminedR2, mode = "a")
      
      summaryDt[filename == "undetermined",
                reads := reads + nrow(undeterminedDt)]
      log.messages(
        Sys.time(),
        paste("...", fq1$`.->.status`[3],
              "read pairs processed"),
        logfile = logfile,
        append = TRUE
      )
    }
    close(fq1)
    close(fq2)
  }
  
  summaryDt[, percent_assigned := 100 * reads /
              summaryDt[filename == "total", reads]]
  
  log.messages(
    Sys.time(),
    paste(
      "... Write",
      i,
      "demultiplex summary to",
      file.path(outDir, i, paste(sampleMetaDt[, unique(project)],
                                 summaryPrefix, i, sep =
                                   "_"))
    ),
    logfile = logfile,
    append = TRUE
  )
  
  data.table::fwrite(summaryDt,
                     file = file.path(
                       outDir,
                       i,
                       paste(sampleMetaDt[, unique(project)],
                             summaryPrefix, i, ".tab", sep =
                               "_")
                     ),
                     sep = "\t")
  
  log.messages(
    Sys.time(),
    paste("... finished demultiplexing sample", i),
    logfile = logfile,
    append = TRUE
  )
  
  return(summaryDt)
}

