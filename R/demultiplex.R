#' Demultiplex cell barcodes and assign cell specific reads
#' 
#' Demultiplex fastq files and write cell specific reads in compressed fastq format to output directory
#' 
#' @param fastq.annot An annotation data table or data frame that contains information about input fastq files. For example, see \code{?exampleannot}.
#' @param bc A vector of pre-determined cell barcodes. For example, see \code{?examplebc}.
#' @param bc.start Integer or vector of integers containing the cell barcode start positions (inclusive, one-based numbering).
#' @param bc.stop Integer or vector of integers containing the cell barcode stop positions (inclusive, one-based numbering).
#' @param bc.edit Maximally allowed edit distance for barcode correction. Barcodes with mismatches equal or fewer than this will be assigned a corrected barcode if the inferred barcode matches uniquely in the provided predetermined barcode list.
#' @param umi.start Integer or vector of integers containing the start positions (inclusive, one-based numbering) of UMI sequences.
#' @param umi.stop Integer or vector of integers containing the stop positions (inclusive, one-based numbering) of UMI sequences.
#' @param keep Read trimming. Read length or number of nucleotides to keep for the read that contains transcript sequence information. Longer reads will be clipped at 3' end. Default is \strong{50}.
#' @param min.qual Minimally acceptable Phred quality score for barcode and umi sequences. Phread quality scores are calculated for each nucleotide in the sequence. Sequences with at least one nucleotide with score lower than this will be filtered out. Default is \strong{10}.
#' @param yield.reads The number of reads to yield when drawing successive subsets from a fastq file, providing the number of successive records to be returned on each yield. Default is \strong{1e6}.
#' @param out.dir Output directory for demultiplexing results. Demultiplexed fastq files will be stored in folders in this directory, respectively. \strong{Make sure the folder is empty.} Default is \code{"../Demultiplex"}.
#' @param summary.prefix Prefix for demultiplex summary file. Default is \code{"demultiplex"}.
#' @param overwrite Overwrite the output directory. Default is \strong{FALSE}.
#' @param cores Number of cores used for parallelization. Default is \code{max(1, parallel::detectCores() / 2)}.
#' @param verbose Print log messages. Useful for debugging. Default to \strong{FALSE}.
#' @param logfile.prefix Prefix for log file. Default is current date and time in the format of \code{format(Sys.time(), "\%Y\%m\%d_\%H\%M\%S")}.
#' @return Demultiplexed annotation \code{data.table}.
#' @import data.table foreach
#' @export
demultiplex <- function(fastq.annot,
                        bc,
                        bc.start,
                        bc.stop,
                        bc.edit = 1,
                        umi.start,
                        umi.stop,
                        keep = 50,
                        min.qual = 27,
                        yield.reads = 1e6,
                        out.dir = "./Demultiplex",
                        summary.prefix = "demultiplex",
                        overwrite = FALSE,
                        cores = max(1, parallel::detectCores() / 2),
                        verbose = FALSE,
                        logfile.prefix = format(Sys.time(), "%Y%m%d_%H%M%S")) {
  
  message(paste(Sys.time(), "Start demultiplexing ..."))
  if (overwrite) {
    message(paste(Sys.time(), "All files in", out.dir,  "will be deleted ..."))
  }
  if (verbose) {
    cat("Input fastq type:", class(fastq.annot), "\n")
    print(fastq.annot)
  }
  
  logfile <- paste0(logfile.prefix, "_demultiplex_log.txt")
  
  fastq.annot.dt <- data.table::data.table(fastq.annot)
  barcode.dt <- data.table::data.table("cell_num" = seq_len(length(bc)),
                                       "barcode" = bc)
  sample.id <- fastq.annot.dt[, unique(sample)]
  
  # disable threading in ShortRead package
  nthreads <- .Call(ShortRead:::.set_omp_threads, 1L)
  on.exit(.Call(ShortRead:::.set_omp_threads, nthreads))
  
  # parallelization
  cl <- if (verbose)
    parallel::makeCluster(cores, outfile = logfile)
  else
    parallel::makeCluster(cores)
  doParallel::registerDoParallel(cl)
  
  res.dt <- foreach::foreach(
    i = sample.id,
    .verbose = verbose,
    .combine = rbind,
    .multicombine = TRUE,
    .packages = c("data.table", "ShortRead")
  ) %dopar% {
    if (verbose) {
      ## Generate a unique log file name based on given prefix and parameters
      logfile <- paste0(logfile.prefix, "_sample_", i, "_log.txt")
      demultiplex.unit(
        i,
        fastq.annot.dt,
        barcode.dt,
        bc.start,
        bc.stop,
        bc.edit,
        umi.start,
        umi.stop,
        keep,
        min.qual,
        yield.reads,
        out.dir,
        summary.prefix,
        overwrite,
        logfile
      )
    } else {
      suppressMessages(
        demultiplex.unit(
          i,
          fastq.annot.dt,
          barcode.dt,
          bc.start,
          bc.stop,
          bc.edit,
          umi.start,
          umi.stop,
          keep,
          min.qual,
          yield.reads,
          out.dir,
          summary.prefix,
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
      file.path(out.dir, paste0(
        format(Sys.time(), "%Y%m%d_%H%M%S"),
        "_",
        summary.prefix,
        ".tab"
      ))
    )
  ))
  
  fwrite(res.dt, file = file.path(out.dir, paste0(
    format(Sys.time(), "%Y%m%d_%H%M%S"),
    "_",
    summary.prefix,
    ".tab"
  )), sep = "\t")
  
  message(paste(Sys.time(), "... Demultiplex done!"))
  return(res.dt)
}


# demultiplex function for one sample (unique id)
demultiplex.unit <- function(i,
                             fastq,
                             barcode.dt,
                             bc.start,
                             bc.stop,
                             bc.edit,
                             umi.start,
                             umi.stop,
                             keep,
                             min.qual,
                             yield.reads,
                             out.dir,
                             summary.prefix,
                             overwrite,
                             logfile) {
    
  log.messages(Sys.time(),
               "... demultiplexing sample",
               i,
               logfile = logfile,
               append = FALSE)
  
  sample.meta.dt <- fastq[sample == i, ]
  lanes <- unique(sample.meta.dt[, lane])
  summary.dt <- data.table::copy(barcode.dt)
  summary.dt[, filename := paste0(sample.meta.dt[,
                                                 paste(unique(project),
                                                       i, sep = "_")],
                                  "_cell_",
                                  sprintf("%04d", cell_num),
                                  ".fastq.gz")]
  summary.dt[, reads := 0]
  summary.dt[, percent_assigned := 0]
  
  # initialize summary data table
  summary.dt <- data.table::rbindlist(
    list(
      summary.dt,
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
  summary.dt[, sample := i]
  
  if (overwrite) {
    # delete existing results
    log.messages(
      Sys.time(),
      "... Delete (if any) existing demultiplex results for sample ",
      i,
      logfile = logfile,
      append = TRUE
    )
    unlink(file.path(out.dir, i), recursive = TRUE)
  } else {
    if (any(file.exists(file.path(out.dir, i,
                                  summary.dt[!(is.na(cell_num)), filename])))) {
      log.messages(
        paste(
          "Stop.",
          summary.dt[!(is.na(cell_num)), ]
          [which(file.exists(file.path(out.dir, i,
                                       summary.dt[!(is.na(cell_num)),
                                                  filename])) == TRUE), filename],
          "already exists in output directory",
          file.path(out.dir, i),
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
    
    f1 <- sample.meta.dt[lane == j, read1_path]
    f2 <- sample.meta.dt[lane == j, read2_path]
    fq1 <- ShortRead::FastqStreamer(f1, n = yield.reads)
    fq2 <- ShortRead::FastqStreamer(f2, n = yield.reads)
    repeat {
      fqy1 <- ShortRead::yield(fq1)
      fqy2 <- ShortRead::yield(fq2)
      if (length(fqy1) != length(fqy2))
      {
        log.messages(
          Sys.time(),
          "Stop. Unequal read lengths between read1 and read2 fastq files:",
          f1,
          f2,
          logfile = logfile,
          append = TRUE
        )
        stop(paste("Unequal read lengths between read1 and read2 fastq files:",
                   f1, f2))
      } else if (length(fqy1) == 0 & length(fqy2) == 0)
        break
      
      summary.dt[filename == "total", reads := reads + length(fqy1)]
      
      #min.base.phred1 <- min(methods::as(Biostrings::PhredQuality(),
      #                                   sapply(lapply(fqy1@quality@quality,
      #                                                 substring,
      #                                                 c(bc.start, umi.start),
      #                                                 c(bc.stop, umi.stop)),
      #                                          paste,
      #                                          collapse = ""),
      #                                   "IntegerList"))
      
      
      umi.bc.qual <- ""
      umi.seq <- ""
      bc.seq <- ""
      for (k in seq_len(length(umi.start))) {
        umi.bc.qual <- paste0(umi.bc.qual,
                              substr(fqy1@quality@quality,
                                     umi.start[k],
                                     umi.stop[k]))
        umi.seq <- paste(umi.seq, substr(fqy1@sread,
                                         umi.start[k],
                                         umi.stop[k]), sep = "_")
        umi.seq <- strip.leading.underscore(umi.seq)
      }
      
      for (k in seq_len(length(bc.start))) {
        umi.bc.qual <- paste0(umi.bc.qual, 
                              substr(fqy1@quality@quality,
                                     bc.start[k],
                                     bc.stop[k]))
        bc.seq <- paste(bc.seq, substr(fqy1@sread,
                                       bc.start[k],
                                       bc.stop[k]), sep = "_")
        bc.seq <- strip.leading.underscore(bc.seq)
      }
      
      min.base.phred1 <- min(methods::as(Biostrings::PhredQuality(umi.bc.qual),
                                         "IntegerList"))
      
      
      #min.base.phred1 <- min(methods::as(Biostrings::PhredQuality(paste0(
      #  substr(fqy1@quality@quality, umi.pos[1], umi.pos[2]),
      #  substr(fqy1@quality@quality, bc.pos[1], bc.pos[2])
      #)), "IntegerList"))
      
      fqy.dt <- data.table::data.table(
        rname1 = data.table::tstrsplit(fqy1@id, " ")[[1]],
        rname2 = data.table::tstrsplit(fqy2@id, " ")[[1]],
        read1 = as.character(fqy1@sread),
        read2 = substr(fqy2@sread, 1, keep),
        qtring1 = as.character(fqy1@quality@quality),
        qtring2 = substr(fqy2@quality@quality, 1, keep),
        min.phred1 = min.base.phred1,
        length1 = S4Vectors::width(fqy1),
        #umi = substr(fqy1@sread, umi.pos[1], umi.pos[2]),
        umi = umi.seq,
        #barcode = substr(fqy1@sread, bc.pos[1], bc.pos[2])
        
        # barcodes are separated by "_"
        barcode = bc.seq
      )
      
      # remove low quality and short R1 reads
      
      #fqy.dt <- fqy.dt[min.phred1 >= min.qual & length1 >=
      #                   (max(umi.pos) - min(umi.pos) + 1) +
      #                   (max(bc.pos) - min(bc.pos) + 1)]
      
      fqy.dt <- fqy.dt[min.phred1 >= min.qual &
                         length1 >= sum(bc.stop - bc.start) + length(bc.start) +
                         sum(umi.stop - umi.start) + length(umi.start), ]
      
      if (bc.edit > 0) {
        # cell barcode correction
        fqy.dt[, bc_correct := sapply(barcode,
                                      bc.correct.mem,
                                      barcode.dt[, barcode],
                                      bc.edit)]
      } else {
        fqy.dt[, bc_correct := barcode]
      }
      
      summary.dt[filename == "low_quality",
                 reads := reads + length(fqy1) - nrow(fqy.dt)]
      
      for (k in barcode.dt[, cell_num]) {
        cell.barcode <- barcode.dt[cell_num == k, barcode]
        cfq.dt <- fqy.dt[bc_correct == cell.barcode, ]
        
        dir.create(file.path(out.dir, i),
                   recursive = TRUE,
                   showWarnings = FALSE)
        # project_sample_"cell"_cellnum.fastq.gz
        out.fname <- summary.dt[cell_num == k, filename]
        out.full <- file.path(out.dir, i, out.fname)
        file.create(out.full, showWarnings = FALSE)
        
        # if barcode exists in fastq reads
        if (nrow(cfq.dt) != 0) {
          fq.out <- ShortRead::ShortReadQ(
            sread = Biostrings::DNAStringSet(cfq.dt[, read2]),
            quality = Biostrings::BStringSet(cfq.dt[, qtring2]),
            id = Biostrings::BStringSet(cfq.dt[, paste0(rname2, ":UMI:",
                                                        umi, ":")])
          )
          # write reads to output file
          ShortRead::writeFastq(fq.out, out.full, mode = "a")
        }
        summary.dt[barcode == cell.barcode, reads := reads + nrow(cfq.dt)]
      }
      
      summary.dt[!(is.na(cell_num)), 
                 fastq_path := file.path(out.dir, sample, filename)]
      
      undetermined.dt <- fqy.dt[!(bc_correct %in% barcode.dt[, barcode]), ]
      undetermined.fq.out.R1 <- ShortRead::ShortReadQ(
        sread = Biostrings::DNAStringSet(undetermined.dt[, read1]),
        quality = Biostrings::BStringSet(undetermined.dt[, qtring1]),
        id = Biostrings::BStringSet(undetermined.dt[, rname1])
      )
      undetermined.fq.out.R2 <- ShortReadQ(
        sread = Biostrings::DNAStringSet(undetermined.dt[, read2]),
        quality = Biostrings::BStringSet(undetermined.dt[, qtring2]),
        id = Biostrings::BStringSet(undetermined.dt[, rname2])
      )
      out.full.undetermined.R1 <- file.path(out.dir,
                                            i,
                                            "Undetermined_R1.fastq.gz")
      out.full.undetermined.R2 <- file.path(out.dir,
                                            i,
                                            "Undetermined_R2.fastq.gz")
      file.create(out.full.undetermined.R1, showWarnings = FALSE)
      file.create(out.full.undetermined.R2, showWarnings = FALSE)
      ShortRead::writeFastq(undetermined.fq.out.R1,
                            out.full.undetermined.R1, mode = "a")
      
      ShortRead::writeFastq(undetermined.fq.out.R2,
                            out.full.undetermined.R2, mode = "a")
      
      summary.dt[filename == "undetermined",
                 reads := reads + nrow(undetermined.dt)]
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
  
  summary.dt[, percent_assigned := 100 * reads /
               summary.dt[filename == "total", reads]]
  
  log.messages(
    Sys.time(),
    paste(
      "... Write",
      i,
      "demultiplex summary to",
      file.path(out.dir, i, paste(sample.meta.dt[, unique(project)],
                                  summary.prefix, i, sep =
                                    "_"))
    ),
    logfile = logfile,
    append = TRUE
  )
  
  data.table::fwrite(summary.dt,
                     file = file.path(
                       out.dir,
                       i,
                       paste(sample.meta.dt[, unique(project)],
                             summary.prefix, i, ".tab", sep =
                               "_")
                     ),
                     sep = "\t")
  
  log.messages(
    Sys.time(),
    paste("... finished demultiplexing sample", i),
    logfile = logfile,
    append = TRUE
  )
  
  return(summary.dt)
}

