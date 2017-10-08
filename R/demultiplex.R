#' Demultiplex cell barcodes and assign cell specific reads
#' 
#' Demultiplex fastq files and write cell specific reads in compressed fastq format to output directory.
#' 
#' @param fastq Can be in one of the following formats: \enumerate{
#'   \item An annotation data table or data frame that contains information about input fastq files. For example, please see \code{?exampleannot}.
#'   \item The directory to fastq files. }
#' @param bc A vector of cell barcodes determined from experimental design. For example, please see \code{?examplebc}.
#' @param bc.pos An integer vector of length 2 consisting of the start and end index of barcodes (one-based numbering). Default is \code{c(6, 11)}.
#' @param umi.pos An integer vector of length 2 consisting of the start and end index of umi sequences (one-based numbering). Default is \code{c(1, 5)}.
#' @param keep Read length or number of nucleotides to keep for read that contains transcript sequence information. Longer reads will be clipped at 3' end. Default is \strong{50}.
#' @param bc.qual Minimal Phred quality score acceptable for barcode and umi sequences. Phread quality scores are calculated for each nucleotide in the sequences. Sequences with at least one score lower than this will be filtered out. Default is \strong{10}.
#' @param out.dir Output directory for demultiplexing results. Demultiplexed fastq files will be stored in folders in this directory, respectively. Default is \code{"../Demultiplex"}.
#' @param summary.prefix Prefix for summary files. Default is \code{"demultiplex"}.
#' @param overwrite Whether to overwrite the output directory or not. Default is \strong{FALSE}.
#' @param cores Number of cores used for parallelization. Default is \code{max(1, parallel::detectCores() - 1)}.
#' @param verbose Print log messages. Useful for debugging. Default to \strong{FALSE}.
#' @param logfile.prefix Prefix for log file. Default is current date and time in the format of \code{format(Sys.time(), "\%Y\%m\%d_\%H\%M\%S")}.
#' @import data.table foreach
#' @export
demultiplex <- function(fastq, bc, bc.pos = c(6, 11), umi.pos = c(1, 5), keep = 50,
                        bc.qual = 10, out.dir = "../Demultiplex", summary.prefix = "demultiplex",
                        overwrite = FALSE, cores = max(1, parallel::detectCores() - 1),
                        verbose = FALSE, logfile.prefix = format(Sys.time(), "%Y%m%d_%H%M%S")) {
  message(paste(Sys.time(), "Start demultiplexing ..."))
  
  if (verbose) {
    cat("Input fastq type:", class(fastq), "\n")
    print(fastq)
  }
  
  logfile <- paste0(logfile.prefix, "_demultiplex_log.txt")
  
  fastq.annot.dt <- parse.fastq(fastq)
  barcode.dt <- data.table::data.table("cell_num" = seq_len(length(bc)), "barcode" = bc)
  
  sample.id <- fastq.annot.dt[, unique(id)]
  
  # parallelization
  cl <- if (verbose) parallel::makeCluster(cores, outfile = logfile) else parallel::makeCluster(cores)
  doParallel::registerDoParallel(cl)
  
  res.dt <- foreach::foreach(i = sample.id, .verbose = verbose,
                               .combine = rbind, .multicombine=TRUE,
                               .packages = c("data.table", "ShortRead")) %dopar% {
    if (verbose) {
      ## Generate a unique log file name based on given prefix and parameters
      logfile = paste0(logfile.prefix, "_sample_", i , "_log.txt")
      demultiplex.unit(i, fastq, barcode.dt, bc.pos, umi.pos, keep, bc.qual,
                         out.dir, summary.prefix, overwrite, logfile)
    } else {
      suppressMessages(demultiplex.unit(i, fastq, barcode.dt, bc.pos, umi.pos, keep, bc.qual,
                         out.dir, summary.prefix, overwrite, logfile = NULL))
    }
  }
  parallel::stopCluster(cl)
  
  print(paste(Sys.time(), paste("... Write demultiplex summary for all samples to", 
                                 file.path(out.dir, paste0(format(Sys.time(), "%Y%m%d_%H%M%S"), "_",
                                                          summary.prefix, ".tab")))))
  fwrite(res.dt, file = file.path(out.dir, paste0(format(Sys.time(), "%Y%m%d_%H%M%S"), "_",
                                                  summary.prefix, ".tab")), sep="\t")
  
  message(paste(Sys.time(), "... Demultiplex done!"))
  return(res.dt)
}


# demultiplex function for one sample (unique id)
demultiplex.unit <- function(i, fastq, barcode.dt, bc.pos, umi.pos, keep, bc.qual,
                               out.dir, summary.prefix, overwrite, logfile) {
  log.messages(Sys.time(), "... demultiplexing sample", i, logfile=logfile, append=FALSE)
  sample.meta.dt <- fastq[id==i,]
  lanes <- unique(sample.meta.dt[,lane])
  summary.dt <- copy(barcode.dt)
  summary.dt[,cell_fname := paste0(sample.meta.dt[, paste(unique(project), i, sep="_")],
                                 "_cell_", sprintf("%04d", cell_num), ".fastq.gz")]
  summary.dt[,reads := 0]
  summary.dt[,percentage := 0]
  summary.dt <- rbindlist(list(summary.dt, list(cell_num=c(NA, NA, NA),
                                            barcode=c(NA, NA, NA),
                                            cell_fname=c("low_quality", "undetermined", "total"),
                                            reads=c(0, 0, 0),
                                            percentage=c(0, 0, 1))),
                        use.names = TRUE, fill = TRUE, idcol = FALSE)
  summary.dt[,id := i]
  
  if (overwrite) {
    # delete results from previous run
    log.messages(Sys.time(), "... Delete demultiplex results from previous run for sample ",
                 i, logfile = logfile, append = TRUE)
    unlink(file.path(out.dir, i), recursive = TRUE)
  } else {
    if (any(file.exists(file.path(out.dir, i, summary.dt[!(is.na(cell_num)), cell_fname])))) {
      log.messages(paste("Abort.", summary.dt[!(is.na(cell_num)),]
                         [which(file.exists(
                           file.path(out.dir, i, summary.dt[!(is.na(cell_num)), cell_fname])) == TRUE),
                           cell_fname], "already exists in output directory", file.path(out.dir, i),
                         "\n"), logfile = logfile, append = TRUE)
      stop("Abort. Try setting overwrite to TRUE\n")
    }
  }
  
  for (j in lanes) {
    log.messages(Sys.time(), "... Processing Lane", j, logfile=logfile, append=TRUE)
    if (is.na(j)) {
      f1 <- sample.meta.dt[read=="R1", dir]
      f2 <- sample.meta.dt[read=="R2", dir]
    } else {
      f1 <- sample.meta.dt[lane==j & read=="R1", dir]
      f2 <- sample.meta.dt[lane==j & read=="R2", dir]
    }
    fq1 <- FastqStreamer(f1)
    fq2 <- FastqStreamer(f2)
    repeat {
      fqy1 <- yield(fq1)
      fqy2 <- yield(fq2)
      if (length(fqy1) != length(fqy2))
      {
        log.messages(Sys.time(), "Abort. Unequal read lengths between read1 and read2 fastq files:",
                     f1, f2, logfile=logfile, append=TRUE)
        stop(paste("Unequal read lengths between read1 and read2 fastq files:", f1, f2))
      } else if (length(fqy1) == 0 & length(fqy2) == 0)
        break
      
      summary.dt[cell_fname == "total", reads := reads + length(fqy1)]
      
      min.base.phred1 <- min(as(PhredQuality(
        paste0(substr(fqy1@quality@quality, umi.pos[1], umi.pos[2]),
               substr(fqy1@quality@quality, bc.pos[1], bc.pos[2]))), "IntegerList"))
      
      fqy.dt <- data.table(rname1=tstrsplit(fqy1@id, " ")[[1]],
                           rname2=tstrsplit(fqy2@id, " ")[[1]],
                           read1=as.character(fqy1@sread),
                           read2=substr(fqy2@sread, 1, keep),
                           qtring1=as.character(fqy1@quality@quality),
                           qtring2=substr(fqy2@quality@quality, 1, keep),
                           min.phred1=min.base.phred1,
                           length1=width(fqy1),
                           umi=substr(fqy1@sread, umi.pos[1], umi.pos[2]),
                           barcode=substr(fqy1@sread, bc.pos[1], bc.pos[2]))
      
      if (!(all(fqy.dt[, rname1] == fqy.dt[, rname2])))
      {
        log.messages(Sys.time(), "Abort. Read1 and read2 have different ids in files:",
                     f1, f2, logfile=logfile, append=TRUE)
        stop(paste("Abort. Read1 and read2 have different ids in files:", f1, "and", f2))
      }
      
      fqy.dt <- fqy.dt[min.phred1 >= bc.qual & length1 >= 
                         (max(umi.pos) - min(umi.pos) + 1) + (max(bc.pos) - min(bc.pos) + 1)]
      
      summary.dt[cell_fname == "low_quality", reads := reads + length(fqy1) - nrow(fqy.dt)]
      
      for (k in barcode.dt[,cell_num]) {
        cell.barcode <- barcode.dt[cell_num == k, barcode]
        cfq.dt <- fqy.dt[barcode == cell.barcode,]
        
        # if barcode exists in fastq reads
        if (nrow(cfq.dt) != 0) {
          fq.out <- ShortReadQ(sread=DNAStringSet(cfq.dt[,read2]),
                               quality=BStringSet(cfq.dt[,qtring2]),
                               id=BStringSet(cfq.dt[,paste0(rname2, ":UMI:", umi, ":")]))
          # project_id_"cell"_cellnum.fastq.gz
          out.fname <- summary.dt[cell_num == k, cell_fname]
          #out.fname <- paste0(sample.meta.dt[, paste(unique(project), i, sep="_")],
          #                    "_cell_", sprintf("%04d", k), ".fastq.gz")
          dir.create(file.path(out.dir, i), recursive = TRUE, showWarnings = FALSE)
          out.full <- file.path(out.dir, i, out.fname)
          if (file.exists(out.full)) {
            writeFastq(fq.out, out.full, mode = "a")
          }
          else {
            writeFastq(fq.out, out.full, mode = "w")
          }
        }
        summary.dt[barcode == cell.barcode, reads := reads + nrow(cfq.dt)]
      }
      
      summary.dt[!(is.na(cell_num)), dir := file.path(out.dir, id, cell_fname)]
      
      undetermined.dt <- fqy.dt[!(barcode %in% barcode.dt[, barcode]), ]
      undetermined.fq.out.R1 <- ShortReadQ(sread=DNAStringSet(undetermined.dt[, read1]),
                                           quality=BStringSet(undetermined.dt[, qtring1]),
                                           id=BStringSet(undetermined.dt[, rname1]))
      undetermined.fq.out.R2 <- ShortReadQ(sread=DNAStringSet(undetermined.dt[, read2]),
                                           quality=BStringSet(undetermined.dt[, qtring2]),
                                           id=BStringSet(undetermined.dt[, rname2]))
      out.full.undetermined.R1 <- file.path(out.dir, i, "Undetermined_R1.fastq.gz")
      out.full.undetermined.R2 <- file.path(out.dir, i, "Undetermined_R2.fastq.gz")
      if (file.exists(out.full.undetermined.R1)) {
        writeFastq(undetermined.fq.out.R1, out.full.undetermined.R1, mode = "a")
      }
      else {
        writeFastq(undetermined.fq.out.R1, out.full.undetermined.R1, mode = "w")
      }
      
      if (file.exists(out.full.undetermined.R2)) {
        writeFastq(undetermined.fq.out.R2, out.full.undetermined.R2, mode = "a")
      }
      else {
        writeFastq(undetermined.fq.out.R2, out.full.undetermined.R2, mode = "w")
      }
      summary.dt[cell_fname == "undetermined", reads := reads + nrow(undetermined.dt)]
      log.messages(Sys.time(), paste("...", fq1$`.->.status`[3], "read pairs processed"), 
                   logfile=logfile, append=TRUE)
    }
    close(fq1)
    close(fq2)
  }
  summary.dt[, percentage := 100*reads/summary.dt[cell_fname == "total", reads]]
  log.messages(Sys.time(), paste("... Write", i, "demultiplex summary to", 
                                 file.path(out.dir, i, paste(sample.meta.dt[, unique(project)],
                                                             summary.prefix, i, sep="_"))),
               logfile=logfile, append=TRUE)
  fwrite(summary.dt, file = file.path(out.dir, i, paste(sample.meta.dt[, unique(project)],
                                                        summary.prefix, i, ".tab", sep="_")), sep="\t")
  log.messages(Sys.time(), paste("... finished demultiplexing sample", i), 
               logfile=logfile, append=TRUE)
  return(summary.dt)
}

