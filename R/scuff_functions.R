
log.messages <- function(...,
                         sep = " ",
                         logfile = NULL,
                         append = FALSE) {
  if (!is.null(logfile)) {
    if (!is.character(logfile) || length(logfile) > 1) {
      stop("The log file parameter needs to be a single character string.")
    }
    cat(paste(..., "\n", sep = sep),
        file = logfile,
        append = append)
    
  } else {
    message(paste(..., sep = sep))
  }
}


# parse fastq filenames (in Illumina Fastq naming convention)
# in this order: project-ID_number_lane_read_001.fastq.gz
# extract project name, sample ID, sample number, lane, and read
# Example fastq names:
# GD-0802-04_S4_L002_R1_001.fastq.gz
# GD-0802-04_S4_L002_R2_001.fastq.gz
parse.fname <- function(fastq_filename) {
  fname <- sub(pattern = "(.*?)\\..*$",
               replacement = "\\1",
               basename(fastq_filename))
  fsplit <- strsplit(fname, "_")[[1]]
  if (length(fsplit) == 5) {
    pr <- strsplit(fsplit[1], "-")[[1]]
    project <- paste(head(pr, length(pr) - 1), collapse = "-")
    id <- tail(pr, 1)
    num <- fsplit[2]
    lane <- fsplit[3]
    read <- fsplit[4]
  } else {
    stop(paste("fastq filename error:", fastq_filename))
  }
  return (
    data.table::data.table(
      project = project,
      id = id,
      num = num,
      lane = lane,
      read = read,
      dir = fastq_filename
    )
  )
}


parse.input.files <- function(input.dir) {
  input_files <- list.files(input.dir, full.names = F)
  # filename does not start with Undetermined and ends with fastq or fastq.gz
  input_files <- grep(
    pattern = "(?=^(?!Undetermined))(?=.*\\.fastq$|.*\\.fastq\\.gz$)",
    input_files,
    ignore.case = T,
    perl = T,
    value = T
  )
  input_files <- file.path(input.dir, input_files)
  return (input_files)
}


parse.fastq <- function(fastq) {
  if ("data.frame" %in% class(fastq))
    return (data.table::data.table(fastq))
  if ("character" %in% class(fastq)) {
    fname <- parse.input.files(fastq)
    meta.dt <- c()
    for (i in seq_len(length(fname))) {
      meta.dt <- data.table::rbindlist(list(meta.dt, parse.fname(fname[i])),
                                       use.names = T,
                                       fill = F)
    }
    meta.dt <- meta.dt[order(id), ]
    return (meta.dt)
  } else {
    stop(
      "Invalid input format for fastq. Need to be of class 'character',
      'data.table' or 'data.frame'."
    )
  }
}


sink.reset <- function() {
  for (i in seq_len(sink.number())) {
    sink(NULL)
  }
}


get.alignment.file.dir <- function(fastq.dir, format, out.dir) {
  filedir <- file.path(out.dir,
                       paste0(
                         sub(
                           pattern = "(.*?)\\..*$",
                           replacement = "\\1",
                           basename(fastq.dir)
                         ),
                         ".",
                         format
                       ))
  return (filedir)
}


# read gtf database and return feature GRangesList by gene ID
gtf.db.read <- function(gtf.file, logfile) {
  gtf.db.file <- paste0(gtf.file, ".sqlite")
  if ((!(file.exists(gtf.file))) & (!(file.exists(gtf.db.file)))) {
    stop(paste("File", gtf.file, "does not exist"))
  }
  
  if (!(file.exists(gtf.db.file))) {
    message(paste(Sys.time(), "... TxDb file", gtf.db.file, "does not exist"))
    message(paste(Sys.time(), "... Creating TxDb object", gtf.db.file))
    gtf.db <- GenomicFeatures::makeTxDbFromGFF(file = gtf.file)
    AnnotationDbi::saveDb(gtf.db, file = gtf.db.file)
    return (GenomicFeatures::exonsBy(gtf.db, by = "gene"))
  }
  
  gtf.db <- tryCatch(
    AnnotationDbi::loadDb(gtf.db.file),
    error = function(e)
      stop(
        paste(
          "Error loading database file. Delete the file",
          gtf.db.file,
          "and try again."
        )
      )
  )
  return (GenomicFeatures::exonsBy(gtf.db, by = "gene"))
}


to.bam <- function(sam,
                   logfile,
                   overwrite = F,
                   index = T) {
  log.messages(
    Sys.time(),
    "... Converting",
    sam,
    "to BAM format (if not exist)",
    logfile = logfile,
    append = TRUE
  )
  tryCatch(
    Rsamtools::asBam(sam, overwrite = overwrite, indexDestination = index),
    error = function(e) {}
  )
  return (sub(
    pattern = "\\.sam$",
    ignore.case = T,
    perl = T,
    replacement = ".BAM",
    x = sam
  ))
}


# collect QC metrics
get.QC.table <- function(de, al, co) {
  de <- data.table::copy(de)
  al <- data.table::copy(al)
  
  colnames(al)[c(1, 3, 4)] <- c("bam_dir", "mapped_reads", "fraction_mapped")
  
  de[, cell := sub(pattern = "(.*?)\\..*$",
                   replacement = "\\1", filename)]
  al[, cell := sub(pattern = "(.*?)\\..*$",
                   replacement = "\\1", basename(bam_dir))]
  
  
  setkey(de, cell)
  setkey(al, cell)
  
  qc.dt <- merge(de[,-"filename"],
                  al[,.(cell, mapped_reads, fraction_mapped)], all.x=TRUE)
  
  # get reads mapped to genes
  reads.mapped.to.genes <- colSums(count.res[,-1])
  reads.mapped.to.genes.dt <- data.table::data.table(
    reads_mapped_to_genes = reads.mapped.to.genes,
    cell = names(reads.mapped.to.genes))
  setkey(reads.mapped.to.genes.dt, cell)
  setkey(qc.dt, cell)
  qc.dt <- merge(qc.dt, reads.mapped.to.genes.dt, all.x = TRUE)
  
  return (qc.dt)
}
