log.messages = function(..., sep = " ", logfile = NULL, append = FALSE) {
  if(!is.null(logfile)) {
    if(!is.character(logfile) || length(logfile) > 1) {
      stop("The log file parameter needs to be a single character string.")
    }    
    cat(paste(..., "\n", sep=sep), file=logfile, append=append)
    
  } else {
    message(paste(..., sep=sep))
  }
}


# parse fastq filenames (in Illumina Fastq naming convention)
# in this order: project-ID_number_lane_read_001.fastq.gz
# extract project name, sample ID, sample number, lane, and read
# Example fastq names:
# GD-0802-04_S4_L002_R1_001.fastq.gz
# GD-0802-04_S4_L002_R2_001.fastq.gz
parse.fname <- function(fastq_filename) {
  fname <- sub(pattern = "(.*?)\\..*$", replacement = "\\1", basename(fastq_filename))
  fsplit <- strsplit(fname, "_")[[1]]
  if (length(fsplit) == 5) {
    pr <- strsplit(fsplit[1], "-")[[1]]
    project <- paste(head(pr, length(pr)-1), collapse="-")
    id <- tail(pr, 1)
    num <- fsplit[2]
    lane <- fsplit[3]
    read <- fsplit[4]
  } else {
    stop(paste("fastq filename error:", fastq_filename))
  }
  return (data.table(project=project, id=id, num=num, lane=lane, read=read, dir=fastq_filename))
}


parse.input.files <- function(input.dir) {
  input_files <- list.files(input.dir, full.names = F)
  # filename does not start with Undetermined and ends with fastq or fastq.gz
  input_files <- grep(pattern= "(?=^(?!Undetermined))(?=.*\\.fastq$|.*\\.fastq\\.gz$)",
                      input_files, ignore.case = T, perl = T, value=T)
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
      meta.dt <- rbindlist(list(meta.dt, parse.fname(fname[i])),
                           use.names=T, fill=F)
    }
    meta.dt <- meta.dt[order(id),]
    return (meta.dt)
  }
  else {
    stop("Invalid input format for fastq. Need to be of class 'character',
         'data.table' or 'data.frame'.")
  }
  }
