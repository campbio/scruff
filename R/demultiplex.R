#!/usr/bin/env Rscript

# demultiplex scRNA-seq fastq.gz with cell assignment
# 
# Zhe Wang
# 20170816


# check minimal length
# if len(read1.qual) < umi+bc:
#  sample_counter['unqualified'] +=1
# if min(quals) >= int(min.bc.quality):
#  ### trim read to cut.length
#  if len(read2)>cut.length:
#  read2 = read2[0:cut.length]
# else 
# unqualified + 1
# if read1 barcode !%in% bc.index.file
# undetermined_R1.fastq + read1
# undetermined_R2.fastq + read2
# undetermined + 1
demultiplex <- function(bc.index.file, input.dir, stats.out = "demultiplex_stats",
                        output.dir = "..", out.folder="Demultiplex", min.bc.quality = 10,
                        umi.length = 0, bc.length = 6, cut.length = 50, fname.delimiter = "_",
                        mc.cores) {
  
  input_files <- parse.input.files(input.dir)
  barcode.dt <- fread(bc.index.file, col.names = c("cell_num", "barcode"))
  
  if (nrow(barcode.dt) != length(barcode.dt[, unique(barcode)])) {
    stop(paste("Abort. Barcode file", bc.index.file, "have duplicate cell barcodes"))
  }
  
  meta.dt <- c()
  for (i in seq_len(length(input_files))) {
    meta.dt <- rbindlist(list(meta.dt, parse.fname(input_files[i], fname.delimiter)),
                      use.names=T, fill=F)
  }
  meta.dt <- meta.dt[order(id),]
  
  # delete results from previous run
  unlink(file.path(output.dir, out.folder), recursive = T)
  i <- meta.dt[,unique(id)]
  
  # parallelization
  mclapply(i, demultiplex.sample, meta.dt, barcode.dt, umi.length, bc.length, cut.length,
           stats.out, output.dir, out.folder, min.bc.quality, fname.delimiter,
           mc.cores = mc.cores)
  
  print(paste(Sys.time(), "Job done!"))
}


demultiplex.sample <- function(i, meta.dt, barcode.dt,  umi.length, bc.length, cut.length,
                               stats.out, output.dir, out.folder, min.bc.quality, fname.delimiter) {
  print(paste(Sys.time(), "Processing sample", i))
  sample.meta.dt <- meta.dt[id==i,]
  lanes <- unique(sample.meta.dt[,lane])
  stats.dt <- copy(barcode.dt)
  stats.dt[,cell_fname := paste0(sample.meta.dt[, paste(unique(project), i, sep=fname.delimiter)],
                                 "_sample_", sprintf("%04d", cell_num), ".fastq.gz")]
  stats.dt[,reads := 0]
  stats.dt[,percentage := 0]
  stats.dt <- rbindlist(list(stats.dt, list(cell_num=c(NA, NA, NA),
                                            barcode=c("low_qual", "other", "total"),
                                            cell_fname=c("unqualified", "undetermined", "total"),
                                            reads=c(0, 0, 0),
                                            percentage=c(0, 0, 1))),
                        use.names=T, fill=T, idcol=F)
  
  for (j in lanes) {
    print(paste("Lane", j))
    if (is.na(j)) {
      f1 <- sample.meta.dt[read=="R1", fname]
      f2 <- sample.meta.dt[read=="R2", fname]
    } else {
      f1 <- sample.meta.dt[lane==j & read=="R1", fname]
      f2 <- sample.meta.dt[lane==j & read=="R2", fname]
    }
    fq1 <- FastqStreamer(f1)
    fq2 <- FastqStreamer(f2)
    repeat {
      fqy1 <- yield(fq1)
      fqy2 <- yield(fq2)
      if ((length(fqy1) == 0 & length(fqy2) != 0) | (length(fqy1) != 0 & length(fqy2) == 0))
        stop(paste("Unequal read lengths between read1 and read2 fastq files:", f1, f2))
      else if (length(fqy1) == 0 & length(fqy2) == 0)
        break
      
      stats.dt[cell_fname == "total", reads := reads + length(fqy1)]
      
      min.base.phred1 <- min(as(PhredQuality(substr(fqy1@quality@quality,
                                                    1, umi.length+bc.length)), "IntegerList"))
      
      fqy.dt <- data.table(rname1=tstrsplit(fqy1@id, " ")[[1]],
                           rname2=tstrsplit(fqy2@id, " ")[[1]],
                           umi=substr(fqy1@sread, 1, umi.length),
                           barcode=substr(fqy1@sread, umi.length+1, umi.length+bc.length),
                           qtring1=as.character(fqy1@quality@quality),
                           min.phred1=min.base.phred1,
                           length1=width(fqy1),
                           read1=as.character(fqy1@sread),
                           read2=substr(fqy2@sread, 1, cut.length),
                           qtring2=substr(fqy2@quality@quality, 1, cut.length))
      
      if (!(all(fqy.dt[, rname1] == fqy.dt[, rname2])))
        stop(paste("Abort. Read1 and read2 have different ids in files:", f1, "and", f2))
      
      fqy.dt <- fqy.dt[min.phred1 >= min.bc.quality & length1 >= umi.length+bc.length, ]
      
      stats.dt[cell_fname == "unqualified", reads := reads + length(fqy1)-nrow(fqy.dt)]
      
      for (k in barcode.dt[,cell_num]) {
        bar <- barcode.dt[cell_num == k, barcode]
        sfq.dt <- fqy.dt[barcode == bar,]
        
        # if barcode exists in fastq
        if (nrow(sfq.dt) != 0) {
          fq.out <- ShortReadQ(sread=DNAStringSet(sfq.dt[,read2]),
                               quality=BStringSet(sfq.dt[,qtring2]),
                               id=BStringSet(sfq.dt[,paste0(rname2, ":UMI:", umi, ":")]))
          # project_id_"sample"_cellnum.fastq.gz
          out.fname <- paste0(sample.meta.dt[, paste(unique(project), i, sep=fname.delimiter)],
                              "_sample_", sprintf("%04d", k), ".fastq.gz")
          dir.create(file.path(output.dir, out.folder, i),
                     recursive = T, showWarnings = F)
          out.full <- file.path(output.dir, out.folder, i, out.fname)
          if (file.exists(out.full)) {
            writeFastq(fq.out, out.full, mode = "a")
          }
          else {
            writeFastq(fq.out, out.full, mode = "w")
          }
        }
        stats.dt[barcode == bar, reads := reads + nrow(sfq.dt)]
      }
      
      undetermined.dt <- fqy.dt[!(barcode %in% barcode.dt[, barcode]), ]
      undetermined.fq.out.R1 <- ShortReadQ(sread=DNAStringSet(undetermined.dt[, read1]),
                                           quality=BStringSet(undetermined.dt[, qtring1]),
                                           id=BStringSet(undetermined.dt[, rname1]))
      undetermined.fq.out.R2 <- ShortReadQ(sread=DNAStringSet(undetermined.dt[, read2]),
                                           quality=BStringSet(undetermined.dt[, qtring2]),
                                           id=BStringSet(undetermined.dt[, rname2]))
      out.full.undetermined.R1 <- file.path(output.dir, out.folder, i,
                                            "Undetermined_R1.fastq.gz")
      out.full.undetermined.R2 <- file.path(output.dir, out.folder, i,
                                            "Undetermined_R2.fastq.gz")
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
      stats.dt[cell_fname == "undetermined", reads := reads + nrow(undetermined.dt)]
      print(paste(fq1$`.->.status`[3],"read pairs processed"))
    }
  }
  stats.dt[, percentage := 100*reads/stats.dt[cell_fname == "total", reads]]
  fwrite(stats.dt, file = file.path(output.dir, out.folder,
                                    i, paste(sample.meta.dt[, unique(project)], i
                                             ,stats.out, ".tab", sep="_")), sep="\t")
}


# parse fastq filenames
# extract project name, sample ID, sample number, lane, read
# fastq files have specific naming convention
# project name, sample ID delimiter: "-" or none
# other fields delimiter: "_"
#
# project-ID_number_lane_read_001.fastq.gz
# sample fastq names:
# GD-0802-04_S4_L002_R1_001.fastq.gz
# GD-0802-04_S4_L002_R2_001.fastq.gz
# TEST3_S3_R1_001.fastq.gz
# TEST3_S3_R2_001.fastq.gz
parse.fname <- function(fastq_filename, fname.delimiter) {
  fname <- sub(pattern = "(.*?)\\..*$", replacement = "\\1", basename(fastq_filename))
  fsplit <- strsplit(fname, fname.delimiter)[[1]]
  if (length(fsplit) == 5) {
    pr <- strsplit(fsplit[1], "-")[[1]]
    project <- paste(head(pr, length(pr)-1), collapse="-")
    id <- tail(pr, 1)
    num <- fsplit[2]
    lane <- fsplit[3]
    read <- fsplit[4]
  } else if (length(fsplit) == 4) {
    pr <- strsplit(gsub("([0-9]+)", "~\\1~", fsplit[1]), "~")[[1]]
    project <- paste(head(pr, length(pr)-1), collapse="-")
    id <- tail(pr, 1)
    num <- fsplit[2]
    lane <- NA
    read <- fsplit[3]
  } else {
    stop(paste("fastq filename error:", fastq_filename))
  }
  return (data.table(project=project, id=id, num=num, lane=lane, read=read, fname=fastq_filename))
}


parse.input.files <- function(input.dir) {
  input_files <- list.files(input.dir, full.names = F)
  # filename does not start with Undetermined and ends with fastq or fastq.gz
  input_files <- grep(pattern= "(?=^(?!Undetermined))(?=.*\\.fastq$|.*\\.fastq\\.gz$)",
                      input_files, ignore.case = T, perl = T, value=T)
  input_files <- file.path(input.dir, input_files)
  return (input_files)
}


demultiplex.wrapper <- function(bc.index.file, input.dir, stats.out, output.dir, out.folder,
                                min.bc.quality, umi.length, bc.length, cut.length, fname.delimiter,
                                mc.cores = 16) {
  suppressPackageStartupMessages(library(parallel))
  suppressPackageStartupMessages(library(data.table))
  suppressPackageStartupMessages(library(ShortRead))
  suppressPackageStartupMessages(library(gtools))
  
  demultiplex(bc.index.file, input.dir, stats.out, output.dir, out.folder, min.bc.quality,
              umi.length, bc.length, cut.length, fname.delimiter, mc.cores)
}

