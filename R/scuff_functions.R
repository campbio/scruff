
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
# in this order: project-cohort_number_lane_read_001.fastq.gz
# extract project name, cohort, sample number, lane, and read
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
    cohort <- tail(pr, 1)
    num <- fsplit[2]
    lane <- fsplit[3]
    read <- fsplit[4]
  } else {
    stop(paste("fastq filename error:", fastq_filename))
  }
  return (
    data.table::data.table(
      project = project,
      cohort = cohort,
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
    meta.dt <- meta.dt[order(cohort), ]
    return (meta.dt)
  } else {
    stop(
      "Invalid input format for fastq. Need to be of class 'character',
      'data.table' or 'data.frame'."
    )
  }
}


strip.leading.underscore <- function (x)  sub("^\\_+", "", x)

remove.last.extension <- function(x) {
  return (sub(pattern = "\\.[^\\.]*$",
              replacement = "\\1",
              basename(x)))
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


# get gene names from biomart
get.gene.annot <- function(co,
                           host = "www.ensembl.org",
                           biomart = "ENSEMBL_MART_ENSEMBL",
                           dataset = "hsapiens_gene_ensembl",
                           GRCh = NULL) {
  ######### Feature annotation
  gene.id <- co[!(gene.id %in% c("reads_mapped_to_genome",
                                 "reads_mapped_to_genes")) &
                  !grepl("ERCC", co[,gene.id]), gene.id]
  # Use Ensembl version 74 (December 2013, hg19)
  #biomart.host <- "dec2013.archive.ensembl.org"
  
  # Initialize featureData slot of ExpressionSet
  #gene.type <- data.frame(row.names = fnames)
  
  # Get feature annotation data
  ensembl <- biomaRt::useEnsembl(biomart = biomart,
                                 dataset = dataset,
                                 host = host,
                                 GRCh = GRCh)
  biomart.result <- data.table::data.table(biomaRt::getBM(
    attributes=c("ensembl_gene_id", "external_gene_name",
                 "gene_biotype", "chromosome_name"),
    filters="ensembl_gene_id",
    values=gene.id,
    mart=ensembl
  ))
  
  biomart.result <- biomart.result[!base::duplicated(biomart.result,
                                                     by = "ensembl_gene_id"), ]
  
  # Reorder/insert rows so that rows of Biomart query result match 
  # rows of the expression matrix
  biomart.result <- data.table::data.table(biomart.result[
    match(co[!(gene.id %in% c("reads_mapped_to_genome",
                              "reads_mapped_to_genes")) &
               !grepl("ERCC", co[, gene.id]), gene.id],
          biomart.result$ensembl_gene_id),])
  
  return (biomart.result)
}


#' Collect QC metrics
#' 
#' Collect QC metrics from demultiplexing, alignment, and counting results. Return a \code{data.table} object containing reads, transcripts, and genes information.
#' 
#' @param de Demultiplex result. Table returned from \code{demultiplex} function.
#' @param al Alignment result. Table returned from \code{align.rsubread} function.
#' @param co Count matrix. Table returned from \code{count.umi} function.
#' @param biomart.result.dt Gene information table generated from running \code{biomaRt} query on gene IDs.
#' @return QC table
#' @export
get.QC.table <- function(de, al, co, biomart.result.dt = NA) {
  de <- data.table::copy(data.table::data.table(de))
  al <- data.table::copy(data.table::data.table(al))
  al <- data.table::copy(data.table::data.table(al))
  
  colnames(al)[c(1, 3)] <- c("bam_dir",
                             "total_mapped_reads")
  
  de[, cell := sub(pattern = "(.*?)\\..*$",
                 replacement = "\\1", filename)]
  al[, cell := remove.last.extension(basename(bam_dir))]
  
  data.table::setkey(de, cell)
  data.table::setkey(al, cell)
  
  qc.dt <- base::merge(de[,-"filename"],
                       al[, .(cell,
                              total_mapped_reads)],
                       all.x=TRUE)
  
  # get reads mapped to genome
  rmtgenome <- co[gene.id == "reads_mapped_to_genome", -"gene.id"]
  rmtgenome <- data.table::data.table(cell = colnames(rmtgenome),
                                      reads_mapped_to_genome = as.numeric(rmtgenome))
  qc.dt <- base::merge(qc.dt, rmtgenome, all.x=TRUE)
  
  # reads mapped to genes
  rmtgene <- co[gene.id == "reads_mapped_to_genes", -"gene.id"]
  rmtgene <- data.table::data.table(cell = colnames(rmtgene),
                                    reads_mapped_to_genes = as.numeric(rmtgene))
  qc.dt <- base::merge(qc.dt, rmtgene, all.x=TRUE)
  
  # UMI filtered transcripts
  transcript <- base::colSums(co[!(gene.id %in% c("reads_mapped_to_genome",
                                                  "reads_mapped_to_genes")) &
                                   !grepl("ERCC", co[,gene.id]), -"gene.id"])
  transcript <- data.table::data.table(
    cell = names(transcript),
    transcripts = as.numeric(transcript))
  qc.dt <- base::merge(qc.dt, transcript, all.x = TRUE)
  
  # MT transcript
  if (!all(is.na(biomart.result.dt))) {
    mt.transcript <- base::colSums(co[gene.id %in% 
                                        biomart.result.dt[chromosome_name == "MT",
                                                          ensembl_gene_id],
                                      -"gene.id"])
    mt.transcript <- data.table::data.table(
      cell = names(mt.transcript),
      mt_transcripts = as.numeric(mt.transcript))
    qc.dt <- base::merge(qc.dt, mt.transcript, all.x = TRUE)
  }
  
  # expressed genes
  cells <- colnames(co[,-"gene.id"])
  gene <- sapply(cells, function(cells) nrow(
    co[!(gene.id %in%
           c("reads_mapped_to_genome",
             "reads_mapped_to_genes")) &
         !grepl("ERCC", co[,gene.id]) &
         eval(parse(text = paste0("`",
                                  cells,
                                  "`"))) != 0,
       cells, with = FALSE]))
  
  gene <- data.table::data.table(
    cell = names(gene),
    gene = as.numeric(gene))
  qc.dt <- base::merge(qc.dt, gene, all.x = TRUE)
  
  # protein coding genes
  if (!all(is.na(biomart.result.dt))) {
    pro.coding.gene <- biomart.result.dt[gene_biotype == "protein_coding",
                                         ensembl_gene_id]
    pro.gene <- sapply(cells, function(cells) nrow(
      co[gene.id %in% pro.coding.gene &
           eval(parse(text = paste0("`", cells, "`"))) != 0,
         cells, with = FALSE]))
    
    pro.gene <- data.table::data.table(
      cell = names(pro.gene),
      protein_coding_genes = as.numeric(pro.gene))
    qc.dt <- base::merge(qc.dt, pro.gene, all.x = TRUE)
    
    # protein coding transcripts
    pro.transcript <- base::colSums(co[gene.id %in% pro.coding.gene,
                                       cells, with = FALSE])
    pro.transcript <- data.table::data.table(
      cell = names(pro.transcript),
      protein_coding_transcripts = as.numeric(pro.transcript))
    qc.dt <- base::merge(qc.dt, pro.transcript, all.x = TRUE)
  }
  
  return (qc.dt)
}
