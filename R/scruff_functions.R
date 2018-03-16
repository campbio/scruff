########################################
####### parsing helper functions #######
########################################

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


get.alignment.file.paths <- function(fastq.paths, format, out.dir) {
  file.paths <- file.path(out.dir,
                          paste0(
                            sub(
                              pattern = "(.*?)\\..*$",
                              replacement = "\\1",
                              basename(fastq.paths)
                            ),
                            ".",
                            format
                          ))
  return (file.paths)
}


########################################
##### procedural helper functions ######
########################################

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


# correct barcode mismatch using memoization
bc.correct.mem <- local({
  res <- list()
  
  f <- function(bc, ref.barcodes, max.edit.dist) {
    if (bc %in% names(res))
      return (res[[bc]])
    if (bc %in% ref.barcodes) {
      res[[bc]] <<- bc
      return (res[[bc]])
    }
    
    sdm <- stringdist::stringdistmatrix(bc,
                                        ref.barcodes,
                                        method = "hamming",
                                        nthread = 1)
    min.dist <- min(sdm)
    if (min.dist <= max.edit.dist) {
      ind <- which(sdm == min.dist)
      if (length(ind) == 1) {
        res[[bc]] <<- ref.barcodes[ind]
        return (res[[bc]])
      }
    }
    res[[bc]] <<- bc
    return (res[[bc]])
  }
})


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
#' @param biomart.annot.dt Gene information table generated from running \code{biomaRt} query on gene IDs. Must contain \emph{ensembl_gene_id} and \emph{gene_biotype}.
#' @return QC metrics table
#' @export
collectqc <- function(de, al, co, biomart.annot.dt = NA) {
  de <- data.table::copy(data.table::data.table(de))
  al <- data.table::copy(data.table::data.table(al))
  
  colnames(al)[c(1, 3)] <- c("bam_dir",
                             "total_mapped_reads_incl_ercc")
  
  de[, cell := sub(pattern = "(.*?)\\..*$",
                   replacement = "\\1", filename)]
  al[, cell := remove.last.extension(basename(bam_dir))]
  
  data.table::setkey(de, cell)
  data.table::setkey(al, cell)
  
  qc.dt <- base::merge(de[,-"filename"],
                       al[, .(cell,
                              total_mapped_reads_incl_ercc)],
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
  if (!all(is.na(biomart.annot.dt))) {
    mt.transcript <- base::colSums(co[gene.id %in% 
                                        biomart.annot.dt[chromosome_name == "MT",
                                                         ensembl_gene_id],
                                      -"gene.id"])
    mt.transcript <- data.table::data.table(
      cell = names(mt.transcript),
      mt_transcripts = as.numeric(mt.transcript))
    qc.dt <- base::merge(qc.dt, mt.transcript, all.x = TRUE)
  }
  
  # expressed genes
  cells <- colnames(co[,-"gene.id"])
  genes <- sapply(cells, function(cells) nrow(
    co[!(gene.id %in%
           c("reads_mapped_to_genome",
             "reads_mapped_to_genes")) &
         !grepl("ERCC", co[,gene.id]) &
         eval(parse(text = paste0("`",
                                  cells,
                                  "`"))) != 0,
       cells, with = FALSE]))
  
  genes <- data.table::data.table(
    cell = names(genes),
    genes = as.numeric(genes))
  qc.dt <- base::merge(qc.dt, genes, all.x = TRUE)
  
  # protein coding genes
  if (!all(is.na(biomart.annot.dt))) {
    pro.coding.gene <- biomart.annot.dt[gene_biotype == "protein_coding",
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


#' Visualize aligned reads
#' 
#' Visualize mapped reads for single cell RNA-sequencing data. Arrow represents orientation of alignment. Reads are colored by their UMI tag.
#' 
#' @param bamGA A GenomicAlignment object
#' @param chr Chromosome. Integer or "X", "Y", "MT".
#' @param start Genomic coordinate of the start position.
#' @param end Genomic coordinate of the end position.
#' @param legend Show legend. Default is FALSE.
#' @return A ggplot object of aligned reads
#' @import ggbio
#' @export
stepping <- function(bamGA,
                     chr = "1",
                     start = 1,
                     end = max(BiocGenerics::end(bamGA)),
                     legend = FALSE) {
  
  reads <- bamGA[BiocGenerics::start(bamGA) >= start &
                   BiocGenerics::end(bamGA) <= end &
                   GenomeInfoDb::seqnames(bamGA) == chr]
  reads.gr <- sort(GenomicRanges::GRanges(reads))
  
  S4Vectors::mcols(reads.gr)$umi <- data.table::last(data.table::tstrsplit(names(reads.gr), ":"))
  
  #g = ggplot2::ggplot(gr) + ggbio::stat_stepping(xlab = "segment",
  #ylab = "stepping",
  #aes(color = umi, fill = umi))
  g = ggplot2::ggplot(reads.gr) +
    ggbio::geom_arrow(ggplot2::aes(color = umi)) +
    theme_Publication()
  
  if (legend == FALSE) {
    g = g + theme(legend.position="none")
  }
  return (g)
}


#' Visualize gene isoforms
#' 
#' Visualize reference genome. Rectangles represent exons. Arrow represents orientation of transcripts.
#' 
#' @param ensemblGenome A 'ensemblGenome' object derived from running \code{ensemblGenome()} function from \emph{refGenome} package.
#' @param chr Chromosome name. Integer or "X", "Y", "MT".
#' @param start Genomic coordinate of the start position.
#' @param end Genomic coordinate of the end position.
#' @param rect_width Exon widths. Default 0.3.
#' @param line_width Line weight. Default 0.5.
#' @param arrow_segments The number of segments lines be divided to. The larger the number, more arrows there are. Default 10.
#' @param arrow_width The angle of the arrow head in degrees (smaller numbers produce narrower, pointier arrows). Essentially describes the width of the arrow head. Passed to the angle parameter of arrow function. Default 30.
#' @param arrow_length The length of the arrow head. Passed to the length argument of arrow function. Default 0.08.
#' @param arrow_type One of "open" or "closed" indicating whether the arrow head should be a closed triangle. Passed to the type argument of arrow function. Default "open".
#' @param text_size Size of text. Passed to the size argument of the geom_text function. Default 4.
#' @return A ggplot object of genomic view
#' @export
gview <- function(ensemblGenome,
                  chr,
                  start,
                  end,
                  rect_width = 0.3,
                  line_width = 0.5,
                  arrow_segments = 10,
                  arrow_width = 30,
                  arrow_length = 0.08,
                  arrow_type = "open",
                  text_size = 4) {
  
  getlevel <- function(txdt) {
    txdt <- txdt[order(start, end), ]
    
    step <- 1
    while (any(txdt[, set != 1])) {
      x1 <- -1
      for (i in seq_len(nrow(txdt))) {
        if (txdt[i, set == 0 & start > x1]) {
          txdt[i, level := step]
          txdt[i, set := 1]
          x1 <- txdt[i, end]
        }
      }
      step <- step + 1
    }
    return (txdt)
  }
  
  
  gettxdt <- function(dt) {
    transcripts <- dt[, unique(transcript_id)]
    txdt <- data.table::data.table()
    for (i in transcripts) {
      exdt <- dt[transcript_id == i, ]
      txdt <- rbind(txdt,
                    data.table::data.table(
                      start = exdt[, min(start)],
                      end = exdt[, max(end)],
                      transcript_id = i,
                      transcript_name = exdt[, unique(transcript_name)],
                      gene_id = exdt[, unique(gene_id)],
                      gene_name = exdt[, unique(gene_name)],
                      strand = exdt[, unique(strand)]))
      
    }
    txdt[, set := 0]
    
    txdt <- getlevel(txdt)
  }
  
  
  trans_rect <- function(dt, txdt) {
    rdt <- data.table::data.table()
    transcripts <- dt[, unique(transcript_id)]
    for (tx in transcripts) {
      dt[transcript_id == tx,
         level := txdt[transcript_id == tx, level]]
    }
    
    exons <- dt[feature == "exon", ]
    
    for (i in seq_len(nrow(exons))) {
      x1 <- exons[i, start]
      x2 <- exons[i, end]
      y1 <- exons[i, level] - rect_width
      y2 <- exons[i, level] + rect_width
      
      rdt <- rbind(rdt,
                   data.table::data.table(x1 = x1,
                                          x2 = x2,
                                          y1 = y1,
                                          y2 = y2,
                                          exon_number = exons[i, exon_number],
                                          transcript_id = exons[i, transcript_id],
                                          gene_id = exons[i, gene_id],
                                          gene_name = exons[i, gene_name]),
                   fill = TRUE)
    }
    return (rdt)
  }
  
  
  trans_arrow <- function(dt) {
    adt <- data.table::data.table()
    for (i in seq_len(nrow(dt))) {
      mi <- dt[i, start]
      ma <- dt[i, end]
      
      if (dt[i, strand] == "+") {
        x1 <- mi + (((ma - mi)/arrow_segments) * seq(0, arrow_segments - 1))
        x2 <- ma - (((ma - mi)/arrow_segments) * seq(arrow_segments - 1, 0))
      } else if (dt[i, strand] == "-") {
        x1 <- ma + (((mi - ma)/arrow_segments) * seq(0, arrow_segments - 1))
        x2 <- mi - (((mi - ma)/arrow_segments) * seq(arrow_segments - 1, 0))
      }
      
      y1 <- rep(dt[i, level], arrow_segments)
      y2 <- y1
      
      adt <- rbind(adt,
                   data.table::data.table(
                     x1 = x1,
                     x2 = x2,
                     y1 = y1,
                     y2 = y2,
                     transcript_id = rep(dt[i, transcript_id],
                                         arrow_segments),
                     transcript_name = rep(dt[i, transcript_name],
                                           arrow_segments),
                     gene_id = rep(dt[i, gene_id],
                                   arrow_segments),
                     gene_name = rep(dt[i, gene_name],
                                     arrow_segments)))
    }
    
    return (adt)
  }
  
  
  trans_text <- function(dt) {
    tdt <- data.table::data.table()
    
    for (i in seq_len(nrow(dt))) {
      mi <- dt[i, start]
      ma <- dt[i, end]
      
      x <- (ma + mi)/2
      y <- dt[i, level] + 0.4
      
      tdt <- rbind(tdt,
                   data.table::data.table(
                     x = x,
                     y = y,
                     transcript_name = dt[i, transcript_name]))
      
    }
    return (tdt)
  }
  
  
  # convert to data.table
  gtf.dt <- data.table::data.table(
    refGenome::getGtf(ensemblGenome)[, c("id",
                                         "seqid",
                                         "feature",
                                         "start",
                                         "end",
                                         "strand",
                                         "gene_biotype",
                                         "gene_name",
                                         "exon_number",
                                         "gene_id",
                                         "transcript_name",
                                         "transcript_id")])
  
  # use new variables to avoid ambiguity
  begin <- start
  stop <- end
  
  # subset features
  gtf.dt <- gtf.dt[end >= begin & start <= stop & seqid == chr, ]
  
  # aggregate transcripts
  txdt <- gettxdt(gtf.dt)
  
  # get tables for plotting
  rectdt <- trans_rect(gtf.dt, txdt)
  arrowdt <- trans_arrow(txdt)
  textdt <- trans_text(txdt)
  
  # plot
  g <- ggplot2::ggplot() +
    ggplot2::geom_rect(data = rectdt,
              mapping = ggplot2::aes(xmin = x1,
                            xmax = x2,
                            ymin = y1,
                            ymax = y2)) +
    ggplot2::geom_segment(data = arrowdt,
                 mapping = ggplot2::aes(x = x1,
                               y = y1,
                               xend = x2,
                               yend = y2),
                 size = line_width,
                 arrow = ggplot2::arrow(angle = arrow_width,
                               length = ggplot2::unit(arrow_length, "inches"),
                               type = arrow_type)) +
    ggplot2::geom_text(data = textdt,
              mapping = ggplot2::aes(x = x,
                            y = y,
                            label = transcript_name),
              size = text_size) +
    theme_Publication() +
    ggplot2::theme(axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank()) +
    ggplot2::xlab(paste0("Chr", chr))
  
  return (g)
}



