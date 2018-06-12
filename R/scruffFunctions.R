########################################
####### parsing helper functions #######
########################################

.logMessages <- function(...,
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


.stripLeadingUnderscore <- function (x)  sub("^\\_+", "", x)


# remove file extension and get basename
.removeLastExtension <- function(x) {
  return (sub(pattern = "\\.[^\\.]*$",
              replacement = "\\1",
              basename(x)))
}


.getAlignmentFilePaths <- function(fastq.paths, format, out.dir) {
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
.gtfReadDb <- function(gtf, logfile) {
  gtf.db.file <- paste0(basename(gtf), ".sqlite")
  if ((!(file.exists(gtf))) & (!(file.exists(gtf.db.file)))) {
    stop(paste("File", gtf, "does not exist"))
  }
  
  if (!(file.exists(gtf.db.file))) {
    message(paste(Sys.time(), "... TxDb file", gtf.db.file, "does not exist"))
    message(paste(Sys.time(), "... Creating TxDb object", gtf.db.file))
    gtf.db <- GenomicFeatures::makeTxDbFromGFF(file = gtf)
    AnnotationDbi::saveDb(gtf.db, file = gtf.db.file)
    return (GenomicFeatures::exonsBy(gtf.db, by = "gene"))
  }
  
  gtf.db <- tryCatch(
    suppressPackageStartupMessages(AnnotationDbi::loadDb(gtf.db.file)),
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
.bcCorrectMem <- local({
  res <- list()
  
  f <- function(bc, refBarcodes, maxEditDist) {
    if (bc %in% names(res))
      return (res[[bc]])
    if (bc %in% refBarcodes) {
      res[[bc]] <<- bc
      return (res[[bc]])
    }
    
    sdm <- stringdist::stringdistmatrix(bc,
                                        refBarcodes,
                                        method = "hamming",
                                        nthread = 1)
    min.dist <- min(sdm)
    if (min.dist <= maxEditDist) {
      ind <- which(sdm == min.dist)
      if (length(ind) == 1) {
        res[[bc]] <<- refBarcodes[ind]
        return (res[[bc]])
      }
    }
    res[[bc]] <<- bc
    return (res[[bc]])
  }
})


.toBam <- function(sam,
                   logfile,
                   overwrite = FALSE,
                   index = FALSE) {
  .logMessages(
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
    ignore.case = TRUE,
    perl = TRUE,
    replacement = ".BAM",
    x = sam
  ))
}


#' Visualize aligned reads
#' 
#' Visualize read alignments for UMI tagged single cell RNA-sequencing data. Read names must contain UMI sequences at the end delimited by "\strong{:}". Arrow represents orientation of alignment. Reads are colored by their UMI and sorted by their start positions and UMI.
#' 
#' @param bamGA A GenomicAlignment object
#' @param chr Chromosome. Integer or "X", "Y", "MT".
#' @param start Genomic coordinate of the start position.
#' @param end Genomic coordinate of the end position.
#' @param legend Show legend. Default is FALSE.
#' @return A ggplot object of aligned reads
#' @examples
#' data(bamExample, package = "scruff")
#' rview(bamExample, chr = "MT", legend = TRUE)
#' @import ggbio
#' @export
rview <- function(bamGA,
                  chr = "1",
                  start = 1,
                  end = max(BiocGenerics::end(bamGA)),
                  legend = FALSE) {
  
  reads <- bamGA[BiocGenerics::start(bamGA) >= start &
                   BiocGenerics::end(bamGA) <= end &
                   GenomeInfoDb::seqnames(bamGA) == chr]
  umi <- data.table::last(data.table::tstrsplit(names(reads), ":"))
  reads <- reads[order(BiocGenerics::start(reads), umi)]
  
  readsGr <- GenomicRanges::GRanges(reads)
  S4Vectors::mcols(readsGr)$umi <- data.table::last(
    data.table::tstrsplit(names(readsGr), ":"))
  
  g <- ggplot2::ggplot(readsGr) +
    ggbio::geom_arrow(ggplot2::aes(color = umi)) +
    .themePublication() + 
    ggplot2::theme(axis.title.y = ggplot2::element_blank())
  
  if (legend == FALSE) {
    g <- g + ggplot2::theme(legend.position="none")
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
#' @param arrow_segments The number of segments lines be divided to. The greater the number, more arrows there are. Default 10.
#' @param arrow_width The angle of the arrow head in degrees (smaller numbers produce narrower, pointier arrows). Essentially describes the width of the arrow head. Passed to the angle parameter of arrow function. Default 30.
#' @param arrow_length The length of the arrow head. Passed to the length argument of arrow function. Default 0.08.
#' @param arrow_type One of "open" or "closed" indicating whether the arrow head should be a closed triangle. Passed to the type argument of arrow function. Default "open".
#' @param text_size Size of text. Passed to the size argument of the geom_text function. Default 4.
#' @return A ggplot object of genomic view
#' @examples
#' gtf <- system.file("extdata", "GRCm38_MT.gtf", package = "scruff")
#' gtfEG = refGenome::ensemblGenome(dirname(gtf))
#' refGenome::read.gtf(gtfEG, filename = basename(gtf))
#' 
#' gview(gtfEG, chr = "MT")
#' @import refGenome
#' @export
gview <- function(ensemblGenome,
                  chr = 1,
                  start = 1,
                  end = max(refGenome::getGtf(ensemblGenome)$end),
                  rect_width = 0.3,
                  line_width = 0.5,
                  arrow_segments = 10,
                  arrow_width = 30,
                  arrow_length = 0.08,
                  arrow_type = "open",
                  text_size = 4) {
  
  .getLevel <- function(txdt) {
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
  
  
  .getTxdt <- function(dt) {
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
    
    txdt <- .getLevel(txdt)
  }
  
  
  .transRect <- function(dt, txdt) {
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
  
  
  .transArrow <- function(dt) {
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
  
  
  .transText <- function(dt) {
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
  gtfDt <- data.table::data.table(
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
  st <- end
  
  # subset features
  gtfDt <- gtfDt[end >= begin & start <= st & seqid == chr, ]
  
  # aggregate transcripts
  txdt <- .getTxdt(gtfDt)
  
  # get tables for plotting
  rectdt <- .transRect(gtfDt, txdt)
  arrowdt <- .transArrow(txdt)
  textdt <- .transText(txdt)
  
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
                                                 length = ggplot2::unit(
                                                   arrow_length, "inches"),
                                                 type = arrow_type)) +
    ggplot2::geom_text(data = textdt,
                       mapping = ggplot2::aes(x = x,
                                              y = y,
                                              label = transcript_name),
                       size = text_size) +
    .themePublication() +
    ggplot2::theme(axis.title.y = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_blank(),
                   axis.ticks.y = ggplot2::element_blank(),
                   axis.line.y = ggplot2::element_blank()) +
    ggplot2::xlab(paste0("Chr", chr))
  
  return (g)
}



