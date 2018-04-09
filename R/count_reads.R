
#' Count the number of UMIs for each transcript/gene and generate the expression matrix
#' 
#' Count unique \emph{UMI:transcript} pairs for single cell RNA-sequencing alignment files. Write resulting table to output directory. Columns are samples (cells) and rows are transcript/gene IDs.
#' 
#' @param alignment A character vector of the paths to input alignment files.
#' @param gtf Path to the gtf reference file. For generation of TxDb objects from gtf files, please refer to \code{makeTxDbFromGFF} function in \code{GenomicFeatures} package.
#' @param format Format of input sequence alignment files. \strong{"BAM"} or \strong{"SAM"}. Default is \strong{"BAM"}.
#' @param out.dir Output directory for UMI counting results. Expression table will be stored in this directory. Default is \code{"../Count"}.
#' @param cores Number of cores used for parallelization. Default is \code{max(1, parallel::detectCores() / 2)}.
#' @param output.prefix Prefix for expression table filename. Default is \code{"countUMI"}.
#' @param verbose Print log messages. Useful for debugging. Default to \strong{FALSE}.
#' @param logfile.prefix Prefix for log file. Default is current date and time in the format of \code{format(Sys.time(), "\%Y\%m\%d_\%H\%M\%S")}.
#' @return A expression matrix \code{data.table} containing the raw counts of unique \emph{UMI:transcript} pairs.
#' @import data.table foreach
#' @export
count.reads <- function(alignment,
                        gtf,
                        format = "BAM",
                        out.dir = "./Count",
                        cores = max(1, parallel::detectCores() / 2),
                        output.prefix = "count",
                        verbose = FALSE,
                        logfile.prefix = format(Sys.time(), "%Y%m%d_%H%M%S")) {
  
  logfile <- paste0(logfile.prefix, "_count_log.txt")
  
  if (verbose) {
    log.messages(Sys.time(),
                 "... Start read counting",
                 logfile = logfile,
                 append = FALSE)
    log.messages(Sys.time(), alignment, logfile = logfile, append = TRUE)
    print("... Input alignment files:")
    print(alignment)
  } else {
    log.messages(Sys.time(),
                 "... Start read counting",
                 logfile = NULL,
                 append = FALSE)
  }
  
  message(paste(Sys.time(),
                "... Creating output directory",
                out.dir))
  dir.create(file.path(out.dir),
             showWarnings = FALSE,
             recursive = TRUE)
  
  message(paste(Sys.time(),
                paste("... Loading TxDb file")))
  
  features.tx <- suppressPackageStartupMessages(
    gtf.db.read(gtf, logfile, grouping = "tx"))
  
  #features.gene <- suppressPackageStartupMessages(
  #  gtf.db.read(features, logfile, grouping = "gene"))
  #tx.to.gene <- GenomicAlignments::findOverlaps(features.gene, features.tx)
  #tx.to.gene <- data.table::data.table(
  #  queryHits = S4Vectors::queryHits(tx.to.gene),
  #  subjectHits = S4Vectors::subjectHits(tx.to.gene),
  #  gene.id = names(features.gene)[S4Vectors::queryHits(tx.to.gene)],
  #  tx.id = names(features.tx)[S4Vectors::subjectHits(tx.to.gene)])
  
  
  # parallelization
  message(paste(Sys.time(),
                paste("... counting")))
  
  cl <- if (verbose)
    parallel::makeCluster(cores, outfile = logfile)
  else
    parallel::makeCluster(cores)
  doParallel::registerDoParallel(cl)
  
  if (format == "SAM") {
    alignment <- foreach::foreach(
      i = alignment,
      .verbose = verbose,
      .combine = c,
      .multicombine = TRUE
    ) %dopar% {
      to.bam(i, logfile, overwrite = FALSE, index = TRUE)
    }
  }
  
  expr.tx <- foreach::foreach(
    i = alignment,
    .verbose = verbose,
    .combine = cbind,
    .multicombine = TRUE,
    .packages = c("BiocGenerics", "S4Vectors",
                  "GenomicFeatures", "GenomicAlignments")
  ) %dopar% {
    count.tx.unit(i, features.tx, format, logfile, verbose)
  }
  
  parallel::stopCluster(cl)
  
  expr.tx <- data.table::data.table(expr.tx, keep.rownames = TRUE)
  colnames(expr.tx)[1] <- "transcript_id"
  
  print(paste(Sys.time(), paste(
    "... Write transcript expression table to",
    file.path(out.dir, paste0(
      format(Sys.time(),
             "%Y%m%d_%H%M%S"), "_",
      output.prefix, "_tx.tab"
    ))
  )))
  
  data.table::fwrite(expr.tx, file.path(out.dir, paste0(
    format(Sys.time(), "%Y%m%d_%H%M%S"), "_",
    output.prefix, "_tx.tab"
  )), sep = "\t")
  
  # get counts for genes
  print(paste(Sys.time(), "... generating gene expression table"))
  
  gtf.eg <- refGenome::ensemblGenome(basedir = dirname(gtf))
  refGenome::read.gtf(gtf.eg, filename = basename(gtf))
  gtf.dt.tx.to.gene <- unique(data.table::data.table(
    refGenome::getGtf(gtf.eg)[, c("transcript_id", "gene_id")]))
  
  expr.gene <- merge(expr.tx, gtf.dt.tx.to.gene,
                     by = "transcript_id")
  expr.gene <- expr.gene[order(gene_id), ]
  
  expr.gene <- data.table::data.table(
    stats::aggregate(.~gene_id, expr.gene[, -"transcript_id"], sum))
  last2rows <- expr.tx[transcript_id %in%
                         c("reads_mapped_to_genome",
                           "reads_mapped_to_transcripts"), ]
  colnames(last2rows)[names(last2rows) == "transcript_id"] <- "gene_id"
  expr.gene <- rbind(expr.gene, last2rows)
  
  print(paste(Sys.time(), paste(
    "... Write gene expression table to",
    file.path(out.dir, paste0(
      format(Sys.time(),
             "%Y%m%d_%H%M%S"), "_",
      output.prefix, "_gene.tab"
    ))
  )))
  
  data.table::fwrite(expr.gene, file.path(out.dir, paste0(
    format(Sys.time(), "%Y%m%d_%H%M%S"), "_",
    output.prefix, "_gene.tab"
  )), sep = "\t")
  
  message(paste(Sys.time(), "... read counting finished!"))
  return(expr.gene)
}


count.tx.unit <- function(i, features.tx, format, logfile, verbose) {
  if (verbose) {
    log.messages(Sys.time(),
                 "... counting reads in sample",
                 i,
                 logfile = logfile,
                 append = TRUE)
  }
  
  # if sequence alignment file is empty
  if (file.size(i) == 0) {
    count.tx.dt <- data.table::data.table(
      tx.id = c(names(features.tx),
                  "reads_mapped_to_genome",
                  "reads_mapped_to_transcripts"))
    cell <- remove.last.extension(i)
    count.tx.dt[[cell]] <- 0
    
    return (data.frame(count.tx.dt,
                       row.names = 1,
                       check.names = FALSE,
                       fix.empty.names = FALSE))
  }
  
  bfl <- Rsamtools::BamFile(i)
  bamGA <- GenomicAlignments::readGAlignments(bfl, use.names = T)
  
  genome.reads <- data.table::data.table(
    name = names(bamGA),
    seqnames = as.vector(GenomicAlignments::seqnames(bamGA)))
  
  if (length(unique(genome.reads[, name])) != nrow(genome.reads)) {
    log.messages(Sys.time(),
                 i,
                 " has duplicate read alignments.",
                 " This is fine if multi-alignment is allowed.",
                 " Otherwise, try re-running demultiplexing",
                 " and alignment functions",
                 " with appropriate number of cores.",
                 logfile = logfile,
                 append = TRUE)
  }
  
  # reads mapped to genome (exclude ERCC spike-in)
  reads.mapped.to.genome <- nrow(
    genome.reads[!grepl("ERCC", genome.reads[, seqnames]), .(name)])
  
  # umi filtering
  ol <- GenomicAlignments::findOverlaps(features.tx, bamGA)
  
  ol.dt <- data.table::data.table(
    tx.id = base::names(features.tx)[S4Vectors::queryHits(ol)],
    name = base::names(bamGA)[S4Vectors::subjectHits(ol)],
    pos = BiocGenerics::start(bamGA)[S4Vectors::subjectHits(ol)]
  )
  
  #dt <- merge(ol.dt, tx.to.gene[, .(tx.id, gene.id)], by = "tx.id")
  
  if (nrow(ol.dt) == 0) {
    reads.mapped.to.tx <- 0
    
    # clean up
    count.tx.dt <- data.table::data.table(
      tx.id = c(names(features.tx),
                "reads_mapped_to_genome",
                "reads_mapped_to_transcripts"))
    cell <- remove.last.extension(i)
    count.tx.dt[[cell]] <- 0
    count.tx.dt[tx.id == "reads_mapped_to_genome",
                 cell] <- reads.mapped.to.genome
    count.tx.dt[tx.id == "reads_mapped_to_transcripts",
                 cell] <- reads.mapped.to.tx
    
    # coerce to data frame to keep rownames for cbind combination
    count.tx.dt <- data.frame(count.tx.dt,
                              row.names = 1,
                              check.names = FALSE,
                              fix.empty.names = FALSE)
    
  } else {
    
    ol.dt[, umi := data.table::last(data.table::tstrsplit(name, ":"))]
    
    # remove ambiguous gene alignments (union mode filtering)
    #ol.dt <- ol.dt[!(
    #  base::duplicated(ol.dt, by = "name") |
    #base::duplicated(ol.dt, by = "name", fromLast = TRUE)
    #), ]
    
    # reads mapped to transcripts
    reads.mapped.to.tx <- nrow(ol.dt[!grepl("ERCC", ol.dt[, tx.id]), ])
    
    # umi filtering
    count.tx <- base::table(unique(ol.dt[, .(tx.id, umi, pos)])[, tx.id])
    
    # clean up
    count.tx.dt <- data.table::data.table(
      tx.id = c(names(features.tx),
                "reads_mapped_to_genome",
                "reads_mapped_to_transcripts"))
    cell <- remove.last.extension(i)
    count.tx.dt[[cell]] <- 0
    count.tx.dt[tx.id == "reads_mapped_to_genome",
                 cell] <- reads.mapped.to.genome
    count.tx.dt[tx.id == "reads_mapped_to_transcripts",
                 cell] <- reads.mapped.to.tx
    count.tx.dt[tx.id %in% names(count.tx),
                 eval(cell) := as.numeric(count.tx[tx.id])]
    
    # coerce to data frame to keep rownames for cbind combination
    count.tx.dt <- data.frame(count.tx.dt,
                              row.names = 1,
                              check.names = FALSE,
                              fix.empty.names = FALSE)
  }
  return (count.tx.dt)
}




# 
# count.umi.unit <- function(i, features.tx, format, logfile, verbose, tx.to.gene) {
#   if (verbose) {
#     log.messages(Sys.time(),
#                  "... UMI counting sample",
#                  i,
#                  logfile = logfile,
#                  append = TRUE)
#   }
#   
#   # if sequence alignment file is empty
#   if (file.size(i) == 0) {
#     count.umi.dt <- data.table::data.table(
#       gene.id = c(tx.to.gene[, unique(gene.id)],
#                   "reads_mapped_to_genome",
#                   "reads_mapped_to_genes"))
#     cell <- remove.last.extension(i)
#     count.umi.dt[[cell]] <- 0
#     
#     return (data.frame(count.umi.dt,
#                row.names = 1,
#                check.names = FALSE,
#                fix.empty.names = FALSE))
#   }
#   
#   bfl <- Rsamtools::BamFile(i)
#   bamGA <- GenomicAlignments::readGAlignments(bfl, use.names = T)
#   
#   genome.reads <- data.table::data.table(
#     name = names(bamGA),
#     seqnames = as.vector(GenomicAlignments::seqnames(bamGA)))
#   
#   if (length(unique(genome.reads[, name])) != nrow(genome.reads)) {
#     print(paste0(i,
#                  " has duplicate read alignments.",
#                  " This is fine if multi-alignment is allowed.",
#                  " Otherwise, try re-running demultiplexing",
#                  " and alignment functions",
#                  " with appropriate number of cores."))
#   }
#   
#   # reads mapped to genome (exclude ERCC spike-in)
#   reads.mapped.to.genome <- nrow(
#     genome.reads[!grepl("ERCC", genome.reads[, seqnames]), .(name)])
#   
#   # umi filtering
#   ol <- GenomicAlignments::findOverlaps(features.tx, bamGA)
#   
#   ol.dt <- data.table::data.table(
#     tx.id = base::names(features.tx)[S4Vectors::queryHits(ol)],
#     name = base::names(bamGA)[S4Vectors::subjectHits(ol)],
#     pos = BiocGenerics::start(bamGA)[S4Vectors::subjectHits(ol)],
#     gene.id = tx.to.gene[, gene.id][S4Vectors::queryHits(ol)]
#   )
#   
#   dt <- merge(ol.dt, tx.to.gene[, .(tx.id, gene.id)], by = "tx.id")
#   
#   if (nrow(ol.dt) == 0) {
#     reads.mapped.to.genes <- 0
#     
#     # clean up
#     count.umi.dt <- data.table::data.table(gene.id = c(names(features),
#                                                        "reads_mapped_to_genome",
#                                                        "reads_mapped_to_genes"))
#     cell <- remove.last.extension(i)
#     count.umi.dt[[cell]] <- 0
#     count.umi.dt[gene.id == "reads_mapped_to_genome",
#                  cell] <- reads.mapped.to.genome
#     count.umi.dt[gene.id == "reads_mapped_to_genes",
#                  cell] <- reads.mapped.to.genes
#     
#     # coerce to data frame to keep rownames for cbind combination
#     count.umi.dt <- data.frame(count.umi.dt,
#                                row.names = 1,
#                                check.names = FALSE,
#                                fix.empty.names = FALSE)
#     
#   } else {
#     
#     ol.dt[, umi := data.table::last(data.table::tstrsplit(name, ":"))]
#     
#     # remove ambiguous gene alignments (union mode filtering)
#     #ol.dt <- ol.dt[!(
#     #  base::duplicated(ol.dt, by = "name") |
#     #base::duplicated(ol.dt, by = "name", fromLast = TRUE)
#     #), ]
#     
#     # reads mapped to genes
#     reads.mapped.to.genes <- nrow(ol.dt[!grepl("ERCC", ol.dt[, gene.id ]), ])
#     
#     # umi filtering
#     count.umi <- base::table(unique(ol.dt[, .(gene.id, umi, pos)])[, gene.id])
#     
#     # clean up
#     count.umi.dt <- data.table::data.table(gene.id = c(names(features),
#                                                        "reads_mapped_to_genome",
#                                                        "reads_mapped_to_genes"))
#     cell <- remove.last.extension(i)
#     count.umi.dt[[cell]] <- 0
#     count.umi.dt[gene.id == "reads_mapped_to_genome",
#                  cell] <- reads.mapped.to.genome
#     count.umi.dt[gene.id == "reads_mapped_to_genes",
#                  cell] <- reads.mapped.to.genes
#     count.umi.dt[gene.id %in% names(count.umi),
#                  eval(cell) := as.numeric(count.umi[gene.id])]
#     
#     # coerce to data frame to keep rownames for cbind combination
#     count.umi.dt <- data.frame(count.umi.dt,
#                                row.names = 1,
#                                check.names = FALSE,
#                                fix.empty.names = FALSE)
#   }
#   return (count.umi.dt)
# }

