
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
                        unique.alignment = TRUE,
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
    gtf.db.read(gtf, feature = "exon", grouping = "transcript"))
  
  features.tx.dt <- data.table::data.table(
    txid = rep(names(features.tx), elementNROWS(features.tx)),
    exon_start = unlist(start(features.tx)),
    exon_end = unlist(end(features.tx)),
    strand = S4Vectors::decode(unlist(strand(features.tx),
                                      use.names = FALSE)),
    length = unlist(end(features.tx)) - unlist(start(features.tx))
  )
  
  features.tx.dt <- features.tx.dt[order(txid, exon_start, exon_end), ]
  
  features.gene <- suppressPackageStartupMessages(
    gtf.db.read(gtf, feature = "transcript", grouping = "gene"))
  
  features.tx.rg <- range(features.tx)
  
  gtf.eg <- refGenome::ensemblGenome(basedir = dirname(gtf))
  refGenome::read.gtf(gtf.eg, basename(gtf))
  gtf.eg.dt <- data.table::data.table(refGenome::getGtf(gtf.eg))
  txtogene <- unique(gtf.eg.dt[,.(gene_id, transcript_id)])
  txtogene <- txtogene[order(transcript_id), ]
  colnames(txtogene)[2] <- "txid"
  
  #features.gene <- suppressPackageStartupMessages(
  #  gtf.db.read(features, logfile, grouping = "gene"))
  #tx.to.gene <- GenomicAlignments::findOverlaps(features.gene, features.tx)
  #tx.to.gene <- data.table::data.table(
  #  queryHits = S4Vectors::queryHits(tx.to.gene),
  #  subjectHits = S4Vectors::subjectHits(tx.to.gene),
  #  gene.id = names(features.gene)[S4Vectors::queryHits(tx.to.gene)],
  #  txid = names(features.tx)[S4Vectors::subjectHits(tx.to.gene)])
  
  
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
  
  expr.gene <- foreach::foreach(
    i = alignment,
    .verbose = verbose,
    .combine = cbind,
    .multicombine = TRUE,
    .packages = c("BiocGenerics", "S4Vectors",
                  "GenomicFeatures", "GenomicAlignments")
  ) %dopar% {
    count.tx.unit(i,
                  features.tx,
                  features.tx.dt,
                  features.gene,
                  features.tx.rg,
                  txtogene,
                  unique.alignment,
                  format,
                  logfile,
                  verbose)
  }
  
  parallel::stopCluster(cl)
  
  expr.gene <- data.table::data.table(expr.gene, keep.rownames = TRUE)
  colnames(expr.gene)[1] <- "gene_id"
  
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
  
  # get counts for genes
#   print(paste(Sys.time(), "... generating gene expression table"))
#   
#   gtf.eg <- refGenome::ensemblGenome(basedir = dirname(gtf))
#   refGenome::read.gtf(gtf.eg, filename = basename(gtf))
#   gtf.dt.tx.to.gene <- unique(data.table::data.table(
#     refGenome::getGtf(gtf.eg)[, c("transcript_id", "gene_id")]))
#   
#   expr.gene <- merge(expr.tx, gtf.dt.tx.to.gene,
#                      by = "transcript_id")
#   expr.gene <- expr.gene[order(gene_id), ]
#   
#   expr.gene <- data.table::data.table(
#     stats::aggregate(.~gene_id, expr.gene[, -"transcript_id"], sum))
#   last2rows <- expr.tx[transcript_id %in%
#                          c("reads_mapped_to_genome",
#                            "reads_mapped_to_transcripts"), ]
#   colnames(last2rows)[names(last2rows) == "transcript_id"] <- "gene_id"
#   expr.gene <- rbind(expr.gene, last2rows)
#   
#   print(paste(Sys.time(), paste(
#     "... Write gene expression table to",
#     file.path(out.dir, paste0(
#       format(Sys.time(),
#              "%Y%m%d_%H%M%S"), "_",
#       output.prefix, "_gene.tab"
#     ))
#   )))
#   
#   data.table::fwrite(expr.gene, file.path(out.dir, paste0(
#     format(Sys.time(), "%Y%m%d_%H%M%S"), "_",
#     output.prefix, "_gene.tab"
#   )), sep = "\t")
  
  message(paste(Sys.time(), "... read counting finished!"))
  return(expr.gene)
}


count.tx.unit <- function(cell,
                          features.tx,
                          features.tx.dt,
                          features.gene,
                          features.tx.rg,
                          txtogene,
                          unique.alignment,
                          format,
                          logfile,
                          verbose) {
  
  if (verbose) {
    log.messages(Sys.time(),
                 "... counting reads in sample",
                 cell,
                 logfile = logfile,
                 append = TRUE)
  }
  
  # if sequence alignment file is empty
  if (file.size(cell) == 0) {
    count.gene.dt <- data.table::data.table(
      gene_id = c(names(features.gene),
                  "reads_mapped_to_genome",
                  "reads_mapped_to_genes"))
    cellname <- remove.last.extension(cell)
    count.gene.dt[[cellname]] <- 0
    
    return (data.frame(count.gene.dt,
                       row.names = 1,
                       check.names = FALSE,
                       fix.empty.names = FALSE))
  }
  
  bfl <- Rsamtools::BamFile(cell)
  bamGA <- GenomicAlignments::readGAlignments(bfl, use.names = T)
  
  genome.reads <- data.table::data.table(
    name = names(bamGA),
    seqnames = as.vector(GenomicAlignments::seqnames(bamGA)))
  
  
  if (unique.alignment) {
    if (length(unique(genome.reads[, name])) != nrow(genome.reads)) {
      stop(Sys.time(),
           " ",
           cell,
           " has duplicate read alignments.",
           " This is fine if multi-alignment is allowed",
           " (change unique.alignment to FALSE).",
           " Otherwise, try re-running demultiplexing",
           " and alignment functions",
           " with appropriate number of cores.")
    }
  }
  
  # reads mapped to genome (exclude ERCC spike-in)
  reads.mapped.to.genome <- nrow(
    genome.reads[!grepl("ERCC", genome.reads[, seqnames]), .(name)])
  
  # count overlaps
  ol <- GenomicAlignments::findOverlaps(features.tx, bamGA)
  
  ol.dt <- data.table::data.table(
    txid = base::names(features.tx)[S4Vectors::queryHits(ol)],
    readname = base::names(bamGA)[S4Vectors::subjectHits(ol)],
    readstart = BiocGenerics::start(bamGA)[S4Vectors::subjectHits(ol)],
    readend = BiocGenerics::end(bamGA)[S4Vectors::subjectHits(ol)],
    strand = S4Vectors::decode(
      BiocGenerics::strand(bamGA))[S4Vectors::subjectHits(ol)],
    txstart = unlist(BiocGenerics::start(features.tx.rg))
    [S4Vectors::queryHits(ol)],
    txend = unlist(BiocGenerics::end(features.tx.rg))[S4Vectors::queryHits(ol)]
  )
  
  #dt <- merge(ol.dt, tx.to.gene[, .(txid, gene.id)], by = "txid")
  
  if (nrow(ol.dt) == 0) {
    reads.mapped.to.gene <- 0
    # clean up
    count.gene.dt <- data.table::data.table(
      gene_id = c(names(features.gene),
                  "reads_mapped_to_genome",
                  "reads_mapped_to_genes"))
    cellname <- remove.last.extension(cell)
    count.gene.dt[[cellname]] <- 0
    count.gene.dt[gene_id == "reads_mapped_to_genome",
                  cellname] <- reads.mapped.to.genome
    count.gene.dt[gene_id == "reads_mapped_to_transcripts",
                  cellname] <- reads.mapped.to.gene
    
    # coerce to data frame to keep rownames for cbind combination
    return (data.frame(count.gene.dt,
                       row.names = 1,
                       check.names = FALSE,
                       fix.empty.names = FALSE))
    
  } else {
    
    get.threeprime.dt <- function(dt) {
      if (dt[, unique(strand)] == "+") {
        threeprime.dt <- dt[txend == min(txend), ]
        # randomly choose one if transcripts have the same min 3' ends
        threeprime.dt <- threeprime.dt[base::sample(nrow(threeprime.dt), 1), ]
      } else if (dt[,unique(strand)] == "-") {
        threeprime.dt <- dt[txstart == max(txstart), ]
        # randomly choose one if transcripts have the same max 5' starts
        threeprime.dt <- threeprime.dt[base::sample(nrow(threeprime.dt), 1), ]
      }
      return (threeprime.dt)
    }
    
    
    # fix unique reads mapped to multiple transcripts
    uni.ali.tx.fix <- function(dt) {
      # if reads mapped to the same gene
      if (length(dt[, unique(gene_id)]) == 1) {
        # select the tx with alignment most close to the 3' end
        return (get.threeprime.dt(dt))
      } else {
        # if reads mapped to multiple genes
        #res.dt <- data.table::data.table()
        #for (i in dt[, unique(gene_id)])
      }
    }
    
    # get relative position of a uniquely aligned read
    getrp <- function(read.dt, ex.dt) {
      if (read.dt[, strand] == "+") {
        exonmatch <- data.table::first(
          ex.dt[txid == read.dt[, txid] &
                  exon_start <= read.dt[, readend] &
                  exon_end >= read.dt[, readstart]])
        tx.ex.dt <- ex.dt[txid == exonmatch[, txid]]
        
        rp <- (min(read.dt[, readend] -
                     exonmatch[, exon_start], exonmatch[, length]) +
                 tx.ex.dt[exon_end < exonmatch[, exon_start],
                          sum(length)]) /
          tx.ex.dt[, sum(length)]
      } else if (read.dt[, strand] == "-") {
        exonmatch <- data.table::last(
          ex.dt[txid == read.dt[, txid] &
                  exon_start <= read.dt[, readend] &
                  exon_end >= read.dt[, readstart]])
        tx.ex.dt <- ex.dt[txid == exonmatch[, txid]]
        
        rp <- (min(exonmatch[, exon_end] -
                     read.dt[, readstart], exonmatch[, length]) +
                 tx.ex.dt[exon_start > exonmatch[, exon_end],
                          sum(length)]) /
          tx.ex.dt[, sum(length)]
      }
      return (rp)
    }
    
    
    # get relative positions of each read to the 3' end of the aligned gene
    getrp.tx <- function(ex.dt, ol.dt.uni) {
      # filter out mt-Rnr2 which shows abnormal positional probabilities
      ol.dt.uni <- ol.dt.uni[gene_id != "ENSMUSG00000064339", ]
      
      rpvc <- vector(mode = "double",
                     length = length(ol.dt.uni[, unique(readname)]))
      j <- 1
      
      for (i in ol.dt.uni[, unique(readname)]) {
        read.dt <- ol.dt.uni[readname == i, ]
        if (nrow(read.dt) > 1) {
          read.dt <- uni.ali.tx.fix(read.dt)
        }
        
        rp <- getrp(read.dt, ex.dt)
        
        rpvc[j] <- rp
        j <- j + 1
        
        #if (j %% 2000 == 0) {
        #  print(paste(j, "reads processed"))
        #}
        
      }
      
      res.dt <- cbind(unique(ol.dt.uni[, .(readname)]),
                      rpvc)
      
      return (res.dt)
    }
    
    
    multimapper.assign <- function(multi.gene.alignments,
                                   ex.dt,
                                   bincounts,
                                   gene.probs.dt) {
      
      sumreads <- gene.probs.dt[, sum(read_pseudo)]
      sumbincounts <- sum(bincounts + 1)
      counter <- 1
      cutpoints <- seq(0, 1, 0.01)
      multi.probs <- log((bincounts + 1) / sumbincounts)
      
      
      for (i in unique(multi.gene.alignments[, readname])) {
        read.alignments <- multi.gene.alignments[readname == i, ]
        
        for (j in seq_len(nrow(read.alignments))) {
          # positional probability
          rp <- getrp(read.alignments[j, ], ex.dt)
          mpb <- cut(rp,
                     cutpoints,
                     labels = FALSE,
                     include.lowest = TRUE)
          
          # k = 100, 
          #x <- rep(0, 100)
          #x[mpb] <- 1
          
          #pp <- dmultinom(x, prob = multi.probs)
          
          pp <- multi.probs[mpb]
          
          read.alignments[j, relative_position := rp]
          read.alignments[j, log_pos_prob := pp]
          
          # abundance probability
          read.alignments[j,
                          log_abund_prob := 
                            gene.probs.dt[gene_id == read.alignments[j, gene_id],
                                          log_abund_prob]]
          
          # joint probability
          read.alignments[j, joint := log_pos_prob + log_abund_prob]
        }
        
        # after calculating all the joints of the read for each alignment
        # assign read to alignment with max(joint)
        # if tie, disregard
        
        if (length(unique(read.alignments[, joint])) > 1) {
          read.alignments[
            data.table::first(order(joint,
                                    decreasing = TRUE)), assign := TRUE]
          
          # abundance probability
          g <- read.alignments[order(joint,
                                     decreasing = TRUE), ][1,
                                                           gene_id]
          #gene.probs.dt[gene_id == g, reads := reads + 1]
          gene.probs.dt[gene_id == g, read_pseudo := read_pseudo + 1]
          sumreads <- sumreads + 1
          
          # update gene.probs.dt
          gene.probs.dt[gene_id == g,
                        log_abund_prob := log(read_pseudo / sumreads)]
          
          # positional probability
          rpmax <- read.alignments[
            order(joint,
                  decreasing = TRUE), ][1,
                                        relative_position]
          updatempb <- cut(rpmax,
                           cutpoints,
                           labels = FALSE,
                           include.lowest = TRUE)
          bincounts[updatempb] <- bincounts[updatempb] + 1
          sumbincounts <- sumbincounts + 1
          
          # update multi.probs
          multi.probs <- log((bincounts + 1) / sumbincounts)
          
          assigned.dt <- read.alignments[assign == TRUE, ]
          
          
          multi.gene.alignments[readname == assigned.dt[, readname] &
                                  txid == assigned.dt[, txid] &
                                  readstart == assigned.dt[, readstart] &
                                  readend == assigned.dt[, readend] &
                                  strand == assigned.dt[, strand] &
                                  gene_id == assigned.dt[, gene_id],
                                assign := TRUE]
          
          #multi.gene.alignments <- merge(multi.gene.alignments,
          #                               read.alignments[assign == TRUE,
          #                                               .(readname,
          #                                                 txid,
          #                                                 readstart,
          #                                                 readend,
          #                                                 strand,
          #                                                 gene_id,
          #                                                 assign)],
          #                               by = c("readname",
          #                                      "txid",
          #                                      "readstart",
          #                                      "readend",
          #                                      "strand",
          #                                      "gene_id"),
          #                               all.x = TRUE)
        }
        
        multi.gene.alignments[readname == i,
                              log_pos_prob :=
                                read.alignments[, log_pos_prob]]
        multi.gene.alignments[readname == i,
                              log_abund_prob :=
                                read.alignments[, log_abund_prob]]
        multi.gene.alignments[readname == i,
                              joint := read.alignments[, joint]]
        
        counter <- counter + 1
        if (counter %% 100 == 0) {
          print(paste(counter, "reads processed"))
        }
        
      }
      return (multi.gene.alignments)
    }
    
    
    # ol.dt[, umi := data.table::last(data.table::tstrsplit(readname, ":"))]
    
    # too slow
    #     gettx <- function(ft) {
    #       return (mcols(ft[,2]))
    #     }
    #     
    #     txtogene <- data.table::as.data.table(
    #       cbind(rep(names(features.gene),
    #                 elementNROWS(features.gene)),
    #             do.call("rbind",
    #                     sapply(features.gene, gettx))))
    #     colnames(txtogene) <- c("gene_id", "txid")
    
    ol.dt <- merge(ol.dt, txtogene)
    
    c1.dt <- unique(ol.dt[, .(readname, readstart)])
    
    # uniquely mapped reads
    c1.dt.uni <- c1.dt[!(
      base::duplicated(c1.dt, by = "readname") |
        base::duplicated(c1.dt, by = "readname", fromLast = TRUE)
    ), ]
    
    c1.dt.uni2 <- merge(c1.dt.uni,
                        unique(ol.dt[, .(readname, gene_id)]),
                        by = "readname")
    
    # remove reads mapped to > 1 genes
    c1.dt.uni3 <- c1.dt.uni2[!(
      base::duplicated(c1.dt.uni2, by = "readname") |
        base::duplicated(c1.dt.uni2, by = "readname", fromLast = TRUE)
    ), ]
    
    
    c1.dt.uni4 <- merge(c1.dt.uni3,
                        unique(ol.dt[, .(txid, readname)]),
                        by = "readname")
    
    # remove reads mapped to > 1 tx
    c1.dt.uni5 <- c1.dt.uni4[!(
      base::duplicated(c1.dt.uni4, by = "readname") |
        base::duplicated(c1.dt.uni4, by = "readname", fromLast = TRUE)
    ), ]
    
    # number of unique mapping, unique gene, unique tx reads
    r1 <- length(unique(c1.dt.uni5[, readname]))
    
    # reads mapped to >= 1 tx
    c1.dt.uni6 <- c1.dt.uni4[(
      base::duplicated(c1.dt.uni4, by = "readname") |
        base::duplicated(c1.dt.uni4, by = "readname", fromLast = TRUE)
    ), ]
    
    # number of unique mapping, unique gene, >1 tx reads
    r2 <- length(unique(c1.dt.uni6[, readname]))
    
    # reads mapped to > 1 gene
    c1.dt.uni7 <- c1.dt.uni2[(
      base::duplicated(c1.dt.uni2, by = "readname") |
        base::duplicated(c1.dt.uni2, by = "readname", fromLast = TRUE)
    ), ]
    
    # number of unique mapping, >1 gene, >=1 tx reads
    r3 <- length(unique(c1.dt.uni7[, readname]))
    
    # reads mapped to > 1 positions
    c1.dt.8 <- c1.dt[(
      base::duplicated(c1.dt, by = "readname") |
        base::duplicated(c1.dt, by = "readname", fromLast = TRUE)
    ), ]
    
    c1.dt.9 <- merge(c1.dt.8,
                     unique(ol.dt[, .(gene_id, readname)]),
                     by = "readname",
                     allow.cartesian = TRUE)
    
    c1.dt.10 <- unique(c1.dt.9[, .(readname, gene_id)])
    
    # remove reads mapped to > 1 genes
    c1.dt.11 <- c1.dt.10[!(
      base::duplicated(c1.dt.10, by = "readname") |
        base::duplicated(c1.dt.10, by = "readname", fromLast = TRUE)
    ), ]
    
    # number of multi mapping, =1 gene, >=1 tx reads
    r4 <- length(unique(c1.dt.11[, readname]))
    
    # reads mapped to >1 gene
    c1.dt.12 <- c1.dt.10[(
      base::duplicated(c1.dt.10, by = "readname") |
        base::duplicated(c1.dt.10, by = "readname", fromLast = TRUE)
    ), ]
    
    # number of multi mapping, >1 gene, >=1 tx reads
    r5 <- length(unique(c1.dt.12[, readname]))
    
    
    # reads mapped to unique genes
    # c1.dt.uni3, c1.dt.11
    unique.gene.dt <- rbind(c1.dt.uni3[, .(readname, gene_id)], c1.dt.11)
    
    # initialize count matrix for reads aligned to unique genes
    # abundance probabilities
    read.counts <- base::table(unique.gene.dt[, .(gene_id)])
    count.reads.dt <- data.table::data.table(gene_id = names(features.gene))
    
    count.reads.dt[["reads"]] <- 0
    count.reads.dt[gene_id %in% names(read.counts),
                   reads := as.numeric(read.counts[gene_id])]
    count.reads.dt[, read_pseudo := reads + 1]
    count.reads.dt[, log_abund_prob := log(read_pseudo / sum(read_pseudo))]
    
    
    # positional probabilities
    
    unique.gene.alignments <- merge(unique.gene.dt[, .(readname)],
                                    ol.dt,
                                    by = "readname")
    
    rp.dt <- getrp.tx(features.tx.dt, unique.gene.alignments)
    
    cutpoints <- seq(0, 1, 0.01)
    multi.probs.bin <- cut(rp.dt[, rpvc],
                           cutpoints,
                           labels = FALSE,
                           include.lowest = TRUE)
    
    bincounts <- table(multi.probs.bin)
    
    # pseudocount
    # a bin might have no read
    # multi.probs <- (bincounts + 1) / sum(bincounts + 1)
    
    # iterate
    # reads mapped to >1 genes
    multi.gene.alignments <- merge(rbind(unique(c1.dt.12[, .(readname)]),
                                         unique(c1.dt.uni7[, .(readname)])),
                                   ol.dt,
                                   by = "readname")
    
    #start_time <- Sys.time()
    multimapper.dt <- multimapper.assign(multi.gene.alignments,
                                         features.tx.dt,
                                         bincounts,
                                         count.reads.dt)
    #end_time <- Sys.time()
    #print(end_time - start_time)
    
    multimapper.dt[,
                   umi := data.table::last(
                     data.table::tstrsplit(readname, ":"))]
    
    unique.gene.dt[, umi := data.table::last(
      data.table::tstrsplit(readname, ":"))]
    count.dt <- rbind(unique.gene.dt,
                      multimapper.dt[assign == TRUE, .(readname,
                                                       gene_id,
                                                       umi)])
    
    # remove ambiguous gene alignments (union mode filtering)
    #ol.dt <- ol.dt[!(
    #  base::duplicated(ol.dt, by = "name") |
    #base::duplicated(ol.dt, by = "name", fromLast = TRUE)
    #), ]
    
    # reads mapped to genes
    reads.mapped.to.gene <- nrow(count.dt[!grepl("ERCC", count.dt[, gene_id]), ])
    
    # umi filtering
    count.gene <- base::table(unique(count.dt[, .(gene_id, umi)])[, gene_id])
    
    # clean up
    count.gene.dt <- data.table::data.table(
      gene_id = c(names(features.gene),
                  "reads_mapped_to_genome",
                  "reads_mapped_to_genes"))
    cellname <- remove.last.extension(cell)
    count.gene.dt[[cellname]] <- 0
    count.gene.dt[gene_id == "reads_mapped_to_genome",
                  cellname] <- reads.mapped.to.genome
    count.gene.dt[gene_id == "reads_mapped_to_genes",
                  cellname] <- reads.mapped.to.gene
    count.gene.dt[gene_id %in% names(count.gene),
                  eval(cellname) := as.numeric(count.gene[gene_id])]
    
    # coerce to data frame to keep rownames for cbind combination
    count.gene.dt <- data.frame(count.gene.dt,
                                row.names = 1,
                                check.names = FALSE,
                                fix.empty.names = FALSE)
  }
  return (count.gene.dt)
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
#     txid = base::names(features.tx)[S4Vectors::queryHits(ol)],
#     name = base::names(bamGA)[S4Vectors::subjectHits(ol)],
#     pos = BiocGenerics::start(bamGA)[S4Vectors::subjectHits(ol)],
#     gene.id = tx.to.gene[, gene.id][S4Vectors::queryHits(ol)]
#   )
#   
#   dt <- merge(ol.dt, tx.to.gene[, .(txid, gene.id)], by = "txid")
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

# no iteration raw joint probabilities
# multimapper.probs <- function(multi.gene.alignments,
#                               ex.dt,
#                               bincounts,
#                               gene.probs.dt) {
#   counter <- 1
#   cutpoints <- seq(0, 1, 0.01)
#   
#   multi.probs <- (bincounts + 1) / sum(bincounts + 1)
#   
#   for (i in unique(multi.gene.alignments[, readname])) {
#     read.alignments <- multi.gene.alignments[readname == i, ]
#     
#     for (j in seq_len(nrow(read.alignments))) {
#       # positional probability
#       rp <- getrp(read.alignments[j, ], ex.dt)
#       mpb <- cut(rp,
#                  cutpoints,
#                  labels = FALSE,
#                  include.lowest = TRUE)
#       
#       # k = 100, 
#       x <- rep(0, 100)
#       x[mpb] <- 1
#       
#       pp <- dmultinom(x, prob = multi.probs)
#       
#       read.alignments[j, relative_position := rp]
#       read.alignments[j, positional_probability := pp]
#       
#       # abundance probability
#       read.alignments[j,
#                       abundance_probability := 
#                         gene.probs.dt[gene_id == read.alignments[j, gene_id],
#                                       abundance_probability]]
#       
#       # joint probability
#       read.alignments[j,
#                       joint := positional_probability *
#                         abundance_probability]
#     }
#     
#     multi.gene.alignments[readname == i,
#                           relative_position := 
#                             read.alignments[, relative_position]]
#     multi.gene.alignments[readname == i,
#                           positional_probability := 
#                             read.alignments[, positional_probability]]
#     multi.gene.alignments[readname == i,
#                           abundance_probability := 
#                             read.alignments[, abundance_probability]]
#     multi.gene.alignments[readname == i,
#                           joint := read.alignments[, joint]]
#     
#     counter <- counter + 1
#     if (counter %% 100 == 0) {
#       print(paste(counter, "reads processed"))
#     }
#     
#   }
#   return (multi.gene.alignments)
# }


