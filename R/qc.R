#' Collect QC metrics
#' 
#' Collect QC metrics from demultiplexing, alignment, and counting results. Return a \code{data.table} object containing reads, transcripts, and genes information.
#' 
#' @param de Demultiplex result. Table returned from \code{demultiplex} function.
#' @param al Alignment result. Table returned from \code{alignRsubread} function.
#' @param co Count matrix. Table returned from \code{countUmi} function.
#' @param geneAnnotation Gene information table generated from parsing GTF file. Must contain \emph{gene_id}, \emph{gene_biotype}, and \emph{seqid} columns.
#' @return QC metrics table
#' @export
collectQC <- function(de, al, co, geneAnnotation = NA) {
  de <- data.table::copy(data.table::data.table(de))
  al <- data.table::copy(data.table::data.table(al))
  
  colnames(al)[c(1, 3)] <- c("bam_dir",
                             "total_mapped_reads_incl_ercc")
  
  de[, cell := sub(pattern = "(.*?)\\..*$",
                   replacement = "\\1", filename)]
  al[, cell := .removeLastExtension(basename(bam_dir))]
  
  data.table::setkey(de, cell)
  data.table::setkey(al, cell)
  
  qcDt <- base::merge(de[, -"filename"],
                       al[, .(cell,
                              total_mapped_reads_incl_ercc)],
                       all.x=TRUE)
  
  # get reads mapped to genome
  rmtgenome <- co[geneid == "reads_mapped_to_genome", -"geneid"]
  rmtgenome <- data.table::data.table(
    cell = colnames(rmtgenome),
    reads_mapped_to_genome = as.numeric(rmtgenome))
  qcDt <- base::merge(qcDt, rmtgenome, all.x=TRUE)
  
  # reads mapped to genes
  rmtgene <- co[geneid == "reads_mapped_to_genes", -"geneid"]
  rmtgene <- data.table::data.table(cell = colnames(rmtgene),
                                    reads_mapped_to_genes = as.numeric(rmtgene))
  qcDt <- base::merge(qcDt, rmtgene, all.x=TRUE)
  
  # UMI filtered transcripts
  transcript <- base::colSums(co[!(geneid %in% c("reads_mapped_to_genome",
                                                  "reads_mapped_to_genes")) &
                                   !grepl("ERCC", co[,geneid]), -"geneid"])
  transcript <- data.table::data.table(
    cell = names(transcript),
    transcripts = as.numeric(transcript))
  qcDt <- base::merge(qcDt, transcript, all.x = TRUE)
  
  # MT transcript
  if (!all(is.na(geneAnnotation))) {
    mtTranscript <- base::colSums(co[geneid %in% 
                                        geneAnnotation[seqid == "MT",
                                                   gene_id],
                                      -"geneid"])
    mtTranscript <- data.table::data.table(
      cell = names(mtTranscript),
      mt_transcripts = as.numeric(mtTranscript))
    qcDt <- base::merge(qcDt, mtTranscript, all.x = TRUE)
  }
  
  # expressed genes
  cells <- colnames(co[, -"geneid"])
  genes <- vapply(cells, function(cells) nrow(
    co[!(geneid %in%
           c("reads_mapped_to_genome",
             "reads_mapped_to_genes")) &
         !grepl("ERCC", co[,geneid]) &
         eval(parse(text = paste0("`",
                                  cells,
                                  "`"))) != 0,
       cells, with = FALSE]), integer(1))
  
  genes <- data.table::data.table(
    cell = names(genes),
    genes = as.numeric(genes))
  qcDt <- base::merge(qcDt, genes, all.x = TRUE)
  
  # protein coding genes
  if (!all(is.na(geneAnnotation))) {
    proteinCodingGene <- geneAnnotation[gene_biotype == "protein_coding",
                                  gene_id]
    proGene <- vapply(cells, function(cells) nrow(
      co[geneid %in% proteinCodingGene &
           eval(parse(text = paste0("`", cells, "`"))) != 0,
         cells, with = FALSE]), integer(1))
    
    proGene <- data.table::data.table(
      cell = names(proGene),
      protein_coding_genes = as.numeric(proGene))
    qcDt <- base::merge(qcDt, proGene, all.x = TRUE)
    
    # protein coding transcripts
    proTranscript <- base::colSums(co[geneid %in% proteinCodingGene,
                                       cells, with = FALSE])
    proTranscript <- data.table::data.table(
      cell = names(proTranscript),
      protein_coding_transcripts = as.numeric(proTranscript))
    qcDt <- base::merge(qcDt, proTranscript, all.x = TRUE)
  }
  
  return (qcDt)
}


#' Visualize data quality
#' 
#' Visualize data quality from QC metrics table and return a list of \code{grobs} objects.
#' 
#' @param qcDt An QC metrics table returnd from \code{collect.qc} function.
#' @return A list of \code{grobs} objects
#' @export
qcplot <- function(qcDt) {
  g1 <- plot.reads.assignment(qcDt)
  g2 <- plot.total.reads(qcDt)
  g3 <- plot.reads.mapped.to.genome(qcDt)
  g4 <- plot.reads.mapped.to.genes(qcDt)
  g5 <- plot.genome.reads.fraction(qcDt)
  g6 <- plot.gene.to.genome.fraction(qcDt)
  g7 <- plot.gene.to.total.fraction(qcDt)
  g8 <- plot.transcripts(qcDt)
  g9 <- plot.MT.transcripts(qcDt)
  g10 <- plot.MT.transcripts.fraction(qcDt)
  g11 <- plot.genes(qcDt)
  g12 <- plot.frac.protein.coding.genes(qcDt)
  g13 <- plot.frac.protein.coding.transcripts(qcDt)
  g14 <- plot.genes.per.million.reads(qcDt)
  return (list(g1,
               gridExtra::marrangeGrob(
                 list(g2,g3,g4,g5,g6,g7,g8,g9,g10,g11,g12,g13,g14),
                 ncol = 1,
                 nrow = 2,
                 top = NULL,
                 bottom = grid::textGrob("sample"))))
}


#' ggplot publication theme
#' Modify theme of ggplots. Make ggplots look better. Adapted from Koundinya Desiraju. https://rpubs.com/Koundy/71792
#' @param base_size Default 12.
#' @param base_family Default \emph{sans}.
#' @export
theme_Publication <- function(base_size = 12,
                              base_family = "sans") {
  (ggthemes::theme_foundation(base_size = base_size,
                              base_family = base_family) +
     ggplot2::theme(plot.title = ggplot2::element_text(
       face = "bold",
       size = ggplot2::rel(1),
       hjust = 0.5),
                    text = ggplot2::element_text(),
                    panel.background = ggplot2::element_rect(color = NA),
                    plot.background = ggplot2::element_rect(color = NA),
                    panel.border = ggplot2::element_rect(color = NA),
                    axis.title = ggplot2::element_text(
                      face = "bold",
                      size = ggplot2::rel(1)),
                    axis.title.y = ggplot2::element_text(angle = 90,
                                                         vjust = 2),
                    axis.title.x = ggplot2::element_text(vjust = -0.2),
                    axis.text = ggplot2::element_text(),
                    axis.line = ggplot2::element_line(color="black"),
                    axis.ticks = ggplot2::element_line(),
                    panel.grid.major = ggplot2::element_line(color = "#f0f0f0"),
                    panel.grid.minor = ggplot2::element_blank(),
                    legend.key = ggplot2::element_rect(color = NA),
                    legend.position = "right",
                    legend.direction = "vertical",
                    legend.key.size= ggplot2::unit(0.2, "cm"),
                    legend.margin = ggplot2::margin(0),
                    legend.title = ggplot2::element_text(face = "bold"),
                    plot.margin = ggplot2::unit(c(10, 5, 5, 5), "mm"),
                    strip.background = ggplot2::element_rect(
                      color = "#f0f0f0", fill = "#f0f0f0"),
                    strip.text = ggplot2::element_text(face = "bold")
   ))
}


# cumulative fraction of reads per cell


# assigned reads
plot.reads.assignment <- function(qcDt) {
  
  plot.reads.assignment.id <- function (i, qc) {
    qc.i <- qc[!is.na(cell_num) & sample == i, ]
    qc.i <- qc.i[order(qc.i$reads, decreasing = TRUE), ]
    ggplot2::ggplot(qc.i) +
      ggplot2::geom_bar(ggplot2::aes(x = seq(nrow(qc.i)),
                                     y = reads),
                        stat="identity") +
      ggplot2::ggtitle(i) +
      ggplot2::scale_y_continuous(labels = scales::comma) +
      theme_Publication() +
      ggplot2::theme(axis.title.y = ggplot2::element_blank(),
                     axis.title.x = ggplot2::element_blank())
  }
  
  qc <- data.table::copy(qcDt)
  
  return (gridExtra::marrangeGrob(
    grobs = lapply(X = qcDt[,unique(sample)],
                   plot.reads.assignment.id,
                   qc = qcDt),
    ncol = 1, nrow = 2,
    top = grid::textGrob("Reads per cell"),
    left = grid::textGrob("Reads", rot = 90),
    bottom = grid::textGrob("Cell indices in descending order")))
}


plot.total.reads <- function(qcDt) {
  g <- ggplot2::ggplot(data = qcDt[!(is.na(cell_num)), ],
                       ggplot2::aes(
                         x = as.factor(sample),
                         y = log10(reads),
                         #y = reads,
                         group = as.factor(sample)
                       )) +
    ggplot2::geom_boxplot(outlier.color = NA, fill = NA) +
    ggplot2::geom_point(
      color = "#424242",
      position = ggplot2::position_jitter(width = 0.3, height = 0),
      size = 0.5
    ) +
    ggplot2::ylab(expression(Log[10]*"reads")) +
    ggplot2::ggtitle("Number of total reads") +
    ggplot2::scale_y_continuous(labels = scales::comma,
                                limits = c(0, NA)) +
    theme_Publication() +
    ggplot2::theme(axis.title.x = ggplot2::element_blank())
  return (g)
}


plot.reads.mapped.to.genome <- function(qcDt) {
  g <- ggplot2::ggplot(data = subset(qcDt, !(cell %in% c("low_quality",
                                                          "total",
                                                          "undetermined"))),
                       ggplot2::aes(
                         x = as.factor(sample),
                         y = log10(reads_mapped_to_genome),
                         group = as.factor(sample)
                       )) +
    ggplot2::geom_boxplot(outlier.color = NA, fill = NA) +
    ggplot2::geom_point(
      color = "#424242",
      position = ggplot2::position_jitter(width = 0.3, height = 0),
      size = 0.5
    ) +
    ggplot2::ylab(expression(Log[10]*"reads")) +
    ggplot2::ggtitle("Number of aligned reads") +
    ggplot2::scale_y_continuous(labels = scales::comma,
                                limits = c(0, NA)) +
    theme_Publication() +
    ggplot2::theme(axis.title.x = ggplot2::element_blank())
    
  return (g)
}


plot.reads.mapped.to.genes <- function(qcDt) {
  g <- ggplot2::ggplot(data = subset(qcDt, !(cell %in% c("low_quality",
                                                          "total",
                                                          "undetermined"))),
                       ggplot2::aes(x = as.factor(sample),
                                    y = log10(reads_mapped_to_genes),
                                    group = as.factor(sample))) +
    ggplot2::geom_boxplot(outlier.color = NA, fill = NA) +
    ggplot2::geom_point(color = "#424242",
               position = ggplot2::position_jitter(width = 0.3, height = 0),
               size = 0.5) +
    ggplot2::ylab(expression(Log[10]*"reads")) +
    ggplot2::ggtitle("Number of reads aligned to an gene") +
    ggplot2::scale_y_continuous(labels = scales::comma,
                                limits = c(0, NA)) +
    theme_Publication() +
    ggplot2::theme(axis.title.x = ggplot2::element_blank())
    
  return (g)
}


plot.genome.reads.fraction <- function(qcDt) {
  g <- ggplot2::ggplot(data = subset(qcDt, !(cell %in% c("low_quality",
                                                          "total",
                                                          "undetermined"))),
                       ggplot2::aes(x = as.factor(sample),
                                    y = reads_mapped_to_genome/reads,
                                    group = as.factor(sample))) +
    ggplot2::geom_boxplot(outlier.color = NA, fill = NA) +
    ggplot2::geom_point(color = "#424242",
                        position = ggplot2::position_jitter(width = 0.3,
                                                            height = 0),
                        size = 0.5) +
    ggplot2::ylim(0, 1) +
    ggplot2::ggtitle("Fraction of aligned reads to total reads") +
    theme_Publication() + 
    ggplot2::theme(axis.title.y = ggplot2::element_blank(),
                   axis.title.x = ggplot2::element_blank())
    
  return (g)
}


plot.gene.to.genome.fraction <- function(qcDt) {
  g <- ggplot2::ggplot(data = subset(qcDt, !(cell %in% c("low_quality",
                                                          "total",
                                                          "undetermined"))),
                       ggplot2::aes(x = as.factor(sample),
                                    y = reads_mapped_to_genes/reads_mapped_to_genome,
                                    group = as.factor(sample))) +
    ggplot2::geom_boxplot(outlier.color = NA, fill = NA) +
    ggplot2::geom_point(color = "#424242",
                        position = ggplot2::position_jitter(width = 0.3,
                                                            height = 0),
                        size = 0.5) +
    ggplot2::ylim(0, 1) +
    ggplot2::ggtitle("Fraction of reads aligned to an gene out of total number of aligned reads") +
    theme_Publication() +
    ggplot2::theme(axis.title.y = ggplot2::element_blank(),
                   axis.title.x = ggplot2::element_blank())
  return (g)
}


plot.gene.to.total.fraction <- function(qcDt) {
  g <- ggplot2::ggplot(data = subset(qcDt, !(cell %in% c("low_quality",
                                                          "total",
                                                          "undetermined"))),
                       ggplot2::aes(x = as.factor(sample),
                                    y = reads_mapped_to_genes/reads,
                                    group = as.factor(sample))) +
    ggplot2::geom_boxplot(outlier.color = NA, fill = NA) +
    ggplot2::geom_point(color = "#424242",
                        position = ggplot2::position_jitter(width = 0.3,
                                                            height = 0),
                        size = 0.5) +
    ggplot2::ylim(0, 1) +
    ggplot2::ggtitle("Fraction of reads aligned to an gene out of total number of reads") +
    theme_Publication() +
    ggplot2::theme(axis.title.y = ggplot2::element_blank(),
                   axis.title.x = ggplot2::element_blank())
  return (g)
}


plot.transcripts <- function(qcDt) {
  g <- ggplot2::ggplot(data = subset(qcDt, !(cell %in% c("low_quality",
                                                          "total",
                                                          "undetermined"))),
                       ggplot2::aes(x = as.factor(sample),
                                    y = log10(transcripts),
                                    group = as.factor(sample))) +
    ggplot2::geom_boxplot(outlier.color = NA, fill = NA) +
    ggplot2::geom_point(color = "#424242",
                        position = ggplot2::position_jitter(width = 0.3,
                                                            height = 0),
                        size = 0.5) +
    ggplot2::ylab(expression(Log[10]*"transcripts")) +
    ggplot2::ggtitle("Total number of transcripts after UMI correction") +
    ggplot2::scale_y_continuous(labels = scales::comma,
                                limits = c(0, NA)) +
    theme_Publication() +
    ggplot2::theme(axis.title.x = ggplot2::element_blank())
  return (g)
}


plot.MT.transcripts <- function(qcDt) {
  g <- ggplot2::ggplot(data = subset(qcDt, !(cell %in% c("low_quality",
                                                          "total",
                                                          "undetermined"))),
                       ggplot2::aes(x = as.factor(sample),
                                    y = log10(mt_transcripts),
                                    group = as.factor(sample))) +
    ggplot2::geom_boxplot(outlier.color = NA, fill = NA) +
    ggplot2::geom_point(color = "#424242",
                        position = ggplot2::position_jitter(width = 0.3,
                                                            height = 0),
                        size = 0.5) +
    ggplot2::ylab(expression(Log[10]*"transcripts")) +
    ggplot2::ggtitle("Number of mitochondrial transcripts after UMI correction") +
    ggplot2::scale_y_continuous(labels = scales::comma,
                                limits = c(0, NA)) +
    theme_Publication() +
    ggplot2::theme(axis.title.x = ggplot2::element_blank())
  return (g)
}


plot.MT.transcripts.fraction <- function(qcDt) {
  g <- ggplot2::ggplot(data = subset(qcDt, !(cell %in% c("low_quality",
                                                          "total",
                                                          "undetermined"))),
                       ggplot2::aes(x = as.factor(sample),
                                    y = mt_transcripts/transcripts,
                                    group = as.factor(sample))) +
    ggplot2::geom_boxplot(outlier.color = NA, fill = NA) +
    ggplot2::geom_point(color = "#424242",
                        position = ggplot2::position_jitter(width = 0.3,
                                                            height = 0),
                        size = 0.5) +
    ggplot2::ylim(0, 1) +
    ggplot2::ggtitle("Fraction of mitochondrial transcripts after UMI correction") +
    theme_Publication() +
    ggplot2::theme(axis.title.y = ggplot2::element_blank(),
                   axis.title.x = ggplot2::element_blank())
  return (g)
}


plot.genes <- function(qcDt) {
  g <- ggplot2::ggplot(data = subset(qcDt, !(cell %in% c("low_quality",
                                                          "total",
                                                          "undetermined"))),
                       ggplot2::aes(x = as.factor(sample),
                                    y = log10(genes),
                                    group = as.factor(sample))) +
    ggplot2::geom_boxplot(outlier.color = NA, fill = NA) +
    ggplot2::geom_point(color = "#424242",
                        position = ggplot2::position_jitter(width = 0.3,
                                                            height = 0),
                        size = 0.5) +
    ggplot2::ylab(expression(Log[10]*"Genes")) +
    ggplot2::ggtitle("Number of genes") +
    ggplot2::scale_y_continuous(labels = scales::comma,
                                limits = c(0, NA)) +
    theme_Publication() +
    ggplot2::theme(axis.title.x = ggplot2::element_blank())
  return (g)
}


plot.frac.protein.coding.genes <- function(qcDt) {
  g <- ggplot2::ggplot(data = subset(qcDt, !(cell %in% c("low_quality",
                                                          "total",
                                                          "undetermined"))),
                       ggplot2::aes(x = as.factor(sample),
                                    y = protein_coding_genes/genes,
                                    group = as.factor(sample))) +
    ggplot2::geom_boxplot(outlier.color = NA, fill = NA) +
    ggplot2::geom_point(color = "#424242",
                        position = ggplot2::position_jitter(width = 0.3,
                                                            height = 0),
                        size = 0.5) +
    ggplot2::ylim(0, 1) +
    ggplot2::ggtitle("Fraction of protein coding genes") +
    theme_Publication() +
    ggplot2::theme(axis.title.y = ggplot2::element_blank(),
                   axis.title.x = ggplot2::element_blank())
  return (g)
}


plot.frac.protein.coding.transcripts <- function(qcDt) {
  g <- ggplot2::ggplot(data = subset(qcDt, !(cell %in% c("low_quality",
                                                          "total",
                                                          "undetermined"))),
                       ggplot2::aes(x = as.factor(sample),
                                    y = protein_coding_transcripts/transcripts,
                                    group = as.factor(sample))) +
    ggplot2::geom_boxplot(outlier.color = NA, fill = NA) +
    ggplot2::geom_point(color = "#424242",
                        position = ggplot2::position_jitter(width = 0.3,
                                                            height = 0),
                        size = 0.5) +
    ggplot2::ylim(0, 1) +
    ggplot2::ggtitle("Fraction of protein coding transcripts after UMI correction") +
    theme_Publication() +
    ggplot2::theme(axis.title.y = ggplot2::element_blank(),
                   axis.title.x = ggplot2::element_blank())
  return (g)
}


plot.genes.per.million.reads <- function(qcDt) {
  g <- ggplot2::ggplot(data = subset(qcDt, !(cell %in% c("low_quality",
                                                          "total",
                                                          "undetermined"))),
                       ggplot2::aes(x = as.factor(sample),
                                    y = log10(genes %x% 1000000/reads),
                                    group = as.factor(sample))) +
    ggplot2::geom_boxplot(outlier.color = NA, fill = NA) +
    ggplot2::geom_point(color = "#424242",
                        position = ggplot2::position_jitter(width = 0.3,
                                                            height = 0),
                        size = 0.5) +
    ggplot2::ylab(expression(paste(Log[10],
                                   "(Genes x 1000000 / total reads)"))) +
    ggplot2::ggtitle("Genes detected divided by total number of reads sequenced per million") +
    ggplot2::scale_y_continuous(labels = scales::comma,
                                limits = c(0, NA)) +
    theme_Publication() +
    ggplot2::theme(axis.title.x = ggplot2::element_blank())
  return (g)
}

