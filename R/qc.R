
#' Visualize data quality
#' 
#' Visualize data quality from QC metrics table and return a list of \code{grobs} objects.
#' 
#' @param qc.dt An QC metrics table returnd from \code{collect.qc} function.
#' @return A list of \code{grobs} objects
#' @export
qcplot <- function(qc.dt) {
  g1 <- plot.reads.assignment(qc.dt)
  g2 <- plot.total.reads(qc.dt)
  g3 <- plot.reads.mapped.to.genome(qc.dt)
  g4 <- plot.reads.mapped.to.genes(qc.dt)
  g5 <- plot.genome.reads.fraction(qc.dt)
  g6 <- plot.gene.to.genome.fraction(qc.dt)
  g7 <- plot.gene.to.total.fraction(qc.dt)
  g8 <- plot.transcripts(qc.dt)
  g9 <- plot.MT.transcripts(qc.dt)
  g10 <- plot.MT.transcripts.fraction(qc.dt)
  g11 <- plot.genes(qc.dt)
  g12 <- plot.frac.protein.coding.genes(qc.dt)
  g13 <- plot.frac.protein.coding.transcripts(qc.dt)
  g14 <- plot.genes.per.million.reads(qc.dt)
  return (list(g1,
               gridExtra::marrangeGrob(
                 list(g2,g3,g4,g5,g6,g7,g8,g9,g10,g11,g12,g13,g14),
                 ncol = 1,
                 nrow = 2,
                 top = NULL,
                 bottom = grid::textGrob("sample"))))
}


# ggplot publication theme
# adapted from Koundinya Desiraju
# https://rpubs.com/Koundy/71792
theme_Publication <- function(base_size=12, base_family="sans") {
  (ggthemes::theme_foundation(base_size=base_size, base_family=base_family)
   + ggplot2::theme(plot.title = ggplot2::element_text(face = "bold",
                                                       size = ggplot2::rel(1), hjust = 0.5),
                    text = ggplot2::element_text(),
                    panel.background = ggplot2::element_rect(colour = NA),
                    plot.background = ggplot2::element_rect(colour = NA),
                    panel.border = ggplot2::element_rect(colour = NA),
                    axis.title = ggplot2::element_text(face = "bold",size = ggplot2::rel(1)),
                    axis.title.y = ggplot2::element_text(angle=90,vjust =2),
                    axis.title.x = ggplot2::element_text(vjust = -0.2),
                    axis.text = ggplot2::element_text(),
                    axis.line = ggplot2::element_line(colour="black"),
                    axis.ticks = ggplot2::element_line(),
                    panel.grid.major = ggplot2::element_line(colour="#f0f0f0"),
                    panel.grid.minor = ggplot2::element_blank(),
                    legend.key = ggplot2::element_rect(colour = NA),
                    legend.position = "right",
                    legend.direction = "vertical",
                    legend.key.size= ggplot2::unit(0.2, "cm"),
                    legend.margin = ggplot2::margin(0),
                    legend.title = ggplot2::element_text(face="bold"),
                    plot.margin = ggplot2::unit(c(10,5,5,5),"mm"),
                    strip.background = ggplot2::element_rect(colour="#f0f0f0",fill="#f0f0f0"),
                    strip.text = ggplot2::element_text(face="bold")
   ))
  
}


# assigned reads
plot.reads.assignment <- function(qc.dt) {
  
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
  
  qc <- data.table::copy(qc.dt)
  
  return (gridExtra::marrangeGrob(
    grobs = lapply(X = qc.dt[,unique(sample)],
                   plot.reads.assignment.id,
                   qc = qc.dt),
    ncol = 1, nrow = 2,
    top = grid::textGrob("Reads per cell"),
    left = grid::textGrob("Reads", rot = 90),
    bottom = grid::textGrob("Cell indices in descending order")))
}


plot.total.reads <- function(qc.dt) {
  g <- ggplot2::ggplot(data = qc.dt[!(is.na(cell_num)), ],
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


plot.reads.mapped.to.genome <- function(qc.dt) {
  g <- ggplot2::ggplot(data = subset(qc.dt, !(cell %in% c("low_quality",
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


plot.reads.mapped.to.genes <- function(qc.dt) {
  g <- ggplot2::ggplot(data = subset(qc.dt, !(cell %in% c("low_quality",
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


plot.genome.reads.fraction <- function(qc.dt) {
  g <- ggplot2::ggplot(data = subset(qc.dt, !(cell %in% c("low_quality",
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


plot.gene.to.genome.fraction <- function(qc.dt) {
  g <- ggplot2::ggplot(data = subset(qc.dt, !(cell %in% c("low_quality",
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


plot.gene.to.total.fraction <- function(qc.dt) {
  g <- ggplot2::ggplot(data = subset(qc.dt, !(cell %in% c("low_quality",
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


plot.transcripts <- function(qc.dt) {
  g <- ggplot2::ggplot(data = subset(qc.dt, !(cell %in% c("low_quality",
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


plot.MT.transcripts <- function(qc.dt) {
  g <- ggplot2::ggplot(data = subset(qc.dt, !(cell %in% c("low_quality",
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


plot.MT.transcripts.fraction <- function(qc.dt) {
  g <- ggplot2::ggplot(data = subset(qc.dt, !(cell %in% c("low_quality",
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


plot.genes <- function(qc.dt) {
  g <- ggplot2::ggplot(data = subset(qc.dt, !(cell %in% c("low_quality",
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


plot.frac.protein.coding.genes <- function(qc.dt) {
  g <- ggplot2::ggplot(data = subset(qc.dt, !(cell %in% c("low_quality",
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


plot.frac.protein.coding.transcripts <- function(qc.dt) {
  g <- ggplot2::ggplot(data = subset(qc.dt, !(cell %in% c("low_quality",
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


plot.genes.per.million.reads <- function(qc.dt) {
  g <- ggplot2::ggplot(data = subset(qc.dt, !(cell %in% c("low_quality",
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

