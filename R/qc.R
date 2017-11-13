

qc.plot <- function(qc.dt) {
  g1 <- plot.reads.assignment(qc.dt)
  g2 <- plot.total.reads(qc.dt)
  g3 <- plot.reads.mapped.to.genome(qc.dt)
  g4 <- plot.reads.mapped.to.genes(qc.dt)
  g5 <- plot.genome.reads.fraction(qc.dt)
  g6 <- plot.gene.to.genome.fraction(qc.dt)
  g7 <- plot.gene.to.total.fraction(qc.dt)
  g8 <- plot.transcripts(qc.dt)
  g9 <- plot.MT.transcripts(qc.dt)
  return (list(g1,
               gridExtra::marrangeGrob(
                 list(g2,g3,g4,g5,g6,g7,g8,g9),
                 ncol = 1,
                 nrow = 2,
                 top = NULL,
                 bottom = grid::textGrob("Cohort"))))
}


# assigned reads
plot.reads.assignment <- function(qc.dt) {
  
  plot.reads.assignment.id <- function (i, qc) {
    qc.i <- qc[!is.na(cell_num) & cohort == i, ]
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
    grobs = lapply(X = qc.dt[,unique(cohort)],
                   plot.reads.assignment.id,
                   qc = qc.dt),
    ncol = 1, nrow = 2,
    top = grid::textGrob("Reads per cell"),
    left = grid::textGrob("Reads", rot = 90),
    bottom = grid::textGrob("Cells in descending order")))
}


plot.total.reads <- function(qc.dt) {
  g <- ggplot2::ggplot(data = qc.dt[!(is.na(cell_num)), ],
                       ggplot2::aes(
                         x = as.factor(cohort),
                         y = log10(reads),
                         #y = reads,
                         group = as.factor(cohort)
                       )) +
    ggplot2::geom_boxplot(outlier.color = NA, fill = NA) +
    ggplot2::geom_point(
      color = "#424242",
      position = ggplot2::position_jitter(width = 0.3, height = 0),
      size = 1
    ) +
    ggplot2::ylab(expression(Log[10]*"reads")) +
    ggplot2::ggtitle("Total reads") +
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
                         x = as.factor(cohort),
                         y = log10(reads_mapped_to_genome),
                         group = as.factor(cohort)
                       )) +
    ggplot2::geom_boxplot(outlier.color = NA, fill = NA) +
    ggplot2::geom_point(
      color = "#424242",
      position = ggplot2::position_jitter(width = 0.3, height = 0),
      size = 1
    ) +
    ggplot2::ylab(expression(Log[10]*"reads")) +
    ggplot2::ggtitle("Reads mapped to genome") +
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
                       ggplot2::aes(x = as.factor(cohort),
                                    y = log10(reads_mapped_to_genes),
                                    group = as.factor(cohort))) +
    ggplot2::geom_boxplot(outlier.color = NA, fill = NA) +
    ggplot2::geom_point(color = "#424242",
               position = ggplot2::position_jitter(width = 0.3, height = 0),
               size = 1) +
    ggplot2::ylab(expression(Log[10]*"reads")) +
    ggplot2::ggtitle("Reads mapped to genes") +
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
                       ggplot2::aes(x = as.factor(cohort),
                                    y = reads_mapped_to_genome/reads,
                                    group = as.factor(cohort))) +
    ggplot2::geom_boxplot(outlier.color = NA, fill = NA) +
    ggplot2::geom_point(color = "#424242",
                        position = ggplot2::position_jitter(width = 0.3,
                                                            height = 0),
                        size = 1) +
    ggplot2::ylim(0, 1) +
    ggplot2::ggtitle("Fraction of mapped reads to total reads") +
    theme_Publication() + 
    ggplot2::theme(axis.title.y = ggplot2::element_blank(),
                   axis.title.x = ggplot2::element_blank())
    
  return (g)
}


plot.gene.to.genome.fraction <- function(qc.dt) {
  g <- ggplot2::ggplot(data = subset(qc.dt, !(cell %in% c("low_quality",
                                                          "total",
                                                          "undetermined"))),
                       ggplot2::aes(x = as.factor(cohort),
                                    y = reads_mapped_to_genes/reads_mapped_to_genome,
                                    group = as.factor(cohort))) +
    ggplot2::geom_boxplot(outlier.color = NA, fill = NA) +
    ggplot2::geom_point(color = "#424242",
                        position = ggplot2::position_jitter(width = 0.3,
                                                            height = 0),
                        size = 1) +
    ggplot2::ylim(0, 1) +
    ggplot2::ggtitle("Fraction of reads mapped to genes to reads mapped to genome") +
    theme_Publication() +
    ggplot2::theme(axis.title.y = ggplot2::element_blank(),
                   axis.title.x = ggplot2::element_blank())
  return (g)
}


plot.gene.to.total.fraction <- function(qc.dt) {
  g <- ggplot2::ggplot(data = subset(qc.dt, !(cell %in% c("low_quality",
                                                          "total",
                                                          "undetermined"))),
                       ggplot2::aes(x = as.factor(cohort),
                                    y = reads_mapped_to_genes/reads,
                                    group = as.factor(cohort))) +
    ggplot2::geom_boxplot(outlier.color = NA, fill = NA) +
    ggplot2::geom_point(color = "#424242",
                        position = ggplot2::position_jitter(width = 0.3,
                                                            height = 0),
                        size = 1) +
    ggplot2::ylim(0, 1) +
    ggplot2::ggtitle("Fraction of reads mapped to genes to total reads") +
    theme_Publication() +
    ggplot2::theme(axis.title.y = ggplot2::element_blank(),
                   axis.title.x = ggplot2::element_blank())
  return (g)
}


plot.transcripts <- function(qc.dt) {
  g <- ggplot2::ggplot(data = subset(qc.dt, !(cell %in% c("low_quality",
                                                          "total",
                                                          "undetermined"))),
                       ggplot2::aes(x = as.factor(cohort),
                                    y = log10(transcripts),
                                    group = as.factor(cohort))) +
    ggplot2::geom_boxplot(outlier.color = NA, fill = NA) +
    ggplot2::geom_point(color = "#424242",
                        position = ggplot2::position_jitter(width = 0.3,
                                                            height = 0),
                        size = 1) +
    ggplot2::ylab(expression(Log[10]*"transcripts")) +
    ggplot2::ggtitle("UMI filtered transcripts") +
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
                       ggplot2::aes(x = as.factor(cohort),
                                    y = log10(mt_transcripts),
                                    group = as.factor(cohort))) +
    ggplot2::geom_boxplot(outlier.color = NA, fill = NA) +
    ggplot2::geom_point(color = "#424242",
                        position = ggplot2::position_jitter(width = 0.3,
                                                            height = 0),
                        size = 1) +
    ggplot2::ylab(expression(Log[10]*"transcripts")) +
    ggplot2::ggtitle("Mitochondrial transcripts") +
    ggplot2::scale_y_continuous(labels = scales::comma,
                                limits = c(0, NA)) +
    theme_Publication() +
    ggplot2::theme(axis.title.x = ggplot2::element_blank())
  return (g)
}


