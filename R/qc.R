
# assigned reads
plot.reads.assignment <- function(qc.dt) {
  
  plot.reads.assignment.id <- function (i, qc) {
    qc.i <- qc[!is.na(cell_num) & id == i, ]
    qc.i <- qc.i[order(qc.i$reads, decreasing = TRUE), ]
    ggplot2::ggplot(qc.i) +
                ggplot2::geom_bar(ggplot2::aes(x = seq(nrow(qc.i)), y = reads),
                                  stat="identity") +
                ggplot2::xlab("cell") +
                ggplot2::ylab("Reads") +
                ggplot2::scale_y_continuous(labels = scales::comma) +
                theme_Publication()
  }
  
  qc <- data.table::copy(qc.dt)
  
  return (gridExtra::arrangeGrob(grobs = lapply(X = qc.dt[,unique(id)],
                                                plot.reads.assignment.id,
                                                qc = qc.dt),
                                 ncol = 2,
                                 top = grid::textGrob("Reads per cell")))
}


plot.total.reads <- function(qc.dt) {
  g <- ggplot2::ggplot(data = qc.dt[!(is.na(cell_num)), ],
                       ggplot2::aes(
                         x = as.factor(id),
                         y = reads,
                         group = as.factor(id)
                       )) +
    ggplot2::geom_boxplot(outlier.color = NA, fill = NA) +
    ggplot2::geom_point(
      color = "#424242",
      position = ggplot2::position_jitter(width = 0.3, height = 0),
      size = 1
    ) +
    ggplot2::xlab("Sample ID") +
    ggplot2::ggtitle("Total reads") +
    ggplot2::scale_y_continuous(labels = scales::comma) +
    theme_Publication() +
    ggplot2::theme(axis.title.y = ggplot2::element_blank())
  return (g)
}


plot.reads.mapped.to.genome <- function(qc.dt) {
  g <- ggplot2::ggplot(data = subset(qc.dt, !(cell %in% c("low_quality",
                                                          "total",
                                                          "undetermined"))),
                       ggplot2::aes(
                         x = as.factor(id),
                         y = reads_mapped_to_genome,
                         group = as.factor(id)
                       )) +
    ggplot2::geom_boxplot(outlier.color = NA, fill = NA) +
    ggplot2::geom_point(
      color = "#424242",
      position = ggplot2::position_jitter(width = 0.3, height = 0),
      size = 1
    ) +
    ggplot2::xlab("Sample ID") +
    ggplot2::ggtitle("Reads mapped to genome") +
    ggplot2::scale_y_continuous(labels = scales::comma) +
    theme_Publication() +
    ggplot2::theme(axis.title.y = ggplot2::element_blank())
    
  return (g)
}


plot.reads.mapped.to.genes <- function(qc.dt) {
  g <- ggplot2::ggplot(data = subset(qc.dt, !(cell %in% c("low_quality",
                                                          "total",
                                                          "undetermined"))),
                       ggplot2::aes(x = as.factor(id), y = reads_mapped_to_genes,
                                    group = as.factor(id))) +
    ggplot2::geom_boxplot(outlier.color = NA, fill = NA) +
    ggplot2::geom_point(color = "#424242",
               position = ggplot2::position_jitter(width = 0.3, height = 0),
               size = 1) +
    ggplot2::xlab("Sample ID") +
    ggplot2::ggtitle("Reads mapped to genes") +
    ggplot2::scale_y_continuous(labels = scales::comma) +
    theme_Publication() +
    ggplot2::theme(axis.title.y = ggplot2::element_blank())
    
  return (g)
}


plot.genome.reads.fraction <- function(qc.dt) {
  g <- ggplot2::ggplot(data = subset(qc.dt, !(cell %in% c("low_quality",
                                                          "total",
                                                          "undetermined"))),
                       ggplot2::aes(x = as.factor(id),
                                    y = reads_mapped_to_genome/reads,
                                    group = as.factor(id))) +
    ggplot2::geom_boxplot(outlier.color = NA, fill = NA) +
    ggplot2::geom_point(color = "#424242",
                        position = ggplot2::position_jitter(width = 0.3,
                                                            height = 0),
                        size = 1) +
    ggplot2::xlab("Sample ID") +
    ggplot2::ylim(0, 1) +
    ggplot2::ggtitle("Fraction of reads mapped to genome") +
    theme_Publication() + 
    ggplot2::theme(axis.title.y = ggplot2::element_blank())
    
  return (g)
}


plot.gene.to.genome.fraction <- function(qc.dt) {
  g <- ggplot2::ggplot(data = subset(qc.dt, !(cell %in% c("low_quality",
                                                          "total",
                                                          "undetermined"))),
                       ggplot2::aes(x = as.factor(id),
                                    y = reads_mapped_to_genes/reads_mapped_to_genome,
                                    group = as.factor(id))) +
    ggplot2::geom_boxplot(outlier.color = NA, fill = NA) +
    ggplot2::geom_point(color = "#424242",
                        position = ggplot2::position_jitter(width = 0.3,
                                                            height = 0),
                        size = 1) +
    ggplot2::xlab("Sample ID") +
    ggplot2::ylim(0, 1) +
    ggplot2::ggtitle("Fraction of gene reads to genome reads") +
    theme_Publication() +
    ggplot2::theme(axis.title.y = ggplot2::element_blank())
  return (g)
}


plot.gene.to.total.fraction <- function(qc.dt) {
  g <- ggplot2::ggplot(data = subset(qc.dt, !(cell %in% c("low_quality",
                                                          "total",
                                                          "undetermined"))),
                       ggplot2::aes(x = as.factor(id),
                                    y = reads_mapped_to_genes/reads,
                                    group = as.factor(id))) +
    ggplot2::geom_boxplot(outlier.color = NA, fill = NA) +
    ggplot2::geom_point(color = "#424242",
                        position = ggplot2::position_jitter(width = 0.3,
                                                            height = 0),
                        size = 1) +
    ggplot2::xlab("Sample ID") +
    ggplot2::ylim(0, 1) +
    ggplot2::ggtitle("Fraction of reads mapped to genes relative to total reads") +
    theme_Publication() +
    ggplot2::theme(axis.title.y = ggplot2::element_blank())
  return (g)
}


plot.transcripts <- function(qc.dt) {
  g <- ggplot2::ggplot(data = subset(qc.dt, !(cell %in% c("low_quality",
                                                          "total",
                                                          "undetermined"))),
                       ggplot2::aes(x = as.factor(id),
                                    y = transcript,
                                    group = as.factor(id))) +
    ggplot2::geom_boxplot(outlier.color = NA, fill = NA) +
    ggplot2::geom_point(color = "#424242",
                        position = ggplot2::position_jitter(width = 0.3,
                                                            height = 0),
                        size = 1) +
    ggplot2::xlab("Sample ID") +
    ggplot2::ggtitle("UMI filtered transcripts") +
    theme_Publication() +
    ggplot2::theme(axis.title.y = ggplot2::element_blank())
  return (g)
}


plot.MT.transcripts <- function() {
  
}


