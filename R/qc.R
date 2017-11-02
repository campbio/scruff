
#percentage of assigned reads
plot.reads.assignment <- function(qc.dt) {
  qc.dt <- data.table::copy(qc.dt)
  for (i in qc.dt$id) {
    qc.dt[id == i & cell != "total",
          cum_per := 100 - cumsum(qc.dt[id == i & cell != "total",
                                        percent_assigned])]
  }

  g <- ggplot2::ggplot() +
    ggplot2::geom_bar(
      data = qc.dt[cell != "total",],
      ggplot2::aes(x = as.factor(id), y = percent_assigned, group = cell),
      stat = "identity",
      fill = "white",
      color = "black"
    ) +
    ggplot2::geom_text(
      data = qc.dt[cell %in% c("undetermined",
                                   "low_quality"),],
      ggplot2::aes(as.factor(id),
                   cum_per,
                   group = cell,
                   label = cell),
      size = 4,
      position = ggplot2::position_nudge(x = 0, y = 2.5)
    ) +
    ggplot2::geom_text(
      data = qc.dt[cell == "total",],
      ggplot2::aes(
        as.factor(id),
        percent_assigned,
        group = cell,
        label = reads
      ),
      size = 4,
      position = ggplot2::position_nudge(x = 0, y = 2.5)
    ) +
    ggplot2::xlab("Sample ID") +
    ggplot2::ylab("Percentage (%)") +
    ggplot2::ggtitle("Percentage of reads per cell for each sample ID") +
    ggplot2::scale_y_continuous(labels = scales::comma) +
    ggplot2::coord_cartesian(ylim = c(0, 100)) +
    theme_Publication()
  return (g)
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
    ggplot2::theme(axis.title.y = element_blank())
  return (g)
}


plot.reads.mapped.to.genome <- function(qc.dt) {
  g <- ggplot2::ggplot(data = subset(qc.dt, !(cell %in% c("low_quality",
                                                          "total",
                                                          "undetermined"))),
                       ggplot2::aes(
                         x = as.factor(id),
                         y = mapped_reads,
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


plot.aligned.reads.fraction <- function(qc.dt) {
  g <- ggplot2::ggplot(data = subset(qc.dt, !(cell %in% c("low_quality",
                                                          "total",
                                                          "undetermined"))),
                       ggplot2::aes(x = as.factor(id), y = fraction_mapped,
                                    group = as.factor(id))) +
    ggplot2::geom_boxplot(outlier.color = NA, fill = NA) +
    ggplot2::geom_point(color = "#424242",
                        position = ggplot2::position_jitter(width = 0.3,
                                                            height = 0),
                        size = 1) +
    ggplot2::xlab("Sample ID") +
    ggplot2::ylim(0, 1) +
    ggplot2::ggtitle("Fraction of aligned reads") +
    theme_Publication() + 
    ggplot2::theme(axis.title.y = ggplot2::element_blank())
    
  return (g)
}


plot.reads.gene.to.aligned <- function(qc.dt) {
  g <- ggplot2::ggplot(data = subset(qc.dt, !(cell %in% c("low_quality",
                                                          "total",
                                                          "undetermined"))),
                       ggplot2::aes(x = as.factor(id),
                                    y = reads_mapped_to_genes/mapped_reads,
                                    group = as.factor(id))) +
    ggplot2::geom_boxplot(outlier.color = NA, fill = NA) +
    ggplot2::geom_point(color = "#424242",
                        position = ggplot2::position_jitter(width = 0.3,
                                                            height = 0),
                        size = 1) +
    ggplot2::xlab("Sample ID") +
    ggplot2::ylim(0, 1) +
    ggplot2::ggtitle("Fraction of reads mapped to genes relative to aligned reads") +
    theme_Publication() +
    ggplot2::theme(axis.title.y = ggplot2::element_blank())
  return (g)
}


plot.reads.gene.to.total <- function(qc.dt) {
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
    ggplot2::ggtitle("Unique transcript:UMI pairs") +
    theme_Publication() +
    ggplot2::theme(axis.title.y = ggplot2::element_blank())
  return (g)
}


