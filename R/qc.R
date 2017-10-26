
#percentage of assigned reads
plot.reads.assignment <- function(demultiplex.res) {
  for (i in demultiplex.res$id) {
    demultiplex.res[id == i,
                    cum_per := 100 - cumsum(demultiplex.res[id == i,
                                                            percentage])]
  }

  g <- ggplot2::ggplot() +
    ggplot2::geom_bar(
      data = demultiplex.res[cell_fname != "total",],
      ggplot2::aes(x = as.factor(id), y = percentage, group = cell_fname),
      stat = "identity",
      fill = "white",
      color = "black"
    ) +
    ggplot2::geom_text(
      data = demultiplex.res[cell_fname %in% c("undetermined",
                                               "low_quality"),],
      ggplot2::aes(as.factor(id),
                   cum_per,
                   group = cell_fname,
                   label = cell_fname),
      size = 4,
      position = ggplot2::position_nudge(x = 0, y = 2.5)
    ) +
    ggplot2::geom_text(
      data = demultiplex.res[cell_fname == "total",],
      ggplot2::aes(
        as.factor(id),
        percentage,
        group = cell_fname,
        label = reads
      ),
      size = 4,
      position = ggplot2::position_nudge(x = 0, y = 2.5)
    ) +
    ggplot2::xlab("Sample ID") +
    ggplot2::ylab("Percentage (%)") +
    ggplot2::ggtitle("Percentage of reads per cell for each sample ID") +
    ggplot2::scale_y_continuous(labels = scales::comma) +
    ggplot2::theme(plot.title = ggplot2::element_text("Total reads"),
                   legend.position = "none") +
    theme_Publication()
  return (g)
}


plot.total.reads <- function(demultiplex.res) {
  g <- ggplot2::ggplot(data = demultiplex.res[!(is.na(cell_num)), ],
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
    ggplot2::ylab("") +
    ggplot2::ggtitle("Total reads") +
    ggplot2::scale_y_continuous(labels = scales::comma) +
    ggplot2::theme(axis.title.y = element_blank()) +
    theme_Publication()
  return (g)
}


plot.reads.mapped.genome <- function(dt) {
  g <- ggplot2::ggplot(data = dt,
                       ggplot2::aes(
                         x = as.factor(id),
                         y = aligned.reads,
                         group = as.factor(id)
                       )) +
    ggplot2::geom_boxplot(outlier.color = NA, fill = NA) +
    ggplot2::geom_point(
      color = "darkgreen",
      position = ggplot2::position_jitter(width = 0.3, height = 0),
      size = 1
    ) +
    ggplot2::xlab("Sample ID") +
    ggplot2::ylab("") +
    ggplot2::ggtitle("Reads mapped to genome") +
    ggplot2::scale_y_continuous(labels = comma) +
    ggplot2::theme(axis.title.y = element_blank()) +
    ggplot2::theme_Publication()
  return (g)
}