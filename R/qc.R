
plot.reads.assignment <- function(demultiplex.res) {
  for (i in demultiplex.res$id) {
    demultiplex.res[id == i,
                    cum_per := 100 - cumsum(demultiplex.res[id == i,
                                                            percentage])]
  }

  g <- ggplot() +
    geom_bar(
      data = demultiplex.res[cell_fname != "total", ],
      aes(x = as.factor(id), y = percentage, group = cell_fname),
      stat = "identity",
      fill = "white",
      color = "black"
    ) +
    geom_text(
      data = demultiplex.res[cell_fname %in% c("undetermined",
                                               "low_quality"), ],
      aes(as.factor(id),
          cum_per,
          group = cell_fname,
          label = cell_fname),
      size = 4,
      position = position_nudge(x = 0, y = 2.5)
    ) +
    geom_text(
      data = demultiplex.res[cell_fname == "total", ],
      aes(
        as.factor(id),
        percentage,
        group = cell_fname,
        label = paste(reads)
      ),
      size = 4,
      position = position_nudge(x = 0, y = 2.5)
    ) +
    xlab("Sample ID") + ylab("Percentage (%)") +
    ggtitle("Percentage of reads per cell for each sample ID") +
    scale_y_continuous(labels = comma) +
    theme(plot.title = element_text("Total reads"),
          legend.position = "none") +
    theme_Publication()
  return (g)
}

