#' Visualize gene isoforms
#'
#' Visualize reference genome. Rectangles represent exons. Arrow represents
#'  orientation of transcripts.
#'
#' @param gtfFile A genome annotation file in GTF format.
#' @param chr Chromosome name. Integer or "X", "Y", "MT".
#' @param start Genomic coordinate of the start position.
#' @param end Genomic coordinate of the end position.
#' @param rect_width Exon widths. Default 0.3.
#' @param line_width Line weight. Default 0.5.
#' @param arrow_segments The number of segments lines be divided to. The
#'  greater the number, more arrows there are. Default 10.
#' @param arrow_width The angle of the arrow head in degrees (smaller numbers
#'  produce narrower, pointier arrows). Essentially describes the width of the
#'  arrow head. Passed to the angle parameter of arrow function. Default 30.
#' @param arrow_length The length of the arrow head. Passed to the length
#'  argument of arrow function. Default 0.08.
#' @param arrow_type One of "open" or "closed" indicating whether the arrow
#'  head should be a closed triangle. Passed to the type argument of arrow
#'  function. Default "open".
#' @param text_size Size of text. Passed to the size argument of the geom_text
#'  function. Default 4.
#' @return A ggplot object of genomic view
#' @examples
#' gtf <- system.file("extdata", "GRCm38_MT.gtf", package = "scruff")
#' g <- gview(gtf, chr = "MT")
#' g
#' @import rtracklayer
#' @export
gview <- function(gtfFile,
    chr = 1,
    start = 1,
    end = max(data.table::as.data.table(
        rtracklayer::import(gtfFile))[seqnames == chr, end]),
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
        return(txdt)
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

        exons <- dt[type == "exon", ]

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
                    exon_number =
                        exons[i, exon_number],
                    transcript_id =
                        exons[i, transcript_id],
                    gene_id = exons[i, gene_id],
                    gene_name =
                        exons[i, gene_name]),
                fill = TRUE)
        }
        return(rdt)
    }


    .transArrow <- function(dt) {
        adt <- data.table::data.table()
        for (i in seq_len(nrow(dt))) {
            mi <- dt[i, start]
            ma <- dt[i, end]

            if (dt[i, strand] == "+") {
                x1 <- mi + (((ma - mi) / arrow_segments) *
                        seq(0, arrow_segments - 1))
                x2 <- ma - (((ma - mi) / arrow_segments) *
                        seq(arrow_segments - 1, 0))
            } else if (dt[i, strand] == "-") {
                x1 <- ma + (((mi - ma) / arrow_segments) *
                        seq(0, arrow_segments - 1))
                x2 <- mi - (((mi - ma) / arrow_segments) *
                        seq(arrow_segments - 1, 0))
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

        return(adt)
    }


    .transText <- function(dt) {
        tdt <- data.table::data.table()

        for (i in seq_len(nrow(dt))) {
            mi <- dt[i, start]
            ma <- dt[i, end]

            x <- (ma + mi) / 2
            y <- dt[i, level] + 0.4

            tdt <- rbind(tdt,
                data.table::data.table(
                    x = x,
                    y = y,
                    transcript_name = dt[i, transcript_name]))

        }
        return(tdt)
    }


    # convert to data.table
    gtfDt <- data.table::as.data.table(
        rtracklayer::import(gtfFile))[type != "gene", c("seqnames",
            "type",
            "start",
            "end",
            "strand",
            "gene_biotype",
            "gene_name",
            "exon_number",
            "gene_id",
            "transcript_name",
            "transcript_id")]

    # use new variables to avoid ambiguity
    begin <- start
    ed <- end

    # subset features
    gtfDt <- gtfDt[end >= begin & start <= ed & seqnames == chr, ]

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
                    arrow_length,
                    "inches"),
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

    return(g)
}
