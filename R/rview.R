#' Visualize aligned reads
#'
#' Visualize read alignments for UMI tagged single cell RNA-sequencing data.
#'  Read names must contain UMI sequences at the end delimited by "\strong{:}".
#'  Arrow represents orientation of alignment. Reads are colored by their UMI
#'  and sorted by their start positions and UMI.
#'
#' @param bamGA A GenomicAlignment object
#' @param chr Chromosome. Integer or "X", "Y", "MT".
#' @param start Genomic coordinate of the start position.
#' @param end Genomic coordinate of the end position.
#' @param legend Show legend. Default is FALSE.
#' @return A ggplot object of aligned reads
#' @examples
#' data(bamExample, package = "scruff")
#' g <- rview(bamExample, chr = "MT", legend = TRUE)
#' g
#' @import ggbio
#' @export
rview <- function(bamGA,
    chr = "1",
    start = 1,
    end = max(BiocGenerics::end(bamGA)),
    legend = FALSE) {

    reads <- bamGA[BiocGenerics::start(bamGA) >= start &
            BiocGenerics::end(bamGA) <= end &
            GenomeInfoDb::seqnames(bamGA) == chr]
    umi <- data.table::last(data.table::tstrsplit(names(reads), ":"))
    reads <- reads[order(BiocGenerics::start(reads), umi)]

    readsGr <- GenomicRanges::GRanges(reads)
    S4Vectors::mcols(readsGr)$umi <- data.table::last(
        data.table::tstrsplit(names(readsGr), ":"))

    g <- ggplot2::ggplot(readsGr) +
        ggbio::geom_arrow(ggplot2::aes(color = umi)) +
        .themePublication() +
        ggplot2::theme(axis.title.y = ggplot2::element_blank())

    if (legend == FALSE) {
        g <- g + ggplot2::theme(legend.position = "none")
    }
    return(g)
}
