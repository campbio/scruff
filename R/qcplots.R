#' Visualize data quality
#' 
#' Visualize data quality from the \code{colData} of the
#'  \code{SingleCellExperiment} object and return a list of figures in
#'  \code{arrangelist} object.
#' 
#' @param sce An \code{SingleCellExperiment} object returned from \code{scruff}
#'  or \code{countUMI} function.
#' @return A list of \code{grobs} objects ready for plotting
#' @examples
#' data(sceExample, package = "scruff")
#' qcplots(sceExample)
#' @import grid
#' @export
qcplots <- function(sce) {
    qcDt <- data.table::as.data.table(SummarizedExperiment::colData(sce))
    
    g2 <- .plotTotalReads(qcDt)
    g3 <- .plotReadsMappedToGenome(qcDt)
    g4 <- .plotReadsMappedToGenes(qcDt)
    g5 <- .plotGenomeReadsFraction(qcDt)
    g6 <- .plotGeneToGenomeFraction(qcDt)
    g7 <- .plotGeneToTotalFraction(qcDt)
    g8 <- .plotCounts(qcDt)
    g9 <- .plotMtCounts(qcDt)
    g10 <- .plotMtCountsFraction(qcDt)
    g11 <- .plotGenes(qcDt)
    g12 <- .plotFracProteinCodingGenes(qcDt)
    g13 <- .plotFracProteinCodingTranscripts(qcDt)
    g14 <- .plotGenesPerMillionReads(qcDt)
    return (gridExtra::marrangeGrob(
        list(g2, g3, g4, g5, g6, g7, g8, g9, g10, g11, g12, g13, g14),
        ncol = 1,
        nrow = 2,
        top = NULL))
}


.plotTotalReads <- function(qcDt) {
    g <- ggplot2::ggplot(data = qcDt,
        ggplot2::aes(
            x = as.factor(experiment),
            y = log10(reads),
            group = as.factor(experiment))) +
        ggplot2::geom_boxplot(outlier.color = NA, fill = NA) +
        ggplot2::geom_point(
            ggplot2::aes(color = as.factor(number_of_cells)),
            position = ggplot2::position_jitter(width = 0.3, height = 0),
            size = 0.5
        ) +
        ggplot2::ylab(expression(bold(Log[10]*"Reads"))) +
        ggplot2::xlab("Experiment") +
        ggplot2::ggtitle("Total reads") +
        ggplot2::labs(color = "Cells/well") +
        ggplot2::scale_y_continuous(labels = scales::comma,
            limits = c(0, NA)) +
        .themePublication()
    return (g)
}


.plotReadsMappedToGenome <- function(qcDt) {
    g <- ggplot2::ggplot(data = qcDt,
        ggplot2::aes(
            x = as.factor(experiment),
            y = log10(reads_mapped_to_genome),
            group = as.factor(experiment))) +
        ggplot2::geom_boxplot(outlier.color = NA, fill = NA) +
        ggplot2::geom_point(
            ggplot2::aes(color = as.factor(number_of_cells)),
            position = ggplot2::position_jitter(width = 0.3, height = 0),
            size = 0.5
        ) +
        ggplot2::ylab(expression(bold(Log[10]*"Reads"))) +
        ggplot2::xlab("Experiment") +
        ggplot2::ggtitle("Reads aligned to reference genome") +
        ggplot2::labs(color = "Cells/well") +
        ggplot2::scale_y_continuous(labels = scales::comma,
            limits = c(0, NA)) +
        .themePublication()
    return (g)
}


.plotReadsMappedToGenes <- function(qcDt) {
    g <- ggplot2::ggplot(data = qcDt,
        ggplot2::aes(
            x = as.factor(experiment),
            y = log10(reads_mapped_to_genes),
            group = as.factor(experiment))) +
        ggplot2::geom_boxplot(outlier.color = NA, fill = NA) +
        ggplot2::geom_point(
            ggplot2::aes(color = as.factor(number_of_cells)),
            position = ggplot2::position_jitter(width = 0.3, height = 0),
            size = 0.5) +
        ggplot2::ylab(expression(bold(Log[10]*"Reads"))) +
        ggplot2::xlab("Experiment") +
        ggplot2::ggtitle("Reads mapped to genes") +
        ggplot2::labs(color = "Cells/well") +
        ggplot2::scale_y_continuous(labels = scales::comma,
            limits = c(0, NA)) +
        .themePublication()
    return (g)
}


.plotGenomeReadsFraction <- function(qcDt) {
    g <- ggplot2::ggplot(data = qcDt,
        ggplot2::aes(
            x = as.factor(experiment),
            y = reads_mapped_to_genome/reads,
            group = as.factor(experiment))) +
        ggplot2::geom_boxplot(outlier.color = NA, fill = NA) +
        ggplot2::geom_point(
            ggplot2::aes(color = as.factor(number_of_cells)),
            position = ggplot2::position_jitter(width = 0.3, height = 0),
            size = 0.5) +
        ggplot2::ylim(0, 1) +
        ggplot2::xlab("Experiment") +
        ggplot2::ggtitle("Fraction of aligned reads to total reads") +
        ggplot2::labs(color = "Cells/well") +
        .themePublication() + 
        ggplot2::theme(axis.title.y = ggplot2::element_blank())
    return (g)
}


.plotGeneToGenomeFraction <- function(qcDt) {
    g <- ggplot2::ggplot(data = qcDt,
        ggplot2::aes(
            x = as.factor(experiment),
            y = reads_mapped_to_genes/reads_mapped_to_genome,
            group = as.factor(experiment))) +
        ggplot2::geom_boxplot(outlier.color = NA, fill = NA) +
        ggplot2::geom_point(
            ggplot2::aes(color = as.factor(number_of_cells)),
            position = ggplot2::position_jitter(width = 0.3, height = 0),
            size = 0.5) +
        ggplot2::ylim(0, 1) +
        ggplot2::xlab("Experiment") +
        ggplot2::ggtitle("Fraction of gene reads out of aligned reads") +
        ggplot2::labs(color = "Cells/well") +
        .themePublication() +
        ggplot2::theme(axis.title.y = ggplot2::element_blank())
    return (g)
}


.plotGeneToTotalFraction <- function(qcDt) {
    g <- ggplot2::ggplot(data = qcDt,
        ggplot2::aes(
            x = as.factor(experiment),
            y = reads_mapped_to_genes/reads,
            group = as.factor(experiment))) +
        ggplot2::geom_boxplot(outlier.color = NA, fill = NA) +
        ggplot2::geom_point(
            ggplot2::aes(color = as.factor(number_of_cells)),
            position = ggplot2::position_jitter(width = 0.3, height = 0),
            size = 0.5) +
        ggplot2::ylim(0, 1) +
        ggplot2::xlab("Experiment") +
        ggplot2::ggtitle("Fraction of gene reads out of total reads") +
        ggplot2::labs(color = "Cells/well") +
        .themePublication() +
        ggplot2::theme(axis.title.y = ggplot2::element_blank())
    return (g)
}


.plotCounts <- function(qcDt) {
    g <- ggplot2::ggplot(data = qcDt,
        ggplot2::aes(
            x = as.factor(experiment),
            y = log10(total_counts),
            group = as.factor(experiment))) +
        ggplot2::geom_boxplot(outlier.color = NA, fill = NA) +
        ggplot2::geom_point(
            ggplot2::aes(color = as.factor(number_of_cells)),
            position = ggplot2::position_jitter(width = 0.3, height = 0),
            size = 0.5) +
        ggplot2::ylab(expression(bold(Log[10]*"Counts"))) +
        ggplot2::xlab("Experiment") +
        ggplot2::ggtitle("Total transcripts") +
        ggplot2::labs(color = "Cells/well") +
        ggplot2::scale_y_continuous(labels = scales::comma,
            limits = c(0, NA)) +
        .themePublication()
    return (g)
}


.plotMtCounts <- function(qcDt) {
    g <- ggplot2::ggplot(data = qcDt,
        ggplot2::aes(
            x = as.factor(experiment),
            y = log10(mt_counts),
            group = as.factor(experiment))) +
        ggplot2::geom_boxplot(outlier.color = NA, fill = NA) +
        ggplot2::geom_point(
            ggplot2::aes(color = as.factor(number_of_cells)),
            position = ggplot2::position_jitter(width = 0.3, height = 0),
            size = 0.5) +
        ggplot2::ylab(expression(bold(Log[10]*"Counts"))) +
        ggplot2::xlab("Experiment") +
        ggplot2::ggtitle("Mitochondrial transcripts") +
        ggplot2::labs(color = "Cells/well") +
        ggplot2::scale_y_continuous(labels = scales::comma,
            limits = c(0, NA)) +
        .themePublication()
    return (g)
}


.plotMtCountsFraction <- function(qcDt) {
    g <- ggplot2::ggplot(data = qcDt,
        ggplot2::aes(
            x = as.factor(experiment),
            y = mt_counts/total_counts,
            group = as.factor(experiment))) +
        ggplot2::geom_boxplot(outlier.color = NA, fill = NA) +
        ggplot2::geom_point(
            ggplot2::aes(color = as.factor(number_of_cells)),
            position = ggplot2::position_jitter(width = 0.3, height = 0),
            size = 0.5) +
        ggplot2::ylim(0, 1) +
        ggplot2::xlab("Experiment") +
        ggplot2::ggtitle("Fraction of mitochondrial transcripts") +
        ggplot2::labs(color = "Cells/well") +
        .themePublication() +
        ggplot2::theme(axis.title.y = ggplot2::element_blank())
    return (g)
}


.plotGenes <- function(qcDt) {
    g <- ggplot2::ggplot(data = qcDt,
        ggplot2::aes(
            x = as.factor(experiment),
            y = log10(genes),
            group = as.factor(experiment))) +
        ggplot2::geom_boxplot(outlier.color = NA, fill = NA) +
        ggplot2::geom_point(
            ggplot2::aes(color = as.factor(number_of_cells)),
            position = ggplot2::position_jitter(width = 0.3, height = 0),
            size = 0.5) +
        ggplot2::ylab(expression(bold(Log[10]*"Genes"))) +
        ggplot2::xlab("Experiment") +
        ggplot2::ggtitle("Transcribed genes") +
        ggplot2::labs(color = "Cells/well") +
        ggplot2::scale_y_continuous(labels = scales::comma,
            limits = c(0, NA)) +
        .themePublication()
    return (g)
}


.plotFracProteinCodingGenes <- function(qcDt) {
    g <- ggplot2::ggplot(data = qcDt,
        ggplot2::aes(
            x = as.factor(experiment),
            y = protein_coding_genes/genes,
            group = as.factor(experiment))) +
        ggplot2::geom_boxplot(outlier.color = NA, fill = NA) +
        ggplot2::geom_point(
            ggplot2::aes(color = as.factor(number_of_cells)),
            position = ggplot2::position_jitter(width = 0.3, height = 0),
            size = 0.5) +
        ggplot2::ylim(0, 1) +
        ggplot2::xlab("Experiment") +
        ggplot2::ggtitle("Fraction of protein coding genes") +
        ggplot2::labs(color = "Cells/well") +
        .themePublication() +
        ggplot2::theme(axis.title.y = ggplot2::element_blank())
    return (g)
}


.plotFracProteinCodingTranscripts <- function(qcDt) {
    g <- ggplot2::ggplot(qcDt,
        ggplot2::aes(
            x = as.factor(experiment),
            y = protein_coding_counts/total_counts,
            group = as.factor(experiment))) +
        ggplot2::geom_boxplot(outlier.color = NA, fill = NA) +
        ggplot2::geom_point(
            ggplot2::aes(color = as.factor(number_of_cells)),
            position = ggplot2::position_jitter(width = 0.3, height = 0),
            size = 0.5) +
        ggplot2::ylim(0, 1) +
        ggplot2::xlab("Experiment") +
        ggplot2::ggtitle("Fraction of protein coding transcripts") +
        ggplot2::labs(color = "Cells/well") +
        .themePublication() +
        ggplot2::theme(axis.title.y = ggplot2::element_blank())
    return (g)
}


.plotGenesPerMillionReads <- function(qcDt) {
    g <- ggplot2::ggplot(data = qcDt,
        ggplot2::aes(
            x = as.factor(experiment),
            y = log10(genes * 1000000/reads),
            group = as.factor(experiment))) +
        ggplot2::geom_boxplot(outlier.color = NA, fill = NA) +
        ggplot2::geom_point(
            ggplot2::aes(color = as.factor(number_of_cells)),
            position = ggplot2::position_jitter(width = 0.3,
                height = 0),
            size = 0.5) +
        ggplot2::ylab(
            expression(bold(paste(Log[10],
                "(Genes x 1000000 / total reads)")))) +
        ggplot2::ggtitle(
            paste("Genes detected divided by total number of reads",
                "sequenced per million")) +
        ggplot2::xlab("Experiment") +
        ggplot2::labs(color = "Cells/well") +
        ggplot2::scale_y_continuous(labels = scales::comma,
            limits = c(0, NA)) +
        .themePublication()
    return (g)
}

