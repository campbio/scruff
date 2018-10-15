#' Visualize data quality
#' 
#' Visualize data quality from the \code{colData} of the
#'  \code{SingleCellExperiment} object and return a list of figures in
#'  \code{arrangelist} object.
#' 
#' @param sce An \code{SingleCellExperiment} object returned from [scruff],
#'  [countUMI], or [tenxBamqc] function.
#' @return A list of \code{grobs} objects ready for plotting
#' @examples
#' data(sceExample, package = "scruff")
#' qcplots(sceExample)
#' @export
qcplots <- function(sce) {
    qcDt <- data.table::as.data.table(SummarizedExperiment::colData(sce))
    
    g1 <- .plotTotalReads(qcDt)
    g2 <- .plotReadsMappedToGenome(qcDt)
    g3 <- .plotReadsMappedToGenes(qcDt)
    g4 <- .plotGenomeReadsFraction(qcDt)
    g5 <- .plotGeneToGenomeFraction(qcDt)
    g6 <- .plotGeneToTotalFraction(qcDt)
    g7 <- .plotCounts(qcDt)
    g8 <- .plotMtCounts(qcDt)
    g9 <- .plotMtCountsFraction(qcDt)
    g10 <- .plotGenes(qcDt)
    g11 <- .plotFracProteinCodingGenes(qcDt)
    g12 <- .plotFracProteinCodingTranscripts(qcDt)
    g13 <- .plotGenesPerMillionReads(qcDt)
    return (list(g1, g2, g3, g4, g5, g6, g7, g8, g9, g10, g11, g12, g13))
}


.plotTotalReads <- function(qcDt) {
    if (!"reads" %in% colnames(qcDt)) {
        return (NULL)
    }
    
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
        ggplot2::xlab("Experiment") +
        ggplot2::ggtitle("Total reads") +
        ggplot2::labs(color = "Cells") +
        ggplot2::scale_y_continuous(name = "Reads",
            limits = c(0, NA),
            labels = scales::math_format(10^.x)) +
        ggplot2::annotation_logticks(sides = "l") +
        .themePublication() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45,
            hjust = 1))
    return (g)
}


.plotReadsMappedToGenome <- function(qcDt) {
    if (!"reads_mapped_to_genome" %in% colnames(qcDt)) {
        return (NULL)
    }
    
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
        ggplot2::xlab("Experiment") +
        ggplot2::ggtitle("Reads aligned to reference genome") +
        ggplot2::labs(color = "Cells") +
        ggplot2::scale_y_continuous(name = "Reads",
            limits = c(0, NA),
            labels = scales::math_format(10^.x)) +
        ggplot2::annotation_logticks(sides = "l") +
        .themePublication() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45,
            hjust = 1))
    return (g)
}


.plotReadsMappedToGenes <- function(qcDt) {
    if (!("reads_mapped_to_genes" %in% colnames(qcDt) &
            "reads" %in% colnames(qcDt))) {
        return (NULL)
    }
    
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
        ggplot2::xlab("Experiment") +
        ggplot2::ggtitle("Reads mapped to genes") +
        ggplot2::labs(color = "Cells") +
        ggplot2::scale_y_continuous(name = "Reads",
            limits = c(0, NA),
            labels = scales::math_format(10^.x)) +
        ggplot2::annotation_logticks(sides = "l") +
        .themePublication() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45,
            hjust = 1))
    return (g)
}


.plotGenomeReadsFraction <- function(qcDt) {
    if (!("reads_mapped_to_genome" %in% colnames(qcDt) &
            "reads" %in% colnames(qcDt))) {
        return (NULL)
    }
    
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
        ggplot2::labs(color = "Cells") +
        .themePublication() + 
        ggplot2::theme(axis.title.y = ggplot2::element_blank(),
            axis.text.x = ggplot2::element_text(angle = 45,
                hjust = 1))
    return (g)
}


.plotGeneToGenomeFraction <- function(qcDt) {
    if (!("reads_mapped_to_genome" %in% colnames(qcDt) &
            "reads_mapped_to_genome" %in% colnames(qcDt))) {
        return (NULL)
    }
    
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
        ggplot2::labs(color = "Cells") +
        .themePublication() +
        ggplot2::theme(axis.title.y = ggplot2::element_blank(),
            axis.text.x = ggplot2::element_text(angle = 45,
                hjust = 1))
    return (g)
}


.plotGeneToTotalFraction <- function(qcDt) {
    if (!("reads_mapped_to_genes" %in% colnames(qcDt) &
            "reads" %in% colnames(qcDt))) {
        return (NULL)
    }
    
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
        ggplot2::labs(color = "Cells") +
        .themePublication() +
        ggplot2::theme(axis.title.y = ggplot2::element_blank(),
            axis.text.x = ggplot2::element_text(angle = 45,
                hjust = 1))
    return (g)
}


.plotCounts <- function(qcDt) {
    if (!"total_counts" %in% colnames(qcDt)) {
        return (NULL)
    }
    
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
        ggplot2::xlab("Experiment") +
        ggplot2::ggtitle("Total transcripts") +
        ggplot2::labs(color = "Cells") +
        ggplot2::scale_y_continuous(name = "Counts",
            limits = c(0, NA),
            labels = scales::math_format(10^.x)) +
        ggplot2::annotation_logticks(sides = "l") +
        .themePublication() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45,
            hjust = 1))
    return (g)
}


.plotMtCounts <- function(qcDt) {
    if (!"mt_counts" %in% colnames(qcDt)) {
        return (NULL)
    }
    
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
        ggplot2::xlab("Experiment") +
        ggplot2::ggtitle("Mitochondrial transcripts") +
        ggplot2::labs(color = "Cells") +
        ggplot2::scale_y_continuous(name = "Counts",
            limits = c(0, NA),
            labels = scales::math_format(10^.x)) +
        ggplot2::annotation_logticks(sides = "l") +
        .themePublication() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45,
            hjust = 1))
    return (g)
}


.plotMtCountsFraction <- function(qcDt) {
    if (!("mt_counts" %in% colnames(qcDt) &
            "total_counts" %in% colnames(qcDt))) {
        return (NULL)
    }
    
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
        ggplot2::labs(color = "Cells") +
        .themePublication() +
        ggplot2::theme(axis.title.y = ggplot2::element_blank(),
            axis.text.x = ggplot2::element_text(angle = 45,
                hjust = 1))
    return (g)
}


.plotGenes <- function(qcDt) {
    if (!"genes" %in% colnames(qcDt)) {
        return (NULL)
    }
    
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
        ggplot2::xlab("Experiment") +
        ggplot2::ggtitle("Transcribed genes") +
        ggplot2::labs(color = "Cells") +
        ggplot2::scale_y_continuous(name = "Genes",
            limits = c(0, NA),
            labels = scales::math_format(10^.x)) +
        ggplot2::annotation_logticks(sides = "l") +
        .themePublication() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45,
            hjust = 1))
    return (g)
}


.plotFracProteinCodingGenes <- function(qcDt) {
    if (!("protein_coding_genes" %in% colnames(qcDt) &
            "genes" %in% colnames(qcDt))) {
        return (NULL)
    }
    
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
        ggplot2::labs(color = "Cells") +
        .themePublication() +
        ggplot2::theme(axis.title.y = ggplot2::element_blank(),
            axis.text.x = ggplot2::element_text(angle = 45,
                hjust = 1))
    return (g)
}


.plotFracProteinCodingTranscripts <- function(qcDt) {
    if (!("protein_coding_counts" %in% colnames(qcDt) &
            "total_counts" %in% colnames(qcDt))) {
        return (NULL)
    }
    
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
        ggplot2::labs(color = "Cells") +
        .themePublication() +
        ggplot2::theme(axis.title.y = ggplot2::element_blank(),
            axis.text.x = ggplot2::element_text(angle = 45,
                hjust = 1))
    return (g)
}


.plotGenesPerMillionReads <- function(qcDt) {
    if (!"reads" %in% colnames(qcDt)) {
        return (NULL)
    }
    
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
        ggplot2::ggtitle(
            paste("Genes detected divided by total number of reads",
                "sequenced per million")) +
        ggplot2::xlab("Experiment") +
        ggplot2::labs(color = "Cells") +
        ggplot2::scale_y_continuous(name = expression(bold(
            "(Genes x 1000000 / total reads)")),
            limits = c(0, NA),
            labels = scales::math_format(10^.x)) +
        ggplot2::annotation_logticks(sides = "l") +
        .themePublication() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45,
            hjust = 1))
    return (g)
}

