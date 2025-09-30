########################################
####### parsing helper functions #######
########################################

.logMessages <- function(...,
    sep = " ",
    logfile = NULL,
    append = FALSE) {
    if (!is.null(logfile)) {
        if (!is.character(logfile) || length(logfile) > 1) {
            stop("The log file parameter needs to be ",
                "a single character string.")
        }
        cat(paste(..., "\n", sep = sep),
            file = logfile,
            append = append)

    } else {
        message(paste(..., sep = sep))
    }
}


.stripLeadingUnderscore <- function(x)  sub("^\\_+", "", x)


# remove file extension and get basename
.removeLastExtension <- function(x) {
    return(sub(pattern = "\\.[^\\.]*$",
        replacement = "\\1",
        basename(x)))
}


.getAlignmentFilePaths <- function(fastq.paths, format, out.dir) {
    fileName <- paste0(sub(pattern = "(.*?)\\..*$",
        replacement = "\\1",
        basename(fastq.paths)), ".", format)

    # replace all punctuations with "." per Rsubread update
    fileName <- gsub("[[:punct:]]+", ".", fileName)
    fileName <- gsub(" ", ".", fileName)

    filePaths <- file.path(out.dir, fileName)
    return(filePaths)
}


########################################
##### procedural helper functions ######
########################################

# read gtf database and return feature GRangesList by gene ID
.gtfReadDb <- function(gtf) {
    message(Sys.time(),
        " ... Reading GTF file ", gtf, " as GRangesList object")
    gtf.db.file <- paste0(basename(gtf), ".sqlite")
    if ((!(file.exists(gtf))) & (!(file.exists(gtf.db.file)))) {
        stop(paste("File", gtf, "does not exist"))
    }

    if (!(file.exists(gtf.db.file))) {
        message(paste(Sys.time(), "... TxDb file",
            gtf.db.file,
            "does not exist"))
        message(paste(Sys.time(), "... Creating TxDb object", gtf.db.file))
        gtf.db <- txdbmaker::makeTxDbFromGFF(file = gtf)
        AnnotationDbi::saveDb(gtf.db, file = gtf.db.file)
        return(GenomicFeatures::exonsBy(gtf.db, by = "gene"))
    }

    gtf.db <- tryCatch(
        suppressPackageStartupMessages(AnnotationDbi::loadDb(gtf.db.file)),
        error = function(e) {
            stop("Error loading database file. Delete the file",
                    gtf.db.file,
                    "and try again.")
        }
    )
    return(GenomicFeatures::exonsBy(gtf.db, by = "gene"))
}


# correct barcode mismatch using memoization
.bcCorrectMem <- local({
    res <- list()

    f <- function(bc, refBarcodes, maxEditDist) {
        if (bc %in% names(res))
            return(res[[bc]])
        if (bc %in% refBarcodes) {
            res[[bc]] <<- bc
            return(res[[bc]])
        }

        sdm <- stringdist::stringdistmatrix(bc,
            refBarcodes,
            method = "hamming",
            nthread = 1)
        min.dist <- min(sdm)
        if (min.dist <= maxEditDist) {
            ind <- which(sdm == min.dist)
            if (length(ind) == 1) {
                res[[bc]] <<- refBarcodes[ind]
                return(res[[bc]])
            }
        }
        res[[bc]] <<- bc
        return(res[[bc]])
    }
})


.toBam <- function(sam,
    logfile,
    overwrite = FALSE,
    index = FALSE) {
    .logMessages(
        Sys.time(),
        "... Converting",
        sam,
        "to BAM format (if not exist)",
        logfile = logfile,
        append = TRUE
    )
    tryCatch(
        Rsamtools::asBam(sam, overwrite = overwrite, indexDestination = index),
        error <- function(e) {
            "Error converting samfiles to bamfiles"
            }
    )
    return(sub(
        pattern = "\\.sam$",
        ignore.case = TRUE,
        perl = TRUE,
        replacement = ".BAM",
        x = sam
    ))
}


.checkCores <- function(cores) {
    if (cores > parallelly::availableCores()) {
        stop("The specified cores ",
            cores,
            " is greater than the number of cores available ",
            parallelly::availableCores())
    }
}


# A function that returns an iterator that reads BAM files
.bamIterator <- function(bamfl, param) {
    return(function() {
        if (Rsamtools::isIncomplete(bamfl)) {
            yld <- GenomicAlignments::readGAlignments(bamfl,
                use.names = TRUE,
                param = param)
            return(yld)
        } else {
            return(NULL)
        }
    })
}


# A function that returns an iterator that reads FASTQ read1 and read2 files
.fastqIterator <- function(fq1, fq2) {
    done <- FALSE
    return(function() {
        if (done) {
            return(NULL)
        }

        yld1 <- ShortRead::yield(fq1)
        yld2 <- ShortRead::yield(fq2)
        if (length(yld1) != length(yld2)) {
            stop("Unequal number of reads",
                " between read1 and read2 fastq files: ",
                fq1,
                fq2)
        } else if (length(yld1) == 0L & length(yld2) == 0L) {
            done <<- TRUE
            return(NULL)
        } else {
            return(list(yld1, yld2))
        }
    })
}


# Check cell barcodes
.checkCellBarcodes <- function(bc, bcStart, bcStop, verbose) {
    if (bcStop < bcStart) {
        stop("bcStop ",
            bcStop,
            " should be equal or greater than bcStart ",
            bcStart)
    }

    if (verbose) {
        message("Cell barcode input vector:\n")
        print(bc)
    }

    if (length(bc) < 10) {
        warning("Length of cell barcode vector is less than 10!")
    }

    if (length(unique(nchar(bc))) > 1) {
        warning("The cell barcode input vector has variable lengths!")
    } else {
        bcL <- bcStop - bcStart + 1
        if (unique(nchar(bc)) != bcL) {
            stop("The number of barcode bases in cell barcode input file is ",
                unique(nchar(bc)),
                " but the number is ",
                bcL,
                " between 'bcStart' and 'bcStop'!")
        }
    }
}


# .getGeneAnnotationRefGenome <- function(reference, features) {
#     gtfEG = refGenome::ensemblGenome(dirname(reference))
#     refGenome::read.gtf(gtfEG, filename = basename(reference))
#     geneAnnotation <- data.table::data.table(
#         unique(refGenome::getGeneTable(gtfEG)[, c("gene_id",
#             "gene_name",
#             "gene_biotype",
#             "seqid")]))
#     geneAnnotation <- geneAnnotation[order(gene_id), ]
#
#     if (length(grep("ERCC", names(features))) > 0) {
#         ercc <- features[grep("ERCC", names(features))]
#         erccDt <- data.table::data.table(
#             gene_id = base::names(ercc),
#             gene_name = base::names(ercc),
#             gene_biotype = "ERCC",
#             seqid = "ERCC"
#         )
#         geneAnnotation <- rbind(geneAnnotation, erccDt)
#     }
#     return(geneAnnotation)
# }


.getGeneAnnotation <- function(reference) {
    message(Sys.time(),
        " ... Reading GTF file ",
        reference, " as data.table object")
    gtf <- rtracklayer::readGFF(reference)

    geneAnnotation <- data.table::as.data.table(gtf)
    seqid <- "seqid"
    if (!seqid %in% colnames(geneAnnotation)) {
      seqid <- "seqnames"
    }

    if ("gene_biotype" %in% colnames(geneAnnotation)) {
        geneAnnotation <- unique(geneAnnotation[type == "gene" |
                source == "ERCC", c("gene_id",
                    "gene_name",
                    "gene_biotype",
                    seqid,
                    "start",
                    "end",
                    #"width",
                    "strand",
                    "source"), with = FALSE])
    } else if ("gene_type" %in% colnames(geneAnnotation)) {
        message("Missing column 'gene_biotype'. Use 'gene_type' instead.")
        geneAnnotation <- unique(geneAnnotation[type == "gene" |
                source == "ERCC", c("gene_id",
                    "gene_name",
                    "gene_type",
                    seqid,
                    "start",
                    "end",
                    #"width",
                    "strand",
                    "source"), with = FALSE])
        colnames(geneAnnotation)[which(colnames(geneAnnotation) ==
                "gene_type")] <- "gene_biotype"
    } else {
        warning("Missing column 'gene_biotype' or 'gene_type'!")
        geneAnnotation <- unique(geneAnnotation[type == "gene" |
                source == "ERCC", c("gene_id",
                    "gene_name",
                    #"gene_biotype",
                    seqid,
                    "start",
                    "end",
                    #"width",
                    "strand",
                    "source"), with = FALSE])
    }

    geneAnnotation <- geneAnnotation[order(gene_id), ]
    geneAnnotation[source == "ERCC", gene_name := gene_id]
    geneAnnotation[source == "ERCC", gene_biotype := source]

    return(geneAnnotation)
}


.checkGTF <- function(gtfDt, cnames) {
    for (i in cnames) {
        if (!(i %in% colnames(gtfDt))) {
            warning("GTF file does not contain ", i, " column")
        }
    }
}


.tenxCheckCellBarcodes <- function(bam, cbfile, validCb, yieldSize = 10000,
    tags = "CB") {

    message(Sys.time(), " Checking cell barcodes")
    utils::data(cbtop10000, package = "scruff", envir = environment())

    bamfl <- Rsamtools::BamFile(bam, yieldSize = yieldSize)
    param <- Rsamtools::ScanBamParam(tag = tags)

    open(bamfl)

    yld <- GenomicAlignments::readGAlignments(bamfl,
        use.names = TRUE,
        param = param)

    close(bamfl)

    cb <- S4Vectors::mcols(yld)$CB
    cb <- data.table::tstrsplit(cb, "-")[[1]]

    if (is.na(cbfile)) {
        l <- nchar(cbtop10000[, v2chemistry][1])
    } else {
        l <- table(length(validCb[[1]]))
    }

    # check barcode length
    if (any(!nchar(cb)[stats::complete.cases(nchar(cb))] %in% l)) {
        cbtb <- table(nchar(cb))
        cbtb2 <- cbtb
        if (length(cbtb2) > 1) {
            cbtb2[as.character(l)] <- NULL
        } else {
            cbtb2 <- 0
        }
        warning("In the first ", length(cb), " alignments in BAM file ", bam,
            ", the lengths of ", sum(cbtb2),
            "/", sum(cbtb, na.rm = TRUE),
            " cell barcodes (", paste(names(cbtb2), collapse = " "),
            ") are not equal to the legnths of cell barcodes",
            " in the provided whitelist (", paste(l, collapse = " "),
            "). Make sure your 'validCb' input is correct!")
    }

    v1s <- sum(cb %in% cbtop10000[, v1chemistry])
    v2s <- sum(cb %in% cbtop10000[, v2chemistry])
    v3s <- sum(cb %in% cbtop10000[, v3chemistry])

    message("In the first ", length(cb), " alignments in BAM file ", bam,
        ", ", v1s, ", ", v2s, ", and ", v3s,
        " cell barcodes are found within the first 10000",
        " cell barcodes in the v1, v2, and v3 chemistry whitelists. Make sure",
        " your 'validCb' input is correct!")
}


########################################
###### Plotting helper functions #######
########################################

# ggplot publication theme
# Make ggplots look better.
# Adapted from Koundinya Desiraju. https://rpubs.com/Koundy/71792
.themePublication <- function(base_size = 12,
    base_family = "sans") {
    (ggthemes::theme_foundation(base_size = base_size,
        base_family = base_family) +
            ggplot2::theme(plot.title = ggplot2::element_text(
                face = "bold",
                size = ggplot2::rel(1),
                hjust = 0.5),
                text = ggplot2::element_text(),
                panel.background = ggplot2::element_rect(fill = NA,
                    color = NA),
                plot.background = ggplot2::element_rect(fill = NA,
                    color = NA),
                panel.border = ggplot2::element_blank(),
                axis.title = ggplot2::element_text(
                    face = "bold",
                    size = ggplot2::rel(1)),
                axis.title.y = ggplot2::element_text(angle = 90,
                    vjust = 2),
                axis.title.x = ggplot2::element_text(vjust = -0.2),
                axis.text = ggplot2::element_text(),
                axis.line = ggplot2::element_line(color = "black"),
                axis.ticks = ggplot2::element_line(),
                panel.grid.major = ggplot2::element_line(color = "#f0f0f0"),
                panel.grid.minor = ggplot2::element_blank(),
                legend.key = ggplot2::element_rect(color = NA),
                legend.position = "right",
                legend.direction = "vertical",
                legend.key.size = ggplot2::unit(0.2, "cm"),
                legend.margin = ggplot2::margin(0),
                legend.title = ggplot2::element_text(face = "bold"),
                plot.margin = ggplot2::unit(c(10, 5, 5, 5), "mm"),
                strip.background = ggplot2::element_rect(
                    color = "#f0f0f0", fill = "#f0f0f0"),
                strip.text = ggplot2::element_text(face = "bold")
            ))
}
