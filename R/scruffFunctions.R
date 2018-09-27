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


.stripLeadingUnderscore <- function (x)  sub("^\\_+", "", x)


# remove file extension and get basename
.removeLastExtension <- function(x) {
    return (sub(pattern = "\\.[^\\.]*$",
        replacement = "\\1",
        basename(x)))
}


.getAlignmentFilePaths <- function(fastq.paths, format, out.dir) {
    file.paths <- file.path(out.dir,
        paste0(
            sub(
                pattern = "(.*?)\\..*$",
                replacement = "\\1",
                basename(fastq.paths)
            ),
            ".",
            format
        ))
    return (file.paths)
}


########################################
##### procedural helper functions ######
########################################

# read gtf database and return feature GRangesList by gene ID
.gtfReadDb <- function(gtf) {
    gtf.db.file <- paste0(basename(gtf), ".sqlite")
    if ((!(file.exists(gtf))) & (!(file.exists(gtf.db.file)))) {
        stop(paste("File", gtf, "does not exist"))
    }
    
    if (!(file.exists(gtf.db.file))) {
        message(paste(Sys.time(), "... TxDb file",
            gtf.db.file,
            "does not exist"))
        message(paste(Sys.time(), "... Creating TxDb object", gtf.db.file))
        gtf.db <- GenomicFeatures::makeTxDbFromGFF(file = gtf)
        AnnotationDbi::saveDb(gtf.db, file = gtf.db.file)
        return (GenomicFeatures::exonsBy(gtf.db, by = "gene"))
    }
    
    gtf.db <- tryCatch(
        suppressPackageStartupMessages(AnnotationDbi::loadDb(gtf.db.file)),
        error = function(e) {
            stop("Error loading database file. Delete the file",
                gtf.db.file,
                "and try again.")
        }
    )
    return (GenomicFeatures::exonsBy(gtf.db, by = "gene"))
}


# correct barcode mismatch using memoization
.bcCorrectMem <- local({
    res <- list()
    
    f <- function(bc, refBarcodes, maxEditDist) {
        if (bc %in% names(res))
            return (res[[bc]])
        if (bc %in% refBarcodes) {
            res[[bc]] <<- bc
            return (res[[bc]])
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
                return (res[[bc]])
            }
        }
        res[[bc]] <<- bc
        return (res[[bc]])
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
        error = function(e) {}
    )
    return (sub(
        pattern = "\\.sam$",
        ignore.case = TRUE,
        perl = TRUE,
        replacement = ".BAM",
        x = sam
    ))
}


.checkCores <- function(cores) {
    if (cores > parallel::detectCores()) {
        stop("The specified cores ",
            cores,
            " is greater than the number of cores available ",
            parallel::detectCores())
    }
}


# A function that returns an iterator that reads BAM files
.bamIterator <- function(bamfl, param) {
    return (function() {
        if (Rsamtools::isIncomplete(bamfl)) {
            yld <- GenomicAlignments::readGAlignments(bamfl,
                use.names = TRUE,
                param = param)
            return(yld)
        } else {
            return (NULL)
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


.getGeneAnnotationRefGenome <- function(reference, features) {
    gtfEG = refGenome::ensemblGenome(dirname(reference))
    refGenome::read.gtf(gtfEG, filename = basename(reference))
    geneAnnotation <- data.table::data.table(
        unique(refGenome::getGeneTable(gtfEG)[, c("gene_id",
            "gene_name",
            "gene_biotype",
            "seqid")]))
    geneAnnotation <- geneAnnotation[order(gene_id), ]
    
    if (length(grep("ERCC", names(features))) > 0) {
        ercc <- features[grep("ERCC", names(features))]
        erccDt <- data.table::data.table(
            gene_id = base::names(ercc),
            gene_name = base::names(ercc),
            gene_biotype = "ERCC",
            seqid = "ERCC"
        )
        geneAnnotation <- rbind(geneAnnotation, erccDt)
    }
    return (geneAnnotation)
}


# merge hash table for 10X BAM file counting
.mergeHashCbGeneUmi <- function(hashCb1, hashCb2) {
    cb1 = ls(hashCb1)
    cb2 = ls(hashCb2)
    
    for (cb in cb2) {
        if (cb %in% cb1) {
            # if cell barcodes intersect
            g1 <- ls(hashCb1[[cb]])
            g2 <- ls(hashCb2[[cb]])
            for (g in g2) {
                if (g %in% g1) {
                    # if genes intersect
                    umi1 <- hashCb1[[cb]][[g]]
                    umi2 <- hashCb2[[cb]][[g]]
                    hashCb1[[cb]][[g]] <- union(umi1, umi2)
                } else {
                    # else append g to hashCb1[[cb]]
                    hashCb1[[cb]][[g]] <- hashCb2[[cb]][[g]]
                }
            }
        } else {
            # else append cb to hashCb1
            hashCb1[[cb]] <- hashCb2[[cb]]
        }
    }
    return (hashCb1)
}


.getGeneFromEnv <- function(env) {
    vapply(env, function(env) {
        genes <- lapply(ls(env), function(x) {get(x, envir = env)})
        return (length(genes))},
        integer(1)
    )
}


.getUMICountsFromEnv <- function(env) {
    vapply(env, function(env) {
        genes <- lapply(ls(env), function(x) {get(x, envir = env)})
        return (sum(vapply(genes, length, integer(1))))
    }, integer(1))
}


.getTopNBarcodesFromEnv <- function(N, env) {
    allBarcodeUMICounts <- sort(.getUMICountsFromEnv(env), decreasing = TRUE)
    firstNBarcodes <- allBarcodeUMICounts[seq_len(N)]
    return (firstNBarcodes)
}


.getUMICountsFromCell <- function(cellEnv) {
    UMIs <- sapply(ls(cellEnv),
        function(x) {get(x, envir = cellEnv)},
        simplify = FALSE,
        USE.NAMES = TRUE)
    return (vapply(UMIs, length, integer(1)))
}


.getFilteredBarcode <- function(N, env, percentile = 0.99, orderRatio = 10) {
    allUMICounts <- .getUMICountsFromEnv(env)
    topNUMICounts <- .getTopNBarcodesFromEnv(N = N, env)
    topNUMISum <- sum(topNUMICounts)
    m <- topNUMISum * percentile
    filteredBarcodes <- allUMICounts[which(allUMICounts > m/orderRatio)]
    return (filteredBarcodes)
}


.subsetEnv <- function(env, filteredBarcodes) {
    newEnv <- new.env()
    for (i in seq_len(length(filteredBarcodes))) {
        newEnv[[filteredBarcodes[i]]] <- env[[filteredBarcodes[i]]]
    }
    return (newEnv)
}


.getFilteredUMICountsTable <- function(fliteredEnv,
    filteredBarcodes,
    geneAnnotation,
    experiment) {
    
    resDt <- data.table::data.table(geneid = geneAnnotation[, gene_id])    
    cells <- paste(experiment, filteredBarcodes, sep = "_")
    filteredCounts <- .getUMICountsFromEnv(fliteredEnv)
    for (i in seq_len(length(filteredBarcodes))) {
        resDt[[cells[i]]] <- 0
        UMIs <- .getUMICountsFromCell(fliteredEnv[[filteredBarcodes[i]]])
        resDt[geneid %in% ls(fliteredEnv[[filteredBarcodes[i]]]),
            eval(cells[i]) := as.numeric(UMIs[geneid])]
    }
    
    return (resDt)
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
                panel.background = ggplot2::element_rect(color = NA),
                plot.background = ggplot2::element_rect(color = NA),
                panel.border = ggplot2::element_rect(color = NA),
                axis.title = ggplot2::element_text(
                    face = "bold",
                    size = ggplot2::rel(1)),
                axis.title.y = ggplot2::element_text(angle = 90,
                    vjust = 2),
                axis.title.x = ggplot2::element_text(vjust = -0.2),
                axis.text = ggplot2::element_text(),
                axis.line = ggplot2::element_line(color="black"),
                axis.ticks = ggplot2::element_line(),
                panel.grid.major = ggplot2::element_line(color = "#f0f0f0"),
                panel.grid.minor = ggplot2::element_blank(),
                legend.key = ggplot2::element_rect(color = NA),
                legend.position = "right",
                legend.direction = "vertical",
                legend.key.size= ggplot2::unit(0.2, "cm"),
                legend.margin = ggplot2::margin(0),
                legend.title = ggplot2::element_text(face = "bold"),
                plot.margin = ggplot2::unit(c(10, 5, 5, 5), "mm"),
                strip.background = ggplot2::element_rect(
                    color = "#f0f0f0", fill = "#f0f0f0"),
                strip.text = ggplot2::element_text(face = "bold")
            ))
}

