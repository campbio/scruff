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
        error = function(e)
            stop(
                paste(
                    "Error loading database file. Delete the file",
                    gtf.db.file,
                    "and try again."
                )
            )
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

