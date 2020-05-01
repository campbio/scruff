#' A wrapper to \code{Rsubread} read alignment function \code{align}
#'
#' This function is \strong{not} available in Windows environment. Align cell
#'  specific reads to reference genome and write sequence alignment results to
#'  output directory. A wrapper to the \code{align} function in \code{Rsubread}
#'  package. For details please refer to \code{Rsubread} manual.
#'
#' @param sce A \code{SingleCellExperiment} object of which the \code{colData}
#'  slot contains the \strong{fastq_path} column with paths to input
#'  cell-specific FASTQ files.
#' @param index Path to the \code{Rsubread} index of the reference genome. For
#'  generation of Rsubread indices, please refer to \code{buildindex} function
#'  in \code{Rsubread} package.
#' @param unique Argument passed to \code{align} function in \code{Rsubread}
#'  package. Boolean indicating if only uniquely mapped reads should be
#'  reported. A uniquely mapped read has one single mapping location that has
#'  less mis-matched bases than any other candidate locations. If set to
#'  \strong{FALSE}, multi-mapping reads will be reported in addition to
#'  uniquely mapped reads. Number of alignments reported for each multi-mapping
#'  read is determined by the nBestLocations parameter.
#'  Default is \strong{FALSE}.
#' @param nBestLocations Argument passed to \code{align} function in
#'  \code{Rsubread} package. Numeric value specifying the maximal number of
#'  equally-best mapping locations that will be reported for a multi-mapping
#'  read. 1 by default. The allowed value is between 1 to 16 (inclusive).
#'  In the mapping output, "NH" tag is used to indicate how many alignments are
#'  reported for the read and "HI" tag is used for numbering the alignments
#'  reported for the same read. This argument is only applicable when unique
#'  option is set to \strong{FALSE}. \code{Scruff} package does not support
#'  counting alignment files with \code{nBestLocations > 1}.
#' @param format File format of sequence alignment results. \strong{"BAM"} or
#'  \strong{"SAM"}. Default is \strong{"BAM"}.
#' @param outDir Output directory for alignment results. Sequence alignment
#'  files will be stored in folders in this directory, respectively.
#'  \strong{Make sure the folder is empty.} Default is \code{"./Alignment"}.
#' @param cores Number of cores used for parallelization. Default is
#'  \code{max(1, parallel::detectCores() - 2)}, i.e. the number of available
#'  cores minus 2.
#' @param threads \strong{Do not change}. Number of threads/CPUs used for
#'  mapping for each core. Refer to \code{align} function in \code{Rsubread}
#'  for details. Default is \strong{1}. It should not be changed in most cases.
#' @param summaryPrefix Prefix for alignment summary filename. Default is
#'  \code{"alignment"}.
#' @param overwrite Boolean indicating whether to overwrite the output
#'  directory. Default is \strong{FALSE}.
#' @param verbose Boolean indicating whether to print log messages. Useful for
#'  debugging. Default to \strong{FALSE}.
#' @param logfilePrefix Prefix for log file. Default is current date and time
#'  in the format of \code{format(Sys.time(), "\%Y\%m\%d_\%H\%M\%S")}.
#' @param ... Additional arguments passed to the \code{align} function in
#'  \code{Rsubread} package.
#' @return A \strong{SingleCellExperiment} object containing the alignment
#'  summary information in the \code{colData} slot. The \code{alignment_path}
#'  column of the annotation table contains the paths to output alignment files.
#' @examples
#' # The SingleCellExperiment object returned by demultiplex function is
#' # required for running alignRsubread function
#'
#' \dontrun{
#' data(barcodeExample, package = "scruff")
#' fastqs <- list.files(system.file("extdata", package = "scruff"),
#'     pattern = "\\.fastq\\.gz", full.names = TRUE)
#'
#' de <- demultiplex(
#'     project = "example",
#'     experiment = c("1h1"),
#'     lane = c("L001"),
#'     read1Path = c(fastqs[1]),
#'     read2Path = c(fastqs[2]),
#'     barcodeExample,
#'     bcStart = 1,
#'     bcStop = 8,
#'     umiStart = 9,
#'     umiStop = 12,
#'     keep = 75,
#'     overwrite = TRUE)
#'
#' # Alignment
#' library(Rsubread)
#' # Create index files for GRCm38_MT.
#' fasta <- system.file("extdata", "GRCm38_MT.fa", package = "scruff")
#' # Specify the basename for Rsubread index
#' indexBase <- "GRCm38_MT"
#' buildindex(basename = indexBase, reference = fasta, indexSplit = FALSE)
#'
#' al <- alignRsubread(de, indexBase, overwrite = TRUE)
#' }
#' @import data.table
#' @importFrom plyr rbind.fill
#' @export
alignRsubread <- function(sce,
    index,
    unique = FALSE,
    nBestLocations = 1,
    format = "BAM",
    outDir = "./Alignment",
    cores = max(1, parallel::detectCores() - 2),
    threads = 1,
    summaryPrefix = "alignment",
    overwrite = FALSE,
    verbose = FALSE,
    logfilePrefix = format(Sys.time(),
        "%Y%m%d_%H%M%S"),
    ...) {

    .checkCores(cores)

    fastqPaths <- SummarizedExperiment::colData(sce)$fastq_path

    if (!all(file.exists(fastqPaths))) {
        stop("Partial or all FASTQ files nonexistent.",
            " Please check paths are correct.\n",
            fastqPaths)
    }

    message(Sys.time(), " Start alignment ...")
    print(match.call(expand.dots = TRUE))

    isWindows <- .Platform$OS.type == "windows"

    logfile <- paste0(logfilePrefix, "_alignment_log.txt")

    if (verbose) {
        .logMessages(Sys.time(),
            "... Start alignment",
            logfile = logfile,
            append = FALSE)
        .logMessages(Sys.time(), fastqPaths, logfile = logfile, append = TRUE)
        message("... Input fastq paths:")
        print(fastqPaths)
    }

    if (overwrite) {
        # delete results from previous run
        message(Sys.time(), " ... Delete (if any) existing alignment results")
        unlink(file.path(outDir), recursive = TRUE)
    } else {
        alignmentPaths <- .getAlignmentFilePaths(fastqPaths, format, outDir)
        if (any(file.exists(alignmentPaths))) {
            message(
                "Abort. ",
                alignmentPaths[which(file.exists(alignmentPaths) == TRUE)],
                " already exists in output directory",
                file.path(outDir)
            )
            stop("Abort. Try re-running the function",
                " by setting overwrite to TRUE")
        }
    }

    message(Sys.time(), " ... Creating output directory ", outDir)
    dir.create(file.path(outDir), showWarnings = FALSE, recursive = TRUE)

    message(Sys.time(), " ... Mapping")
    # parallelization BiocParallel
    if (isWindows) {
        # Windows
        if (verbose) {
            resL <- BiocParallel::bplapply(
                X = fastqPaths,
                FUN = .alignRsubreadUnit,
                BPPARAM = BiocParallel::SnowParam(
                    workers = cores),
                index,
                unique,
                nBestLocations,
                format,
                outDir,
                threads,
                logfile,
                ...)
        } else {
            invisible(capture.output(resL <-
                    BiocParallel::bplapply(X = fastqPaths,
                        FUN = .alignRsubreadUnit,
                        BPPARAM = BiocParallel::SnowParam(
                            workers = cores),
                        index,
                        unique,
                        nBestLocations,
                        format,
                        outDir,
                        threads,
                        logfile = NULL,
                        ...), type = "message"))
        }
    } else {
        # Linux or macOS
        if (verbose) {
            resL <- BiocParallel::bplapply(
                X = fastqPaths,
                FUN = .alignRsubreadUnit,
                BPPARAM = BiocParallel::MulticoreParam(workers = cores),
                index,
                unique,
                nBestLocations,
                format,
                outDir,
                threads,
                logfile,
                ...)
        } else {
            invisible(capture.output(resL <-
                    BiocParallel::bplapply(X = fastqPaths,
                        FUN = .alignRsubreadUnit,
                        BPPARAM = BiocParallel::MulticoreParam(
                            workers = cores),
                        index,
                        unique,
                        nBestLocations,
                        format,
                        outDir,
                        threads,
                        logfile = NULL,
                        ...), type = "message"))
        }
    }

    resDt <- data.table::as.data.table(plyr::rbind.fill(resL))

    message(Sys.time(),
        " ... Write alignment summary to ",
        file.path(outDir, paste0(
            format(Sys.time(), "%Y%m%d_%H%M%S"),
            "_",
            summaryPrefix,
            ".tsv"
        ))
    )

    fwrite(resDt, file.path(outDir, paste0(
        format(Sys.time(), "%Y%m%d_%H%M%S"),
        "_",
        summaryPrefix,
        ".tsv"
    )), sep = "\t")

    colnames(resDt) <- c("alignment_path",
        "reads",
        "mapped_reads_incl_ercc",
        "uniquely_mapped_reads_incl_ercc",
        "multi_mapping_reads_incl_ercc",
        "Unmapped_reads",
        "indels",
        "fraction_aligned")

    message(Sys.time(), " ... Add alignment information to SCE colData.")
    SummarizedExperiment::colData(sce) <-
        cbind(SummarizedExperiment::colData(sce), resDt[, -"reads"])

    message(Sys.time(), " ... Alignment done!")
    return(sce)
}


.alignRsubreadUnit <- function(i,
    index,
    unique,
    nBestLocations,
    format,
    outDir,
    threads,
    logfile,
    ...) {

    filePath <- .getAlignmentFilePaths(i, format, outDir)

    if (file.size(i) == 0) {
        file.create(filePath, showWarnings = FALSE)
        rest <- data.frame(alignmentFilePaths = filePath,
            Total_reads = 0,
            Mapped_reads = 0,
            Uniquely_mapped_reads = 0,
            Multi_mapping_reads = 0,
            Unmapped_reads = 0,
            Indels = 0,
            PropMapped = 0,
            row.names = basename(filePath))
        return(rest)
    } else {
        .logMessages(Sys.time(),
            "... mapping sample",
            i,
            logfile = logfile,
            append = TRUE)

        res <- Rsubread::align(
            index = index,
            readfile1 = i,
            unique = unique,
            nBestLocations = nBestLocations,
            nthreads = threads,
            output_format = format,
            output_file = filePath,
            ...)

        rest <- as.data.frame(t(res))
        rest$PropMapped <- rest$Mapped_reads / rest$Total_reads
        rest$alignmentFilePaths <- filePath
        rest <- rest[c(colnames(rest)[ncol(rest)],
            colnames(rest)[-ncol(rest)])]
        return(rest)
    }
}
