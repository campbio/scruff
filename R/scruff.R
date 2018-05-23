#' Run scruff pipeline
#'
#' Run the \code{scruff} pipeline. This function performs all \code{demultiplex}, \code{alignRsubread}, and \code{countUMI} functions. Write demultiplex statistics, alignment statistics, and UMI filtered count matrix in output directories. Return a SingleCellExperiment object containing the count matrix, cell and gene annotations, and all QC metrics.
#'
#' @param project The project name. Default is \code{paste0("project_", Sys.Date())}.
#' @param sample A character vector of sample names. Represents the group label for each FASTQ file, e.g. "patient1, patient2, ...".
#' @param lane A character or character vector of flow cell lane numbers. If FASTQ files from multiple lanes are concatenated, any placeholder would be sufficient, e.g. "L001".
#' @param read1Path A character vector of file paths to the read1 FASTQ files. These are the read files with UMI and cell barcode information.
#' @param read2Path A character vector of file paths to the read2 FASTQ files. These read files contain genomic sequences.
#' @param bc A vector of pre-determined cell barcodes. For example, see \code{?barcodeExample}.
#' @param index Path to the \code{Rsubread} index of the reference genome. For generation of Rsubread indices, please refer to \code{buildindex} function in \code{Rsubread} package.
#' @param reference Path to the reference GTF file. The TxDb object of the GTF file will be generated and saved in the current working directory with ".sqlite" suffix.
#' @param bcStart Integer or vector of integers containing the cell barcode start positions (inclusive, one-based numbering).
#' @param bcStop Integer or vector of integers containing the cell barcode stop positions (inclusive, one-based numbering).
#' @param bcEdit Maximally allowed edit distance for barcode correction. Barcodes with mismatches equal or fewer than this will be assigned a corrected barcode if the inferred barcode matches uniquely in the provided predetermined barcode list. Default is 0, meaning no cell barcode correction is performed.
#' @param umiStart Integer or vector of integers containing the start positions (inclusive, one-based numbering) of UMI sequences.
#' @param umiStop Integer or vector of integers containing the stop positions (inclusive, one-based numbering) of UMI sequences.
#' @param keep Read trimming. Read length or number of nucleotides to keep for read 2 (the read that contains transcript sequence information). Longer reads will be clipped at 3' end. Shorter reads will not be affected. This number should be determined based on the sequencing kit that was used in library preparation step.
#' @param cellPerWell Number of cells per well. Can be an integer (e.g. 1) indicating the number of cells in each well or an vector with length equal to the total number of cells in the input alignment files specifying the number of cells in each file. Default is 1.
#' @param unique Argument passed to \code{align} function in \code{Rsubread} package. Boolean indicating if only uniquely mapped reads should be reported. A uniquely mapped read has one single mapping location that has less mis-matched bases than any other candidate locations. If set to \strong{FALSE}, multi-mapping reads will be reported in addition to uniquely mapped reads. Number of alignments reported for each multi-mapping read is determined by the nBestLocations parameter. Default is \strong{TRUE}.
#' @param nBestLocations Argument passed to \code{align} function in \code{Rsubread} package. Numeric value specifying the maximal number of equally-best mapping locations that will be reported for a multi-mapping read. 1 by default. The allowed value is between 1 to 16 (inclusive). In the mapping output, "NH" tag is used to indicate how many alignments are reported for the read and "HI" tag is used for numbering the alignments reported for the same read. This argument is only applicable when unique option is set to \strong{FALSE}.
#' @param minQual Minimally acceptable Phred quality score for cell barcode and UMI sequences. Phread quality scores are calculated for each nucleotide in these tags. Tags with at least one nucleotide with score lower than this will be filtered out. Default is \strong{10}.
#' @param yieldReads The number of reads to yield when drawing successive subsets from a fastq file, providing the number of successive records to be returned on each yield. This parameter is passed to the \code{n} argument of the \code{FastqStreamer} function in \emph{ShortRead} package. Default is \strong{1e06}.
#' @param alignmentFileFormat File format of sequence alignment results. \strong{"BAM"} or \strong{"SAM"}. Default is \strong{"BAM"}.
#' @param demultiplexOutDir Output folder path for demultiplex results. Demultiplexed cell specifc FASTQ files will be stored in folders in this path, respectively. \strong{Make sure the folder is empty.} Default is \code{"./Demultiplex"}.
#' @param alignmentOutDir Output directory for alignment results. Sequence alignment maps will be stored in folders in this directory, respectively. \strong{Make sure the folder is empty.} Default is \code{"./Alignment"}.
#' @param countUmiOutDir Output directory for UMI counting results. UMI filtered count matrix will be stored in this directory. Default is \code{"./Count"}.
#' @param demultiplexSummaryPrefix Prefix for demultiplex summary filename. Default is \code{"demultiplex"}.
#' @param alignmentSummaryPrefix Prefix for alignment summary filename. Default is \code{"alignment"}.
#' @param countPrefix Prefix for UMI filtered count matrix filename. Default is \code{"countUMI"}.
#' @param logfilePrefix Prefix for log file. Default is current date and time in the format of \code{format(Sys.time(), "\%Y\%m\%d_\%H\%M\%S")}.
#' @param overwrite Boolean indicating whether to overwrite the output directory. Default is \strong{FALSE}.
#' @param verbose Boolean indicating whether to print log messages. Useful for debugging. Default to \strong{FALSE}.
#' @param cores Number of cores to use for parallelization. Default is \code{max(1, parallel::detectCores() / 2)}, i.e. the number of available cores divided by 2.
#' @param threads \strong{Do not change}. Number of threads/CPUs used for mapping for each core. Refer to \code{align} function in \code{Rsubread} for details. Default is \strong{1}. It should not be changed in most cases.
#' @param ... Additional arguments passed to the \code{align} function in \code{Rsubread} package.
#' @return A \code{SingleCellExperiment} object.
#' @export
scruff <- function(project = paste0("project_", Sys.Date()),
                   sample,
                   lane,
                   read1Path,
                   read2Path,
                   bc,
                   index,
                   reference,
                   bcStart,
                   bcStop,
                   bcEdit = 0,
                   umiStart,
                   umiStop,
                   keep,
                   cellPerWell = 1,
                   unique = TRUE,
                   nBestLocations = 1,
                   minQual = 10,
                   yieldReads = 1e06,
                   alignmentFileFormat = "BAM",
                   demultiplexOutDir = "./Demultiplex",
                   alignmentOutDir = "./Alignment",
                   countUmiOutDir = "./Count",
                   demultiplexSummaryPrefix = "demultiplex",
                   alignmentSummaryPrefix = "alignment",
                   countPrefix = "countUMI",
                   logfilePrefix = format(Sys.time(), "%Y%m%d_%H%M%S"),
                   overwrite = FALSE,
                   verbose = FALSE,
                   cores = max(1, parallel::detectCores() / 2),
                   threads = 1,
                   ...) {

  # run pipeline
  message(paste(Sys.time(), "Start running scruff ..."))
  print(match.call(expand.dots = TRUE))

  de <- demultiplex(
    project = project,
    sample = sample,
    lane = lane,
    read1Path = read1Path,
    read2Path = read2Path,
    bc = bc,
    bcStart = bcStart,
    bcStop = bcStop,
    bcEdit = bcEdit,
    umiStart = umiStart,
    umiStop = umiStop,
    keep = keep,
    minQual = minQual,
    yieldReads = yieldReads,
    outDir = demultiplexOutDir,
    summaryPrefix = demultiplexSummaryPrefix,
    overwrite = overwrite,
    verbose = verbose,
    cores = cores,
    logfilePrefix = logfilePrefix
  )

  al <- alignRsubread(
    sce = de,
    index = index,
    unique = unique,
    nBestLocations = nBestLocations,
    format = alignmentFileFormat,
    outDir = alignmentOutDir,
    cores = cores,
    threads = threads,
    summaryPrefix = alignmentSummaryPrefix,
    overwrite = overwrite,
    verbose = verbose,
    logfilePrefix = logfilePrefix
  )

  co <- countUMI(
    sce = al,
    reference = reference,
    format = alignmentFileFormat,
    outDir = countUmiOutDir,
    cellPerWell = cellPerWell,
    cores = cores,
    outputPrefix = countPrefix,
    verbose = verbose,
    logfilePrefix = logfilePrefix
  )

  message(paste(Sys.time(), "Finished running scruff ..."))

  return (co)
}



