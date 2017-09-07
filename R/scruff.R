#!/usr/bin/env Rscript

## scRNA-seq pipeline
## fastq.gz to count matrix

#' Run Scruff pipeline
#' 
#' This fuction runs the whole Scruff pipeline
#' 
#' @param dir directory to input fastq files.
#' @export
#' scruff
scruff <- function() {
  source("demultiplex.R")
  source("align.R")
  source("count_umi.R")
  
  # required packages
  suppressPackageStartupMessages(library(parallel))
  suppressPackageStartupMessages(library(data.table))
  suppressPackageStartupMessages(library(ShortRead))
  suppressPackageStartupMessages(library(gtools))
  suppressPackageStartupMessages(library(Rsubread))
  suppressPackageStartupMessages(library(data.table))
  suppressPackageStartupMessages(library(GenomicAlignments))
  suppressPackageStartupMessages(library(GenomicFeatures))
  suppressPackageStartupMessages(library(Rsamtools))
  
  # required fields
  bc.index.file <- "barcodes.txt"
  input.dir <- "/restricted/projectnb/pulmseq/fastq/170802_NB500996_0083_AH5FC2BGX3"
  
  cut.length <- 50
  min.bc.quality <- 10
  output.dir <- ".."
  umi.length <- 5
  bc.length <- 6
  fname.delimiter <- "_"
  out.folder <- "Demultiplex"
  stats.out <- "demultiplex_stats"
  out.format <- "BAM"
  align.dir <- "../Alignment"
  count.out <- "../Count"
  
  # alignment
  GRCh38.index <- "/restricted/projectnb/cbmhive/references/Homo_Sapiens/GRCh38/Rsubread_index/GRCh38"
  gtf.file <- "/restricted/projectnb/cbmhive/references/Homo_Sapiens/GRCh38/gtf/Homo_sapiens.GRCh38.89.chr_ercc.gtf"
  nthreads <- 16
  
  # cores
  mc.cores <- 16
  
  # run pipeline
  demultiplex.wrapper(bc.index.file, input.dir, stats.out, output.dir, out.folder, min.bc.quality,
                      umi.length, bc.length, cut.length, fname.delimiter, mc.cores)
  fastq.dir <- file.path(output.dir, out.folder)
  align.wrapper(fastq.dir, GRCh38.index, out.format, align.dir, nthreads, mc.cores)
  count.wrapper(alignment.dir = align.dir, gtf.file, if.bam = T, count.out, mc.cores)
}

# run pipeline
#scruff()

