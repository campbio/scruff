#' Run scuff pipeline
#' 
#' scuff runs the whole pipeline including demultiplex, align, and count_umi. Write tables of count matrices in output directory.
#' 
#' @param fastqs An annotation data table or data frame of input fastq files.
#' @param barcodes A vector of cell barcodes known from experimental design.
#' @param index Rsubread index for reference sequences.
#' @param annot Directory of the gene annotation file.
#' @param length Read length. Longer reads will be clipped.
#' @param umi.pos Start and end index number of umi sequences.
#' @param bc.pos Start and end index number of barcodes.
#' @param bc.qual Minimal acceptable quality score for barcode sequence.
#' @param out Output directory.
#' @param overwrite Overwrite the output directory. Default is TRUE.
#' @param ncore Number of cores to use.
#' @param nthreads Number of threads to run for each core.
#' @export
scuff <- function(fastqs, barcodes, index, annot, length, umi.pos, bc.pos, bc.qual=10, out="./",
                  overwrite=TRUE, ncore=1, nthreads=1) {
  # run pipeline
  demultiplex.wrapper(bc.index.file, input.dir, stats.out, output.dir, out.folder, min.bc.quality,
                      umi.length, bc.length, cut.length, fname.delimiter, mc.cores)
  fastq.dir <- file.path(output.dir, out.folder)
  align.wrapper(fastq.dir, GRCh38.index, out.format, align.dir, nthreads, mc.cores)
  count.wrapper(alignment.dir = align.dir, gtf.file, if.bam = T, count.out, mc.cores)
}



