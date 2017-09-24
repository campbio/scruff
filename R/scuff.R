#' Run scuff pipeline
#' 
#' scuff runs the whole pipeline including demultiplex, align, and count_umi. Write tables of count matrices in output directory.
#' 
#' @param fastq An annotation data table or data frame of input fastq files.
#' @param index Rsubread index for reference sequences.
#' @param bc A vector of cell barcodes known from experimental design.
#' @param annot Directory of the gene annotation file.
#' @param clip.at Read length or number of nucleotides to keep. Longer reads will be clipped at 3' end.
#' @param umi.pos Start and end index number of umi sequences.
#' @param bc.pos Start and end index number of barcodes.
#' @param bc.qual Minimal acceptable quality score for barcode sequence.
#' @param out Output directory.
#' @param overwrite Overwrite the output directory. Default is TRUE.
#' @param ncore Number of cores to use.
#' @param threads Number of threads to run for each core. Default is \strong{16}.
#' @export
scuff <- function(fastq, barcodes, index, annot, length, umi.pos, bc.pos, bc.qual=10, out="./",
                  overwrite=TRUE, ncore=1, nthreads=1) {
  # run pipeline
  message("Starting scuff ...")
  
  demultiplex.wrapper(bc.index.file, input.dir, stats.out, output.dir, out.folder, min.bc.quality,
                      umi.length, bc.length, cut.length, fname.delimiter, mc.cores)
  fastq.dir <- file.path(output.dir, out.folder)
  align.wrapper(fastq.dir, GRCh38.index, out.format, align.dir, nthreads, mc.cores)
  count.wrapper(alignment.dir = align.dir, gtf.file, if.bam = T, count.out, mc.cores)
}



