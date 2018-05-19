fdir <- list.files("data-raw/Mus_musculus/vandenBrink_fastq", full.names=TRUE)

fastqExample <- vector("list", length(fdir))

for (i in seq_len(length(fdir))) {
  fastqExample[[i]] <- ShortRead::readFastq(fdir[i])
}

names(fastqExample) <- basename(fdir)

devtools::use_data(fastqExample, overwrite = TRUE, compress="xz")

