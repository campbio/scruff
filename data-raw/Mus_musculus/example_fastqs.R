fdir <- list.files("data-raw/Mus_musculus/vandenBrink_fastq", full.names=TRUE)

examplefastq <- vector("list", length(fdir))

for (i in seq_len(length(fdir))) {
  examplefastq[i] <- ShortRead::readFastq(fdir[i])
}

names(examplefastq) <- basename(fdir)

devtools::use_data(examplefastq, overwrite = TRUE, compress="xz")

