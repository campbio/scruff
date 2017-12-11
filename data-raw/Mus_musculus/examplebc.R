
examplebc <- as.character(read.table(
  "data-raw/Mus_musculus/GSE85755_Cel-seq_barcodes_96.txt")[,2])

devtools::use_data(examplebc, overwrite = TRUE)

