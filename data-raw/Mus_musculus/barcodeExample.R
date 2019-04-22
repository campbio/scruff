
barcodeExample <- as.character(read.table(
    "data-raw/Mus_musculus/GSE85755_CEL-Seq_barcodes_96.txt")[,2])

devtools::use_data(barcodeExample, overwrite = TRUE)

