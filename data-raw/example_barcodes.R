library(data.table)

examplebc <- fread("data-raw/example_barcodes.txt", col.names=c("index", "barcode"))

devtools::use_data(examplebc, overwrite = TRUE)

