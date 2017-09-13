
examplebc <- data.table::fread("data-raw/example_barcodes.txt", col.names=c("index", "barcode"))
examplebc <- examplebc$barcode

devtools::use_data(examplebc, overwrite = TRUE)

