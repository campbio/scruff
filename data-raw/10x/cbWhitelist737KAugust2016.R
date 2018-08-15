
validCb <- read.table("data-raw/10x/737K-august-2016.txt")

devtools::use_data(validCb, overwrite = TRUE, compress = "xz")

