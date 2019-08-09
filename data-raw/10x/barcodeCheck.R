
validCbv1 <- read.table("data-raw/10x/737K-april-2014_rc.txt", as.is = TRUE)
validCbv2 <- read.table("data-raw/10x/737K-august-2016.txt", as.is = TRUE)
validCbv3 <- read.table("data-raw/10x/3M-february-2018.txt", as.is = TRUE)

cb1top10000 <- validCbv1[seq(10000), ]
cb2top10000 <- validCbv2[seq(10000), ]
cb3top10000 <- validCbv3[seq(10000), ]

cbtop10000 <- data.table::data.table(v1chemistry = cb1top10000,
    v2chemistry = cb2top10000,
    v3chemistry = cb3top10000)

devtools::use_data(cbtop10000, overwrite = TRUE, compress = "xz")
