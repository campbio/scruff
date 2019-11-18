
library(data.table)

fe <- fread("../scruff_upgrade/data/Solo.out_PBMC_1k_v3/Gene/filtered/features.tsv", header = FALSE)
cb <- fread("../scruff_upgrade/data/Solo.out_PBMC_1k_v3/Gene/filtered/barcodes.tsv", header = FALSE)
ma <- Matrix::readMM(gzfile("../scruff_upgrade/data/Solo.out_PBMC_1k_v3/Gene/filtered/matrix.mtx"))


mat <- as.matrix(ma)
colnames(mat) <- cb[["V1"]]
rownames(mat) <- fe[["V1"]]

# get top 20 expressed genes
rs <- rowSums(mat)
top20 <- sort(rs, decreasing = TRUE)[1:20]

fesub <- fe[V1 %in% names(top20), ]
cbsub <- cb[1:20, ]

matsub <- mat[fesub[["V1"]], 1:20]
matsubma <- as(matsub, "dgTMatrix")

fwrite(fesub, file = "./data-raw/10x/PBMC_1k_v3_20x20/Gene/filtered/features.tsv",
    sep = "\t", col.names = FALSE)
fwrite(cbsub, file = "./data-raw/10x/PBMC_1k_v3_20x20/Gene/filtered/barcodes.tsv",
    sep = "\t", col.names = FALSE)
Matrix::writeMM(matsubma,
    file = "./data-raw/10x/PBMC_1k_v3_20x20/Gene/filtered/matrix.mtx")
