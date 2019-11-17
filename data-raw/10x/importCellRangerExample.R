
library(data.table)

fe <- fread("./hgmm_1k_v3_20x20/outs/filtered_feature_bc_matrix/features.tsv.gz",
    header = FALSE)
cb <- fread("./hgmm_1k_v3_20x20/outs/filtered_feature_bc_matrix/barcodes.tsv.gz",
    header = FALSE)
ma <- Matrix::readMM(gzfile("./hgmm_1k_v3_20x20/outs/filtered_feature_bc_matrix/matrix.mtx.gz"))


mat <- as.matrix(ma)
colnames(mat) <- cb[["V1"]]
rownames(mat) <- fe[["V1"]]

# get top 10 expressed human genes and mouse genes
hg19ind <- grep("^hg19_", fe[["V1"]])
mm10ind <- grep("^mm10_", fe[["V1"]])

hg19rs <- rowSums(mat[hg19ind, ])
mm10rs <- rowSums(mat[mm10ind, ])

hg19rstop10 <- sort(hg19rs, decreasing = TRUE)[1:10]
mm10rstop10 <- sort(mm10rs, decreasing = TRUE)[1:10]

fesub <- fe[V1 %in% names(hg19rstop10) | V1 %in% names(mm10rstop10), ]
cbsub <- cb[1:20, ]

matsub <- mat[fesub[["V1"]], 1:20]
matsubma <- as(matsub, "dgTMatrix")

fwrite(fesub, file = "./data-raw/10x/hgmm_1k_v3_20x20/outs/filtered_feature_bc_matrix/features.tsv",
    sep = "\t", col.names = FALSE)
fwrite(cbsub, file = "./data-raw/10x/hgmm_1k_v3_20x20/outs/filtered_feature_bc_matrix/barcodes.tsv",
    sep = "\t", col.names = FALSE)
Matrix::writeMM(matsubma,
    file = "./data-raw/10x/hgmm_1k_v3_20x20/outs/filtered_feature_bc_matrix/matrix.mtx")

R.utils::gzip("./data-raw/10x/hgmm_1k_v3_20x20/outs/filtered_feature_bc_matrix/features.tsv")
R.utils::gzip("./data-raw/10x/hgmm_1k_v3_20x20/outs/filtered_feature_bc_matrix/barcodes.tsv")
R.utils::gzip("./data-raw/10x/hgmm_1k_v3_20x20/outs/filtered_feature_bc_matrix/matrix.mtx")
