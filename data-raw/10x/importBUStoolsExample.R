
library(data.table)
library(Matrix)

fe <- fread("../scruff_upgrade/data/BUStools_PBMC_1k_v3/genecount/genes.genes.txt", header = FALSE)
cb <- fread("../scruff_upgrade/data/BUStools_PBMC_1k_v3/genecount/genes.barcodes.txt", header = FALSE)
# BUStools (unfiltered) sparse matrix result is transposed
ma <- Matrix::readMM(gzfile("../scruff_upgrade/data/BUStools_PBMC_1k_v3/genecount/genes.mtx"))

# Only keep top 1000 cells

rorder <- order(rowSums(ma), decreasing = T)[1:1000]

mat <- as.matrix(t(ma[rorder, ]))
colnames(mat) <- cb[["V1"]][rorder]
rownames(mat) <- fe[["V1"]]

# get top 20 expressed genes
rs <- rowSums(mat)
top20 <- sort(rs, decreasing = TRUE)[1:20]

fesub <- fe[V1 %in% names(top20), ]
cbsub <- as.data.table(cb[["V1"]][rorder][1:20])

matsub <- mat[fesub[["V1"]], 1:20]
matsubma <- as(t(matsub), "dgTMatrix")

fwrite(fesub, file = "./data-raw/10x/BUStools_PBMC_1k_v3_20x20/genecount/genes.genes.txt",
    sep = "\t", col.names = FALSE)
fwrite(cbsub, file = "./data-raw/10x/BUStools_PBMC_1k_v3_20x20/genecount/genes.barcodes.txt",
    sep = "\t", col.names = FALSE)
Matrix::writeMM(matsubma,
    file = "./data-raw/10x/BUStools_PBMC_1k_v3_20x20/genecount/genes.mtx")
