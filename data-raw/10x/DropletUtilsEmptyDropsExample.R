
library(data.table)
library(Matrix)

fe <- fread("../scruff_upgrade/data/hg_1k_v3_raw_feature_bc_matrix/raw_feature_bc_matrix/features.tsv.gz",
    header = FALSE)
cb <- fread("../scruff_upgrade/data/hg_1k_v3_raw_feature_bc_matrix/raw_feature_bc_matrix/barcodes.tsv.gz",
    header = FALSE)
ma <- Matrix::readMM(gzfile("../scruff_upgrade/data/hg_1k_v3_raw_feature_bc_matrix/raw_feature_bc_matrix/matrix.mtx.gz"))


#mat <- as.matrix(ma)
colnames(ma) <- cb[["V1"]]
rownames(ma) <- fe[["V1"]]

hg19cs <- colSums(ma)
#mm10rs <- rowSums(mat[mm10ind, ])

hg19cssort <- sort(hg19cs, decreasing = TRUE)

# get the top 10 cells with most counts
hg19csfirst10 <- hg19cssort[1:10]
# get the last 10 cells with non-zero counts
hg19cssortnonzero <- hg19cssort[hg19cssort != 0]
hg19cslast10 <- hg19cssortnonzero[seq(length(hg19cssortnonzero) - 9,
    length(hg19cssortnonzero))]

fesub <- fe
cbsub <- data.table(V1 = c(names(hg19csfirst10), names(hg19cslast10)))
matsub <- ma[, c(names(hg19csfirst10), names(hg19cslast10))]

fwrite(fesub, file = "./data-raw/10x/hg_1k_v3_33538x20/outs/raw_feature_bc_matrix/features.tsv",
    sep = "\t", col.names = FALSE)
fwrite(cbsub, file = "./data-raw/10x/hg_1k_v3_33538x20/outs/raw_feature_bc_matrix/barcodes.tsv",
    sep = "\t", col.names = FALSE)
Matrix::writeMM(matsub,
    file = "./data-raw/10x/hg_1k_v3_33538x20/outs/raw_feature_bc_matrix/matrix.mtx")

R.utils::gzip("./data-raw/10x/hg_1k_v3_33538x20/outs/raw_feature_bc_matrix/features.tsv")
R.utils::gzip("./data-raw/10x/hg_1k_v3_33538x20/outs/raw_feature_bc_matrix/barcodes.tsv")
R.utils::gzip("./data-raw/10x/hg_1k_v3_33538x20/outs/raw_feature_bc_matrix/matrix.mtx")


emptyDropsSceExample <- scruff::importCellRanger(
    cellRangerDirs = "./data-raw/10x/",
    cellRangerOuts = "outs/raw_feature_bc_matrix/",
    samples = "hg_1k_v3_33538x20",
    class = "Matrix")

usethis::use_data(emptyDropsSceExample, overwrite = TRUE, compress = "xz")


