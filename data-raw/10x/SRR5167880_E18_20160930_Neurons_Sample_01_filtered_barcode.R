
library(data.table)
library(TENxBrainData)
tenx <- TENxBrainData()

bc <- colData(tenx)[1:9970, ][["Barcode"]]
write.table(bc,
    "SRR5167880_E18_20160930_Neurons_Sample_01_filtered_barcode.tsv",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE)
