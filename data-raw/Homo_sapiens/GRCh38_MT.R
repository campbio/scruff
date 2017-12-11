library(data.table)
ref = "data-raw/gtf/Homo_sapiens.GRCh38.89.chr_MT.gtf"

GRCh38_MT_gtf = fread(ref, sep = "\t")

devtools::use_data(GRCh38_MT_gtf, overwrite = TRUE)
