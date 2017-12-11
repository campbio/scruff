library(data.table)
ref = "data-raw/Mus_musculus/Mus_musculus.GRCm38.90.MT.gtf"

GRCm38_MT_gtf = fread(ref, sep = "\t")

devtools::use_data(GRCm38_MT_gtf, overwrite = TRUE)
