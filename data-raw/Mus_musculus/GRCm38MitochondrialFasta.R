ref = "data-raw/Mus_musculus/GCA_000001635.8_GRCm38.p6_genomic_modified_MT.fna"

GRCm38MitochondrialFasta = data.table::fread(ref, sep = "\t")

devtools::use_data(GRCm38MitochondrialFasta, overwrite = TRUE)
