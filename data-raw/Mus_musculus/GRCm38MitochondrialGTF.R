ref = "data-raw/Mus_musculus/Mus_musculus.GRCm38.90.MT.gtf"

GRCm38MitochondrialGTF = data.table::fread(ref, sep = "\t")

devtools::use_data(GRCm38MitochondrialGTF, overwrite = TRUE)
