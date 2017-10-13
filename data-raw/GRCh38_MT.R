ref = "data-raw/rsubread_index/Homo_sapiens.GRCh38.dna.chromosome.MT.fa"

# buildindex(basename="GRCh38_MT", reference=ref, indexSplit=F, memory=16000)

GRCh38_MT = scan(ref, what = "character", sep = "\t")

devtools::use_data(GRCh38_MT, overwrite = TRUE)

