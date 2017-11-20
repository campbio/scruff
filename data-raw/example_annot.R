
exampleannot <- data.table::data.table(project = "Example",
                    sample = c("01", "02"),
                    lane = c(rep("L001", 2)),
                    read1_path = c("Example-01_S1_L001_R1_001.fastq.gz",
                                   "Example-02_S2_L001_R1_001.fastq.gz"),
                    read2_path = c("Example-01_S1_L001_R2_001.fastq.gz",
                                   "Example-02_S2_L001_R2_001.fastq.gz")
                    )

devtools::use_data(exampleannot, overwrite = TRUE)

