
exampleannot <- data.table::data.table(project="Example",
                    cohort=c(rep("01", 2), rep("02", 2)),
                    num=c(rep("S1", 2), rep("S2", 2)),
                    lane=c(rep("L001", 4)),
                    read=rep(paste0("R", 1:2), 2),
                    dir=c("Example-01_S1_L001_R1_001.fastq.gz",
                           "Example-01_S1_L001_R2_001.fastq.gz",
                           "Example-02_S2_L001_R1_001.fastq.gz",
                           "Example-02_S2_L001_R2_001.fastq.gz")
                    )

devtools::use_data(exampleannot, overwrite = TRUE)

