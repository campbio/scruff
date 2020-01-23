
library(ShortRead)


v1h1R1 <- "../scruff/inst/extdata/vandenBrink_1h1_L001_R1_001.fastq.gz"
v1h1R2 <- "../scruff/inst/extdata/vandenBrink_1h1_L001_R2_001.fastq.gz"
vb1R1 <- "../scruff/inst/extdata/vandenBrink_b1_L001_R1_001.fastq.gz"
vb1R2 <- "../scruff/inst/extdata/vandenBrink_b1_L001_R2_001.fastq.gz"

v1r1fastq <- readFastq(v1h1R1)
v1r2fastq <- readFastq(v1h1R2)
vb1fastq <- readFastq(vb1R1)
vb2fastq <- readFastq(vb1R2)


ShortRead::writeFastq(v1r1fastq[seq(2000)], "vandenBrink_1h1_L001_R1_001.fastq.gz")
ShortRead::writeFastq(v1r2fastq[seq(2000)], "vandenBrink_1h1_L001_R2_001.fastq.gz")
ShortRead::writeFastq(vb1fastq[seq(2000)], "vandenBrink_b1_L001_R1_001.fastq.gz")
ShortRead::writeFastq(vb2fastq[seq(2000)], "vandenBrink_b1_L001_R2_001.fastq.gz")
