library(testthat)
library(scruff)

# extdata exists
test_that(
    desc = "Making sure files in inst/extdata folder exist", {
        expect_true(
            file.exists(system.file("extdata",
                "vandenBrink_1h1_L001_R1_001.fastq.gz",
                package = "scruff")))
        expect_true(
            file.exists(system.file("extdata",
                "vandenBrink_1h1_L001_R2_001.fastq.gz",
                package = "scruff")))
        expect_true(
            file.exists(system.file("extdata",
                "GRCm38_MT.fa", package = "scruff")))
        expect_true(
            file.exists(system.file("extdata",
                "GRCm38_MT.gtf", package = "scruff")))
    })
