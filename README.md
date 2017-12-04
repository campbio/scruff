# scruff: Single Cell RNA-Seq UMI Filtering Facilitator (under development)

**scruff** (**S**ingle **C**ell **R**NA-Seq **U**MI **F**iltering **F**acilitator) is a package for processing single cell RNA-seq raw fastq files. Demultiplex single cell RNA-seq fastq files, align reads to reference genome using Rsubread, and generate transcript count matrices with UMI filtering.

## Installation Instructions

To install `scruff` using `devtools`:
```
library(devtools)
install_github("compbiomed/scruff")
```

## Package user guide

An introduction to `scruff` package is available within the package as a vignette which can be opened via running 
```
vignette("scruff")
```
## Example QC plots

The following figures are generated using data from [Van den Brink et al. 2017](https://www.nature.com/articles/nmeth.4437).

![Alt text](/data-raw/figure/to/20171204_qc_reads_assigned_excl_bulk_Page_8.png?raw=true "Optional Title")



