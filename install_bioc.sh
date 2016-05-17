#!/bin/bash

Rscript -e 'install.packages("devtools"); library(devtools);  source("https://bioconductor.org/biocLite.R"); biocLite("BiocInstaller"); biocLite("BSgenome.Hsapiens.UCSC.hg19"); biocLite("GenomicRanges")'
