#!/bin/bash

echo "Rscript -e source(https://bioconductor.org/biocLite.R); biocLite(BiocInstaller)"
Rscript -e 'source("https://bioconductor.org/biocLite.R"); biocLite("BiocInstaller"); library(BSgenome.Hsapiens.UCSC.hg19); library(GenomicRanges);'
