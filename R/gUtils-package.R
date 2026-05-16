#' gUtils - Extending GenomicRanges
#' 
#' @import GenomicRanges
#' @import data.table
#' @importFrom GenomeInfoDb Seqinfo keepSeqlevels seqlengths seqlengths<- genome<- seqinfo seqinfo<- seqnames seqnames<- seqlevels seqlevels<- seqlevelsStyle seqlevelsStyle<- isCircular orderSeqlevels rankSeqlevels
#' @importFrom utils read.delim relist
#' @importFrom stats setNames
#' @importFrom methods as is setMethod
#' @importFrom IRanges IRanges
#' @importFrom BiocGenerics unlist
#' @importFrom data.table data.table rbindlist is.data.table := setkeyv as.data.table
#' @importFrom IRanges findOverlaps overlapsRanges
#' @importFrom parallel mclapply
#' @importFrom S4Vectors queryHits subjectHits
#' @importMethodsFrom S4Vectors elementNROWS Rle mcols mcols<- values values<- elementMetadata elementMetadata<- from to with within
#' @importMethodsFrom BiocGenerics width
#' @importMethodsFrom IRanges relist
"_PACKAGE"