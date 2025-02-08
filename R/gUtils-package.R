#' gUtils - Extending GenomicRanges
#' 
#' @import GenomicRanges
#' @importFrom GenomeInfoDb Seqinfo seqinfo keepSeqlevels seqlevels seqlengths seqlevels<- seqlengths<- genome<- seqnames seqlengths seqlengths<- seqinfo seqinfo<- seqnames seqnames<- seqlevels seqlevels<- seqlevelsStyle seqlevelsStyle<-
#' @importFrom utils read.delim relist
#' @importFrom stats setNames
#' @importFrom methods as is setMethod
#' @importFrom IRanges IRanges
#' @importFrom BiocGenerics unlist
#' @importFrom data.table data.table rbindlist is.data.table := setkeyv as.data.table
#' @importFrom IRanges findOverlaps overlapsRanges
#' @importFrom parallel mclapply
#' @importFrom S4Vectors queryHits subjectHits
#' @importMethodsFrom S4Vectors elementNROWS Rle mcols mcols<- values values<- elementMetadata elementMetadata<- from to 
#' @importMethodsFrom BiocGenerics width
"_PACKAGE"