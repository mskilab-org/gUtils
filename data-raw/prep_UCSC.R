library(gUtils)

## prepare a Seqinfo object
sl <- hg_seqlengths(BSgenome.Hsapiens.UCSC.hg19::Hsapiens)
si <- Seqinfo(seqnames=names(sl), seqlengths=sl, genome="hg19")
isCircular(si) <- rep(FALSE, length(sl))

## download data from UCSC Table Browser to include as example data in gUtils
library (rtracklayer)
mySession = browserSession("UCSC")
genome(mySession) <- "hg19"

## DNAaseI hypersensitivity sites
track.name <- "wgEncodeUwDgf"
table.name <- "wgEncodeUwDgfK562Hotspots"
tbl.DNAase <- getTable( ucscTableQuery(mySession, track=track.name, table=table.name))
example_dnase <- with(tbl.DNAase, GRanges(gsub("chr(.*?)", "\\1", chrom), IRanges(chromStart, chromEnd), seqinfo=si))
mcols(example_dnase) <- tbl.DNAase[,-which(colnames(tbl.DNAase) %in% c("chrom", "chromStart","chromEnd","strand","name","qValue","bin","score"))]

## RefSeq genes
track.name <- "RefSeq Genes"
table.name <- "refGene"
tbl.genes <- getTable( ucscTableQuery(mySession, track=track.name, table=table.name))
example_genes <- with(tbl.genes[ix <- !grepl("(hap)|(gl)", tbl.genes$chrom),], GRanges(gsub("chr(.*?)", "\\1", chrom), IRanges(cdsStart, cdsEnd), strand=strand, seqinfo=si))
mcols(example_genes) <- tbl.genes[ix,-which(colnames(tbl.genes) %in% c("chrom","strand","exonStarts","exonEnds","exonFrames"))]
mcols(example_genes)$name <- as.character(mcols(example_genes)$name2)
mcols(example_genes)[, c("cdsEnd","cdsStart","cdsEndStat","cdsStartStat","score","txStart","txEnd","bin","name2")] <- NULL
example_genes <- sort(example_genes[!duplicated(example_genes$name) & width(example_genes) < 1e6 & width(example_genes) > 10])
example_genes <- gr.fix(example_genes, si)

## import some HiC data
library(LiebermanAidenHiC2009)
data("HiC_GM_chr14")
chain <- import.chain("/Users/jwala/hg18ToHg19.over.chain")
gr1 <- with(HiC_GM_chr14,GRanges(paste("chr",c(chromosome1,chromosome2),sep=""), IRanges(c(position1, position2),width=1), strand=ifelse(c(strand1, strand2)==0,"+","-"), id=c(seq_along(strand1), seq_along(strand1))))
gr1.hg19 <- gr.fix(gr.nochr(unlist(liftOver(gr1, chain))), si)
tab <- table(gr1.hg19$id)
gr1.hg19 <- gr1.hg19[as.character(gr1.hg19$id) %in% names(tab[tab==2])]
grl.hiC <- GenomicRanges::split(gr1.hg19, gr1.hg19$id)

## make some fake rearrangement data
set.seed(137)
grg <- example_genes
mcols(grg) <- NULL
strand(grg) <- ifelse(runif(length(grg)) > 0.5, "+", "-")
mcols(grg)$bin <- sample(rep(seq(length(grg)/2), each=2))
grg <- grg[grg$bin %in% sample(unique(grg$bin), 500)]
grl <- S4Vectors::split(gr.start(grg), grg$bin)
sam <- sample(length(grl), length(grl), replace=FALSE)
grl1 <- grl[1:(length(grl)/2)]
grl2 <- grl[(length(grl)/2):length(grl)]

## save it
example_genes <- sort(example_genes)
example_dnase <- sort(example_dnase[sample(length(example_dnase), 10000)])
grl.hiC <- grl.hiC[sample(length(grl.hiC), 10000)]
save(example_genes, example_dnase, si, grl1, grl2, grl.hiC, file="data/grdata.rda", compress='xz')
