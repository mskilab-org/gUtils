library(gUtils)
library(BSgenome.Hsapiens.UCSC.hg19)
context("Range ops")

gr  <- GRanges(1, IRanges(c(3,7,13), c(5,9,16)), strand=c('+','-','-'), seqinfo=Seqinfo("1", 25), name=c("A","B","C"))
gr2 <- GRanges(1, IRanges(c(1,9), c(6,14)), strand=c('+','-'), seqinfo=Seqinfo("1", 25), field=c(1,2))
dt <- data.table(seqnames=1, start=c(2,5,10), end=c(3,8,15))

#test_that("hg_seqlengths", {
#  expect_error(hg_seqlengths())
#  ee <- structure(names="1", 249250621L)
#  expect_identical(hg_seqlengths(Hsapiens)[1],ee)
#  expect_equal(names(hg_seqlengths(Hsapiens, chr=TRUE)[1]), "chr1")
#  expect_equal(length(hg_seqlengths(Hsapiens, include.junk = TRUE)), 93)
#})

test_that("gr.start ", {

  expect_identical(start(gr.start(gr)), c(3L,7L,13L))
  expect_identical(end(gr.start(gr)),   c(3L,7L,13L))
  expect_identical(end(gr.start(gr, width=100)),   c(25L,25L,25L))
  expect_identical(suppressWarnings(end(gr.start(gr, width=100, force=TRUE))),   c(102L,106L,112L))
  expect_identical(end(gr.start(gr, width=100, ignore.strand=FALSE)),   c(25L,9L,16L))
  expect_identical(suppressWarnings(start(gr.start(gr, width=100, ignore.strand=FALSE, force=TRUE))),   c(3L,-90L,-83L))
  expect_identical(suppressWarnings(end(gr.start(gr, width=100, force=TRUE))),   c(102L,106L,112L))
  expect_identical(end(gr.start(gr, width=10000, clip=TRUE)), c(5L, 9L, 16L))

})

test_that("gr.end", {

  expect_identical(start(gr.end(gr)), c(5L,9L,16L))
  expect_identical(start(gr.end(gr, width=100)),   c(1L,1L,1L))
  expect_identical(suppressWarnings(start(gr.end(gr, width=100, force=TRUE))),   c(-94L,-90L,-83L))
  expect_identical(start(gr.end(gr, width=100, ignore.strand=FALSE)),   c(1L,7L,13L))
  expect_identical(suppressWarnings(start(gr.end(gr, width=100, ignore.strand=FALSE, force=TRUE))),   c(-94L,7L,13L))

})

test_that("gr.mid", {

  expect_identical(start(gr.mid(gr)), c(4L,8L,14L))

})

#test_that("gr.dist", {

  #m <- gr.dist(gr, gr2, ignore.strand=TRUE)
  #expect_equal(m[1,2], 3)

#})

test_that("gr.rand", {

  set.seed(137)
  gg <- gr.rand(c(3,5), si)
  print(start(gg))
  expect_equal(start(gg)[1], 59221325L)

})

test_that("si2gr", {

  gg <- si2gr(si, strip.empty = TRUE)
  expect_equal(start(gg)[3], 1)
  expect_equal(end(gg)[1], 249250621)
  expect_equal(as.character(strand(gg)[1]), "+")

})

test_that("grbind", {

  expect_that(length(grbind(example_genes, example_dnase)) > 0, is_true())
  expect_equal(colnames(mcols(grbind(example_genes, example_dnase)))[1], "name")

})

test_that("gr.dice", {

  expect_equal(length(gr.dice(gr)[[3]]),4)

})

test_that("gr.findoverlaps", {

  fo <- gr.findoverlaps(example_genes, example_dnase)
  expect_equal(ncol(mcols(fo)), 2)
  expect_that(length(fo) > 0, is_true())

  ## null input
  expect_equal(length(gr.findoverlaps(example_genes, GRanges())), 0)
  expect_equal(length(gr.findoverlaps(GRanges(), GRanges())), 0)
  expect_equal(length(gr.findoverlaps(GRanges(), example_dnase)), 0)

})

test_that("gr.findoverlaps, return as data.table", {

  expect_error(gr.findoverlaps(example_genes, example_dnase, return.type = "data.frame"))

  fo <- gr.findoverlaps(example_genes, example_dnase, return.type = 'data.table')
  expect_identical(colnames(fo), c("start", "end", "query.id", "subject.id", "seqnames", "strand"))

})

test_that("rrbind", {

  expect_that(ncol(rrbind(mcols(example_genes), mcols(example_dnase))) > 2, is_true())
  expect_equal(ncol(rrbind(mcols(example_genes), mcols(example_dnase), union=FALSE)), 0)

})

test_that("gr.match", {

  ## accepts data.table
  expect_that(sum(!is.na(gr.match(gr2dt(example_genes), example_genes))) > 0, is_true())

  ## gives back overlapping matches
  gr11 <- GRanges(1, IRanges(c(10,20), width=5), strand=c("+", "-"))
  gr12 <- GRanges(1, IRanges(c(8,18, 100), width=5), strand=c("-", "+", "+"))

  expect_identical(gr.match(gr11, gr12), c(1L,2L))

  ## ignore strand is successfully passed
  expect_identical(gr.match(gr11, gr12, ignore.strand = FALSE), c(NA,NA))

})

test_that("gr.findoverlaps chunk", {

  fo  <- gr.findoverlaps(example_genes, example_dnase)
  fo2 <- gr.findoverlaps(example_genes, example_dnase, max.chunk = 1e7, verbose=TRUE)
  expect_identical(fo, fo2)

})

test_that("gr.findoverlaps, input data.table", {

  expect_equal(class(gr.findoverlaps(gr2dt(example_genes), example_genes, return.type='GRanges'))[1], "GRanges")
  expect_equal(class(gr.findoverlaps(gr2dt(example_genes), example_genes))[1], "data.table")
  expect_equal(class(gr.findoverlaps(gr2dt(example_genes), example_genes, max.chunk = 1e7))[1], "data.table")

})

test_that("gr.findoverlaps ignore.strand", {

  ## make a stranded DNAase track (for testing only)
  example_dnase2 = example_dnase
  set.seed(137)
  strand(example_dnase2) <- ifelse(runif(length(example_dnase)) > 0.5, '+', '-')

  ## get the overlaps with the original unstranded, and with ignore.strand
  fo1 <- gr.findoverlaps(example_dnase, example_genes)
  fo2 <- gr.findoverlaps(example_dnase2, example_genes, ignore.strand=TRUE)
  expect_identical(fo1, fo2)

  ## make sure no strands overlap
  fo1 <- gr.findoverlaps(example_dnase2, example_genes, ignore.strand=FALSE)
  expect_that(!any(strand(example_dnase2)[fo1$query.id] != strand(example_genes)[fo1$subject.id]), is_true())

})

# test_that("gr.findoverlaps foverlaps", {
#   ## make stranded DNAase track (for testing only)
#   example_dnase2 = example_dnase
#   set.seed(137)
#   strand(example_dnase2) <- ifelse(runif(length(example_dnase)) > 0.5, '+', '-')
#
#   ## assure that things are same with/without foverlaps
#   fo1 <- gr.findoverlaps(example_dnase, example_genes, foverlaps=TRUE)
#   fo1b <- gr.findoverlaps(example_genes, example_dnase, foverlaps=TRUE)
#   expect_identical(fo1$query.id, fo1b$subject.id)
#   expect_identical(fo1$subject.id, fo1b$query.id)
#   expect_identical(start(fo1), start(fo1b))
#
#   fo2 <- gr.findoverlaps(example_dnase, example_genes, foverlaps=FALSE)
#   fo2b <- gr.findoverlaps(example_genes, example_dnase, foverlaps=FALSE)
#   expect_identical(fo1$query.id, fo1b$subject.id)
#   expect_identical(fo1$subject.id, fo1b$query.id)
#   expect_identical(start(fo1), start(fo1b))
#
#   expect_identical(start(fo1), start(fo2))
#   expect_identical(end(fo1), end(fo2))
#   ## sort subject id because order for multi-overlaps is different, but it's OK
#   expect_identical(sort(fo1$subject.id), sort(fo2$subject.id))
# })

test_that("gr.findoverlap by", {

  e1 <- example_genes
  e2 <- example_dnase
  e1$bin <- e2$bin <- 1
  expect_error(gr.findoverlaps(example_genes, example_dnase, by = "dummy"))
  expect_that(length(gr.findoverlaps(e1, e2, by = "bin")) > 0, is_true())

})

#test_that("gr2dt works as expected", {

  #expect_identical(colnames(gr2dt(gr)), c("seqnames","start",'end','width','strand','name'))
  #expect_equal(nrow(gr2dt(gr)), length(gr))

  #subjectdt <- gr2dt(example_genes)
  #expect_that(!any(subjectdt$start!=start(example_genes)), is_true)
  #expect_that(!any(subjectdt$end!=end(example_genes)), is_true)
  #expect_that(!any(subjectdt$width!=width(example_genes)), is_true)
  #expect_that(!any(subjectdt$strand!=strand(example_genes)), is_true)
  #expect_that(!any(subjectdt$seqnames!=seqnames(example_genes)), is_true)
#})

test_that("dt2gr", {

  dt <- data.table(seqnames=1, start=1, end=10, strand='+', name="A")
  expect_equal(as.character(strand(dt2gr(dt))), '+')
  expect_equal(start(dt2gr(dt)), 1)
  expect_equal(dt2gr(dt)$name, "A")
  
  expect_equal(start(dt2gr(as.data.frame(dt)))[1], 1)
  expect_error(dt2gr(1))
  
  dt <- data.table(seqnames1=1, start=1, end=10, strand='+', name="A")
  expect_error(dt2gr(dt))
})

test_that("gr.val", {
  
  gr <- GRanges(1, IRanges(1e6,2e6))
  
  expect_equal(colnames(mcols(gr.val(gr, example_genes, val = 'name'))),"name")
    
})

test_that("gr.duplicated", {
  gr <- GRanges(c(1,1,1), IRanges(c(2,5,5), width=1), val=c(1,2,3))
  
  expect_identical(gr.duplicated(gr), c(FALSE, FALSE, TRUE))
  expect_identical(gr.duplicated(gr, by="val"), c(FALSE, FALSE, FALSE))
})

test_that("gr.sample", {

  set.seed(137)
  gg <- gr.sample(reduce(example_genes), 10, len=1)
  expect_equal(unique(width(gg)), 1)

  ## query width less than output
  expect_error(gr.sample(gr.start(example_genes), c(1:3), len=5)) 
  
  gg <- gr.sample(example_genes[1:5], c(2,2,3,4,5), len=2)
  expect_equal(length(gg), 16)
  expect_equal(sum(width(gg)), 32)

})


test_that("gr.sample without replace", {

  gg <- gr.sample(reduce(example_genes), 10, len=1, replace=FALSE)
  expect_equal(unique(width(gg)), 1)
  
  expect_error(gr.sample(example_genes[1:3], 1e7, len=1, replace=FALSE))

  gg <- gr.sample(example_genes[1:5], c(2,2,3,4,5), len=2, replace=FALSE)
  expect_equal(length(gg), 16)
  expect_equal(sum(width(gg)), 32)

})

test_that("gr.chr", {

  expect_equal(as.character(seqnames(gr.chr(GRanges(c(1,"chrX"), IRanges(c(1,2), 1))))), c("chr1","chrX"))

})

test_that("gr.nochr",{

  expect_identical(gr.nochr(gr.chr(example_genes)), example_genes)

})

test_that("gr.string", {

  expect_that(grepl(":",gr.string(example_genes)[1]), is_true())
  expect_that(grepl("-",gr.string(example_genes)[1]), is_true())
  expect_that(grepl("(+|-)",gr.string(example_genes)[1]), is_true())

})

test_that("gr.fix", {

  gg <- GRanges(c("X",1), IRanges(c(1,2), width=1))
  expect_equal(length(seqlengths(gr.fix(gg, si))), 25)
  expect_equal(length(seqlengths(gr.fix(gg, BSgenome.Hsapiens.UCSC.hg19::Hsapiens))), 95)

})

test_that("gr.flipstrand", {

  expect_identical(as.character(strand(gr.flipstrand(gr))), c("-","+","+"))

})

#test_that("gr.flatten", {

  #df <- gr.flatten(gr)
  #expect_equal(as.character(class(df)), "data.frame")
  #expect_identical(colnames(df), c("start", 'end','name'))
  #expect_identical(df$start, c(1,4,7))

  #df <- gr.flatten(gr, gap=5)
  #expect_equal(df$start[2] - df$end[1], 5 + 1)

#})

test_that("grlbind", {

  grl.hiC2 <- grl.hiC[1:20]
  mcols(grl.hiC2)$test = 1
  suppressWarnings(gg <- grlbind(grl.hiC2, grl.hiC[1:30]))
  expect_equal(length(gg), 50)
  expect_equal(colnames(mcols(gg)), "test")

  names(grl.hiC) <- NULL
  out <- grlbind(grl.hiC)
  names(out) <- NULL
  expect_identical(out, grl.hiC)

  ## expect error
  expect_error(grlbind('d'))

})

test_that("streduce", {

  gg <- streduce(grl.hiC, pad=10)
  expect_equal(length(gg), length(reduce(gg)))

  gg <- streduce(example_genes, pad=10)
  expect_equal(length(gg), length(reduce(gg)))

})

test_that("parse.gr", {

  parse.gr(c('1:1e6-5e6+', '2:2e6-5e6-'))

})

test_that("parse.grl", {

  parse.grl(c('1:1e6-5e6+;5:10-2000', '2:2e6-5e6-;10:100231321-100231399')) 
              
})


test_that("grl.pivot", {

  gg <- grl.pivot(grl.hiC)
  expect_equal(as.character(class(gg)), "GRangesList")
  expect_equal(length(gg),2)
  expect_equal(length(gg[[1]]), 10000)

})

test_that("gr.tile", {

  expect_identical(start(gr.tile(gr, w=3)), c(3L, 7L, 13L, 16L))
  expect_equal(length(gr.tile(GRanges())), 0)

})

test_that("grl.string", {

  expect_that(nchar(names(grl.string(grl.hiC[1:5])[1])) > 0, is_true())
  expect_that(grepl(",",grl.string(grl.hiC[1:5])[1]), is_true())

})

test_that("grl.unlist", {

  gg <- grl.unlist(grl.hiC)
  expect_equal(length(gg), length(grl.hiC)*2)
  expect_equal(max(mcols(gg)$grl.iix), 2)
  expect_equal(max(mcols(gg)$grl.ix), length(grl.hiC))

})

test_that("grl.in", {

  gg <- grl.in(grl.hiC[1:100], example_genes)
  expect_equal(length(gg), 100)

})

test_that("gr.fix with null genome", {

  gg <- GRanges(c("X",1), IRanges(c(1,2), width=1))
  es <- structure(c(1L,2L), names=c("X","1"))
  expect_identical(seqlengths(seqinfo(gr.fix(gg))), es)

  es <- structure(c("gen","gen"), names=c("X","1"))
  expect_identical(genome(seqinfo(gr.fix(gg, gname='gen'))), es)

})

#test_that("grfo", {

  #fo <- example_genes %*% example_dnase
  #expect_that(ncol(mcols(fo)) > 2, is_true())
  #expect_that(length(fo) > 0, is_true())

#})

#test_that("gr.simplify", {

  #gg <- gr.simplify(gr, pad=4)
  #expect_identical(end(gg), c(5L, 16L))
  #expect_equal(length(gg), 2)
  
  #gr$field <- c("A","B","B")
  #expect_equal(length(gr.simplify(gr, pad=4, field="name")), 3)
  #expect_equal(length(gr.simplify(gr, pad=4, field="field")), 2)

  #expect_equal(class(gr.simplify(gr, pad=4, field="name", split = TRUE))[1], "GRangesList")
  #expect_equal(ncol(mcols((gr.simplify(gr, pad=4, field="name", include.val = FALSE)))), 0)

#})

test_that("gr.tile.map", {

  gr1 <- gr.tile(GRanges(1, IRanges(1,100)), w=10)
  gr2 <- gr.tile(GRanges(1, IRanges(1,100)), w=5)
  gg <- gr.tile.map(gr1, gr2, verbose=TRUE)
  expect_equal(length(gg), 10)
  expect_equal(length(unlist(gg)), 20)

})

#test_that('ra.overlaps throws error', {

 #ss<- split(example_genes, example_genes$name)
 #expect_error(ra.overlaps(ss,ss))

#})

#test_that("ra.overlaps handles empty",{

  ## test empty inputs and no overlaps inputs
  #expect_equal(ra.overlaps(GRangesList(), grl1)[1], NA)
  #expect_equal(ra.overlaps(grl2[2:3], grl1)[1], NA)

#})

#test_that("ra.overlaps handles wrong signs", {

  ## make one that overlaps, but wrong signs
  #grl3 <- grl1[115]
  #strand(grl3[[1]]) <- c("+", "-")
  #expect_equal(ra.overlaps(grl3, grl2)[1], NA)

#})

#test_that('ra.overlaps', {

  #grl1 <- grlbind(grl1, grl2)
  #ro <- ra.overlaps(grl1, grl2)
  #expect_equal(class(ro), "matrix")
  #expect_equal(nrow(ro), 252)
  #expect_equal(ncol(ro), 2)
  #expect_equal(nrow(ra.overlaps(grl2, grl2)), length(grl2))

#})

test_that('%_% works', {
  
  gr1 <- GRanges(1, IRanges(10,20), strand="+")
  gr2 <- GRanges(1, IRanges(15,25), strand="-")
  gr3 <- "1:1-15"
  expect_equal(width(gr1 %_% gr2), 5)
  expect_equal(width(gr1 %_% gr3), 5)
})

test_that('%Q% works', {

  testset <- GRanges(seqnames = Rle(1,c(5)) , ranges = IRanges(1:5 , end = 2000:2004) , strand = Rle(strand(c("+")) , c(5)) , mean = c(1, -2 , -5 , 5 ,6))
 # expect_equal(length(testset %Q% (mean > 0)) , 3)    

})

test_that('%^% works', {

  testset <- GRanges(seqnames = Rle(1,c(5)) , ranges = IRanges(1:5 , end = 2000:2004) , strand = Rle(strand(c("+")) , c(5)) , mean = c(1, -2 , -5 , 5 ,6)) 
  expect_equal(length(testset %^% testset) , 5)

})
