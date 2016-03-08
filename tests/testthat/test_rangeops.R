library(gUtils)
library(BSgenome.Hsapiens.UCSC.hg19)
context("Range ops")

gr  <- GRanges(1, IRanges(c(3,7,13), c(5,9,16)), strand=c('+','-','-'), seqinfo=Seqinfo("1", 25), name=c("A","B","C"))
gr2 <- GRanges(1, IRanges(c(1,9), c(6,14)), strand=c('+','-'), seqinfo=Seqinfo("1", 25), field=c(1,2))
dt <- data.table(seqnames=1, start=c(2,5,10), end=c(3,8,15))

test_that("gr.start ", {

  expect_identical(start(gr.start(gr)), c(3L,7L,13L))
  expect_identical(end(gr.start(gr)),   c(3L,7L,13L))
  expect_identical(end(gr.start(gr, width=100)),   c(25L,25L,25L))
  expect_identical(suppressWarnings(end(gr.start(gr, width=100, force=TRUE))),   c(102L,106L,112L))
  expect_identical(end(gr.start(gr, width=100, ignore.strand=FALSE)),   c(25L,9L,16L))
  expect_identical(suppressWarnings(start(gr.start(gr, width=100, ignore.strand=FALSE, force=TRUE))),   c(3L,-90L,-83L))
  expect_identical(suppressWarnings(end(gr.start(gr, width=100, force=TRUE))),   c(102L,106L,112L))
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

test_that("gr.dist", {
  m <- gr.dist(gr, gr2, ignore.strand=TRUE)
  expect_equal(m[1,2], 3)
})

test_that("gr.rand", {
  set.seed(137)
  gg <- gr.rand(c(3,5), si)
  print(start(gg))
  expect_equal(start(gg)[1], 59221325L)
})

test_that("si2gr", {
  gg <- si2gr(si)
  expect_equal(start(gg)[3], 1)
  expect_equal(end(gg)[1], 249250621)
  expect_equal(as.character(strand(gg)[1]), "+")
})

test_that("grbind", {
  expect_equal(length(grbind(gr.genes, gr.DNAase)), 310779)
  expect_equal(colnames(mcols(grbind(gr.genes, gr.DNAase)))[10], "cdsStartStat")
})

test_that("gr.dice", {
  expect_equal(length(gr.dice(gr)[[3]]),4)
})

test_that("gr.findoverlaps", {
  fo <- gr.findoverlaps(gr.genes, gr.DNAase)
  expect_equal(start(fo)[3], 67042300)
  expect_equal(length(fo), 256058)
})

test_that("gr.in", {
  expect_equal(sum(gr.in(gr.genes, gr.DNAase)),35570)
})

test_that("gr2dt works as expected", {
  expect_identical(colnames(gr2dt(gr)), c("seqnames","start",'end','width','strand','name'))
  expect_equal(nrow(gr2dt(gr)), length(gr))
})

test_that("dt2gr", {
  dt <- data.table(seqnames=1, start=1, end=10, strand='+', name="A")
  expect_equal(as.character(strand(dt2gr(dt))), '+')
  expect_equal(start(dt2gr(dt)), 1)
  expect_equal(dt2gr(dt)$name, "A")

  dt <- data.table(seqns=1, start=1, end=10, strand='+', name="A")
  expect_error(dt2gr(dt))

})

test_that("gr.sample", {
  set.seed(137)
  gg <- gr.sample(reduce(gr.genes), 10, len=1)
  expect_equal(start(gg)[1], 77055451)
  expect_equal(unique(width(gg)), 1)

  expect_error(gr.sample(reduce(gr.genes), c(1:3), len=5))
  set.seed(137)
  gg <- gr.sample(gr.genes[1:5], c(2,2,3,4,5), len=2)
  expect_equal(length(gg), 16)
  expect_equal(sum(width(gg)), 32)
})

test_that("gr.chr", {
  expect_equal(as.character(seqnames(gr.chr(GRanges(c(1,"chrX"), IRanges(c(1,2), 1))))), c("chr1","chrX"))
})

test_that("gr.nochr",{
  expect_identical(gr.nochr(gr.chr(gr.genes)), gr.genes)
})

test_that("gr.string", {
  expect_equal(gr.string(gr.genes)[1], "1:67000041-67208778+")
})

test_that("gr.fix", {
  library(BSgenome.Hsapiens.UCSC.hg19)
  gg <- GRanges(c("X",1), IRanges(c(1,2), width=1))
  expect_equal(length(seqlengths(gr.fix(gg, si))), 25)
  expect_equal(length(seqlengths(gr.fix(gg, Hsapiens))), 95)
})

test_that("gr.flipstrand", {
  expect_identical(as.character(strand(gr.flipstrand(gr))), c("-","+","+"))
})

test_that("gr.flatten", {
  df <- gr.flatten(gr)
  expect_equal(as.character(class(df)), "data.frame")
  expect_identical(colnames(df), c("start", 'end','name'))
  expect_identical(df$start, c(1L,4L,7L))

  df <- gr.flatten(gr, gap=5)
  expect_equal(df$start[2] - df$end[1], 5 + 1)
})

test_that("grlbind", {
  grl.hiC2 <- grl.hiC[1:20]
  mcols(grl.hiC2)$test = 1
  suppressWarnings(gg <- grlbind(grl.hiC2, grl.hiC[1:30]))
  expect_equal(length(gg), 50)
  expect_equal(colnames(mcols(gg)), "test")
})

test_that("streduce", {
  gg <- streduce(grl.hiC, pad=10)
  expect_equal(length(gg), length(reduce(gg)))

  gg <- streduce(gr.genes, pad=10)
  expect_equal(length(gg), length(reduce(gg)))
})

test_that("grl.pivot", {
  gg <- grl.pivot(grl.hiC)
  expect_equal(as.character(class(gg)), "GRangesList")
  expect_equal(length(gg),2)
  expect_equal(length(gg[[1]]), 537341)
})

test_that("gr.tile", {
  expect_identical(start(gr.tile(gr, w=3)), c(3L, 7L, 13L, 16L))
  expect_equal(length(gr.tile(GRanges())), 0)
})

#test_that("grl.string", {
#  expect_equal(grl.string(grl.hiC[1:5])[1]=="14:29864023-29864023+,14:19001056-19001056+", TRUE)
#})

test_that("grl.unlist", {
  gg <- grl.unlist(grl.hiC)
  expect_equal(length(gg), length(grl.hiC)*2)
  expect_equal(max(mcols(gg)$grl.iix), 2)
  expect_equal(max(mcols(gg)$grl.ix), length(grl.hiC))

})

test_that("grl.in", {
  gg <- grl.in(grl.hiC[1:100], gr.genes)
  expect_equal(length(gg), 100)
})
