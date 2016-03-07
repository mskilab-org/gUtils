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
  gg <- gr.rand(c(3,5), si, 137)
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
