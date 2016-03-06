library(gUtils)
library(BSgenome.Hsapiens.UCSC.hg19)
context("Range ops")

gr  <- GRanges(1, IRanges(c(3,7,13), c(5,9,16)), strand=c('+','-','-'), seqinfo=Seqinfo("1", 25), name=c("A","B","C"))
gr2 <- GRanges(1, IRanges(c(1,9), c(6,14)), strand=c('+','-'), seqinfo=Seqinfo("1", 25), field=c(1,2))
dt <- data.table(seqnames=1, start=c(2,5,10), end=c(3,8,15))

test_that("gr.start works as expected", {

  expect_identical(start(gr.start(gr)), c(3L,7L,13L))
  expect_identical(end(gr.start(gr)),   c(3L,7L,13L))
  
})

test_that("gr.dist", {
  
  m <- gr.dist(gr, gr2, ignore.strand=TRUE)
  expect_equal(m[1,2], 3)
  
})

#test_that("gr2dt works as expected", {
#  gr <- GRanges(1, IRanges(10,20), seqinfo=Seqinfo("1", 200))
#  expect_identical(colnames(gr2dt(gr)), c("seqnames","start",'end','width','strand'))
#  expect_equal(nrow(gr2dt(gr)), length(gr))
#})
