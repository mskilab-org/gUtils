library(gUtils)
context("Range ops")

test_that("gr.start works as expected", {

  gr <- GRanges(1, IRanges(10,20), seqinfo=Seqinfo("1", 200))
  expect_equal(start(gr.start(gr)), 10)
  expect_equal(end(gr.start(gr)), 10)

})

#test_that("gr2dt works as expected", {
#  gr <- GRanges(1, IRanges(10,20), seqinfo=Seqinfo("1", 200))
#  expect_identical(colnames(gr2dt(gr)), c("seqnames","start",'end','width','strand'))
#  expect_equal(nrow(gr2dt(gr)), length(gr))
#})
