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
  gg <- si2gr(si, strip.empty = TRUE)
  expect_equal(start(gg)[3], 1)
  expect_equal(end(gg)[1], 249250621)
  expect_equal(as.character(strand(gg)[1]), "+")
})

test_that("grbind", {
  expect_equal(length(grbind(gr.genes, gr.DNAase)), 20000)
  expect_equal(colnames(mcols(grbind(gr.genes, gr.DNAase)))[10], "cdsStartStat")
})

test_that("gr.dice", {
  expect_equal(length(gr.dice(gr)[[3]]),4)
})

test_that("gr.findoverlaps", {
  fo <- gr.findoverlaps(gr.genes, gr.DNAase)
  expect_equal(start(fo)[3], 5946813)
  expect_equal(length(fo), 1856)

  ## null input
  expect_equal(length(gr.findoverlaps(gr.genes, GRanges())), 0)
  expect_equal(length(gr.findoverlaps(GRanges(), GRanges())), 0)
  expect_equal(length(gr.findoverlaps(GRanges(), gr.DNAase)), 0)
})

test_that("gr.findoverlaps, return as data.table", {
  expect_error(gr.findoverlaps(gr.genes, gr.DNAase, return.type = "data.frame"))

  fo <- gr.findoverlaps(gr.genes, gr.DNAase, return.type = 'data.table')
  expect_identical(colnames(fo), c("start", "end", "query.id", "subject.id", "seqnames", "strand"))
})

test_that("gr.findoverlaps ignore.strand", {

  ## make a stranded DNAase track (for testing only)
  gr.DNAase2 = gr.DNAase
  set.seed(137)
  strand(gr.DNAase2) <- ifelse(runif(length(gr.DNAase)) > 0.5, '+', '-')

  ## get the overlaps with the original unstranded, and with ignore.strand
  fo1 <- gr.findoverlaps(gr.DNAase, gr.genes)
  fo2 <- gr.findoverlaps(gr.DNAase2, gr.genes, ignore.strand=TRUE)
  expect_identical(fo1, fo2)

  ## make sure no strands overlap
  fo1 <- gr.findoverlaps(gr.DNAase2, gr.genes, ignore.strand=FALSE)
  expect_that(!any(strand(gr.DNAase2)[fo1$query.id] != strand(gr.genes)[fo1$subject.id]), is_true())
})

# test_that("gr.findoverlaps foverlaps", {
#   ## make stranded DNAase track (for testing only)
#   gr.DNAase2 = gr.DNAase
#   set.seed(137)
#   strand(gr.DNAase2) <- ifelse(runif(length(gr.DNAase)) > 0.5, '+', '-')
#
#   ## assure that things are same with/without foverlaps
#   fo1 <- gr.findoverlaps(gr.DNAase, gr.genes, foverlaps=TRUE)
#   fo1b <- gr.findoverlaps(gr.genes, gr.DNAase, foverlaps=TRUE)
#   expect_identical(fo1$query.id, fo1b$subject.id)
#   expect_identical(fo1$subject.id, fo1b$query.id)
#   expect_identical(start(fo1), start(fo1b))
#
#   fo2 <- gr.findoverlaps(gr.DNAase, gr.genes, foverlaps=FALSE)
#   fo2b <- gr.findoverlaps(gr.genes, gr.DNAase, foverlaps=FALSE)
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
  expect_error(gr.findoverlaps(gr.genes, gr.DNAase, by = "dummy"))
  gr.findoverlaps(gr.genes, gr.DNAase, by = "bin")
})

test_that("gr2dt works as expected", {
  expect_identical(colnames(gr2dt(gr)), c("seqnames","start",'end','width','strand','name'))
  expect_equal(nrow(gr2dt(gr)), length(gr))

  #subjectdt <- gr2dt(gr.genes)
  #expect_that(!any(subjectdt$start!=start(gr.genes)), is_true)
  #expect_that(!any(subjectdt$end!=end(gr.genes)), is_true)
  #expect_that(!any(subjectdt$width!=width(gr.genes)), is_true)
  #expect_that(!any(subjectdt$strand!=strand(gr.genes)), is_true)
  #expect_that(!any(subjectdt$seqnames!=seqnames(gr.genes)), is_true)
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
  expect_equal(start(gg)[1], 131788860)
  expect_equal(unique(width(gg)), 1)

  expect_error(gr.sample(reduce(gr.genes), c(1:3), len=5))
  set.seed(137)
  gg <- gr.sample(gr.genes[1:5], c(2,2,3,4,5), len=2)
  expect_equal(length(gg), 16)
  expect_equal(sum(width(gg)), 32)

})


test_that("gr.sample without replace", {
  set.seed(137)
  gg <- gr.sample(reduce(gr.genes), 10, len=1, replace=FALSE)
  expect_equal(start(gg)[1], 131788861)
  expect_equal(unique(width(gg)), 1)

  expect_error(gr.sample(reduce(gr.genes)[1:3], 10000, len=1, replace=FALSE))

  set.seed(137)
  gg <- gr.sample(gr.genes[1:5], c(2,2,3,4,5), len=2, replace=FALSE)
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
  expect_equal(gr.string(gr.genes)[1], "12:10772742-10787217-")
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
  expect_identical(df$start, c(1,4,7))

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
  expect_equal(length(gg[[1]]), 10000)
})

test_that("gr.tile", {
  expect_identical(start(gr.tile(gr, w=3)), c(3L, 7L, 13L, 16L))
  expect_equal(length(gr.tile(GRanges())), 0)
})

test_that("grl.string", {
  exp_result = structure(names="450448", "14:66569495-66569495-,14:66716403-66716403+")
  expect_identical(grl.string(grl.hiC[1:5])[1], exp_result)
})

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

test_that("gr.fix with null genome", {
  gg <- GRanges(c("X",1), IRanges(c(1,2), width=1))
  es <- structure(c(1L,2L), names=c("X","1"))
  expect_identical(seqlengths(seqinfo(gr.fix(gg))), es)

  es <- structure(c("gen","gen"), names=c("X","1"))
  expect_identical(genome(seqinfo(gr.fix(gg, gname='gen'))), es)

})

test_that("grfo", {
  fo <- gr.genes %*% gr.DNAase
  expect_equal(ncol(mcols(fo)), 17)
  expect_equal(length(fo), 1856)
})

test_that("gr.simplify", {
  gg <- gr.simplify(gr, pad=4)
  expect_identical(end(gg), c(5L, 16L))
  expect_equal(length(gg), 2)

  gr$field <- c("A","B","B")
  expect_equal(length(gr.simplify(gr, pad=4, field="name")), 3)
  expect_equal(length(gr.simplify(gr, pad=4, field="field")), 2)

  expect_equal(class(gr.simplify(gr, pad=4, field="name", split = TRUE))[1], "GRangesList")
  expect_equal(ncol(mcols((gr.simplify(gr, pad=4, field="name", include.val = FALSE)))), 0)
})

test_that("gr.tile.map", {

  gr1 <- gr.tile(GRanges(1, IRanges(1,100)), w=10)
  gr2 <- gr.tile(GRanges(1, IRanges(1,100)), w=5)
  gg <- gr.tile.map(gr1, gr2, verbose=TRUE)
  expect_equal(length(gg), 10)
  expect_equal(length(unlist(gg)), 20)

})

test_that('ra.overlaps', {
  ro <- ra.overlaps(grl1, grl2)
  expect_equal(class(ro), "matrix")
  expect_equal(nrow(ro), 1)
  expect_equal(ncol(ro), 2)
})
