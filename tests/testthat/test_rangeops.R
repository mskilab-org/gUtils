library(gUtils)
library(BSgenome.Hsapiens.UCSC.hg19)
Sys.setenv(DEFAULT_BSGENOME = "BSgenome.Hsapiens.UCSC.hg19::Hsapiens")

context("Range ops")




gr  <- GRanges(1, IRanges(c(3,7,13), c(5,9,16)), strand=c('+','-','-'), seqinfo=Seqinfo("1", 25), name=c("A","B","C"))
gr2 <- GRanges(1, IRanges(c(1,9), c(6,14)), strand=c('+','-'), seqinfo=Seqinfo("1", 25), field=c(1,2))
dt <- data.table(seqnames=1, start=c(2,5,10), end=c(3,8,15))





test_that("hg_seqlengths()", {
    
    gr  <- GRanges(1, IRanges(c(3,7,13), c(5,9,16)), strand=c('+','-','-'), seqinfo=Seqinfo("1", 25), name=c("A","B","C"))
    expect_identical(as.numeric(length(hg_seqlengths())), 25)
    ee = structure(names="1", 249250621L)
    expect_identical(hg_seqlengths(Hsapiens)[1],ee)
    expect_equal(names(hg_seqlengths(Hsapiens, chr=TRUE)[1]), "chr1")
    expect_equal(length(hg_seqlengths(Hsapiens, include.junk = TRUE)), 93)

})


test_that("gr2dt", {

    gr = GRanges(1, IRanges(c(3,7,13), c(5,9,16)), strand=c('+','-','-'), seqinfo=Seqinfo("1", 25), name=c("A","B","C"))
    example_genes = GRanges(2, IRanges(c(233101, 233101, 231023, 231023, 229966), c(233229, 233229, 231191, 231191, 230044)), strand = c("-"), type = c("exon", "CDS", "exon", "CDS", "exon"))
    expect_identical(colnames(gr2dt(gr)), c("seqnames", "start", "end", "strand", "width", "name"))
    expect_equal(nrow(gr2dt(gr)), length(gr))

    subjectdt <- gr2dt(example_genes)
    expect_equal(!any(subjectdt$start!=start(example_genes)), TRUE)
    expect_equal(!any(subjectdt$end!=end(example_genes)), TRUE)
    expect_equal(!any(subjectdt$width!=width(example_genes)), TRUE)
    expect_equal(all(as.character(subjectdt$strand) == as.character(strand(example_genes))), TRUE)
    expect_equal(all(subjectdt$seqnames == as.vector(seqnames(example_genes))), TRUE)
    
})


test_that("gr.start ", {
    
    gr = GRanges(1, IRanges(c(3,7,13), c(5,9,16)), strand=c('+','-','-'), seqinfo=Seqinfo("1", 25), name=c("A","B","C"))
    expect_identical(start(gr.start(gr)), c(3L,7L,13L))
    expect_identical(end(gr.start(gr)),   c(3L,7L,13L))
    expect_identical(end(gr.start(gr, width=100)),   c(5L,9L,16L))
    expect_identical(end(gr.start(gr, width=100, clip=FALSE)),   c(25L,25L,25L))
    expect_identical(suppressWarnings(end(gr.start(gr, width=100, force=TRUE))),   c(5L,9L,16L))
    expect_identical(suppressWarnings(end(gr.start(gr, width=100, force=TRUE, clip=FALSE))),   c(102L,106L,112L))
    expect_identical(end(gr.start(gr, width=100, ignore.strand=FALSE)),   c(5L,9L,16L))
    expect_identical(end(gr.start(gr, width=100, ignore.strand=FALSE, clip=FALSE)),   c(25L,9L,16L))
    expect_identical(suppressWarnings(start(gr.start(gr, width=100, ignore.strand=FALSE, force=TRUE))),   c(3L,7L,13L))
    expect_identical(suppressWarnings(start(gr.start(gr, width=100, ignore.strand=FALSE, force=TRUE, clip=FALSE))),   c(3L,-90L,-83L))
    expect_identical(suppressWarnings(end(gr.start(gr, width=100, force=TRUE))),   c(5L, 9L,16L))
    expect_identical(suppressWarnings(end(gr.start(gr, width=100, force=TRUE, clip=FALSE))),   c(102L,106L,112L))
    expect_identical(end(gr.start(gr, width=10000)), c(5L, 9L, 16L))  ## default clip=TRUE
    
})


test_that("dt2gr", {

    Sys.setenv(DEFAULT_BSGENOME = "BSgenome.Hsapiens.UCSC.hg19::Hsapiens")
    dt <- data.table(seqnames=1, start=1, end=10, strand='+', name="A")
    expect_equal(as.character(strand(dt2gr(dt))), '+')
    expect_equal(start(dt2gr(dt)), 1)
    expect_equal(dt2gr(dt)$name, "A")

    expect_equal(start(dt2gr(as.data.frame(dt)))[1], 1)
    expect_error(suppressWarnings(dt2gr(1)))

    dt <- data.table(sdf=1, start=1, end=10, strand='+', name="A")
    expect_error(suppressWarnings(dt2gr(dt)))    ### warning within error---warning: coercing to GRanges via non-standard columns
})


test_that("gr.end", {

    gr = GRanges(1, IRanges(c(3,7,13), c(5,9,16)), strand=c('+','-','-'), seqinfo=Seqinfo("1", 25), name=c("A","B","C"))
    expect_identical(start(gr.end(gr)), c(5L,9L,16L))
    expect_identical(start(gr.end(gr, width=100)),   c(3L,7L,13L))
    expect_identical(start(gr.end(gr, width=100, clip=FALSE)),   c(1L,1L,1L))
    expect_identical(suppressWarnings(start(gr.end(gr, width=100, force=TRUE))),   c(3L,7L,13L))
    expect_identical(suppressWarnings(start(gr.end(gr, width=100, force=TRUE, clip=FALSE))),   c(-94L,-90L,-83L))
    expect_identical(start(gr.end(gr, width=100, ignore.strand=FALSE)),   c(3L,7L,13L))
    expect_identical(start(gr.end(gr, width=100, ignore.strand=FALSE, clip=FALSE)),   c(1L,7L,13L))
    expect_identical(suppressWarnings(start(gr.end(gr, width=100, ignore.strand=FALSE, force=TRUE, clip=FALSE))),   c(-94L,7L,13L))

})


test_that("gr.mid", {
    
    gr = GRanges(1, IRanges(c(3,7,13), c(5,9,16)), strand=c('+','-','-'), seqinfo=Seqinfo("1", 25), name=c("A","B","C"))
    expect_identical(start(gr.mid(gr)), c(4L,8L,14L))

})


test_that('gr.rand', {

    set.seed(137)
    gg <- gr.rand(c(3,5), si)
    print(start(gg))
    expect_equal(start(gg)[1], 59221325L)

})


test_that('gr.trim', {

    ## from example
    ## trim the first 20 and last 50 bases
    ## gr.trim(GRanges(1, IRanges(1e6, width=1000)), starts=20, ends=950)
    ## return value: GRanges on 1:1,000,019-1,000,949
    expect_equal(width(gr.trim(GRanges(1, IRanges(1e6, width=1000)), starts=20, ends=950)), 931)

})


test_that("gr.sample", {

    ## ALERT: change of argument name, "k" instead of "len"
    set.seed(42)
    gg <- gr.sample(reduce(example_genes), 10, k=1)
    expect_equal(unique(width(gg)), 10)

    ## query width less than output
    ## expect_error(gr.sample(gr.start(example_genes), c(1:3), k=5))

    gg <- suppressWarnings(gr.sample(example_genes[1:5], c(2,2,3,4,5), k=2))    ### expect warning: longer object length is not a multiple of shorter object length
    expect_equal(length(gg), 5)
    expect_equal(sum(width(gg)), 16)

})


test_that("gr.sample without replace", {
    
    gene1 = GRanges("3:3540000-234329000")
    gene2 = GRanges("2:24440-30000")
    gene3 = GRanges("2:278444-321000")
    gene4 = GRanges("3:1000-1500")
    example_genes = suppressWarnings(c(gene1, gene2, gene3, gene4))  ##  The 2 combined objects have no sequence levels in common. (Use suppressWarnings() to suppress this warning
    gg <- gr.sample(reduce(example_genes), 10, k=1, replace=FALSE)
    expect_equal(unique(width(gg)), 10)

    expect_equal(width(gr.sample(example_genes[1:3], 1e7, k=1, replace=FALSE)),  10000000)

    gg <- gr.sample(example_genes, c(2,2,3,4), k=2, replace=FALSE)
    expect_equal(length(gg), 4)
    expect_equal(sum(width(gg)), 11)

})


test_that("si2gr", {

    gg <- si2gr(si, strip.empty = TRUE)
    expect_equal(start(gg)[3], 1)
    expect_equal(end(gg)[1], 249250621)
    expect_equal(as.character(strand(gg)[1]), "+")

})


test_that("gr.bind", {

    example_dnase = GRanges(1, IRanges(c(562757, 564442, 564442), c(563203, 564813, 564813)), strand = c("-", "+", "+"))
    example_genes = GRanges(2, IRanges(c(233101, 233101, 231023, 231023, 229966), c(233229, 233229, 231191, 231191, 230044)), strand = c("-"), type = c("exon", "CDS", "exon", "CDS", "exon"))
    expect_equal(length(suppressWarnings(gr.bind(example_genes, example_dnase))) > 0, TRUE)

})


test_that("grl.bind", {

    grl.hiC2 <- grl.hiC[1:20]
    mcols(grl.hiC2)$test = 1
    suppressWarnings(gg <- grl.bind(grl.hiC2, grl.hiC[1:30]))
    expect_equal(length(gg), 50)
    expect_equal(colnames(mcols(gg)), "test")

    ## names(grl.hiC) <- NULL
    ## out <- grl.bind(grl.hiC)
    ## names(out) <- NULL
    ## expect_identical(out, grl.hiC)

    ## expect error
    expect_error(grl.bind('d'))

})


test_that("gr.chr", {

    expect_equal(as.character(seqnames(gr.chr(GRanges(c(1,"chrX"), IRanges(c(1,2), 1))))), c("chr1", "chrX"))

})


test_that("streduce", {

    gg = streduce(grl.hiC, pad=10)
    expect_equal(length(gg), length(reduce(gg)))

    gg2 = streduce(example_genes, pad=10)
    expect_equal(length(gg2), length(reduce(gg2)))

})


test_that("gr.string", {

    expect_that(grepl(":", gr.string(example_genes)[1]), is_true())
    expect_that(grepl("-", gr.string(example_genes)[1]), is_true())
    expect_that(grepl("(+|-)", gr.string(example_genes)[1]), is_true())

})


test_that('grl.reduce', {

    gr = GRanges(1, IRanges(c(3,7), c(5,9)), strand=c('+','-'))
    gr1 = GRanges(1, IRanges(c(10,20), width=5), strand=c("+", "-"))
    grl1 = GRangesList("gr"=gr, "gr1"=gr1)
    expect_equal(start(suppressWarnings(ranges(grl.reduce(grl1, 50)[[1]]))), -47)
    expect_equal(end(suppressWarnings(ranges(grl.reduce(grl1, 50)[[1]]))), 59)
    expect_equal(width(suppressWarnings(ranges(grl.reduce(grl1, 50)[[1]]))), 107)
    expect_equal(suppressWarnings(width(grl.reduce(grl1, 2000, clip=TRUE)[[1]])), 10)
    expect_equal(suppressWarnings(width(grl.reduce(grl1, 2000, clip=TRUE)[[2]])), 25)

})


test_that("grl.string", {

    expect_that(nchar(names(grl.string(grl.hiC[1:5])[1])) > 0, is_true())
    expect_that(grepl(",", grl.string(grl.hiC[1:5])[1]), is_true())

})



test_that("gr.fix", {

    gg <- GRanges(c("X",1), IRanges(c(1,2), width=1))
    expect_equal(length(seqlengths(gr.fix(gg, si))), 25)
    expect_equal(length(seqlengths(gr.fix(gg, BSgenome.Hsapiens.UCSC.hg19::Hsapiens))), 95)

})


test_that("gr.fix with null genome", {

    gg <- GRanges(c("X",1), IRanges(c(1,2), width=1))
    es <- structure(c(1L,2L), names=c("X","1"))
    expect_identical(seqlengths(seqinfo(gr.fix(gg))), es)
    
    es <- structure(c("gen","gen"), names=c("X","1"))
    expect_identical(genome(seqinfo(gr.fix(gg, gname='gen'))), es)

})



test_that("gr.flatten", {

    gr = GRanges(1, IRanges(c(3,7,13), c(5,9,16)), strand=c('+','-','-'), seqinfo=Seqinfo("1", 25), name=c("A","B","C"))
    df = gr.flatten(gr)
    expect_equal(as.character(class(df)), "data.frame")
    expect_identical(colnames(df), c("start", "end", "name"))
    expect_identical(df$start, c(1,4,7))

    df = gr.flatten(gr, gap=5)
    expect_equal(df$start[2] - df$end[1], 6)

})


test_that('gr.stripstrand', {

    expect_identical(as.character(strand(gr.stripstrand(gr))), c('*', '*', '*'))

})


test_that('gr.pairflip', {

    expect_identical(as.character(strand(gr.pairflip(gr)[[1]])), c('+', '-'))
    expect_identical(as.character(strand(gr.pairflip(gr)[[2]])), c('-', '+'))
    expect_identical(as.character(strand(gr.pairflip(gr)[[3]])), c('-', '+'))

})


test_that('gr.strandflip', {

    expect_identical(as.character(strand(gr.strandflip(gr))), c("-","+","+"))

})


test_that("gr.tile", {

    expect_identical(start(gr.tile(gr, w=3)), c(3L, 7L, 13L, 16L))
    expect_equal(length(gr.tile(GRanges())), 0)

})


test_that("gr.tile.map", {

    gr1 <- gr.tile(GRanges(1, IRanges(1,100)), width=10)
    gr2 <- gr.tile(GRanges(1, IRanges(1,100)), width=5)
    gg <- gr.tile.map(gr1, gr2, verbose=TRUE)
    expect_equal(length(gg), 10)
    expect_equal(length(unlist(gg)), 20)

})


test_that("gr.val", {

    gr <- GRanges(1, IRanges(1e6, 2e6))
    expect_equal(colnames(mcols(gr.val(gr, example_genes, val = 'name'))), "name")

})


test_that('gr.duplicated', {

    gr = GRanges(c(1,1,1), IRanges(c(2,5,5), width=1), val=c(1,2,3))
    
    expect_identical(gr.duplicated(gr), c(FALSE, FALSE, TRUE))
    expect_identical(gr.duplicated(gr, by='val'), c(FALSE, FALSE, FALSE))
    
})


test_that('gr.dice', {
    
    gr  = GRanges(1, IRanges(c(3,7,13), c(5,9,16)), strand=c('+','-','-'), seqinfo=Seqinfo("1", 25), name=c("A","B","C"))
    expect_equal(length(gr.dice(gr)[[3]]), 4)

})


test_that('gr.dist', {

    gr  = GRanges(1, IRanges(c(3,7,13), c(5,9,16)), strand=c('+','-','-'), seqinfo=Seqinfo("1", 25), name=c("A","B","C"))
    gr2 = GRanges(1, IRanges(c(1,9), c(6,14)), strand=c('+','-'), seqinfo=Seqinfo("1", 25), field=c(1,2))
    m = gr.dist(gr, gr2, ignore.strand=TRUE)
    expect_equal(m[1,2], 3)

})


## grl.stripnames; not exported in dev


## rle.query
test_that('rle.query', {

    gg = grl.in(grl.hiC[1:100], example_genes)
    expect_equal(length(gg), 100)

})


test_that('grl.in', {

    gg = grl.in(grl.hiC[1:100], example_genes)
    expect_equal(length(gg), 100)

})


test_that("grl.unlist", {

    gg <- grl.unlist(grl.hiC)
    expect_equal(length(gg), length(grl.hiC)*2)
    expect_equal(max(mcols(gg)$grl.iix), 2)
    expect_equal(max(mcols(gg)$grl.ix), length(grl.hiC))

})


test_that("grl.pivot", {

    gg <- grl.pivot(grl.hiC)
    expect_equal(as.character(class(gg)), "GRangesList")
    expect_equal(length(gg),2)
    expect_equal(length(gg[[1]]), 10000)

})


test_that("rrbind", {
    
    gr  = GRanges(1, IRanges(c(3,7,13), c(5,9,16)), strand=c('+','-','-'), seqinfo=Seqinfo("1", 25), name=c("A","B","C"))
    gr2 = GRanges(1, IRanges(c(1,9), c(6,14)), strand=c('+','-'), seqinfo=Seqinfo("1", 25), field=c(1,2))
    expect_that(ncol(rrbind(mcols(gr), mcols(gr2))) > 0, is_true())
    expect_equal(ncol(rrbind(mcols(gr), mcols(gr2), union=FALSE)), 0)

})


## gr.sub


## seg2gr


## standardize_segs


test_that("gr.nochr",{

    expect_identical(gr.nochr(gr.chr(example_genes)), example_genes)

})

## gr.findoverlaps() tests
test_that("gr.findoverlaps", {

    example_dnase = GRanges(1, IRanges(c(562757, 564442, 564442), c(563203, 564813, 564813)), strand = c("-", "+", "+"))
    example_genes = GRanges(2, IRanges(c(233101, 233101, 231023, 231023, 229966), c(233229, 233229, 231191, 231191, 230044)), strand = c("-"), type = c("exon", "CDS", "exon", "CDS", "exon"))
    fo = suppressWarnings(gr.findoverlaps(example_genes, example_dnase))   ## The 2 combined objects have no sequence levels in common. (Use suppressWarnings() to suppress this warning.)
    expect_equal(ncol(mcols(fo)), 0)
    expect_that(length(fo) == 0, is_true())

    expect_equal(length(gr.findoverlaps(example_genes, GRanges())), 0)
    expect_equal(length(gr.findoverlaps(GRanges(), GRanges())), 0)
    expect_equal(length(gr.findoverlaps(GRanges(), example_dnase)), 0)

})

test_that("gr.findoverlaps, return as data.table", {

    example_dnase = GRanges(1, IRanges(c(562757, 564442, 564442), c(563203, 564813, 564813)), strand = c("-", "+", "+"))
    example_genes = GRanges(2, IRanges(c(233101, 233101, 231023, 231023, 229966), c(233229, 233229, 231191, 231191, 230044)), strand = c("-"), type = c("exon", "CDS", "exon", "CDS", "exon"))
    expect_error(gr.findoverlaps(example_genes, example_dnase, return.type = "data.frame"))

    fo = gr.findoverlaps(example_dnase[1:3], example_dnase, return.type = 'data.table')
    expect_identical(colnames(fo), c("start", "end", "query.id", "subject.id", "seqnames", "strand"))

})

test_that("gr.findoverlaps chunk", {

    example_dnase = GRanges(1, IRanges(c(562757, 564442, 564442), c(563203, 564813, 564813)), strand = c("-", "+", "+"))
    example_genes = GRanges(2, IRanges(c(233101, 233101, 231023, 231023, 229966), c(233229, 233229, 231191, 231191, 230044)), strand = c("-"), type = c("exon", "CDS", "exon", "CDS", "exon"))

    fo  = suppressWarnings(gr.findoverlaps(example_genes, example_dnase))
    fo2 = suppressWarnings(gr.findoverlaps(example_genes, example_dnase, max.chunk = 1e7, verbose=TRUE))
    expect_identical(fo, fo2)

})

test_that("gr.findoverlaps, input data.table", {

    example_genes = GRanges(2, IRanges(c(233101, 233101, 231023, 231023, 229966), c(233229, 233229, 231191, 231191, 230044)), strand = c("-"), type = c("exon", "CDS", "exon", "CDS", "exon"))
    expect_equal(class(suppressWarnings(gr.findoverlaps(gr2dt(example_genes), example_genes, return.type='GRanges')))[1], "GRanges")  
    expect_equal(class(suppressWarnings(gr.findoverlaps(gr2dt(example_genes), example_genes)))[1], "data.table")
    expect_equal(class(suppressWarnings(gr.findoverlaps(gr2dt(example_genes), example_genes, max.chunk = 1e7)))[1], "data.table")

})

test_that("gr.findoverlaps ignore.strand", {

    ## make a stranded DNAase track (for testing only)
    example_dnase2 = example_dnase
    set.seed(137)
    strand(example_dnase2) <- ifelse(runif(length(example_dnase)) > 0.5, '+', '-')

    ## get the overlaps with the original unstranded, and with ignore.strand
    fo1 <- suppressWarnings(gr.findoverlaps(example_dnase, example_genes))
    fo2 <- suppressWarnings(gr.findoverlaps(example_dnase2, example_genes, ignore.strand=TRUE))
    expect_identical(fo1, fo2)

    ## make sure no strands overlap
    fo1 <- suppressWarnings(gr.findoverlaps(example_dnase2, example_genes, ignore.strand=FALSE))
    expect_that(!any(strand(example_dnase2)[fo1$query.id] != strand(example_genes)[fo1$subject.id]), is_true())

})

test_that("gr.findoverlap by", {
    
    example_genes = GRanges("2:23000-255555") 
    e1 <- example_genes
    e2 <- example_dnase
    e1$bin <- e2$bin <- 1
    expect_error(gr.findoverlaps(example_genes, example_dnase, by = "dummy"))
    expect_that(length(gr.findoverlaps(e1, e2, by = "bin")) == 0, is_true())

})


## grl.eval


## gr.merge


## gr.disjoin


## gr.in 


## gr.sum


## gr.collapse


test_that("gr.match", {
    
    ## gives back overlapping matches
    gr1 = GRanges(1, IRanges(c(10,20), width=5), strand=c("+", "-"))
    gr2 = GRanges(1, IRanges(c(8,18, 100), width=5), strand=c("-", "+", "+"))

    expect_identical(suppressWarnings(gr.match(gr1, gr2)), c(1L,2L))

    ## ignore strand is successfully passed
    expect_error(gr.match(gr1, gr2, ignore.strand = FALSE))

})


## %+%


## %-%


## %&%


## %&&%


## %O%


## %OO%


## %o%


## %oo%


## %N%


## %NN%


## %_%
test_that('%_% works', {

    gr1 <- GRanges(1, IRanges(10,20), strand="+")
    gr2 <- GRanges(1, IRanges(15,25), strand="-")
    gr3 <- GRanges("1:1-15")
    expect_equal(width(gr1 %_% gr2), 5)
    expect_equal(width(gr1 %_% gr3), 5)
    
})


## %Q%
test_that('%Q% works', {

    testset <- GRanges(seqnames = Rle(1,c(5)) , ranges = IRanges(1:5 , end = 2000:2004) , strand = Rle(strand(c("+")) , c(5)) , mean = c(1, -2 , -5 , 5 ,6))
    expect_equal(length(testset %Q% (mean > 0)) , 3)

})

test_that('second test that %Q% works', {
    
    foo = GRanges(seqnames = Rle(1,c(5)) , ranges = IRanges(1:5 , end = 2000:2004) , strand = Rle(strand(c("+")) , c(5)) , mean = c(50, 200, 300, 400, 500))
    expect_equal(length(foo %Q% (mean == 300)) , 1)
    
})

test_that('third test that %Q% works', {
    
    foo = GRanges(seqnames = Rle(1,c(5)) , ranges = IRanges(1:5 , end = c(5, 10, 15, 20, 25)) , strand = Rle(strand(c("+", "-", "-", "-", "-"))) , mean = c(1, -2 , -5 , 5 ,6))
    expect_equal(length(foo %Q% (strand == "-") ) , 4)
    
})


## %^%
test_that('%^% works', {

    testset <- GRanges(seqnames = Rle(1,c(5)) , ranges = IRanges(1:5 , end = 2000:2004) , strand = Rle(strand(c("+")) , c(5)) , mean = c(1, -2 , -5 , 5 ,6))
    expect_equal(length(testset %^% testset), 5)

})


## %$%


## %*%


## %**%


## %^^%


## %$$%


## gr.setdiff



## XT plans to re-write, Jan 18
## test_that('ra.overlaps', {
##
##    gr = GRanges(1, IRanges(c(3,7), c(5,9)), strand=c('+','-'))
##    gr1 = GRanges(1, IRanges(c(10,20), width=5), strand=c("+", "-"))
##    gr2 = GRanges(1, IRanges(c(1,9), c(6,14)), strand=c('+','-'))
##    grl1 = GRangesList("gr"=gr, "gr1"=gr1)
##    grl2 = GRangesList("gr1"=gr1, "gr2"=gr2)
##    foobar = suppressWarnings(grl.bind(grl1, grl2))
##    ro = ra.overlaps(grl1, grl2)
##    expect_equal(class(ro), "matrix")
##    expect_equal(nrow(ro), 2)
##    expect_equal(ncol(ro), 2)
##    expect_equal(nrow(ra.overlaps(grl2, grl2)), length(grl2))
##
## })
## test_that("ra.overlaps handles empty",{
##     
##     ## test empty inputs and no overlaps inputs
##     gr = GRanges(1, IRanges(c(10,20), width=5), strand=c("+", "-"))
##     grl1 = GRangesList("gr1" = gr)
##     expect_equal(ra.overlaps(GRangesList(), grl1)[1], NA)
##     expect_equal(ra.overlaps(grl2[2:3], grl1)[1], NA)
##     
## })  
## test_that("ra.overlaps handles wrong signs", {
## 
##     ## make one that overlaps, but wrong signs
##     gr = GRanges(1, IRanges(c(3,7), c(5,9)), strand=c('+','-'))
##     gr1 = GRanges(1, IRanges(c(10,20), width=5), strand=c("+", "-"))
##     gr2 = GRanges(1, IRanges(c(1,9), c(6,14)), strand=c('+','-'))
##     grl1 = GRangesList("gr" = gr, "gr1" = gr1, "gr2" = gr2)
##     grl3 <- grl1[2]
##     strand(grl3[[1]]) <- c("+", "-")
##     expect_equal(ra.overlaps(grl3, grl2)[1], NA)
## 
##  })



test_that('ra.merge' {

    gr1 = GRanges(1, IRanges(1:10, width = 1), strand = rep(c('+', '-'), 5))
    gr2 = GRanges(1, IRanges(4 + 1:10, width = 1), strand = rep(c('+', '-'), 5))
    ra1 = split(gr1, rep(1:5, each = 2))
    ra2 = split(gr2, rep(1:5, each = 2))

    ram = ra.merge(ra1, ra2)
    ## 'values(ram)' shows the metadata with TRUE / FALSE flags
    expect_equal(values(ram)[, 1], c(TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE))


  
    ##ram2 = ra.merge(ra1, ra2, pad = 5) # more inexact matching results in more merging
    ##values(ram2)

    ##ram3 = ra.merge(ra1, ra2, ind = TRUE) #indices instead of flags
    ##values(ram3)

})



## XT, ra.dedup
## Jan 18, correspondence from XT 'leave that internal for now'


## XT, ra.duplicated
## Jan 18, correspondence from XT 'leave that internal for now'


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


test_that("parse.gr", {

    gr_example = parse.gr(c('1:1e6-5e6+', '2:2e6-5e6-'))
    expect_equal(width(gr_example[1]), 4000001)
    expect_equal(width(gr_example[2]), 3000001)

})


test_that("parse.grl", {

    grl_example = parse.grl(c('1:1e6-5e6+;5:10-2000', '2:2e6-5e6-;10:100231321-100231399'))
    expect_equal(width(grl_example[[1]][1]), 4000001)
    expect_equal(width(grl_example[[1]][2]), 1991)
    expect_equal(width(grl_example[[2]][1]), 3000001)
    expect_equal(width(grl_example[[2]][2]), 79)

})


## anchorlift


## XT Yao function
#### gr.breaks










