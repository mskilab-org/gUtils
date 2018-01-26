library(gUtils)
library(BSgenome.Hsapiens.UCSC.hg19)

library(testthat)

Sys.setenv(DEFAULT_BSGENOME = "BSgenome.Hsapiens.UCSC.hg19::Hsapiens")


context("unit testing gUtils operations")

gr = GRanges(1, IRanges(c(3,7,13), c(5,9,16)), strand=c('+','-','-'), seqinfo=Seqinfo("1", 25), name=c("A","B","C"))
gr2 = GRanges(1, IRanges(c(1,9), c(6,14)), strand=c('+','-'), seqinfo=Seqinfo("1", 25), field=c(1,2))
dt = data.table(seqnames=1, start=c(2,5,10), end=c(3,8,15))


test_that("hg_seqlengths()", {
    
    Sys.setenv(DEFAULT_BSGENOME = "")
    expect_warning(hg_seqlengths(genome=NULL))
    ## throw error with incorrect DEFAULT_BSGENOME 
    Sys.setenv(DEFAULT_BSGENOME = "incorrect")
    expect_error(hg_seqlengths())
    ## set DEFAULT_BSGENOME as hg19
    Sys.setenv(DEFAULT_BSGENOME = "BSgenome.Hsapiens.UCSC.hg19::Hsapiens")
    expect_identical(as.numeric(length(hg_seqlengths())), 25)
    ee = structure(names="1", 249250621L)
    expect_identical(hg_seqlengths(Hsapiens)[1], ee)
    expect_equal(names(hg_seqlengths(Hsapiens, chr=TRUE)[1]), "chr1")
    expect_equal(length(hg_seqlengths(Hsapiens, include.junk = TRUE)), 93)

})


test_that("gr2dt", {

    example_genes = GRanges(2, IRanges(c(233101, 233101, 231023, 231023, 229966), c(233229, 233229, 231191, 231191, 230044)), strand = c("-"), type = c("exon", "CDS", "exon", "CDS", "exon"))
    expect_identical(colnames(gr2dt(gr)), c("seqnames", "start", "end", "strand", "width", "name"))
    expect_equal(nrow(gr2dt(gr)), length(gr))
    ## using subjectdt
    subjectdt <- gr2dt(example_genes)
    expect_equal(!any(subjectdt$start!=start(example_genes)), TRUE)
    expect_equal(!any(subjectdt$end!=end(example_genes)), TRUE)
    expect_equal(!any(subjectdt$width!=width(example_genes)), TRUE)
    expect_equal(all(as.character(subjectdt$strand) == as.character(strand(example_genes))), TRUE)
    expect_equal(all(subjectdt$seqnames == as.vector(seqnames(example_genes))), TRUE)
    ## check 'if (any(duplicated(names(x))))'
    ## check 'if (is.null(out)){'
    expect_equal(gr2dt(NULL), data.table())
    
})



test_that("gr.start", {
    
    expect_equal(gr.start(NULL), NULL)  ### check 'if (length(x)==0){ return(x) }' 
    seqlengths(gr) = 0
    expect_warning(gr.start(gr)) ## check 'if (any(seqlengths(x)==0) | any(is.na(seqlengths(x))))' warning :: warning('Warning: Check or fix seqlengths, some are equal 0 or NA, may lead to negative widths')
    gr = GRanges(1, IRanges(c(3,7,13), c(5,9,16)), strand=c('+','-','-'), seqinfo=Seqinfo("1", 25), name=c("A","B","C"))
    expect_identical(start(gr.start(gr)), c(3L,7L,13L))
    expect_identical(end(gr.start(gr)), c(3L,7L,13L))
    expect_identical(end(gr.start(gr, width=100)),  c(5L,9L,16L))
    expect_identical(end(gr.start(gr, width=100, clip=FALSE)),  c(25L,25L,25L))
    expect_identical(suppressWarnings(end(gr.start(gr, width=100, force=TRUE))),  c(5L,9L,16L))
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
    ## expect_error(suppressWarnings(dt2gr(dt)))    ### warning within error---warning: coercing to GRanges via non-standard columns
    ## as.integer(seqnames(seqinfo(dt2gr(dt, seqlengths=NULL, seqinfo=NULL)))
    ## check stop("Error: Needs to be data.table or data.frame")
    expect_error(dt2gr(matrix()))
    expect_error(dt2gr(GRanges()))

})


test_that("gr.end", {

    expect_equal(gr.end(NULL), NULL)  ### check 'if (length(x)==0){ return(x) }' 
    seqlengths(gr) = 0
    expect_warning(gr.end(gr)) ## check 'if (any(seqlengths(x)==0) | any(is.na(seqlengths(x))))' warning :: warning('Warning: Check or fix seqlengths, some are equal 0 or NA, may lead to negative widths')
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
    expect_error(gr.rand(c(3,500000000), si)) ## Error: Allocation failed. Supplied widths are likely too large
    ## check 'if (!is(genome, 'Seqinfo')){'
    ref = si2gr(BSgenome.Hsapiens.UCSC.hg19::Hsapiens)
    expect_equal(length(gr.rand(c(3,5), ref)), 2)

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
    ### checks 'if (!inherits(gr, 'GRanges')){ gr = si2gr(gr) }'
    expect_equal(length(gr.sample(si, 10, 1)), 10)   
    expect_equal(length(gr.sample(reduce(example_genes), 10, k=3)), 3)

    ## check 'if (length(k)==1)''
    ###
    ## query width less than output
    ## expect_error(gr.sample(gr.start(example_genes), c(1:3), k=5))

    gg <- suppressWarnings(gr.sample(example_genes[1:5], c(2,2,3,4,5), k=2))    ### expect warning: longer object length is not a multiple of shorter object length
    expect_equal(length(gg), 5)
    expect_equal(sum(width(gg)), 16)
    ## check ' stop('Error: Input territory has zero regions of sufficient width')'
    expect_error(gr.sample(GRanges(), 10, k=1))

})


test_that("gr.sample without replace", {
    
    gene1 = GRanges("3:3540000-234329000")
    gene2 = GRanges("2:24440-30000")
    gene3 = GRanges("2:278444-321000")
    gene4 = GRanges("3:1000-1500")
    example_genes = suppressWarnings(c(gene1, gene2, gene3, gene4))  ##  The 2 combined objects have no sequence levels in common. (Use suppressWarnings() to suppress this warning
    gg = gr.sample(reduce(example_genes), 10, k=1, replace=FALSE)
    expect_equal(unique(width(gg)), 10)

    expect_equal(width(gr.sample(example_genes[1:3], 1e7, k=1, replace=FALSE)),  10000000)

    gg <- gr.sample(example_genes, c(2,2,3,4), k=2, replace=FALSE)
    expect_equal(length(gg), 4)
    expect_equal(sum(width(gg)), 11)

})


test_that("si2gr", {

    gg = si2gr(si, strip.empty = TRUE)
    expect_equal(start(gg)[3], 1)
    expect_equal(end(gg)[1], 249250621)
    expect_equal(as.character(strand(gg)[1]), "+")
    ## check 'if (is(si, 'BSgenome'))'
    expect_equal(length(si2gr(BSgenome.Hsapiens.UCSC.hg19::Hsapiens)), 93)   
    ## check ' else if (!is(si, 'Seqinfo'))'
    expect_equal(length(si2gr(GRanges(si))), 25)   

})


test_that("grbind", {

    example_dnase = GRanges(1, IRanges(c(562757, 564442, 564442), c(563203, 564813, 564813)), strand = c("-", "+", "+"))
    example_genes = GRanges(2, IRanges(c(233101, 233101, 231023, 231023, 229966), c(233229, 233229, 231191, 231191, 230044)), strand = c("-"), type = c("exon", "CDS", "exon", "CDS", "exon"))
    expect_equal(length(suppressWarnings(grbind(example_genes, example_dnase))) > 0, TRUE)
    expect_equal(grbind(0, 1, 2, 3), NULL)
    ## check ' grs <- c(x, list(...))'
    expect_equal(width(grbind(GRanges(), GRanges('1:1-10'),GRanges(),GRanges(),GRanges('1:100-200'),GRanges())[2]), 101)
    ## check ' if (is.null(tmp)){'
    expect_equal(grbind(data.frame(),data.frame(),data.frame()), NULL)

})


test_that("grl.bind", {

    grl.hiC2 = grl.hiC[1:20]
    mcols(grl.hiC2)$test = 1
    ##gg <- grl.bind(grl.hiC2, grl.hiC[1:30])
    expect_equal(length(grl.bind(grl.hiC2, grl.hiC[1:30])), 50)
    expect_equal(colnames(mcols(grl.bind(grl.hiC2, grl.hiC[1:30]))), "test")

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

    ## check 'if (any(is.na(seqlengths(gr)))){'
    seqlengths(grl.hiC) = NA
    expect_equal(length(streduce(grl.hiC, pad=10)), 19701)

})




test_that("gr.string", {

    expect_that(grepl(":", gr.string(example_genes)[1]), is_true())
    expect_that(grepl("-", gr.string(example_genes)[1]), is_true())
    expect_that(grepl("(+|-)", gr.string(example_genes)[1]), is_true())
    ## check 'if (length(gr)==0){'
    ## add.chr
    expect_equal(gr.string(example_genes, add.chr=TRUE)[1], 'chr1:69090-70008+')
    ## mb
    expect_equal(gr.string(example_genes[1], mb = TRUE), '1:0.069-0.07+')
    ## round
    expect_equal(gr.string(example_genes[1], round = 1, mb=TRUE), '1:0.1-0.1+')
    ## other.cols
    expect_equal(gr.string(example_genes[1], other.cols='name'), '1:69090-70008+  OR4F5')
    ## pretty
    expect_equal(gr.string(example_genes, pretty=TRUE)[1], '1:69,090-70,008+')

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
    ## check 'if (class(grl) == "GRanges"){'
    expect_equal(grl.string(example_genes[1]), '1:69090-70008+')
    ## check 'error handling'
    expect_error(grl.string('foo')) ##  Error: Input must be GRangesList (or GRanges, which is sent to gr.string)
    ## check 'else{ nm = 1:length(grl) ! 1138 }'
    names(grl1) = NULL
    expect_equal(as.character(grl.string(grl1[2])), '9:140100229-140100229-,19:24309057-24309057-')


})


test_that("gr.fix", {

    gg = GRanges(c("X",1), IRanges(c(1,2), width=1))
    expect_equal(length(seqlengths(gr.fix(gg, si))), 25)
    expect_equal(length(seqlengths(gr.fix(gg, BSgenome.Hsapiens.UCSC.hg19::Hsapiens))), 95)
    ## check 'if (is(gr, 'GRangesList')){'
    expect_equal(as.integer(seqnames(gr.fix(grl1[1])[1])), 5)

})


test_that("gr.fix with null genome", {

    gg <- GRanges(c("X",1), IRanges(c(1,2), width=1))
    es <- structure(c(1L,2L), names=c("X","1"))
    expect_identical(seqlengths(seqinfo(gr.fix(gg))), es)
    
    es <- structure(c("gen","gen"), names=c("X","1"))
    expect_identical(genome(seqinfo(gr.fix(gg, gname='gen'))), es)

})


test_that("gr.flatten", {

    ## check 'if (length(gr) == 0)'
    foo = gr.flatten(GRanges())
    expect_true(is(foo, 'data.frame'))
    expect_equal(length(foo), 0)
    ##
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


test_that('gr.flipstrand', {

    expect_identical(as.character(strand(gr.flipstrand(gr))), c("-","+","+"))
    expect_error(gr.flipstrand(data.frame()))
    expect_equal(length(gr.flipstrand(GRanges())), 0)

})


test_that("gr.tile", {

    expect_identical(start(gr.tile(gr, w=3)), c(3L, 7L, 13L, 16L))
    expect_equal(length(gr.tile(GRanges())), 0)
    ## check 'if (is(gr, 'data.table'))'
    expect_equal(length(gr.tile(dt, w=3)), 5)
    ## check 'else if (!is(gr, 'GRanges')){'
    expect_equal(length(gr.tile(si['M'], w=3)), 5524)


})


test_that("gr.tile.map", {

    gr1 = gr.tile(GRanges(1, IRanges(1,100)), width=10)
    gr2 = gr.tile(GRanges(1, IRanges(1,100)), width=5)
    gg = gr.tile.map(gr1, gr2, verbose=TRUE)
    expect_equal(length(gg), 10)
    expect_equal(length(unlist(gg)), 20)

})


test_that("gr.val", {

    gr = GRanges(1, IRanges(1e6, 2e6))
    expect_equal(colnames(mcols(gr.val(gr, example_genes, val = 'name'))), "name")
    ## check val = NULL
    expect_equal(width(gr.val(gr, example_genes)), 1000001)
    ## check mean 
    expect_equal(gr.val(gr, example_genes, mean =TRUE)$value, 1)
    expect_equal(gr.val(gr, example_genes, mean =FALSE)$value, 41)
    ## check 'max.slice'
    ## check's  'if (length(query)>max.slice)'
    expect_equal(length(gr.val(example_dnase, example_genes, max.slice = 50)), 10000)
    ## check 'if (inherits(target, 'GRangesList'))'
    expect_equal(as.integer(seqnames(gr.val(grl1, grl.unlist(grl2))[[2]])), c(9, 19))
    ## check mc.cores
    expect_equal(width(gr.val(gr, example_genes, mc.cores=2)), 1000001)
    ## check FUN
    expect_equal(as.data.frame(gr.val(gr, example_genes, FUN = function(x, w, na.rm = FALSE){ return(w*(x**2))}))$value, 68671)

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
    ## check 'if (is.null(gr2)){ gr2 = gr1 }'
    expect_equal(dim(gr.dist(gr)), c(3, 3))

})


## grl.stripnames; not exported in dev


## rle.query
test_that('rle.query', {

    chr1_Rle = Rle(10:1, 1:10)
    chr2_Rle = Rle(10:1, 1:10)
    example_rlelist = RleList( chr1=chr1_Rle, chr2=chr2_Rle)
    expect_equal(length(rle.query(example_rlelist, gr)), 10)
    expect_equal(length(rle.query(example_rlelist, gr2)), 12)  
    ## check 'if (is(query.gr, 'GRangesList'))'
    expect_equal(length(rle.query(example_rlelist, grl2)), 251)
    ## check 'chunksize'
    expect_equal(length(rle.query(example_rlelist, grl2[1:15], chunksize=5)), 15)


})


test_that('grl.in', {

    gg = grl.in(grl.hiC[1:100], example_genes)
    expect_equal(length(gg), 100)
    ## check  'if (length(grl)==0){''
    expect_equal(grl.in(GRangesList(),  example_genes), logical(0))
    ## check  'if (length(windows)==0){''
    expect_false(all(grl.in(grl.hiC[1:100], GRanges())))

})


test_that("grl.unlist", {

    gg <- grl.unlist(grl.hiC)
    expect_equal(length(gg), length(grl.hiC)*2)
    expect_equal(max(mcols(gg)$grl.iix), 2)
    expect_equal(max(mcols(gg)$grl.ix), length(grl.hiC))
    ## check 'if (length(grl) == 0)'
    expect_equal(length(grl.unlist(GRangesList())), 0)
    ## check 'if (is(grl, 'GRanges'))'
    expect_equal(grl.unlist(gr2)[1]$grl.ix, 1)


})


test_that("grl.pivot", {

    gg <- grl.pivot(grl.hiC)
    expect_equal(as.character(class(gg)), "GRangesList")
    expect_equal(length(gg),2)
    expect_equal(length(gg[[1]]), 10000)
    ## check 'if (length(x) == 0)'
    expect_equal(length(grl.pivot(GRangesList())), 2)
    expect_equal(length(grl.pivot(GRangesList())[[1]]), 0)
    expect_equal(length(grl.pivot(GRangesList())[[2]]), 0)

})


test_that("rrbind", {
    
    gr  = GRanges(1, IRanges(c(3,7,13), c(5,9,16)), strand=c('+','-','-'), seqinfo=Seqinfo("1", 25), name=c("A","B","C"))
    gr2 = GRanges(1, IRanges(c(1,9), c(6,14)), strand=c('+','-'), seqinfo=Seqinfo("1", 25), field=c(1,2))
    expect_that(ncol(rrbind(mcols(gr), mcols(gr2))) > 0, is_true())
    expect_equal(ncol(rrbind(mcols(gr), mcols(gr2), union=FALSE)), 0)
    ## check 'if (any(mix <- sapply(dfs, class) == 'matrix')){'
    expect_equal(dim(rrbind( mcols(gr), matrix(gr2dt(gr2))))[1], 6)
    expect_equal(dim(rrbind( mcols(gr), matrix(gr2dt(gr2))))[2], 1)
    ## check 'if (is.null(rout)){'
    expect_equal(rrbind(mcols(GRanges()), mcols(GRanges()), mcols(GRanges())), data.frame())

})


## gr.sub
test_that('gr.sub', {
    
    gr1  = GRanges('chr1', IRanges(c(3,7,13), c(5,9,16)), strand=c('+','-','-'), seqinfo=Seqinfo("chr1", 25), name=c("A","B","C"))
    expect_error(gr.sub(gr1), NA)  ## check works
    ## check 'if (is.null(tmp.gr))'
    expect_equal(length(gr.sub(GRanges())), 0)

})

## seg2gr
test_that('seg2gr', {

    expect_equal(width(seg2gr(dt)), c(2, 4, 6))    

})

## standardize_segs
test_that('standardize_segs', {

    expect_error(standardize_segs(gr2))
    
})


test_that("gr.nochr", {

    expect_identical(gr.nochr(gr.chr(example_genes)), example_genes)

})

## gr.findoverlaps() tests
test_that("gr.findoverlaps", {

    example1 = GRanges(1, IRanges(c(562757, 564442, 564442), c(563203, 564813, 564813)), strand = c("-", "+", "+"))
    example2 = GRanges(2, IRanges(c(233101, 233101, 231023, 231023, 229966), c(233229, 233229, 231191, 231191, 230044)), strand = c("-"), type = c("exon", "CDS", "exon", "CDS", "exon"))
    fo = suppressWarnings(gr.findoverlaps(example1, example2))   ## The 2 combined objects have no sequence levels in common. (Use suppressWarnings() to suppress this warning.)
    expect_equal(ncol(mcols(fo)), 0)
    expect_that(length(fo) == 0, is_true())
    ## qcol
    expect_equal(length(gr.findoverlaps(example_genes, example_dnase, qcol = 'exonCount')), 4184)
    ## scol
    expect_equal(round(mean(gr.findoverlaps(example_genes, example_dnase, scol='pValue')$pValue), 2), 66.12)
    expect_equal(length(gr.findoverlaps(example_genes, GRanges())), 0)
    expect_equal(length(gr.findoverlaps(GRanges(), GRanges())), 0)
    expect_equal(length(gr.findoverlaps(GRanges(), example_dnase)), 0)
    ## check input error catching
    expect_error(gr.findoverlaps(1, 2)) ##   Error: Both subject and query have to be GRanges or data.table
    ## check 'if ((as.numeric(length(query)) * as.numeric(length(subject))) > max.chunk)'
    ## here == 188120000
    expect_equal(length(gr.findoverlaps(example_genes, example_dnase, max.chunk = 1e6)), 4184)
    ## check 'if (is.null(h))'
    expect_equal(length(gr.findoverlaps(GRanges(), GRanges())), 0)


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
test_that('grl.eval', {

    expect_equal(grl.eval(grl1, bin**2)[1:5], c(100, 100, 676, 676, 2401))
    expect_equal(grl.eval(grl1, width-5)[1:5], rep(-4, 5))

})


## gr.merge
test_that('gr.merge', {

    sv1 = grl.unlist(grl1)
    sv2 = grl.unlist(grl2)
    ## default
    expect_equal(length(suppressWarnings(gr.merge(sv1, sv2))), 3)
    expect_equal(suppressWarnings(gr.merge(sv1, sv2)$query.id), c(367, 499, 500))
    expect_equal(suppressWarnings(gr.merge(sv1, sv2)$subject.id), c(443, 1, 2))
    ## all = TRUE
    expect_equal(length(suppressWarnings(gr.merge(sv1, sv2, all=TRUE))), 999)
    ## by
    expect_equal(length((gr.merge(sv1, sv2, all=TRUE, by='bin'))), 2)
})


## gr.disjoint
test_that('gr.disjoin', {

    sv1 = grl.unlist(grl1)
    sv2 = grl.unlist(grl2)
    ## cf. http://bioconductor.org/packages/devel/bioc/vignettes/GenomicRanges/inst/doc/GenomicRangesIntroduction.pdf
    ## example: 'disjoin(g)'
    gr = GRanges(
        seqnames = Rle(c("chr1", "chr2", "chr1", "chr3"), c(1, 3, 2, 4)),
        ranges = IRanges(101:110, end = 111:120, names = head(letters, 10)),
        strand = Rle(strand(c("-", "+", "*", "+", "-")), c(1, 2, 2, 3, 2)),
        score = 1:10,
        GC = seq(1, 0, length=10))
    singles = split(gr, names(gr))
    g = gr[1:3]
    g = append(g, singles[[10]])
    expect_equal(suppressWarnings(gr.disjoin(g)$score), c(1, 2, 2, 3, 10))
    expect_equal(round(gr.disjoin(g)$GC, 3), c(1.000, 0.889, 0.889, 0.778, 0.000))

})


## gr.in 
test_that('gr.in', {

    sv1 = grl.unlist(grl1)
    sv2 = grl.unlist(grl2)
    table(gr.in(sv1, sv2))
    expect_equal(as.numeric(table(gr.in(sv1, sv2))[1]), 497)
    expect_equal(as.numeric(table(gr.in(sv1, sv2))[2]), 3)

})


## gr.sum
test_that('gr.sum', {
    
    ## cf. http://bioconductor.org/packages/devel/bioc/vignettes/GenomicRanges/inst/doc/GenomicRangesIntroduction.pdf
    gr = GRanges(
        seqnames = Rle(c("chr1", "chr2", "chr1", "chr3"), c(1, 3, 2, 4)),
        ranges = IRanges(101:110, end = 111:120, names = head(letters, 10)),
        strand = Rle(strand(c("-", "+", "*", "+", "-")), c(1, 2, 2, 3, 2)),
        score = 1:10,
        GC = seq(1, 0, length=10))
    ## default
    expect_equal(length(gr.sum(gr)), 20)
    ## 'field' argument
    expect_error(length(gr.sum(gr, field = 'CG')))
    expect_equal(length(gr.sum(gr, field = 'GC')), 19)
    expect_equal(round(gr.sum(gr, field = 'GC')$GC[1:6], 3), c(0.000, 1.000, 1.556, 2.000, 1.000, 0.444))
    ## 'field' argument and 'mean' = TRUE
    expect_equal(round(gr.sum(gr, field = 'GC', mean=TRUE)$GC[1:6], 3), c(NaN, 1.000, 0.778, 0.667, 0.500, 0.444))
    singles = split(gr, names(gr))
    g = gr[1:3]
    g = append(g, singles[[10]])
    expect_equal(width(gr.sum(g)), c(100, 11, 101, 1, 10, 1, 109, 11))
    sv1 = grl.unlist(grl1)
    expect_equal(length(gr.sum(sv1)), 1025)

})


## gr.collapse  
test_that('gr.collapse', {
    
    ## cf. http://bioconductor.org/packages/devel/bioc/vignettes/GenomicRanges/inst/doc/GenomicRangesIntroduction.pdf
    ## example: 'reduce(g)'
    gr = GRanges(
        seqnames = Rle(c("chr1", "chr2", "chr1", "chr3"), c(1, 3, 2, 4)),
        ranges = IRanges(101:110, end = 111:120, names = head(letters, 10)),
        strand = Rle(strand(c("-", "+", "*", "+", "-")), c(1, 2, 2, 3, 2)),
        score = 1:10,
        GC = seq(1, 0, length=10))
    singles = split(gr, names(gr))
    g = gr[1:3]
    g = append(g, singles[[10]])
    ## expect_equal(length(gr.collapse(gr)), 3) ## must fix
    expect_error(length(gr.collapse(gr)))
    ## expect_equal(length(gr.collapse(g)), 1)  ## must fix
    expect_error(length(gr.collapse(g)))
    sv1 = grl.unlist(grl1)
    ## expect_equal(length(gr.collapse(sv1)), 0)  ## must fix
    expect_error(length(gr.collapse(sv1)))

})


test_that("gr.match", {
    
    ## gives back overlapping matches
    gr1 = GRanges(1, IRanges(c(10,20), width=5), strand=c("+", "-"))
    gr2 = GRanges(1, IRanges(c(8,18, 100), width=5), strand=c("-", "+", "+"))

    expect_identical(suppressWarnings(gr.match(gr1, gr2)), c(1L,2L))

    ## ignore strand is successfully passed
    expect_error(gr.match(gr1, gr2, ignore.strand = FALSE))
    ## check 'if (length(query)>max.slice)'
    expect_equal(length(gr.match(gr1, gr2, max.slice=1)), 2)


})


test_that('%+% works', {
    
    expect_warning(gr %+% 20000) ## GRanges object contains 3 out-of-bound ranges located on sequence 1.
    shifted_gr = gr %+% 20000
    expect_equal(width(shifted_gr), c(3, 3, 4))
    expect_equal(start(shifted_gr), c(20003, 20007, 20013))
    expect_equal(end(shifted_gr), c(20005, 20009, 20016))
    gr_zero = gr %+% 0
    expect_equal(gr_zero, gr)
    ### expect_error(gr %+% -200) ## error: invalid class “IRanges” object: 'width(x)' cannot contain negative integers  ## MUST FIX. If length(shift) < width(), error
    expect_warning(gr %+% -3) ## warning: GRanges object contains 1 out-of-bound range located on sequence 1.
    gr_shifted_neg = gr %+% -3
    expect_equal(width(gr_shifted_neg), c(3, 3, 4))
    expect_equal(start(gr_shifted_neg), c(0, 4, 10))
    expect_equal(end(gr_shifted_neg), c(2, 6, 13))

})


test_that('%-% works', {
    
    expect_warning(gr %-% 20000) ## GRanges object contains 3 out-of-bound ranges located on sequence 1.
    shifted_gr = gr %-% 20000
    expect_equal(width(shifted_gr), c(3, 3, 4))
    expect_equal(start(shifted_gr), c(-19997, -19993, -19987))
    expect_equal(end(shifted_gr), c(-19995, -19991, -19984))
    gr_zero = gr %-% 0
    expect_equal(gr_zero, gr)
    ##expect_error(gr %-% -200) ## error: invalid class “IRanges” object: 'width(x)' cannot contain negative integers ## MUST FIX
    ## expect_error(gr %-% -4) ## error: invalid class “IRanges” object: 'width(x)' cannot contain negative integers ## MUST FIX
    gr_neg = gr %-% -3
    expect_equal(width(gr_neg), c(3, 3, 4))
    expect_equal(start(gr_neg), c(6, 10, 16))
    expect_equal(end(gr_neg), c(8, 12, 19))

})


## a %&% b # strand agnostic
## Return the subset of ranges in a that overlap with at least one range in b.
test_that('%&% works', {
    
    expect_equal(gr %&% gr2, gr)
    expect_equal(gr2 %&% gr, gr2)
    gr_A = GRanges(2, IRanges(c(233101, 233105, 231023, 231028, 229966), c(233229, 233229, 231191, 231191, 230044)), strand = c("-"))
    gr_B = GRanges(2, IRanges(c(233101, 333105, 331023, 331028, 329966), c(333229, 333229, 331191, 331191, 330044)), strand = c("+"))
    ## gr_A %&% gr_B
    expect_equal(length(gr_A %&% gr_B), 2)
    expect_equal(start(gr_A %&% gr_B), c(233101, 233105))
    expect_equal(end(gr_A %&% gr_B), c(233229, 233229))
    expect_equal(width(gr_A %&% gr_B), c(129, 125))
    expect_equal(as.vector(strand(gr_A %&% gr_B)), c('-', '-'))
    ## gr_B %&% gr_A
    expect_equal(length(gr_B %&% gr_A), 1)
    expect_equal(start(gr_B %&% gr_A), 233101)
    expect_equal(end(gr_B %&% gr_A), 333229)
    expect_equal(width(gr_B %&% gr_A), 100129)
    expect_equal(as.vector(strand(gr_B %&% gr_A)), '+')

})


## a %&&% b # strand specific
## Return the subset of ranges in a that overlap with at least one range in b.
test_that('%&&% works', {
    
    expect_equal(gr %&&% gr2, gr)
    expect_equal(gr2 %&&% gr, gr2)
    gr_A = GRanges(2, IRanges(c(233101, 233105, 231023, 231028, 229966), c(233229, 233229, 231191, 231191, 230044)), strand = c("-"))
    gr_B = GRanges(2, IRanges(c(233101, 333105, 331023, 331028, 329966), c(333229, 333229, 331191, 331191, 330044)), strand = c("+"))
    ## gr_A %&&% gr_B
    ## gr_A %&&% gr_B !=  GRanges(). The former has seqinfo() Seqinfo object with 1 sequence from an unspecified genome; seqnames 2
    expect_equal(length(gr_B %&&% gr_A), 0)
    ## gr_B %&&% gr_A
    expect_equal(length(gr_A %&&% gr_B), 0)

})


## %O% ## strand-agnostic
## Returns a length(a) numeric vector whose item i is the number of bases in a[i] that overlaps at least one range in b.
test_that('%O% works', {

    expect_equal(gr %O% gr, c(1, 1, 1))
    expect_equal(gr2 %O% gr2, c(1, 1))
    ## 
    ## expect_equal(gr %O% gr2, c(1.0000000, 0.3333333, 0.5000000))  
    ## expect_equal(gr2 %O% gr, c(0.5, 0.5))
    gr_A = GRanges(2, IRanges(c(233101, 233105, 231023, 231028, 229966), c(233229, 233229, 231191, 231191, 230044)), strand = c("-"))
    gr_B = GRanges(2, IRanges(c(233101, 333105, 331023, 331028, 329966), c(333229, 333229, 331191, 331191, 330044)), strand = c("+"))
    ## gr_A %O% gr_B
    ## 1 1 0 0 0
    ##  gr_B %O% gr_A
    ## 0.001288338 0.000000000 0.000000000 0.000000000 0.000000000


})


## %OO% ## strand-specific
## Returns a length(a) numeric vector whose item i is the number of bases in a[i] that overlaps at least one range in b.
test_that('%OO% works', {

    expect_equal(gr %OO% gr, c(1, 1, 1))
    expect_equal(gr2 %OO% gr2, c(1, 1))
    ## 
    ##expect_equal(gr %OO% gr2, c(1.0000000, 0.3333333, 0.5000000))  
    ##expect_equal(gr2 %OO% gr, c(0.5, 0.5))
    gr_A = GRanges(2, IRanges(c(233101, 233105, 231023, 231028, 229966), c(233229, 233229, 231191, 231191, 230044)), strand = c("-"))
    gr_B = GRanges(2, IRanges(c(233101, 333105, 331023, 331028, 329966), c(333229, 333229, 331191, 331191, 330044)), strand = c("+"))
    ## gr_A %O% gr_B
    expect_equal(gr_A %OO% gr_B, c(0, 0, 0, 0, 0))  
    ## gr_B %O% gr_A
    expect_equal(gr_B %OO% gr_A, c(0, 0, 0, 0, 0))  

})


## %o%
## strand-agnostic
## Returns a length(a) numeric vector whose item i is the fraction of the width of a[i] that overlaps at least one range in b.
test_that('%o% works', {

    expect_equal(gr %o% gr, c(3, 3, 4))
    expect_equal(gr2 %o% gr2, c(6, 6))
    expect_equal(gr %o% gr2, c(3, 1, 2))
    expect_equal(gr2 %o% gr, c(3, 3))
    gr_A = GRanges(2, IRanges(c(233101, 233105, 231023, 231028, 229966), c(233229, 233229, 231191, 231191, 230044)), strand = c("-"))
    gr_B = GRanges(2, IRanges(c(233101, 333105, 331023, 331028, 329966), c(333229, 333229, 331191, 331191, 330044)), strand = c("+"))
    expect_equal(gr_A %o% gr_B, c(129, 125, 0, 0, 0))  
    expect_equal(gr_B %o% gr_A, c(129, 0, 0, 0, 0))  

})


## %oo%
## strand-specific
## Returns a length(a) numeric vector whose item i is the fraction of the width of a[i] that overlaps at least one range in b.
test_that('%oo% works', {
    
    expect_equal(gr %oo% gr, c(3, 3, 4))
    expect_equal(gr2 %oo% gr2, c(6, 6))
    expect_equal(gr %oo% gr2, c(3, 1, 2))
    expect_equal(gr2 %oo% gr, c(3, 3))
    gr_A = GRanges(2, IRanges(c(233101, 233105, 231023, 231028, 229966), c(233229, 233229, 231191, 231191, 230044)), strand = c("-"))
    gr_B = GRanges(2, IRanges(c(233101, 333105, 331023, 331028, 329966), c(333229, 333229, 331191, 331191, 330044)), strand = c("+"))
    expect_equal(gr_A %oo% gr_B, c(0, 0, 0, 0, 0))  
    expect_equal(gr_B %oo% gr_A, c(0, 0, 0, 0, 0))  

})


## %N%
## strand-agnostic
## Returns a length(a) numeric vector whose item i is the total number of ranges in b that overlap with a[i].
test_that('%N% works', {
    
    expect_equal(gr %N% gr, c(1, 1, 1))
    expect_equal(gr2 %N% gr2, c(1, 1))
    expect_equal(gr %N% gr2, c(1, 1, 1))
    expect_equal(gr2 %N% gr, c(1, 2))
    gr_A = GRanges(2, IRanges(c(233101, 233105, 231023, 231028, 229966), c(233229, 233229, 231191, 231191, 230044)), strand = c("-"))
    gr_B = GRanges(2, IRanges(c(233101, 333105, 331023, 331028, 329966), c(333229, 333229, 331191, 331191, 330044)), strand = c("+"))
    expect_equal(gr_A %N% gr_B, c(1, 1, 0, 0, 0))  
    expect_equal(gr_B %N% gr_A, c(2, 0, 0, 0, 0))  

})


## %NN%
## strand-specific
## Returns a length(a) numeric vector whose item i is the total number of ranges in b that overlap with a[i].
test_that('%NN% works', {
    
    expect_equal(gr %NN% gr, c(1, 1, 1))
    expect_equal(gr2 %NN% gr2, c(1, 1))
    expect_equal(gr %NN% gr2, c(1, 1, 1))
    expect_equal(gr2 %NN% gr, c(1, 2))
    gr_A = GRanges(2, IRanges(c(233101, 233105, 231023, 231028, 229966), c(233229, 233229, 231191, 231191, 230044)), strand = c("-"))
    gr_B = GRanges(2, IRanges(c(233101, 333105, 331023, 331028, 329966), c(333229, 333229, 331191, 331191, 330044)), strand = c("+"))
    expect_equal(gr_A %NN% gr_B, c(0, 0, 0, 0, 0))  
    expect_equal(gr_B %NN% gr_A, c(0, 0, 0, 0, 0))  

})


## %_%
## Shortcut for GenomicRanges::setdiff
test_that('%_% works', {

    gr1 = GRanges(1, IRanges(10,20), strand="+")
    gr2 = GRanges(1, IRanges(15,25), strand="-")
    gr3 = GRanges("1:1-15")
    expect_equal(width(gr1 %_% gr2), 5)
    expect_equal(start(gr1 %_% gr2), 10)
    expect_equal(end(gr1 %_% gr2), 14)
    expect_equal(width(gr1 %_% gr3), 5)
    expect_equal(start(gr1 %_% gr3), 16)
    expect_equal(end(gr1 %_% gr3), 20)
    expect_equal(width(gr3 %_% gr1), 9)
    expect_equal(start(gr3 %_% gr1), 1)
    expect_equal(end(gr3 %_% gr1), 9)
    expect_equal(width(gr2 %_% gr1), 5)
    expect_equal(start(gr2 %_% gr1), 21)
    expect_equal(end(gr2 %_% gr1), 25)

})


## %Q%
## Subsets or re-orders a based on a logical or integer valued expression that operates on the GRanges metadata columns of a.
test_that('%Q% works', {

    testset = GRanges(seqnames = Rle(1,c(5)) , ranges = IRanges(1:5 , end = 2000:2004) , strand = Rle(strand(c("+")) , c(5)) , mean = c(1, -2 , -5 , 5 ,6))
    expect_equal(length(testset %Q% (mean > 0)) , 3)
    foo = GRanges(seqnames = Rle(1,c(5)) , ranges = IRanges(1:5 , end = 2000:2004) , strand = Rle(strand(c("+")) , c(5)) , mean = c(50, 200, 300, 400, 500))
    expect_equal(length(foo %Q% (mean == 300)) , 1)
    testgr = GRanges(seqnames = Rle(1,c(5)) , ranges = IRanges(1:5 , end = c(5, 10, 15, 20, 25)) , strand = Rle(strand(c("+", "-", "-", "-", "-"))) , mean = c(1, -2 , -5 , 5 ,6), type = c("exon", "CDS", "exon", "CDS", "exon"))
    ## string matching strand
    expect_equal(length(testgr %Q% (strand == "-")), 4)
    expect_equal(width(testgr %Q% (strand == '+')), 5)
    expect_equal(length(testgr %Q% (strand == '*')), 0)
    expect_error(testgr %Q% (strand > 100))
    ## strng matching type
    expect_equal(length(testgr %Q% ((type != 'exon') | (type != 'CDS'))), 5)
    expect_equal(length(testgr %Q% ((type != 'exon') & (type != 'CDS'))), 0)
    expect_equal(width(testgr %Q% ((mean > 1) & (type == 'CDS'))), 17)  ## only GRanges 1 [4, 20] - | 5 CDS
    expect_equal(length(testgr %Q% (mean > 'foo')), 0)  ## nonsense %Q% subsets don't throw error, but do return empty GRanges

})


## %^%
## strand-agnostic
## Returns a length(a) logical vector whose item i TRUE if the a[i] overlaps at least on range in b (similar to %over% just less fussy about Seqinfo).
test_that('%^% works', {

    testset = GRanges(seqnames = Rle(1,c(5)) , ranges = IRanges(1:5 , end = 2000:2004) , strand = Rle(strand(c("+")) , c(5)) , mean = c(1, -2 , -5 , 5 ,6))
    expect_equal(length(testset %^% testset), 5)
    expect_equal(gr %^% testset, c(TRUE, TRUE, TRUE))
    expect_equal(testset %^% gr, c(TRUE, TRUE, TRUE, TRUE, TRUE))
    gr_A = GRanges(2, IRanges(c(233101, 233105, 231023, 231028, 229966), c(233229, 233229, 231191, 231191, 230044)), strand = c("-"))
    gr_B = GRanges(2, IRanges(c(233101, 333105, 331023, 331028, 329966), c(333229, 333229, 331191, 331191, 330044)), strand = c("+"))
    expect_equal(as.vector(gr_A %^% gr_B), c(TRUE, TRUE, FALSE, FALSE, FALSE))
    expect_equal(as.vector(gr_B %^% gr_A), c(TRUE, FALSE, FALSE, FALSE, FALSE))

})


## %^^%
## strand-specific
test_that('%^^% works', {

    testset = GRanges(seqnames = Rle(1,c(5)) , ranges = IRanges(1:5 , end = 2000:2004) , strand = Rle(strand(c("+")) , c(5)) , mean = c(1, -2 , -5 , 5 ,6))
    expect_equal(length(testset %^^% testset), 5)
    expect_equal(gr %^^% testset, c(TRUE, FALSE, FALSE))
    expect_equal(testset %^^% gr, c(TRUE, TRUE, TRUE, TRUE, TRUE))
    gr_A = GRanges(2, IRanges(c(233101, 233105, 231023, 231028, 229966), c(233229, 233229, 231191, 231191, 230044)), strand = c("-"))
    gr_B = GRanges(2, IRanges(c(233101, 333105, 331023, 331028, 329966), c(333229, 333229, 331191, 331191, 330044)), strand = c("+"))
    expect_equal(as.vector(gr_A %^^% gr_B), c(FALSE, FALSE, FALSE, FALSE, FALSE))
    expect_equal(as.vector(gr_B %^^% gr_A), c(FALSE, FALSE, FALSE, FALSE, FALSE))

})


## %$%  ## strand-specific
## Aggregates the metadata in b across the territory of each range in a. 
test_that('%$% works', {

    gr1 = GRanges(1, IRanges(c(3,7,13), c(5,9,16)), strand=c('+','-','-'), seqinfo=Seqinfo("1", 25), name=c('A', 'B', 'C'))
    gr1$gene_name = c('BRCA1', 'BRCA2', 'KRAS')
    gr2 =  GRanges(1, IRanges(c(3,7,13), c(50, 90, 160)), strand='+')
    gr2$aux = c(42.42, 'Arp2/3', NA)
    expect_equal(colnames(gr2dt(gr1 %$% gr2))[6:8], c('name', 'gene_name', 'aux'))
    expect_equal(colnames(gr2dt(gr2 %$% gr1))[6:8], c('aux', 'name', 'gene_name'))
    expect_equal((gr1 %$% gr2)$name, c('A', 'B', 'C'))
    expect_equal((gr1 %$% gr2)$gene_name, c('BRCA1', 'BRCA2', 'KRAS'))
    expect_equal((gr1 %$% gr2)$aux, c('42.42', '42.42, Arp2/3', '42.42, Arp2/3'))
    expect_equal((gr2 %$% gr1)$aux, c('42.42', 'Arp2/3', NA))
    expect_equal((gr2 %$% gr1)$name, c('A, B, C', 'B, C', 'C'))
    expect_equal((gr2 %$% gr1)$gene_name, c('BRCA1, BRCA2, KRAS', 'BRCA2, KRAS', 'KRAS'))

})


## %$$%  ## strand-agnostic
## Aggregates the metadata in b across the territory of each range in a. 
test_that('%$$% works', {

    gr1 = GRanges(1, IRanges(c(3,7,13), c(5,9,16)), strand=c('+','-','-'), seqinfo=Seqinfo("1", 25), name=c('A', 'B', 'C'))
    gr1$gene_name = c('BRCA1', 'BRCA2', 'KRAS')
    gr2 =  GRanges(1, IRanges(c(3,7,13), c(50, 90, 160)), strand='+')
    gr2$aux = c(42.42, 'Arp2/3', NA)
    expect_equal(colnames(gr2dt(gr1 %$$% gr2))[6:8], c('name', 'gene_name', 'aux'))
    expect_equal(colnames(gr2dt(gr2 %$$% gr1))[6:8], c('aux', 'name', 'gene_name'))
    expect_equal((gr1 %$$% gr2)$name, c('A', 'B', 'C'))
    expect_equal((gr1 %$$% gr2)$gene_name, c('BRCA1', 'BRCA2', 'KRAS'))
    expect_equal((gr1 %$$% gr2)$aux, c('42.42', '', ''))
    expect_equal((gr2 %$$% gr1)$aux, c('42.42', 'Arp2/3', NA))
    expect_equal((gr2 %$$% gr1)$name, c('A', '', ''))
    expect_equal((gr2 %$$% gr1)$gene_name, c('BRCA1', '', ''))

})


## %*%
## strand-specific
test_that('%*% works', {

    gr1 = GRanges(1, IRanges(c(3,7,13), c(5,9,16)), strand=c('+','-','-'), seqinfo=Seqinfo("1", 25), name=c('A', 'B', 'C'))
    gr1$gene_name = c('BRCA1', 'BRCA2', 'KRAS')
    gr2 =  GRanges(1, IRanges(c(3,7,13), c(50, 90, 160)), strand='+')
    gr2$aux = c(42.42, 'Arp2/3', NA)
    expect_equal(colnames(gr2dt(gr2 %*% gr1))[8:10], c('aux', 'name', 'gene_name'))
    expect_equal(colnames(gr2dt(gr1 %*% gr2))[8:10], c('name', 'gene_name', 'aux'))
    expect_equal((gr1 %*% gr2)$query.id, c(1, 2, 2, 3, 3, 3))
    expect_equal((gr2 %*% gr1)$query.id, c(1, 1, 2, 1, 2, 3))

})


## %**%
## strand-agnostic
test_that('%**% works', {

    gr1 = GRanges(1, IRanges(c(3,7,13), c(5,9,16)), strand=c('+','-','-'), seqinfo=Seqinfo("1", 25), name=c('A', 'B', 'C'))
    gr1$gene_name = c('BRCA1', 'BRCA2', 'KRAS')
    gr2 =  GRanges(1, IRanges(c(3,7,13), c(50, 90, 160)), strand='+')
    gr2$aux = c(42.42, 'Arp2/3', NA)
    expect_equal(colnames(gr2dt(gr2 %**% gr1))[8:10], c('aux', 'name', 'gene_name'))
    expect_equal(colnames(gr2dt(gr1 %**% gr2))[8:10], c('name', 'gene_name', 'aux'))
    expect_equal(length(gr1 %**% gr2), 1)
    expect_equal(length(gr2 %**% gr1), 1)

})


## gr.setdiff
## gr.setdiff = function(query, subject, ignore.strand = TRUE, by = NULL)
test_that('gr.setdiff', {

    gr1 = GRanges(1, IRanges(c(3,7,13), c(5,9,16)), strand=c('+','-','-'), seqinfo=Seqinfo("1", 25), name=c('A', 'B', 'C'))
    gr1$gene_name = c('BRCA1', 'BRCA2', 'KRAS')
    gr2 =  GRanges(1, IRanges(c(3,7,13), c(50, 90, 160)), strand='+')
    gr2$aux = c(42.42, 'Arp2/3', NA)
    expect_equal(length(gr.setdiff(gr1, gr2)), 0)
    expect_equal(width(gr.setdiff(gr2, gr1)), c(1, 3, 3, 9, 9, 9))  
    expect_equal(gr.setdiff(gr2, gr1)$query.id, c(1, 1, 2, 1, 2, 3))
    expect_equal(gr.setdiff(gr2, gr1)$subject.id, c(2, 3, 3, 4, 4, 4))  
    expect_equal(gr.setdiff(gr2, gr1)$aux, c('42.42', '42.42', 'Arp2/3', '42.42', 'Arp2/3', NA))  
    expect_equal(length(gr1 %**% gr2), 1)
    expect_equal(length(gr2 %**% gr1), 1)
    expect_equal(width(gr.setdiff(gr2, gr1, ignore.strand = FALSE)), c(23, 20, 19, 19, 13, 13))
    expect_equal(gr.setdiff(gr2, gr1, ignore.strand = FALSE)$query.id, c(1, 1, 2, 2, 3, 3))
    expect_equal(gr.setdiff(gr2, gr1, ignore.strand = FALSE)$subject.id, c(6, 2, 2, 6, 2, 6))
    expect_equal(gr.setdiff(gr2, gr1, ignore.strand = FALSE)$aux, c('42.42', '42.42', 'Arp2/3', 'Arp2/3', NA, NA))
   
})

#### 
#### 
#### test_that('ra.merge', {
#### 
####     ## beginning with fake rearrangment data grl1 and grl2
####     expect_equal(length(ra.merge(grl1, grl2)), 500)
####     expect_true(unique(values(ra.merge(grl1, grl2))[, 1][1:250]))
####     expect_false(unique(values(ra.merge(grl1, grl2))[, 1][251:500]))
####     expect_false(unique(values(ra.merge(grl1, grl2))[, 2][1:249]))
####     expect_true(unique(values(ra.merge(grl1, grl2))[, 2][250:500]))
####     ## ignore.strand == TRUE makes no difference...
####     expect_equal(ra.merge(grl1, grl2, ignore.strand=TRUE), ra.merge(grl1, grl2))
####     ## example in function 'ra.merge()'
####     gr1 = GRanges(1, IRanges(1:10, width = 1), strand = rep(c('+', '-'), 5))
####     gr2 = GRanges(1, IRanges(4 + 1:10, width = 1), strand = rep(c('+', '-'), 5))
####     expect_error(ra.merge(gr1, gr2)) ##   Error: All inputs must be a GRangesList
####     ## create GRangesLists
####     ra1 = split(gr1, rep(1:5, each = 2))
####     ra2 = split(gr2, rep(1:5, each = 2))
####     ram = ra.merge(ra1, ra2)
####     expect_warning(ra.merge(ra1, ra2))  ## warning: GRanges object contains 10 out-of-bound ranges located on sequence 1.
####     expect_equal(length(ram), 7)
####     ## 'values(ram)' shows the metadata with TRUE / FALSE flags
####     expect_equal(values(ram)[, 1], c(TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE)) 
####     expect_equal(values(ram)[, 2], c(FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE)) 
####     ## ram2 = ra.merge(ra1, ra2, pad = 5) # more inexact matching results in more merging
####     ram2 = ra.merge(ra1, ra2, pad = 500) ## more inexact matching results in more merging
####     ##values(ram2)
####     expect_equal(length(ram2), 5)    
####     expect_equal(values(ram2)[, 1], c(TRUE, TRUE, TRUE, TRUE, TRUE)) 
####     expect_equal(values(ram2)[, 2], c(TRUE, TRUE, TRUE, TRUE, TRUE)) 
####     expect_error(ra.merge(ra1, ra2, pad =-1)) ## adjustment would result in ranges with negative widths
####     ## ram3
####     ram3 = ra.merge(ra1, ra2, ind = TRUE) ## indices instead of flags
####     ##values(ram3)
####     expect_equal(length(ram3), 7)
####     expect_equal(values(ram3)[, 1], c(1, 2, 3, 4, 5, NA, NA)) 
####     expect_equal(values(ram3)[, 2], c(NA, NA, 3, 4, 5, 4, 5)) 
####     ## test both 'pad', 'ind'
####     ram4 = ra.merge(ra1, ra2, pad = 500, ind = TRUE) 
####     expect_equal(values(ram4)[, 1], c(1, 2, 3, 4, 5))
####     expect_equal(values(ram4)[, 2], c(1, 2, 3, 4, 5))
####     ## ignore.strand == TRUE
####     expect_error(ra.merge(ra1, ra2, ignore.stand = TRUE)) ##  unable to find an inherited method for function ‘values’ for signature ‘"logical"’
####     ### all args
####     expect_error(ra.merge(ra1, ra2, pad = 500, ind = TRUE, ignore.stand = TRUE)) ## unable to find an inherited method for function ‘values’ for signature ‘"logical"’
#### 
#### })
#### 



test_that("gr.simplify", {

    gg = gr.simplify(gr, pad=4)
    expect_identical(end(gg), c(5L, 16L))
    expect_equal(length(gg), 2)

    gr$field = c("A","B","B")
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




test_that('anchorlift', {

    ### check returns NULL
    expect_equal(anchorlift(GRanges(), GRanges()), NULL)
    ## check 'if (length(ov) == 0){ return(NULL) }'
    expect_equal(anchorlift(GRanges('2:2000-3000'), GRanges('1:10-100'), window=100), NULL)
    ## unlist rearrangement datasets grl1 and grl2
    sv1 = grl.unlist(grl1)
    sv2 = grl.unlist(grl2)
    ## default
    expect_equal(length(suppressWarnings(anchorlift(sv1, sv2))), 13158)
    ## 'by' argument
    anchor1 = anchorlift(sv1, sv2, by='bin')
    expect_equal(width(anchor1), c(1, 1))
    expect_equal(anchor1[1]$subject.id, 1)
    expect_equal(anchor1[2]$subject.id, 2)
    expect_equal(anchor1[1]$query.id, 499)
    expect_equal(anchor1[2]$query.id, 500)   
    expect_equal(anchor1[1]$bin, 4678)   
    expect_equal(anchor1[2]$bin, 4678)   
    ## 'window' argument
    ## error if over 1e9
    expect_error(anchorlift(gr, gr2, window=1.1e9))
    ## include.values
    expect_equal(dim(gr2dt(suppressWarnings(anchorlift(sv1, sv2, by='bin', include.values = FALSE))))[2], 7)  ## check only 7 columns

})


## tests for data functions



## hg_seqlengths()
test_that('hg_seqlengths()', {

    expect_equal(length(hg_seqlengths()), 25)
    expect_equal(as.vector(hg_seqlengths()[1]), 249250621)
    expect_equal(as.vector(hg_seqlengths()[2]), 243199373)
    expect_equal(as.vector(hg_seqlengths()[3]), 198022430)
    expect_equal(as.vector(hg_seqlengths()[4]), 191154276)
    expect_equal(as.vector(hg_seqlengths()[5]), 180915260)
    expect_equal(as.vector(hg_seqlengths()[6]), 171115067)
    expect_equal(as.vector(hg_seqlengths()[7]), 159138663)
    expect_equal(as.vector(hg_seqlengths()[8]), 146364022)
    expect_equal(as.vector(hg_seqlengths()[9]), 141213431)
    expect_equal(as.vector(hg_seqlengths()[10]), 135534747)
    expect_equal(as.vector(hg_seqlengths()[11]), 135006516)
    expect_equal(as.vector(hg_seqlengths()[12]), 133851895)
    expect_equal(as.vector(hg_seqlengths()[13]), 115169878)
    expect_equal(as.vector(hg_seqlengths()[14]), 107349540)
    expect_equal(as.vector(hg_seqlengths()[15]), 102531392)
    expect_equal(as.vector(hg_seqlengths()[16]), 90354753)
    expect_equal(as.vector(hg_seqlengths()[17]), 81195210)
    expect_equal(as.vector(hg_seqlengths()[18]), 78077248)
    expect_equal(as.vector(hg_seqlengths()[19]), 59128983)
    expect_equal(as.vector(hg_seqlengths()[20]), 63025520)
    expect_equal(as.vector(hg_seqlengths()[21]), 48129895)
    expect_equal(as.vector(hg_seqlengths()[22]), 51304566)
    expect_equal(as.vector(hg_seqlengths()[23]), 155270560)
    expect_equal(as.vector(hg_seqlengths()[24]), 59373566)
    expect_equal(as.vector(hg_seqlengths()[25]), 16571)

})



## example_genes
test_that('example_genes', {

    expect_equal(length(example_genes), 18812)
    expect_true(is(example_genes, 'GRanges'))
    expect_equal(length(unique(example_genes$exonCount)), 102)

})



## example_dnase
test_that('example_dnase', {

    expect_equal(length(example_dnase), 10000)
    expect_true(is(example_dnase, 'GRanges'))
    expect_equal(max(example_dnase$signalValue), 714.731)
    expect_equal(max(example_dnase$pValue), 324)

})



## grl1
test_that('grl1, SV rearrangements', {

    expect_equal(length(grl1), 250)
    expect_true(is(grl1, 'GRangesList'))
    expect_equal(length(unlist(grl1)), 500)
    expect_equal(grl1[[1]]$bin, c(10, 10))
    expect_equal(width(grl1[[1]]), c(1, 1))

})



## grl2
test_that('grl2, SV rearrangements', {

    expect_equal(length(grl2), 251)
    expect_true(is(grl2, 'GRangesList'))
    expect_equal(length(unlist(grl2)), 502)
    expect_equal(grl2[[1]]$bin, c(4678, 4678))
    expect_equal(width(grl2[[1]]), c(1, 1))

})


## 
test_that('si', {
    
    expect_equal(length(si), 25)
    expect_false(all(as.data.frame(si)$isCircular))
    expect_match(unique(as.data.frame(si)$genome), 'hg19')

})

test_that('grl.hiC', {

    expect_equal(length(grl.hiC), 10000)
    expect_equal(length(unlist(grl.hiC)), 20000)


})










## XT Yao function
#### gr.breaks


## XT, ra.dedup
## Jan 18, correspondence from XT 'leave that internal for now'


## XT, ra.duplicated
## Jan 18, correspondence from XT 'leave that internal for now'



## XT plans to re-write, Jan 18
test_that('ra.overlaps', {

    gr = GRanges(1, IRanges(c(3,7), c(5,9)), strand=c('+','-'))
    gr1 = GRanges(1, IRanges(c(10,20), width=5), strand=c("+", "-"))
    gr2 = GRanges(1, IRanges(c(1,9), c(6,14)), strand=c('+','-'))
    grl1 = GRangesList("gr"=gr, "gr1"=gr1)
    grl2 = GRangesList("gr1"=gr1, "gr2"=gr2)
    foobar = suppressWarnings(grl.bind(grl1, grl2))
    ro = ra.overlaps(grl1, grl2)
    expect_equal(class(ro), "matrix")
    expect_equal(nrow(ro), 2)
    expect_equal(ncol(ro), 2)
    expect_equal(nrow(ra.overlaps(grl2, grl2)), length(grl2))

 })

test_that("ra.overlaps handles empty",{
     
    ## test empty inputs and no overlaps inputs
    gr = GRanges(1, IRanges(c(10,20), width=5), strand=c("+", "-"))
    grl1 = GRangesList("gr1" = gr)
    expect_equal(ra.overlaps(GRangesList(), grl1)[1], NA)
    expect_equal(ra.overlaps(grl2[2:3], grl1)[1], NA)
    
})  

test_that("ra.overlaps handles wrong signs", {
 
    ## make one that overlaps, but wrong signs
    gr = GRanges(1, IRanges(c(3,7), c(5,9)), strand=c('+','-'))
    gr1 = GRanges(1, IRanges(c(10,20), width=5), strand=c("+", "-"))
    gr2 = GRanges(1, IRanges(c(1,9), c(6,14)), strand=c('+','-'))
    grl1 = GRangesList("gr" = gr, "gr1" = gr1, "gr2" = gr2)
    grl3 <- grl1[2]
    strand(grl3[[1]]) <- c("+", "-")
    expect_equal(ra.overlaps(grl3, grl2)[1], NA)
 
})







