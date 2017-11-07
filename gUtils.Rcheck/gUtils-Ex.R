pkgname <- "gUtils"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('gUtils')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
cleanEx()
nameEx("dt2gr")
### * dt2gr

flush(stderr()); flush(stdout())

### Name: dt2gr
### Title: Convert data.table to GRanges
### Aliases: dt2gr

### ** Examples

gr <- dt2gr(data.table(start=c(1,2), seqnames=c("X", "1"), end=c(10,20), strand = c('+', '-')))



cleanEx()
nameEx("gr.chr")
### * gr.chr

flush(stderr()); flush(stdout())

### Name: gr.chr
### Title: Prepend "chr" to 'GRanges seqlevels'
### Aliases: gr.chr

### ** Examples

gr <-  gr.chr(GRanges(c(1,"chrX"), IRanges(c(1,2), 1)))
seqnames(gr)



cleanEx()
nameEx("gr.dice")
### * gr.dice

flush(stderr()); flush(stdout())

### Name: gr.dice
### Title: Dice up 'GRanges' into 'width = 1' 'GRanges' spanning the input
###   (warning can produce a very large object)
### Aliases: gr.dice

### ** Examples

gr.dice(GRanges(c(1,4), IRanges(c(10,10),20)))



cleanEx()
nameEx("gr.duplicated")
### * gr.duplicated

flush(stderr()); flush(stdout())

### Name: gr.duplicated
### Title: Allows to restrict duplicates using "by" columns and allows in
###   exact matching
### Aliases: gr.duplicated

### ** Examples

gr.duplicated(GRanges(c(1,1,1), IRanges(c(2,5,5), width=1)))

gr.duplicated(GRanges(c(1,1,1), IRanges(c(2,5,5), width=1)))




cleanEx()
nameEx("gr.end")
### * gr.end

flush(stderr()); flush(stdout())

### Name: gr.end
### Title: Get the right ends of a 'GRanges'
### Aliases: gr.end

### ** Examples

gr.end(example_dnase, width=200, clip=TRUE)



cleanEx()
nameEx("gr.flipstrand")
### * gr.flipstrand

flush(stderr()); flush(stdout())

### Name: gr.flipstrand
### Title: Flip strand on 'GRanges'
### Aliases: gr.flipstrand

### ** Examples

gr.flipstrand(GRanges(1, IRanges(c(10,10,10),20), strand=c("+","*","-")))



cleanEx()
nameEx("gr.mid")
### * gr.mid

flush(stderr()); flush(stdout())

### Name: gr.mid
### Title: Get the midpoints of 'GRanges' ranges
### Aliases: gr.mid

### ** Examples

gr.mid(GRanges(1, IRanges(1000,2000), seqinfo=Seqinfo("1", 2000)))



cleanEx()
nameEx("gr.rand")
### * gr.rand

flush(stderr()); flush(stdout())

### Name: gr.rand
### Title: Generate random 'GRanges' on genome
### Aliases: gr.rand

### ** Examples

## Generate 5 non-overlapping regions of width 10 on hg19
gr.rand(rep(10,5), BSgenome.Hsapiens.UCSC.hg19::Hsapiens)



cleanEx()
nameEx("gr.start")
### * gr.start

flush(stderr()); flush(stdout())

### Name: gr.start
### Title: Get GRanges corresponding to beginning of range
### Aliases: gr.start

### ** Examples

gr.start(example_dnase, width=200)
gr.start(example_dnase, width=200, clip=TRUE)



cleanEx()
nameEx("gr.string")
### * gr.string

flush(stderr()); flush(stdout())

### Name: gr.string
### Title: Return UCSC style interval string corresponding to 'GRanges'
###   pile (ie chr:start-end)
### Aliases: gr.string

### ** Examples

gr.string(example_genes, other.cols = c("name", "name2"))



cleanEx()
nameEx("gr.tile")
### * gr.tile

flush(stderr()); flush(stdout())

### Name: gr.tile
### Title: Tile ranges across 'GRanges'
### Aliases: gr.tile

### ** Examples

## 10 tiles of width 10
gr1 <- gr.tile(GRanges(1, IRanges(1,100)), w=10)
## make them overlap each other by 5
gr1 + 5



cleanEx()
nameEx("gr.trim")
### * gr.trim

flush(stderr()); flush(stdout())

### Name: gr.trim
### Title: Trims pile of 'GRanges' relative to the specified <local>
###   coordinates of each range
### Aliases: gr.trim

### ** Examples

## trim the first 20 and last 50 bases
gr.trim(GRanges(1, IRanges(1e6, width=1000)), starts=20, ends=950)
## return value: GRanges on 1:1,000,019-1,000,949



cleanEx()
nameEx("grfo")
### * grfo

flush(stderr()); flush(stdout())

### Name: %*%
### Title: Metadata join with coordinates as keys (wrapper to
###   'gr.findoverlaps')
### Aliases: %*% %*%,GRanges-method

### ** Examples

example_genes %*% example_dnase



cleanEx()
nameEx("grl.pivot")
### * grl.pivot

flush(stderr()); flush(stdout())

### Name: grl.pivot
### Title: Pivot a 'GRangesList', inverting "x" and "y"
### Aliases: grl.pivot

### ** Examples

grl.pivot(grl.hiC)



cleanEx()
nameEx("grl.reduce")
### * grl.reduce

flush(stderr()); flush(stdout())

### Name: grl.reduce
### Title: grl.reduce
### Aliases: grl.reduce

### ** Examples


grl.reduce(grl, 1000)

unlist(grl.reduce(split(reads+10000, reads$BX)))




cleanEx()
nameEx("grl.string")
### * grl.string

flush(stderr()); flush(stdout())

### Name: grl.string
### Title: Create string representation of 'GRangesList'
### Aliases: grl.string

### ** Examples

grl.string(grl.hiC, mb=TRUE)



cleanEx()
nameEx("grl.unlist")
### * grl.unlist

flush(stderr()); flush(stdout())

### Name: grl.unlist
### Title: Robust unlisting of 'GRangesList' that keeps track of origin
### Aliases: grl.unlist

### ** Examples

grl.unlist(grl.hiC)



cleanEx()
nameEx("grlbind")
### * grlbind

flush(stderr()); flush(stdout())

### Name: grlbind
### Title: Concatenate 'GRangesList' objects.
### Aliases: grlbind

### ** Examples

## Concatenate
grl.hiC2 <- grl.hiC[1:20]
mcols(grl.hiC2)$test = 1
grlbind(grl.hiC2, grl.hiC[1:30])



cleanEx()
nameEx("ra.merge")
### * ra.merge

flush(stderr()); flush(stdout())

### Name: ra.merge
### Title: Merges rearrangements represented by 'GRangesList' objects
### Aliases: ra.merge

### ** Examples


# generate some junctions
gr1 <- GRanges(1, IRanges(1:10, width = 1), strand = rep(c('+', '-'), 5))
gr2 <- GRanges(1, IRanges(4 + 1:10, width = 1), strand = rep(c('+', '-'), 5))
ra1 = split(gr1, rep(1:5, each = 2))
ra2 = split(gr2, rep(1:5, each = 2))

ram = ra.merge(ra1, ra2)
values(ram) # shows the metadata with TRUE / FALSE flags

ram2 = ra.merge(ra1, ra2, pad = 5) # more inexact matching results in more merging
values(ram2)

ram3 = ra.merge(ra1, ra2, ind = TRUE) #indices instead of flags
values(ram3)



cleanEx()
nameEx("si2gr")
### * si2gr

flush(stderr()); flush(stdout())

### Name: si2gr
### Title: Create 'GRanges' from 'Seqinfo' or 'BSgenome'
### Aliases: si2gr

### ** Examples

si2gr(BSgenome.Hsapiens.UCSC.hg19::Hsapiens)



cleanEx()
nameEx("streduce")
### * streduce

flush(stderr()); flush(stdout())

### Name: streduce
### Title: Reduce 'GRanges' and 'GRangesList' to miminal footprint
### Aliases: streduce

### ** Examples

streduce(grl.hiC, pad=10)
streduce(example_genes, pad=1000)



### * <FOOTER>
###
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
