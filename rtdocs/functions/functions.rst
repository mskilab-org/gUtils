Function Use
============

si2gr
~~~~~

## Convert a numeric vector to a GRanges object.


.. sourcecode:: r
    

    ## created a named numeric vector. 
    fake.genome = c('1'=1e4, '2'=1e3, '3'=5e3)
    si2gr(fake.genome)


::

    ## GRanges object with 3 ranges and 0 metadata columns:
    ##     seqnames     ranges strand
    ##        <Rle>  <IRanges>  <Rle>
    ##   1        1 [1, 10000]      +
    ##   2        2 [1,  1000]      +
    ##   3        3 [1,  5000]      +
    ##   -------
    ##   seqinfo: 3 sequences from an unspecified genome



gr.fix
~~~~~~


.. sourcecode:: r
    

    ## fix a GRanges object that has incorrectly formated 'seqnames'
    ## first create a GRanges object that uses character vectors in the 'seqnames' field.
     gr <- GRanges(seqnames = Rle(c("chr1" , "chr2" , "chr1" , "chr3"),c(1,3,2,4)), ranges = IRanges(c(1,3,5,7,9,11,13,15,17,19) ,end = c(2,4,6,8,10,12,14,16,18,20),names = head(letters,10)),GC=seq(1,10,length=10),name=seq(5,10,length=10))
     gr


::

    ## GRanges object with 10 ranges and 2 metadata columns:
    ##     seqnames    ranges strand |        GC             name
    ##        <Rle> <IRanges>  <Rle> | <numeric>        <numeric>
    ##   a     chr1  [ 1,  2]      * |         1                5
    ##   b     chr2  [ 3,  4]      * |         2 5.55555555555556
    ##   c     chr2  [ 5,  6]      * |         3 6.11111111111111
    ##   d     chr2  [ 7,  8]      * |         4 6.66666666666667
    ##   e     chr1  [ 9, 10]      * |         5 7.22222222222222
    ##   f     chr1  [11, 12]      * |         6 7.77777777777778
    ##   g     chr3  [13, 14]      * |         7 8.33333333333333
    ##   h     chr3  [15, 16]      * |         8 8.88888888888889
    ##   i     chr3  [17, 18]      * |         9 9.44444444444444
    ##   j     chr3  [19, 20]      * |        10               10
    ##   -------
    ##   seqinfo: 3 sequences from an unspecified genome; no seqlengths


.. sourcecode:: r
    

     gr.fix(gr)


::

    ## GRanges object with 10 ranges and 2 metadata columns:
    ##     seqnames    ranges strand |        GC             name
    ##        <Rle> <IRanges>  <Rle> | <numeric>        <numeric>
    ##   a     chr1  [ 1,  2]      * |         1                5
    ##   b     chr2  [ 3,  4]      * |         2 5.55555555555556
    ##   c     chr2  [ 5,  6]      * |         3 6.11111111111111
    ##   d     chr2  [ 7,  8]      * |         4 6.66666666666667
    ##   e     chr1  [ 9, 10]      * |         5 7.22222222222222
    ##   f     chr1  [11, 12]      * |         6 7.77777777777778
    ##   g     chr3  [13, 14]      * |         7 8.33333333333333
    ##   h     chr3  [15, 16]      * |         8 8.88888888888889
    ##   i     chr3  [17, 18]      * |         9 9.44444444444444
    ##   j     chr3  [19, 20]      * |        10               10
    ##   -------
    ##   seqinfo: 3 sequences from an unspecified genome


