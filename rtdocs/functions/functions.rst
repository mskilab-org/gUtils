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



