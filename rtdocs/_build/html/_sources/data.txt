Data manipulation
-----------------

R provides a number of data structures for storing genomic data, each with its advantages and drawbacks. 

The most useful structures for this purpose are:

GRanges
  Store ranges along with metadata, sequences and the coordaintes of the reference genome.

GRangesList
  Store groups of ranges, with additional metadata belonging to the group.

data.table
  Fast and efficient general-purpose container similar to data.frame, but with significant performance improvements.


In gUtils functions, we often manipulate the data to move between these data structures where one is more useful than another. A 
key example is in ``gr.findoverlaps``, which converts the input ``GRanges`` into ``data.table`` objects to take advantage of the
blazing fast ``foverlaps`` util. For the most part, these conversions should be invisible to the user. 

However, often there are data structures conversions that may be useful to the end user. This includes unlisting GRangesList objects
into GRanges, making data.table objects from GRanges, and binding together multiple GRanges or GRangesList objects, among others. This
section will describe and demonstrate the functionality ``gUtils`` provides for manipulating these data structures.

.. code-block:: bash

   ref19 <- readRDS(system.file("extdata","refGene.hg19.gr.rds", package="gUtils"))
   gr <- GRanges(1, IRanges(c(2,5,10), c(4,9,16)), seqinfo=Seqinfo("1", 20)) 
   dt <- data.table(seqnames=1, start=c(2,5,10), end=c(3,8,15))

grbind

grlbind

dtgr

grdt

si2gr

gr2gatk

gr.flatten

gr.flatmap

grl.split

grl.stripnames

grl.unlist

grl.span

grl.pivot

rrbind
