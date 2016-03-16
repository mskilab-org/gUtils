Range Operations
----------------

This section will describe additional GRanges operations provided by gUtils.

.. code-block:: bash

   ## make some example data sets
   ref19 <- readRDS(system.file("extdata","refGene.hg19.gr.rds", package="gUtils"))
   gr  <- GRanges(1, IRanges(c(2,5,10), c(4,9,16)), seqinfo=Seqinfo("1", 20))
   gr2 <- c(gr, GRanges(1, IRanges(c(1,9), c(6,14)), seqinfo=Seqinfo("1", 20)))
   dt <- data.table(seqnames=1, start=c(2,5,10), end=c(3,8,15))

.. figure:: figures/shift.png
   :alt:
   :scale: 125 %

.. figure:: figures/flank.png
   :alt:
   :scale: 125 %

.. figure:: figures/flank_start.png
   :alt:
   :scale: 125 %

.. figure:: figures/gr.start.png
   :alt:
   :scale: 125 %

.. figure:: figures/gr.start_wstrand.png
   :alt:
   :scale: 125 %

.. figure:: figures/gr.end.png
   :alt:
   :scale: 125 %

.. figure:: figures/gr.end_wstrand.png
   :alt:
   :scale: 125 %

.. figure:: figures/gr.mid.png
   :alt:
   :scale: 125 %

.. figure:: figures/gr.flatten.png
   :alt:
   :scale: 125 %

.. figure:: figures/gr.flipstrand.png
   :alt:
   :scale: 125 %

.. figure:: figures/gr.tile.png
   :alt:
   :scale: 125 %

.. figure:: figures/streduce.png
   :alt:
   :scale: 125 %

``grbind``

.. code-block:: bash

   ## add metadata to one field
   mcols(gr)$score = 3
   ## try to concatenate
   c(gr,gr2)  ## ERROR
   ## with grbind
   grbind(gr, gr2) ## SUCCESS. Adds NA for missing fields
   ## GenomicRanges::c does this already for GRangesList

``gr.sample(gr2, 2, len=2, replace=TRUE)``

.. code-block:: bash

   ## output GRanges
   GRanges object with 3 ranges and 1 metadata column:
      seqnames    ranges strand |  query.id
         <Rle> <IRanges>  <Rle> | <integer>
   [1]        1  [ 8,  9]      * |         2
   [2]        1  [ 5,  6]      * |         2
   [3]        1  [11, 12]      * |         3

.. figure:: figures/gr.sample.png
   :alt:
   :scale: 125 %

``gr.rand(w=c(2,5,3), seqinfo(gr))``

.. figure:: figures/gr.rand.png
   :alt:
   :scale: 125 %

``gr.simplify``

.. figure:: figures/gr.simplify.png
   :alt:
   :scale: 125 %

``gr.tile(GRanges(1, IRanges(1,9)), w=3) + 1``

.. figure:: figures/gr.tile.png
   :alt:
   :scale: 125 %

``gr.tile.map``
