Range Operations
----------------

This section will describe additional GRanges operations provided by gUtils.

Load an example data set. Here, we'll use the refererence gene annotations included with gUtils:

.. code-block:: bash

   ref19 <- readRDS(system.file("extdata","refGene.hg19.gr.rds", package="gUtils"))
   gr <- GRanges(1, IRanges(c(2,5,10), c(4,9,16)), seqinfo=Seqinfo("1", 20)) 
   dt <- data.table(seqnames=1, start=c(2,5,10), end=c(3,8,15))

	
``shift(gr, 2)``

.. figure:: figures/shift.png
   :alt: 
   :scale: 50 %	

``flank(gr, width=2)``

.. figure:: figures/flank.png
   :alt: 
   :scale: 50 %

``gr.start(gr, width=3)``

.. figure:: figures/gr.start.png
   :alt: 
   :scale: 50 %

``gr.end(gr, width=3)``

.. figure:: figures/gr.end.png
   :alt: 
   :scale: 50 %

``gr.mid(gr) + 2``

.. figure:: figures/gr.mid.png
   :alt: 
   :scale: 50 %

``gr.round``

``gr.rand``

``gr.sample``

``grbind`` and ``grlbind``

``streduce``

``gr.sample``

``gr.rand``

``gr.simplify``

``gr.refactor``

``gr.tile``

``gr.tile.amp``

``gr.collapse``
