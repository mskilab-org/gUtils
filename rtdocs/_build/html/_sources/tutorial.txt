Tutorial
--------

``gr1 %^% gr2``        

gives a length(gr1) logical vector with TRUE if gr1 intersects some interval in gr2

``gr1 %*% gr2``

gives a length n ( where n <= length(gr1)*length(gr2)) of all pairwise overlaps of gr1 and gr2 with merged meta.data

``gr1 %$% gr2``       

returns length(gr1) GRanges aggregating all metadata values of gr2 within gr1, taking mean if the meta data item is numeric and concatenating the value if it is a string

``gr1 %Q% (expression)``

subsets gr1 using indices or logical values resulting from expression

``gr1 %Q% (gene== “EGFR”)``

will return the subset of gr1 entries for whose metadata column $gene has “EGFR”

.. code-block:: bash 

   library(gTrack)
   library(gUtils)   

   ## tutorial on using the above operations
   
   ## set current working directory to the folder where the OSMIC TSV was saved in.
   setwd("~/Downloads")
   
   ## load the load 
   OSMICgenes <- read.delim("osmicdata.tsv")

   ## save the gene symbols into a factor
   geneSymbols <- OSMICgenes[,1]

   ## loading gene definitions which we can use to plot windows (GRanges)
   genes = readRDS('files/genes.rds')

   ## save each gene's GRange object in a variable referenced by its name.
   for (i in geneSymbols) { assign(i , genes %Q% (gene_name==i)) }    
   
   TLX3

  '''
   GRanges object with 1 range and 20 metadata columns:
      seqnames                 ranges strand |   source     type     score
         <Rle>              <IRanges>  <Rle> | <factor> <factor> <numeric>
  [1]        5 [170736288, 170739138]      + |   HAVANA     gene      <NA>
          phase           gene_id     transcript_id      gene_type gene_status
      <integer>       <character>       <character>    <character> <character>
  [1]      <NA> ENSG00000164438.5 ENSG00000164438.5 protein_coding       KNOWN
        gene_name transcript_type transcript_status transcript_name     level
      <character>     <character>       <character>     <character> <numeric>
  [1]        TLX3  protein_coding             KNOWN            TLX3         1
               havana_gene         tag havana_transcript exon_number
               <character> <character>       <character>   <numeric>
  [1] OTTHUMG00000163207.3        <NA>              <NA>        <NA>
          exon_id         ont      ccdsid
      <character> <character> <character>
  [1]        <NA>        <NA>        <NA>
  '''
  
  ## this loads coverage data from cancer cell line (GRanges object, have to make it  into a gTrack)
  cov = readRDS('files/coverage.rds')

  genes %$% cov 

  ## loading the GENCODE gene model gTrack
  gt.ge = track.gencode()

  ###############

  genezz <- genes %Q%(gene_name==geneSymbols[1])  

  for (i in 2:length(geneSymbols)) {genezz <- c( genes %Q% (gene_name==geneSymbols[i]) , genezz)}

  genezz <- genezz %$% cov
  
  gt.cov = gTrack(genezz , y.field = 'mean' , circles = TRUE , col = 'blue' , name = 'Cov')

  window = genes[genes$gene_name == "RB1"] + 1.5e6

  plot(c(gt.ge , gt.cov) , window)