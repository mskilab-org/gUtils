[![Build Status](https://travis-ci.org/mskilab/gUtils.svg?branch=master)](https://travis-ci.org/mskilab/gUtils)
[![Documentation Status](https://readthedocs.org/projects/gutils/badge/?version=latest)](https://readthedocs.org/projects/gutils/?badge=latest)


gUtils
=======

Set of utility functions for use with `GenomicRanges`



Installation
------------

1. Install dependent packages and latest Bioconductor (if you haven't already)

```{r}
source("https://bioconductor.org/biocLite.R")
biocLite("GenomicRanges")
```

2. Install devtools from CRAN (if you don't have it already)

```{r}
install.packages('devtools')
```

3. Install gUtils

```{r}
devtools::install_github('mskilab/gUtils')
```

gUtils cheat sheet 
------------------

One feature of gUtils is syntactic sugar on top of basic `GenomicRanges` functionality, enabling easy piping of interval operations as part of interactive "genomic data science" exploration in R. In all these examples `a` and `b` are `GRanges` (e.g `a` are gene territories and `b` might be copy number segments or Chip-Seq peaks). 

### `%Q%`

Subsets or re-orders `a` based on a logical or integer valued expression that operates on the `GRanges` metadata columns of `a`.  
```{r}
  a %Q% (expr)
  a %Q% (col1  == "value" & col2 > 0 & col3 < 100)
  a %Q% (order(col1))  
```


### `%*%`

Performs "natural join" or merge of metadata columns of `a` and `b` using interval overlap as a "primary key", outputs a new `GRanges` whose maximum length is `length(a)*length(b)`. (See `gr.findoverlaps` for more complex queries, including `by` argument that merging based on a hybrid primary key combining both metadata and interval territories).
```{r}	 	  
  a %*% b # strand agnostic merging
  a %**% b # strand specific merging

  ## fancier merges
  gr.findoverlaps(a, b, 
      by = 'column_in_both_a_and_b', qcol = c('acolumn1', 'acolumn2'), scol = c('bcolumn1', 'bcolumn2'))	     
```


### `%$%`
Aggregates the metadata in `b` across the territory of each range in `a`.  This returns `a` appended with additional metadata columns of `b` with values aggregated over the `a` and `b` overlap. For character or factor-valued metadata columns of `b`, aggregation will return a comma collapsed character value of all `b` values (e.g. gene names) that overlap `a[i]`, For numeric columns of `b` it will return the width-weighted mean value (e.g. peak intensity) of that column across the  `a[i]` and `b` overlap.  For custom aggregations please see `gr.val` function. 
```{r}	   
  a %$% b # strand agnostic aggregation
  a %$$% b # strand specific aggregation

  # for additional customization
  # gr.val aggregates and casts data using levels of column "sample_id"				   
  # and a custom function (e.g. max, mode, median) that takes three values as input,
  # where width refers to the width of the overlaps between a[i] and b[jj]
  gr.val(a, b, val = c('field1', 'field2'),
        by = 'sample_id', FUN = function(value, width, is.na) my_cool_fn(value, width, is.na))
	

```

### `%&%`
Return the subset of ranges in `a` that overlap with at least one range in `b`.
```{r}
  a %&% b # strand agnostic
  a %&&% b # strand specific
```

### `%O%`
Returns a `length(a)` numeric vector whose item `i` is the number of bases in `a[i]` that overlaps at least one range in `b`.
```{r}
  a %O% b # strand agnostic
  a %OO% b # strand specific
```


### `%o%`
Returns a `length(a)` numeric vector whose item `i` is the fraction of the width of `a[i]` that overlaps at least one range in `b`.
```{r}
  a %o% b # strand agnostic
  a %oo% b # strand specific
```


### `%N%`
Returns a `length(a)` numeric vector whose item `i` is the total number of ranges in `b` that overlap with `a[i]`.
```{r}
  a %N% b # strand agnostic
  a %NN% b # strand specific
```

### `%^%`
Returns a `length(a)` logical vector whose item `i` TRUE if the  `a[i]` overlaps at least on range in `b` (similar to `%over%` just less fussy about `Seqinfo`).
```{r}
  a %^% b # strand agnostic
  a %^^% b # strand specific
```

### `gr.match`
Returns a `length(a)` integer vector whose item `i` contains the *first* index in `b` overlapping `a[i]` (this function is the match cousin to `%over%` and `%^%`).
```{r}
  gr.match(a, b) # strand agnostic
  gr.match(a, b, ignore.strand = FALSE) # strand specific	
  gr.match(a, b, by = 'sample_id') # match on metadata column "sample_id" as well as interval
```

### `%+%`
Shifts intervals right by `k` bases.
```{r}
  a %+% k
```

### `%-%`
Shifts intervals left by `k` bases.
```{r}
  a %-% k 
```

### `gr.tile`
Tiles `a` or the genome in which `a` resides (as defined by `seqlengths(a)`) with non-overlapping bins of width `w`.
```{r}	  
  gr.tile(a, w) ## outputs non-overlapping tiles of a
  gr.tilexs(seqlengths(a), w) ## outputs non-overlapping tiles of a's genome
  gr.tile(seqlengths(a), 100)+450 # tiles a's genome with 1kbp bins having 900bp overlap
```

### `gr.start`
Returns a `GRanges` of the first coordinate (or first k coordinates) in each interval (in a strand agnostic or specific manner)
```{r}	  
  gr.start(a) # returns the an interval corresponding to the left coordinate
  gr.start(a, k) # returns the first k bases on the left end of a
  
  # returns an interval corresponding to the left coordinate in '+' and '*' ranges and the right coordinate in '-' ranges
  gr.start(a, ignore.strand = FALSE) 
```

### `gr.end`
Returns a `GRanges` of the last coordinate (or last k coordinates) in each interval (in a strand agnostic or specific manner)
```{r}	  
  gr.end(a) # returns the an interval corresponding to the right coordinate
  gr.end(a, k) # returns the last k bases on the right end of a

  # returns an interval corresponding to the right coordinate in '+' and '*' ranges and the left coordinate in '-' ranges
  gr.end(a, ignore.strand = FALSE) 
```


Full documentation with examples is available here: [Documentation][docs]

Attributions
------------
> Marcin Imielinski - Assistant Professor, Weill Cornell Medicine; Core Member, New York Genome Center

> Jeremiah Wala - Harvard MD-PhD candidate, Bioinformatics and Integrative Genomics, Rameen Beroukhim Lab, Dana Farber Cancer Institute

[license]: https://github.com/mskilab/gUtils/blob/master/LICENSE
[docs]: http://gutils.readthedocs.org/en/latest/index.html

