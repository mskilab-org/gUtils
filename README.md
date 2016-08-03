[![Build Status](https://travis-ci.org/mskilab/gUtils.svg?branch=master)](https://travis-ci.org/mskilab/gUtils)
[![Documentation Status](https://readthedocs.org/projects/gutils/badge/?version=latest)](https://readthedocs.org/projects/gutils/?badge=latest)
[![Coverage Status](https://coveralls.io/repos/github/mskilab/gUtils/badge.svg?branch=master)](https://coveralls.io/github/mskilab/gUtils?branch=master)


gUtils
=======

Set of utility functions for use with GenomicRanges



Installation
------------

1. Install dependent packages and latest Bioconductor (if you haven't already)
  ```source("http://bioconductor.org/biocLite.R") 
     bioclite()
     biocLite(c("GenomicRanges"))
  ```

2. Install devtools from CRAN (if you don't have it already)

  ```
  install.packages('devtools')
  ```

3. Load devtools

  ```
  library(devtools)
  ````

4. Install gUtils

  ```
  install_github('mskilab/gUtils')
  ````

gUtils cheat sheet 
-----

One feature of gUtils is syntactic sugar on top of basic GenomicRanges, enabling easy piping of interval operations as part of interactive "genomic data science" exploration in R. In all these examples `a` and `b` are GRanges (e.g `a` are gene territories and `b` might be copy number segments or Chip-Seq peaks). 

### `%Q%`

Subsets or re-orders `a` based on a logical or integer valued expression that operates on the GRanges metadata columns of `a` 
```{r}
  a %Q% (expr)
  a %Q% (col1  == "value" & col2 > 0 & col3<100)
  a %Q% (order(col1))  
```


### `%*%`

Performs "natural join" or merge of metadata columns of `a` and `b` using interval overlap as a "primary key", outputs a new GRanges whose maximum length is `length(a)*length(b)`. (See `gr.findoverlaps` for more complex queries, including `by` argument that merging based on a hybrid primary key combining both metadata and interval territories)
```{r}	 	  
  a %*% b # strand agnostic merging
  a %**% b # strand specific merging

  ## fancier merges
  gr.findoverlaps(a, b, 
      by = 'key_common_to_a_and_b', qcol = c('acolumn1', 'acolumn2'), scol = c('bcolumn1', 'bcolumn2'))	     
```


### `%$%`
Aggregates the metadata in `b` across the territory of each range in `a`, returning `a` with additional columns of `b` populated with values. For character or factor-valued metadata columns of `b`, aggregation will return a comma collapsed character value of all `b` values (e.g. gene names) that overlap `a[i]`, For numeric columns of `b` it will return the width-weighted mean value (e.g. peak intensity) of that column across the territory `a[i]` and `b` overlap.  For custom aggregations please see `gr.val` function. 
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
Return the subset of ranges in `a` that intersect with `b'
```{r}
  a %&% b # strand agnostic
  a %&&% b # strand specific
```

### `%o%`
Returns a `length(a)` numeric vector whose item `i` is the fraction of `a[i]` that overlaps at least one range in `b' 
```{r}
  a %o% b # strand agnostic
  a %oo% b # strand specific
```

### `%O%`
Returns a `length(a)` numeric vector whose item `i` is the fraction of `a[i]` that overlaps at least one range in `b' 
```{r}
  a %O% b # strand agnostic
  a %OO% b # strand specific
```

### `%N%`
Returns a `length(a)` numeric vector whose item `i` is the total number of ranges in `b` that overlap with `a[i]` 
```{r}
  a %N% b # strand agnostic
  a %NN% b # strand specific
```

### `%^%`
Returns a `length(a)` logical vector whose item `i` TRUE if the  `a[i]` overlaps at least on range in `b' (similar to "%over%" just less fussy)
```{r}
  a %^% b # strand agnostic
  a %^^% b # strand specific
```

### `gr.match`
Returns a `length(a)` integer vector whose item `i` contains the *first* index in `b` overlapping `a[i]` (this function is the match cousin to "%over%" and %^%)
```{r}
  gr.match(a, b) # strand agnostic
  gr.match(a, b, ignore.strand = FALSE) # strand specific	
  gr.match(a, b, by = 'sample_id') # match on metadata column "sample_id" as well as interval
```


Full documentation with examples is available here: [Documentation][docs]

Attributions
------------
> Marcin Imielinski - Assistant Professor, Weill Cornell Medicine; Core Member, New York Genome Center

> Jeremiah Wala - Harvard MD-PhD candidate, Bioinformatics and Integrative Genomics, Rameen Beroukhim Lab, Dana Farber Cancer Institute

[license]: https://github.com/jwalabroad/gTrack/blob/master/LICENSE
[docs]: http://gutils.readthedocs.org/en/latest/index.html

