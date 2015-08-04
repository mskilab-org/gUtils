#############################################################################
## Marcin Imielinski
## Jeremiah Wala
## The Broad Institute of MIT and Harvard / Cancer program.
## marcin@broadinstitute.org
## jwala@broadinstitute.org

## This program is free software: you can redistribute it and/or modify it
## under the terms of the GNU Lesser General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.

## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
## GNU General Public License for more details.

## You should have received a copy of the GNU Lesser General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.
###############################################################################

########################################################################
#
#
##  $$$$$$\  $$$$$$$\        $$\   $$\   $$\     $$\ $$\ 
## $$  __$$\ $$  __$$\       $$ |  $$ |  $$ |    \__|$$ |
## $$ /  \__|$$ |  $$ |      $$ |  $$ |$$$$$$\   $$\ $$ |
## $$ |$$$$\ $$$$$$$  |      $$ |  $$ |\_$$  _|  $$ |$$ |
## $$ |\_$$ |$$  __$$<       $$ |  $$ |  $$ |    $$ |$$ |
## $$ |  $$ |$$ |  $$ |      $$ |  $$ |  $$ |$$\ $$ |$$ |
## \$$$$$$  |$$ |  $$ |      \$$$$$$  |  \$$$$  |$$ |$$ |
##  \______/ \__|  \__|       \______/    \____/ \__|\__|
#
#
# GR Util
#########################################################################

## WORLD'S WORST HACK TO GET NOTES TO DISAPPEAR FOR DATA TABLE CALLS
## WHEN BUILDING FOR CRAN
globalVariables(c("V1", "V2", "V3", "len", "chr", "bin", "count", "rowid", "bin1", "bin2", "newcount",
                  "qname", "reads", "last.line", "uix", "id", "i.start", "i.end", "sn", "subject.id", "query.id",
                  "CIRCOS.DIR"))


#' Get GRanges corresponding to beginning of range
#'
#' Alternative to \code{flank} that will provide start positions *within* intervals
#' 
#' @param x \code{GRanges} object to operate on
#' @param width [default = 1] Specify subranges of greater width including the start of the range.
#' @param force [default = F] Allows returned \code{GRanges} to have ranges outside of its \code{Seqinfo} bounds.
#' @param clip [default = F] Trims returned \code{GRanges} so that it does not extend beyond bounds of the input \code{GRanges}
#' @param ignore.strand [default = T] If set to \code{FALSE}, will extend '-' strands from the other direction.
#' @return GRanges object of width 1 ranges representing start of each genomic range in the input.
#' @import GenomicRanges
#' @examples
#'   \dontrun{st <- c(1,2,3)
#'   si <- Seqinfo("1",3)
#'   gr.start(GRanges(1, IRanges(st, width=101), seqinfo=si), width=20)}
#' @export
gr.start = function(x, width = 1, force = FALSE, ignore.strand = TRUE, clip = FALSE)
  {
    if (length(x)==0)
      return(x)

    width = pmax(width, 1)

    if (any(seqlengths(x)==0) | any(is.na(seqlengths(x))))
        warning('Check or fix seqlengths, some are equal 0 or NA, may lead to negative widths')

    if (force)
      {
        if (ignore.strand)
          {
            st = as.vector(start(x))
            en = as.vector(start(x))+width-1
          }
        else
          {
            st = ifelse(as.logical(strand(x)=='+'),              
              as.vector(start(x)),
              as.vector(end(x))-width+1)
                          
            en = ifelse(as.logical(strand(x)=='+'),
              as.vector(start(x))+width-1,
              as.vector(end(x))
              )
          }       
      }
    else
      {
        if (ignore.strand)
          {
            st = start(x)
            en = pmin(as.vector(start(x))+width-1, seqlengths(x)[as.character(seqnames(x))], na.rm = T)
          }
        else
          {
            st = ifelse(as.logical(strand(x)=='+'),
              as.vector(start(x)),
              pmax(as.vector(end(x))-width+1, 1)
              )

            en = ifelse(as.logical(strand(x)=='+'),
              pmin(as.vector(start(x))+width-1, seqlengths(x)[as.character(seqnames(x))], na.rm = T),
              as.vector(end(x)))
          }        
      }

    if (clip)
      {
        en = pmin(en, end(x))
        st = pmax(st, start(x))
      }
    
    out = GRanges(seqnames(x), IRanges(st, en), seqlengths = seqlengths(x), strand = strand(x))
    values(out) = values(x)
    return(out)    
  }

#' Convert data.table to GRanges
#'
#' Takes as input a data.table which must have the fields: start, end, strand, seqnames.
#' All of the remaining fields are added as meta data to the GRanges
#' @param dt data.table to convert to GRanges
#' @return GRanges object of length = nrow(dt)
#' @import data.table
#' @import IRanges
#' @import GenomicRanges
#' @examples
#' \dontrun{r <- dtgr(data.table(start=1, seqnames="X", end=2, strand='+'))}
#' @export
dtgr <- function(dt) {
  
  rr <- IRanges(dt$start, dt$end)
  if (!'strand' %in% colnames(dt))
    dt$strand <- '*'
  sf <- factor(dt$strand, levels=c('+', '-', '*'))
  ff <- factor(dt$seqnames, levels=unique(dt$seqnames))
  out <- GRanges(seqnames=ff, ranges=rr, strand=sf)
  if (inherits(dt, 'data.table'))
    mc <- as.data.frame(dt[, setdiff(colnames(dt), c('start', 'end', 'seqnames', 'strand')), with=FALSE])
  else if (inherits(dt, 'data.frame'))
    mc <- as.data.frame(dt[, setdiff(colnames(dt), c('start', 'end', 'seqnames', 'strand'))])
  else
    warning("Needs to be data.table or data.frame")
  if (nrow(mc))
    mcols(out) <- mc
  return(out)
}

####################
#' Get GRanges corresponding to beginning of end
#'
#' Alternative to \code{flank} that will provide end positions *within* intervals
#' 
#' @param x \code{GRanges} object to operate on
#' @param width [default = 1] Specify subranges of greater width including the start of the range.
#' @param force [default = FALSE] Allows returned \code{GRanges} to have ranges outside of its \code{Seqinfo} bounds.
#' @param clip [default = FALSE] Trims returned \code{GRanges} so that it does not extend beyond bounds of the input \code{GRanges}
#' @param ignore.strand [default = T] If set to \code{FALSE}, will extend '-' strands from the other direction.
#' @return GRanges object of width 1 ranges representing end of each genomic range in the input.
#' @examples
#'   \dontrun{st <- c(1,2,3)
#'   si <- Seqinfo("1",3)
#'   gr.end(GRanges(1, IRanges(st, width=101), seqinfo=si), width=200, clip=TRUE)}
#' @export
gr.end = function(x, width = 1, force = FALSE, ignore.strand = TRUE, clip = TRUE)
  {
    if (length(x)==0)
      return(x)
    
    if (any(seqlengths(x)==0) | any(is.na(seqlengths(x))))
      warning('Check or fix seqlengths, some are equal 0 or NA, may lead to negative widths')
    
    width = pmax(width, 1)

    if (force)
      {
        if (ignore.strand)
          {
            st = as.vector(end(x))-width+1
            en = as.vector(end(x))
          }
        else
          {
            st = ifelse(as.logical(strand(x)=='+'),
              as.vector(end(x))-width+1,
              as.vector(start(x)))
            
            en = ifelse(as.logical(strand(x)=='+'),
              as.vector(end(x)), 
              as.vector(start(x))+width-1)
          }        
        out = GRanges(seqnames(x), IRanges(st, en), seqlengths = seqlengths(x), strand = strand(x))
      }
    else
      {
        if (ignore.strand)
          {
            st = pmax(as.vector(end(x))-width+1, 1)
            en = as.vector(end(x))
          }
        else
          {
            st = ifelse(as.logical(strand(x)=='+'),
              pmax(as.vector(end(x))-width+1, 1),
              as.vector(start(x)))

            en = ifelse(as.logical(strand(x)=='+'),
              as.vector(end(x)), 
              pmin(as.vector(start(x))+width-1, seqlengths(x)[as.character(seqnames(x))], na.rm = T)
              )
          }
        
        out = GRanges(seqnames(x), IRanges(st, en), seqlengths = seqlengths(x), strand = strand(x))
      }          

    if (clip)
      {
        en = pmin(en, end(x))
        st = pmax(st, start(x))
      }
    
    values(out) = values(x)
    return(out)    
  }

#' Get the midpoint of range
#'
#' @param x \code{GRanges} object to operate on
#' @param pad Amount to pad the output width to (default 0, so width of 1)
#' @return \code{GRanges} of the midpoint, calculated from \code{floor(width(x)/2)}
#' @examples
#' \dontrun{gr.mid(GRanges(1, IRanges(1000,2000), seqinfo=Seqinfo("1", 2000)))}
gr.mid = function(x)
  {
      start(x) = end(x) = rowMeans(cbind(start(x), end(x)))
      return(x)
  }

#' Round a set of GRanges to another set 
#
#' "rounds" a set of query ranges Q to a given interval set S using the following rule:
#' 1) If q in Q is partially / fully within S then return intersection of q with S.
#' 2) If q intersects multiple ranges in S and \code{up = F} then return the "first" range, otherwise the last range
#' 3) If q in Q is fully outside of S (ie fully inside not S) then return the \code{start-1} (if \code{up = T}) or \code{end+1} (if \code{up = F})
#' of the matching range in not S
#'
#' @param Q Query \code{GRanges} (strand is ignored)
#' @param S Subject \code{GRanges} (strand is ignored)
#' @param up [default TRUE] See description.
#' @param parallel [default FALSE] If \code{TRUE}, assumes Q and S are same length and this analysis is only performed between the corresponding Q and S pairs. 
#' @return Rounded \code{GRanges}
#' @examples
#' \dontrun{query   <- GRanges(1, IRanges(c(100,110),width=201), seqinfo=Seqinfo("1", 500))
#' subject <- GRanges(1, IRanges(c(160,170),width=201), seqinfo=Seqinfo("1", 500))
#' gr.round(query, subject)}
#' @export
gr.round = function(Q, S, up = TRUE, parallel = FALSE)
  {
    str = strand(Q)
    Q = gr.stripstrand(Q)
    S = gr.stripstrand(S)
    nS = gaps(S);
    QS = gr.findoverlaps(Q, S)
    tmp = gr.findoverlaps(Q, nS)
    QnotS = nS[tmp$subject.id]
    QnotS$query.id = tmp$query.id
    
    if (parallel)
      {
        QS = QS[QS$query.id==QS$subject.id]
        QnotS = QnotS[QnotS$query.id==QnotS$subject.id]
      }

    if (up)
      suppressWarnings(end(QnotS) <- start(QnotS) <- end(QnotS)+1)
    else
      suppressWarnings(start(QnotS) <- end(QnotS) <- start(QnotS)-1)

    suppressWarnings(out <- sort(grbind(QS, QnotS)))
    
    if (up)
      {
        out = rev(out)
        out = out[!duplicated(out$query.id)]
      }
    else
      out = out[!duplicated(out$query.id)]
    
    out = out[order(out$query.id)]
    values(out) = values(Q)
    names(out) = names(Q)
    strand(out) = str;
    return(out)
  }


#' Generate random GRanges on genome
#'
#' Randomly generates non-overlapping GRanges with supplied widths on supplied genome
#'
#' @param w Vector of widths (length of w determines length of output)
#' @param genome Genome which can be a \code{GRanges}, \code{GRangesList}, or \code{Seqinfo} object. Default is "hg19" from the \code{\link{BSgenome}} package.
#' @param seed [default NA] Optionally specify a seed for the RNG. Defualt behavior is random seed.
#' @return \code{GRanges} with random intervals on the specifed "chromosomes"
#' @note This function is currently quite slow, needs optimization
#' @examples
#' ## Generate a single random interval of width 10, on "chr" of length 1000
#' \dontrun{gr.rand(10, Seqinfo("1", 1000))
#' ## Generate 5 non-overlapping regions of width 10 on hg19
#' gr.rand(rep(10,5))}
#' @export
gr.rand = function(w, genome = Seqinfo(names(hg_seqlengths()), hg_seqlengths()), seed=NA)
  {

    if (!is.na(seed))
      set.seed(seed)
    
    if (!is(genome, 'Seqinfo'))
      genome = seqinfo(genome)
    
    sl = seqlengths(genome);
    available = seqinfo2gr(genome);
    
    out = GRanges(rep(names(sl)[1], length(w)), IRanges(rep(1, length(w)), width = 1), seqlengths = seqlengths(genome));
    for (i in 1:length(w))
      {
        if (i == 1)
          available = seqinfo2gr(genome)
        else
          {
            available = gaps(out[1:(i-1)])
            available = available[strand(available)=='*']
          }
        
        available = available[width(available)>w[i]]

        if (length(available)>0)
          {
            end(available) = end(available)-w[i]
            starts = c(1, cumsum(as.numeric(width(available))+1))
            rstart = ceiling(runif(1)*starts[length(starts)])-starts
            rind = max(which(rstart>0))
            new.chr = seqnames(available[rind])
            new.ir = IRanges(rstart[rind]+start(available[rind])-1, width = w[i])

            ## FIX: this is the slowest part
            seqnames(out)[i] = new.chr;
            ranges(out)[i] = new.ir;
          }
        else
          stop('Allocation failed.  Supplied widths are likely too large')                
      }

    return(out)
  }


#' Sample GRanges intervals
#'
#' Samples k intervals of length "len" from a pile of gr's.
#' If k is a scalar then will (uniformly) select k intervals from the summed territory of gr's
#' If k is a vector of length(gr) then will uniformly select k intervals from each. 
#' from a tiling of the set (and may do fewer than k samples if width(gr[i])<= k[i] *len)
#' If k[i] = NA, will return tiling of that interval, if k = NA will return tiling of the entire
#' gr's (with length len tiles).
#'
#' @param gr \code{GRanges} object to operate on
#' @param k See function description
#' @param len Length param. Default 100
#' @param replace If TRUE, will bootstrap, otherwise will sample without replacement. Default TRUE
#' @return GRanges of max length sum(k) [if k is vector) or k*length(gr) (if k is scalar) with labels indicating the originating range.
#'
#' @note This is different from overloaded sample() function implemented in GenomicRanges class, which just samples from a pile of GRanges
#' @export
gr.sample = function(gr, k, len = 100, replace = TRUE)
{
  if (!inherits(gr, 'GRanges'))
    gr = seqinfo2gr(gr)
  
  if (length(k)==1)
    {
      gr.f = gr.flatten(trim(gr, starts = 1, ends = width(gr)-len), gap = 0);
      terr = sum(gr.f$end-gr.f$start)
      st = gr.f$start;

      if (!replace)
        {
          if (!is.na(k))
            s = len*sample(floor(terr/len), k, replace = FALSE)
          else
            s = seq(1, terr, len)
        }
      else
        s = terr*runif(k)      
     
      # map back to original gr's
#      si = sapply(s, function(x) {i = 1; while (st[i]<x & i<=length(st)) {i = i+1}; return(i)})
      si = rep(NA, length(s))
      for (i in 1:nrow(gr.f))
        si[s>=gr.f$start[i] & s<=gr.f$end[i]] = i
      
      sg = s-st[si]+start(gr)[si];
      
      return(GRanges(seqnames(gr)[si], IRanges(sg, sg+len-1), strand = strand(gr)[si], seqlengths = seqlengths(gr), query.id = si))
    }
  else
    {
      gr.df = data.frame(chr = as.character(seqnames(gr)), start = start(gr), end = end(gr))
      gr.df$k = k;
      gr.df$length = len
      gr.df$replace = replace
      tmp = lapply(1:length(gr), function(i)
                    {
                      if (!gr.df$replace[i])
                        {
                          if (!is.na(k[i]))
                            {
                              w = floor(width(gr)[i]/len)
                              k[i] = min(k[i], w)
                              if (k[i]>0)
                                s = len*sample(w, k[i], replace = FALSE) + gr.df$start[i]
                              else
                                return(NULL)
                            }
                          else
                            s = seq(gr.df$start[i], gr.df$end[i], len)
                        }
                      else
                        s = (gr.df$end[i]-gr.df$start[i]-gr.df$len[i])*runif(k[i])+gr.df$start[i]
                      
                      return(data.frame(chr = gr.df$chr[i], start=s, end =s+len-1, strand = as.character(strand(gr)[i]), query.id = i))
                    })
      return(gr.fix(seg2gr(do.call('rbind', tmp)), gr))
    }
}

footprint = function(gr)
    cat(prettyNum(sum(as.numeric(width(reduce(gr)))), big.mark = ','), '\n')


#' Create GRanges from Seqinfo
#'
#' Creates a genomic ranges from seqinfo object
#' ie a pile of ranges spanning the genome
#' @param si Seqinfo object
#' @param strip.empty Don't know. Default FALSE
#' @examples
#' \dontrun{si <- Seqinfo(names(hg_seqlength(), hg_seqlengths()))
#' si2gr(si)}
#' @export
si2gr <- seqinfo2gr <- function(si, strip.empty = FALSE)
  {
    if (is(si, 'vector')) ## treat si as seqlengths if vector
      si = Seqinfo(seqlengths = si, seqnames = names(si))
    else if (!is(si, 'Seqinfo'))
      si = seqinfo(si)

    sl = seqlengths(si)
    sn = seqnames(si);
    sl[is.na(sl)] = 0;

    if (strip.empty)
      {
        sn = sn[sl!=0];
        sl = sl[sl!=0];
      }

    sigr = GRanges(sn, IRanges(rep(1, length(sl)), width = sl), seqlengths = seqlengths(si), strand = rep('+', length(sl)))
    names(sigr) = sn;
    
    return(sigr)
  }

#' rle.query
#'
#' Queries an \code{\link{RleList}} representing genomic data (ie a list whose names represent
#' seqnames ie chromosomes, and lengths represent seqlengths)
#' via \code{GRanges} object
#'
#' @param subject.rle Subject
#' @param query.gr Query
#' @param verbose Default FALSE
#' @param mc.cores number of cores
#' @param chunksize Default 1e9
#' @return Rle representing the (concatenated) vector of data (reversing order in case of negative strand input)
#' @note Throws warning if seqlengths(gr) do not correspond to the lengths of the \code{RleList} components
#' @export
####################
rle.query = function(subject.rle, query.gr, verbose = FALSE, mc.cores = 1, chunksize = 1e9) ## mc.cores only relevant if there are over 1e9 bases to be queried from subject.rle
  {
    was.grl = FALSE
    
    if (is(query.gr, 'GRangesList'))
    {
      was.grl = TRUE
      query.gr = grl.unlist(query.gr)
    }
    
##     if (!identical(names(subject.rle), seqlevels(query.gr)))
##       warning('seqlevels of subject and query are not the same')

##     com.seq = intersect(seqlevels(query.gr), names(subject.rle))
    
##     if (!identical(sapply(subject.rle[com.seq], length), seqlengths(query.gr)[com.seq]))
##       warning('seqlengths of subject and query are not the same')

    chunksize = pmin(1e9, chunksize)
    if ((sum(as.numeric(width(query.gr))))>chunksize) ## otherwise integer overflow
      {
        tmp = rle(ceiling(cumsum(as.numeric(width(query.gr)))/chunksize))
        chunks = cbind(cumsum(c(1, tmp$lengths[-length(tmp$lengths)])), cumsum(c(tmp$lengths)))
        if (verbose)
          cat(sprintf('chunking up into %s chunks \n', nrow(chunks)))
        out = do.call('c', mclapply(1:nrow(chunks), function(x) rle.query(subject.rle, query.gr[chunks[x,1]:chunks[x,2]]), mc.cores = mc.cores))
      }
    else
      {    
        out = Rle(NA, sum(as.numeric(width(query.gr))));
        
        if (length(query.gr)>1)
          {
            st.ix = cumsum(c(1, width(query.gr)[1:(length(query.gr)-1)]))
          }else
        {
          st.ix = 1    
        }
        out.ix = IRanges(st.ix, st.ix + width(query.gr)-1) ## ranges in out corresponding to query

        for (chr in intersect(names(subject.rle), unique(as.character(seqnames(query.gr)))))
          {
            ix = which(as.character(seqnames(query.gr)) == chr)
            rix = ranges(query.gr)[ix]
            m = max(end(rix))
            if (length(subject.rle[[chr]]) < m) ## pad subject rle if not long enough
              {
                subject.rle[[chr]] = c(subject.rle[[chr]], Rle(NA, m - length(subject.rle[[chr]])))
              }
            out[unlist(as.integer(out.ix[ix]))] = subject.rle[[chr]][rix]
          }

        if("-" %in% as.character(strand(query.gr)))
          {
              tmp = data.table(ix = 1:sum(width(out.ix)), id = rep(1:length(out.ix), width(out.ix)), strand = rep(as.character(strand(query.gr)), width(out.ix)), key = 'ix')
              out = out[tmp[, rev(ix), by = id][, V1]]                                                             
            ## strand.col = c(1, 2)
            ## names(strand.col) = c("+", "-")
            ## cumsums = cumsum(width(out.ix))
            ## exon.max = rep(cumsums, times = width(out.ix))
            ## exon.min = rep(c(0, cumsums[1:(length(cumsums) - 1)]), times = width(out.ix))
            ## negs = exon.max - 0:(length(out)-1) + exon.min 

            ## pos.neg.mat = cbind(1:length(out), negs)
            ## index = pos.neg.mat[cbind(1:length(out), strand.col[strand])]
            ## out = out[index]
          }
      }
    
    if (was.grl)
      out = split(out, Rle(query.gr$grl.ix, width(query.gr)))
            
    return(out)        
  }


#' Concatenate GRanges
#'
#' Concatenates \code{GRanges} objects, taking the union of their features if they have non-overlapping features
#' @param x First \code{GRanges}
#' @param ... Additional \code{GRanges} to concatenate to
#' @note Wraps a call to \code{\link{rrbind}}
#' @note Does not fill in the \code{Seqinfo} for the output \code{GRanges}
#' @return Concatenated \code{GRanges} 
#' gr1 <- GRanges(1, IRanges(100,1000), my.gene='X', seqinfo=Seqinfo("1", 1000))
#' gr2 <- GRanges(1, IRanges(200,2000), my.annot='Y', seqinfo=Seqinfo("1",2000))
#' gr3 <- GRanges(2, IRanges(1,3000), my.annot='Z', seqinfo=Seqinfo("2",3000))
#' grbind(gr1, gr2, gr3)
#' @export
grbind = function(x, ...)
  {
    if (missing('x'))
      grs = list(...)
    else if (class(x) != 'list')
      grs <- c(list(x), list(...))
    else
      grs <- c(x, list(...))

    force.rrbind = FALSE
    if ('force.rrbind' %in% names(grs))
      {
        if (is.logical(grs[['force.rrbind']]))
          force.rrbind = grs[['force.rrbind']]
        grs = grs[-match('force.rrbind', names(grs))]
      }            
          
    #grs = list(...);  # gets list of gr's
    keep = sapply(grs, length)>0 & sapply(grs, function(x) inherits(x, 'GRanges'))
    grs = grs[keep]
    
    if (length(grs)==0)
      return(NULL)

    if (length(grs)==1)
      return(grs[[1]])        
    
    vals = lapply(grs, function(x) values(x))
    
    ## DataFrame from IRanges package can hold XStringSets. Convert first

    isDataFrame <- sapply(vals, class) == 'DataFrame'
    if (any(isDataFrame))
        vals[isDataFrame] = lapply(vals[isDataFrame], grdt)
    
#      if (any(sapply(vals[isDataFrame][[1]], function(y) inherits(y, 'XStringSet'))))
#        stop('grbind: DataFrame from IRanges contains XStringSet objects. Convert to character, then grbind, then convert back')
   ### TOFIX: get this to work with XStringSet
   #   vals <- lapply(vals, function(x) {
   #     isXStringSet <- sapply(x, function(y) inherits(y, 'XStringSet'))
   #     x[,isXStringSet] <- DataFrame(sapply(x[,isXStringSet], as.character))
   #     return(x)
   #   })

    ### FIXING seqlengths when not exactly matching    
    sls = lapply(grs, seqlengths)
    names(sls) = NULL
    levs = unique(names(unlist(sls)))
    sl.new = structure(rep(0, length(levs)), names = levs)
    for (sl in sls)
        sl.new[names(sl)] = pmax(sl.new[names(sl)], sl, na.rm = T)

    bare.grs = lapply(grs, function(x) gr.fix(x[,c()], sl.new))
    out = tryCatch(do.call('c', bare.grs), error = function(e) NULL)
        
    if (is.null(out) | is.list(out)) ## something failed with concatenation, likely some weird ghost with concatenating GRanges with 'c', below is a temp fix
        {
            getsridofghostsomehow =  c(bare.grs[[1]], bare.grs[[2]])
            out = tryCatch(do.call('c', bare.grs), error = function(e) NULL)

            if (is.null(out) | is.list(out)) ## now we are really reaching                
                out = seg2gr(do.call('rrbind', lapply(bare.grs, as.data.frame)), seqlengths = sl.new)[, c()]
        }


    # hack to deal with empty value columns
    ix <- (sapply(vals, ncol)==0)
    if (any(ix))
        vals[ix] = lapply(which(ix), function(x) data.frame(col.4214124124124 = rep(NA, length(grs[[x]]))))

    if (!force.rrbind)
      tmp = tryCatch(do.call('rrbind2', vals), error = function(e) NULL)
    else
      tmp = NULL
    
    if (is.null(tmp) | force.rrbind) ## sometimes rrbind2 gets picky because of type checking (rbindlist) .. so just run rrbind then
       values(out) = do.call('rrbind', vals)
     else
       values(out) = tmp

    if (any(ix))
      out$col.4214124124124 = NULL
    return(out)     
  }

#' Concatenate GRangesList 
#'
#' Concatenates \code{GRangesList} objects taking the union of their \code{mcols} features if they have non-overlapping features
#'
#' @param ... list of \code{GRangesList} object to bind
#' @return Concatenated GRangesList 
#' @examples
#' ## Create some dummy data
#' \dontrun{gr1 <- GRanges(1, IRanges(100,1000), my.gene='X', seqinfo=Seqinfo("1", 1000))
#' gr2 <- GRanges(1, IRanges(200,2000), my.annot='Y', seqinfo=Seqinfo("1",2000))
#' gr3 <- GRanges(2, IRanges(1,3000), my.annot='Z', seqinfo=Seqinfo("2",3000))
#' grl1 <- GRangesList(grbind(gr1, gr2))
#' grl2 <- GRangesList(grbind(gr1, gr3))
#' ## Add unique annotation to just one \code{mcols}.
#' mcols(gr1)$my.new.annot=1
#' ## Concatenate
#' grlbind(grl1, grl2)}
#' @export
grlbind = function(...)
  {
    ## TODO: make this work for when underlying grs do not have matching features
    ## currently will loose gr level features
    grls = list(...)

    ## annoying acrobatics to reconcile gr and grl level features for heterogenous input gr / grls
    grls.ul = lapply(grls, grl.unlist)
    grls.ul.rb = do.call('grbind', grls.ul)
    sp = unlist(lapply(1:length(grls), function(x) rep(x, length(grls.ul[[x]]))))
    gix = split(grls.ul.rb$grl.ix, sp)
    gjx = split(1:length(grls.ul.rb), sp)
    grls.ul.rb$grl.iix = grls.ul.rb$grl.ix = NULL
        
    grls.vals = lapply(grls, function(x)
      { if (ncol(values(x))>0)  return(as.data.frame(values(x))) else return(data.frame(dummy241421 = rep(NA, length(x))))})

    grls.new = mapply(function(x,y) split(grls.ul.rb[x],y), gjx, gix)
    
    out = do.call('c', grls.new)

    if (is.list(out))
      {
        if (length(grls.new)>1)
          {
            bla = c(grls.new[[1]], grls.new[[2]]) ## fix R ghost
            out = do.call('c', grls.new)
            if (is.list(out)) ## if still is list then do manual 'c'
                {
                   out = grls.new[[1]]
                   for (i in 2:length(grls.new))
                       out = c(out, grls.new[[i]])
                }
          }
        else
          out = grls.new[[1]]
      }        
    
    out.val = do.call('rrbind', grls.vals)
    out.val$dummy241421 = NULL
    values(out) = out.val

    return(out)      
  }

#' Add "chr" to GRanges seqlevels
#'
#' Adds "chr" to seqlevels of gr / grl object to make compatible with UCSC genome
#' @param gr GRanges object to append 'chr' to
#' @return Identical GRanges, but with 'chr' appended to each seqlevel
#' @export
gr.chr = function(gr)
  {
   # if (!grepl('^chr', seqlevels(gr)[1]))  # Jeremiah

    if (any(ix <- !grepl('^chr', seqlevels(gr))))
      seqlevels(gr)[ix] = paste('chr', seqlevels(gr)[ix], sep = "")
    return(gr)
  }

#' Dump GRanges to GATK file
#'
#' Dumps gr object into gatk intervals in file path "file"
#' @param gr GRanges
#' @param file file
#' @param add.chr Flag to add "chr" to seqnames. Default FALSE
#' @return returns 0 if completed
#' @export
gr2gatk = function(gr, file, add.chr = FALSE)
  {
    sn = as.character(seqnames(gr));
    if (add.chr)
      sn = paste('chr', sn, sep = '');

    writeLines(paste(sn, ':', start(gr), '-', end(gr), sep = ''), con = file)
    return(0)
  }

#' Shortcut for \code{reduce(sort(gr.stripstrand(unlist(x))))}
#'
#' @param gr takes in gr or grl
#' @param pad asdf. Default 0
#' @param sort Flag to sort the output. Default TR#' @return GRanges
#' @export
streduce = function(gr, pad = 0, sort = TRUE)
  {

    if (inherits(gr, 'GRangesList'))
      gr = unlist(gr)

    if (any(is.na(seqlengths(gr))))
      gr = gr.fix(gr)
    
    #out = suppressWarnings(sort(reduce(gr.stripstrand(gr+pad))))
        out = suppressWarnings(sort(reduce(gr.stripstrand(gr.pad(gr, pad)))))
    suppressWarnings(start(out) <-pmax(1, start(out)))
#    out <- gr.tfix(out)
    end(out) = pmin(end(out), seqlengths(out)[as.character(seqnames(out))])


    return(out)
  }

#' Simplify granges by collapsing all non-overlapping adjacent ranges that share a given "field" value
#' (adjacent == adjacent in the input GRanges object)
#'
#' @param gr takes in gr or grl
#' @param field character scalar, corresponding to value field of gr
#' @param val Default NULL
#' @param include.val scalar logical, will include in out gr values field of first matching record in input gr
#' @param split Default FALSE
#' @param pad Default 1
#' @return Simplified GRanges with "field" populated with uniquely contiguous values
#' @export
gr.simplify = function(gr, field = NULL, val = NULL, include.val = TRUE, split = FALSE, pad = 1)
  {
    tmp = as.logical(suppressWarnings(width(pintersect(ranges(gr[-length(gr)]), ranges(gr[-1]+pad), resolve.empty = 'max.start'))>0) &
      seqnames(gr[-length(gr)]) == seqnames(gr[-1]) & strand(gr[-length(gr)]) == strand(gr[-1]))
    
    tmp = as.vector(c(0, cumsum(!tmp)))
    
    if (!is.null(field))      
      tmp = paste(tmp, values(gr)[, field])

    if (!is.null(val))
      tmp = paste(tmp, val)
    
    r = rle(tmp)
    
    lix = unlist(lapply(1:length(r$lengths), function(x) rep(x, r$lengths[x])))
    sn.gr = split(as.character(seqnames(gr)), lix)
    st.gr = split(start(gr), lix)
    en.gr = split(end(gr), lix)
    str.gr = split(as.character(strand(gr)), lix)

    st.gr.min = sapply(st.gr, min)
    en.gr.max = sapply(en.gr, max)
    
    out = GRanges(sapply(sn.gr, function(x) x[1]), IRanges(st.gr.min, en.gr.max),
            strand = sapply(str.gr, function(x) x[1]), seqinfo = seqinfo(gr))

    ix = match(1:length(out), lix)
    
    if (include.val)
      values(out) = values(gr)[ix, ]
        
    if (split)
      if (!is.null(field))
        out = split(out, values(gr)[ix, field])
      else
        out = GRangesList(out)
            
    return(out)
  }  

#' gr.peaks
#'
#' Finds "peaks" in an input GRanges with value field y.
#' first piles up ranges according to field score (default = 1 for each range)
#' then finds peaks.  If peel > 0, then recursively peels segments 
#' contributing to top peak, and recomputes nextpeak "peel" times
#' if peel>0, bootstrap controls whether to bootstrap peak interval nbootstrap times
#' if id.field is not NULL will peel off with respect to unique (sample) id of segment and not purely according to width
#' if FUN preovided then will complex aggregating function of piled up values in dijoint intervals prior to computing "coverage"
#' (FUN must take in a single argument and return a scalar)
#' if id.field is not NULL, AGG.FUN is a second fun to aggregate values from id.field to output interval
gr.peaks = function(gr, field = 'score', minima = F, peel = 0, id.field = NULL, bootstrap = TRUE, na.rm = TRUE, pbootstrap = 0.95, nbootstrap = 1e4, FUN = NULL, AGG.FUN = sum,
    peel.gr = NULL, ## when peeling will use these segs instead of gr (which can just be a standard granges of scores)
    score.only = FALSE, 
    verbose = peel>0)
    
    {


      if (!is(gr, 'GRanges'))
          gr = seg2gr(gr)
          
      if (is.null(field))
          field = 'score'
      
      if (!(field %in% names(values(gr))))
          values(gr)[, field] = 1

      if (is.logical(values(gr)[, field]))
          values(gr)[, field] = as.numeric(values(gr)[, field])
      
      if (peel>0 & !score.only)
          {
              if (verbose)
                  cat('Peeling\n')
              out = GRanges()

              if (bootstrap)
                  pbootstrap = pmax(0, pmin(1, pmax(pbootstrap, 1-pbootstrap)))

              ## peel.gr are an over-ride if we have pre-computed the score and only want to match peaks to their supporting segments
              if (is.null(peel.gr)) 
                  peel.gr = gr
              
              for (p in 1:peel)
                  {
                      if (verbose)
                          cat('Peel', p, '\n')
                      if (p == 1)
                          last = gr.peaks(gr, field, minima, peel = 0, FUN = FUN, AGG.FUN = AGG.FUN, id.field = id.field)
                      else
                          {
                              ## only need to recompute peak in region containing any in.peak intervals
                              tmp = gr.peaks(gr[gr.in(gr, peak.hood), ], field, minima, peel = 0, FUN = FUN, AGG.FUN = AGG.FUN, id.field = id.field)
                              last = c(last[!gr.in(last, peak.hood)], tmp)
                          }

                      ## these are the regions with the maximum peak value
                      mix = which(values(last)[, field] == max(values(last)[, field]))

                      ## there can be more than one peaks with the same value
                      ## and some are related since they are supported by the same gr
                      ## we group these peaks and define a tmp.peak to span all the peaks that are related 
                      ## to the top peak
                      ## the peak is the span beteween the first and last interval with the maximum 
                      ## peak value that are connected through at least one segment to the peak value

                      ##
                      tmp.peak = last[mix]

                      if (length(tmp.peak)>1)
                          {
                              tmp.peak.gr = gr[gr.in(gr, tmp.peak)]                              
                              ov = gr.findoverlaps(tmp.peak, tmp.peak.gr)
                              ed = rbind(ov$query.id, ov$subject.id+length(tmp.peak))[1:(length(ov)*2)]
                              cl = clusters(graph(ed), 'weak')$membership                              
                              tmp = tmp.peak[cl[1:length(tmp.peak)] %in% cl[1]]
                              peak = GRanges(seqnames(tmp)[1], IRanges(min(start(tmp)), max(end(tmp))))
                              values(peak)[, field] = values(tmp.peak)[, field][1]
                          }
                      else
                          peak = tmp.peak
                      ## tmp.peak is the interval spanning all the top values in this region

                      in.peak1 =  gr.in(peel.gr, gr.start(peak))
                      in.peak2 = gr.in(peel.gr, gr.end(peak))
                      in.peak = in.peak1 | in.peak2

                      ## peak.gr are the gr supporting the peak
                      peak.gr = peel.gr[in.peak1 & in.peak2] ## want to be more strict with segments used for peeling
                      peak.hood = reduce(peak.gr) ## actual peak will be a subset of this, and we can this in further iterations to limit peak revision

                      if (bootstrap)
                          {
                              ## asking across bootstrap smaples how does the intersection fluctuate
                              ## among segments contributing to the peak

                              if (!is.null(id.field))
                                  {                                      
                                      peak.gr = seg2gr(grdt(peak.gr)[, list(seqnames = seqnames[1], start = min(start),
                                          eval(parse(text = paste(field, '= sum(', field, '*(end-start))/sum(end-start)'))),end = max(end)),
                                          by = eval(id.field)])
                                      names(values(peak.gr))[3] = field ## not sure why I need to do this line, should be done above
                                  }
                                                                                                              
                              B = matrix(sample(1:length(peak.gr), nbootstrap * length(peak.gr), prob = abs(values(peak.gr)[, field]), replace = TRUE), ncol = length(peak.gr))
                               ## bootstrap segment samples
                              ## the intersection is tha max start and min end among the segments in each
                              st = apply(matrix(start(peak.gr)[B], ncol = length(peak.gr)), 1, max)
                              en = apply(matrix(end(peak.gr)[B], ncol = length(peak.gr)), 1, min)
                              
                              ## take the left tail of the start position as the left peak boundary
                              start(peak) = quantile(st, (1-pbootstrap)/2)
                              
                              ## and the right tail of the end position as the right peak boundary
                              end(peak) = quantile(en, pbootstrap + (1-pbootstrap)/2)
                              
                              in.peak =  gr.in(gr, peak)
                          }
                      gr = gr[!in.peak]
                      peak$peeled = TRUE
                      out = c(out, peak)
                      if (length(gr)==0)
                          return(out)
                  }
              last$peeled = FALSE
              return(c(out, last[-mix]))              
          }

      if (na.rm)
          if (any(na <- is.na(values(gr)[, field])))
              gr = gr[!na]
      
      if (!is.null(FUN))
          {
              agr = disjoin(gr)
              values(agr)[, field] = NA
              tmp.mat = cbind(as.matrix(values(gr.val(agr[, c()], gr, field, weighted = FALSE, verbose = verbose, by = id.field, FUN = FUN, default.val = 0))))
              values(agr)[, field] = apply(tmp.mat, 1, AGG.FUN)
              gr = agr
          }

      cov = as(coverage(gr, weight = values(gr)[, field]), 'GRanges')

      if (score.only)
          return(cov)

      dcov = diff(cov$score)
      dchrom = diff(as.integer(seqnames(cov)))
      
      if (minima)
          peak.ix = (c(0, dcov) < 0 & c(0, dchrom)==0) & (c(dcov, 0) > 0 & c(dchrom, 0)==0)
      else
          peak.ix = (c(0, dcov) > 0 & c(0, dchrom)==0) & (c(dcov, 0) < 0 & c(dchrom, 0)==0)

      out = cov[which(peak.ix)]
     
      if (minima)
          out = out[order(out$score)]
      else
          out = out[order(-out$score)]
      
      names(values(out)) = field
      
      return(out)
  }

#' gr.string
#'
#' return ucsc style interval string corresponding to gr pile (ie chr:start-end)
#'
#' if mb will return as MB and round to "round"
#'
#' @param gr dummy
#' @param add.chr dummy
#' @param mb dummy
#' @param round dummy
#' @param other.cols dummyz
#' @name gr.string
#' @export
gr.string = function(gr, add.chr = FALSE, mb = FALSE, round = 3, other.cols = c())
  {
      if (length(gr)==0)
          return(as.character(NULL))
    sn = as.character(seqnames(gr));
    if (add.chr)
      sn = paste('chr', sn, sep = '');

    other.cols = intersect(names(values(gr)), other.cols)
    if (length(other.cols)>0)
      other.str = paste(' ', do.call('paste', c(lapply(other.cols, function(x) values(gr)[, x]), list(sep = ' '))))
    else
      other.str = ''

    str = ifelse(as.logical(strand(gr)!='*'), as.character(strand(gr)), '')

    if (mb)      
      return(paste(sn, ':', round(start(gr)/1e6, round), '-', round(end(gr)/1e6, round), str, other.str, sep = ''))
    else
      return(paste(sn, ':', start(gr), '-', end(gr), str, other.str, sep = ''))
  }


#' grl.string
#'
#' return ucsc style interval string corresponding to each gr in each grl, one line per grl item, with gr's in 
#' grl separated by sep
#'
#' if mb will return as MB and round to "round"
#' @param grl GenomicRangesList to get strings from
#' @param sep [Default ","]
#' @param ... dummy
#' @name grl.string
#' @export
grl.string = function(grl, sep = ',', ...)
  {    
    gr = grl.unlist(grl)
    if (!is.null(names(grl)))
      nm = names(grl)
    else
      nm = 1:length(grl)
    
    grs = gr.string(gr, ...)
    out = sapply(split(grs, gr$grl.ix), paste, collapse = sep)
    names(out) = nm[as.numeric(names(out))]
    return(out)
  }

#' parse.grl
#'
#' quick function to parse \code{GRangesList} from character vector IGV / UCSC style strings of format gr1;gr2;gr3 where each gr is of format chr:start-end[+/-]
#'
#' @param x String to parse
#' @param seqlengths [Default \link{hg_seqlengths}
#' @param ... dummy
#' @name parse.grl
#' @export
parse.grl = function(x, seqlengths = hg_seqlengths())
  {
    nm = names(x)
    tmp = strsplit(x, '[;]')
    tmp.u = unlist(tmp)
    tmp.u = gsub('\\,', '', tmp.u)
    tmp.id = rep(1:length(tmp), sapply(tmp, length))
    str = gsub('.*([\\+\\-])$','\\1', tmp.u)
    spl = strsplit(tmp.u, '[\\:\\-\\+]', perl = T)
    if (any(ix <- sapply(spl, length)!=3))
      spl[ix] = strsplit(gr.string(seqinfo2gr(seqlengths)[sapply(spl[ix], function(x) x[[1]])], mb = F), '[\\:\\-\\+]', perl = T)

    if (any(ix <- !str %in% c('+', '-')))
      str[ix] = '*'    
    df = cbind(as.data.frame(matrix(unlist(spl), ncol = 3, byrow = T), stringsAsFactors = F), str)
    names(df) = c('chr', 'start', 'end', 'strand')
    rownames(df) = NULL
    gr = seg2gr(df, seqlengths = seqlengths)[, c()]
    grl = split(gr, tmp.id)
    names(grl) = nm
    return(grl)
  }

#' parse.gr
#'
#' quick function to parse gr from character vector IGV / UCSC style strings of format gr1;gr2;gr3 where each gr is of format chr:start-end[+/-]
#'
#' @name parse.grl
#' @export
gstring = parse.gr = function(...)
  {
    return(unlist(parse.grl(...)))
  }


#' ra.overlaps
#'
#' Determines overlaps between two piles of rearrangement junctions ra1 and ra2 (each GRangesLists of signed locus pairs)
#' against each other, returning a sparseMatrix that is T at entry ij if junction i overlaps junction j.
#'
#' if argument pad = 0 (default) then only perfect overlap will validate, otherwise if pad>0 is given, then
#' padded overlap is allowed
#'
#' strand matters, though we test overlap of both ra1[i] vs ra2[j] and gr.flip(ra2[j])
#'
#' @param ra1 \code{GRangesList} with rearrangement set 1
#' @param ra2 \code{GRangesList} with rearrangement set 2
#' @param pad Amount to pad the overlaps by. Larger is more permissive. Default is exact (0)
#' @param arr.ind Default TRUE
#' @param ignore.strand Ignore rearrangement orientation when doing overlaps. Default FALSE
#' @param ... params to be sent to \code{\link{gr.findoverlaps}}
#' @name ra.overlaps
#' @importFrom Matrix sparseMatrix
#' @export
ra.overlaps = function(ra1, ra2, pad = 0, arr.ind = T, ignore.strand=FALSE, ...)
  {    
    bp1 = grl.unlist(ra1) + pad
    bp2 = grl.unlist(ra2) + pad 
    ix = gr.findoverlaps(bp1, bp2, ignore.strand = ignore.strand, ...)
#    ix.rev = gr.findoverlaps(bp1, gr.flip(bp2), ignore.strand = F) ## match even if flipped

    .make_matches = function(ix, bp1, bp2)
      {
        if (length(ix) == 0)
          return(NULL)
        tmp.match = cbind(bp1$grl.ix[ix$query.id], bp1$grl.iix[ix$query.id], bp2$grl.ix[ix$subject.id], bp2$grl.iix[ix$subject.id])
        tmp.match.l = lapply(split(1:nrow(tmp.match), paste(tmp.match[,1], tmp.match[,3])), function(x) tmp.match[x, , drop = F])
        
        ## match only occurs if each range in a ra1 junction matches a different range in the ra2 junction
        matched.l = sapply(tmp.match.l, function(x) all(c('11','22') %in% paste(x[,2], x[,4], sep = '')) | all(c('12','21') %in% paste(x[,2], x[,4], sep = '')))
        
        return(do.call('rbind', lapply(tmp.match.l[matched.l], function(x) cbind(x[,1], x[,3])[!duplicated(paste(x[,1], x[,3])), , drop = F])))        
      }


#    tmp = rbind(.make_matches(ix, bp1, bp2), .make_matches(ix.rev, bp1, bp2))
#    rownames(tmp) = NULL

    tmp = .make_matches(ix, bp1, bp2)
    
    if (is.null(tmp))
      {
        if (arr.ind)
          return(matrix())
        else
          return(sparseMatrix(length(ra1), length(ra2), x = 0))
      }
    
    rownames(tmp) = NULL
    
    colnames(tmp) = c('ra1.ix', 'ra2.ix')

    if (arr.ind) {
      ro <- tmp[order(tmp[,1], tmp[,2]), ]
      if (class(ro)=='integer')
        ro <- matrix(ro, ncol=2, nrow=1, dimnames=list(c(), c('ra1.ix', 'ra2.ix'))) 
      return(ro)
    } else {
      ro <- sparseMatrix(tmp[,1], tmp[,2], x = 1, dims = c(length(ra1), length(ra2)))
      return(ro)
    }
  }

#' gr.sub
#'
#' will apply gsub to seqlevels of gr, by default removing 'chr', and "0.1" suffixes, and replacing "MT" with "M"
#' @name gr.sub
#' @export
gr.sub = function(gr, a = c('(^chr)(\\.1$)', 'MT'), b= c('', 'M'))
  {    
    tmp = mapply(function(x, y) seqlevels(gr) <<- gsub(x, y, seqlevels(gr)), a, b)
    return(gr)
  }

#' gr.fix
#'
#' "Fixes" seqlengths / seqlevels 
#' If "genome" not specified will replace NA seq lengths in GR to reflect largest coordinate per seqlevel
#' and removes all NA seqlevels after this fix. 
#'
#' if "genome" defined (ie as Seqinfo object, or a BSgenome, GRanges, GRnagesList object with populated seqlengths) then will replace
#' seqlengths in gr with those for that genome (and if drop = T, drop all ranges without
#' seqlevels in that genome)
#' @param gr \code{GRanges} to be fixed
#' @param genome Genome to fix to
#' @param gname dummy
#' @param drop Remove all ranges without seqlevelts in the genome. Default FALSE
#' @name gr.fix
#' @import GenomicRanges
#' @import data.table
#' @export
gr.fix = function(gr, genome = NULL, gname = NULL,  drop = FALSE)
  {
    #### marcin: now it is 
    ## if (inherits(gr, "GRangesList"))
    ##   stop("gr.fix not setup to take GRangesList")
    
    if (!is.null(genome))
      {
        if (is.vector(genome))  ## assume seqlengths was provided (ie named vector of lengths)
          genome = Seqinfo(names(genome), seqlengths = genome)
        else if (!(is(genome, 'character') | inherits(genome, 'GRanges') | inherits(genome, 'BSgenome') | inherits(genome, 'GRangesList') | inherits(genome, 'Seqinfo')) & !is.vector(genome))
          genome = seqinfo(genome)
        
        if (!is.vector(genome))
          if (!drop)
            {
              levs = union(seqlevels(genome), as.character(seqnames(seqinfo(gr))))
              lens = structure(rep(NA, length(levs)), names = levs)
              lens[seqlevels(genome)] = seqlengths(genome);
              lens[seqlevels(gr)] = pmax(seqlengths(gr), lens[seqlevels(gr)], na.rm = T)
              
            }
          else
            lens = structure(seqlengths(genome), names = seqlevels(genome))              
        else
          {
            if (is.character(genome))
              genome = structure(rep(NA, length(genome)), names = genome)
            
            lens = structure(NA, names = union(names(genome), seqlevels(gr)));
            
            lens[seqlevels(gr)] = seqlength(gr);
            lens[names(genome)[!is.na(genome)]] = pmax(lens[names(genome)[!is.na(genome)]], genome[!is.na(genome)], na.rm = T)

            if (drop)
              lens = lens[names(genome)]
          }
                              
        seqlevels(gr, force = T) = names(lens)
        seqlengths(gr) = lens;
      }
    else
      {          
        if (length(gr)>0)
            {
                if (is(gr, 'GRangesList'))
                    tmp.gr = unlist(gr)
                else
                    tmp.gr = gr
                
                tmp.sl = data.table(sn = as.character(seqnames(tmp.gr)), end = end(tmp.gr))[, max(end, na.rm = T), by = sn][,  structure(V1, names = sn)][seqlevels(tmp.gr)]
                names(tmp.sl) = seqlevels(tmp.gr)
                seqlengths(tmp.gr)[!is.na(tmp.sl)] = suppressWarnings(pmax(tmp.sl[!is.na(tmp.sl)], seqlengths(tmp.gr)[!is.na(tmp.sl)], na.rm = T))
                gr = tmp.gr
            }
        if (drop)
          gr = keepSeqlevels(gr, seqlevels(gr)[!is.na(seqlengths(gr))])

      }

    ## hack to get rid of annoying "genome"
    if (!is.null(gname))
        {
            si = seqinfo(gr)
            genome(si) = gname
            gr@seqinfo = si
        }
    
    return(gr)
  }


#' gr.flatten
#'
#' Takes pile of GRanges and returns into a data.frame with nrow = length(gr) with each
#' representing the corresponding input range superimposed onto a single "flattened"
#' chromosome.
#' @param gr \code{GRanges} to flatten
#' @param gap Default 0
#' @name gr.flatten
#' @export
###########################################
gr.flatten = function(gr, gap = 0)
  {
    if (length(gr) == 0)
      return(data.frame())
    else if (length(gr) == 1)
      return(data.frame(start = 1, end = width(gr)))
    else
      {
        starts = cumsum(c(1, width(gr[1:(length(gr)-1)])+gap))
        ends = starts+width(gr)-1
        #return(data.frame(start = starts, end = ends)) ## MARCIN
		return(cbind(data.frame(start=starts, end=ends), as.data.frame(mcols(gr)))) ## JEREMIAH
      }
}

#' gr.refactor
#'
#' Takes a pile of ranges gr and new seqnames "sn" (either of length 1 or
#' of length(gr)) and returns a gr object with the new seqnames and same
#' widths and new start coordinates.  These coordinates are determined by placing
#' each gr on the corresponding chromosome at 1 + gap after previous gr (or at 1)
#' @param gr \code{GRanges} to refactor
#' @param sn character vector of new seqnames
#' @param gap Default 0
#' @param rev Default FALSE
#' @name gr.refactor
#' @export
gr.refactor = function(gr, sn, gap = 0, rev = FALSE)
  {    
    if (is.factor(sn))
      slev = levels(sn)
    else
      slev = unique(sn);
    
    sn = cbind(as.character(start(gr)), as.character(sn))[,2]
    w = width(gr)
    gap = pmax(cbind(gap, w)[,1], 0);
#    gap = pmax(gap, 0);
#    starts = levapply(width(gr), sn, function(x) cumsum(gap+c(1, x[1:(length(x)-1)]))[1:length(x)]-gap)

    starts = levapply(1:length(w), sn, function(x) cumsum(gap[x] + c(1, w[x[1:(length(x)-1)]])[1:length(x)])-gap[x])
    ir = IRanges(starts, width = width(gr))

    # figure out seqlevels so that the order matches seqlevels of gr with
    # derivative chromosomes "next" to their original
    sl = aggregate(end(ir)+gap, by = list(sn), FUN = max); sl = structure(sl[,2], names = sl[,1])

    # reorder and add any missing levels
    oth.names = setdiff(slev, names(sl))
    if (length(oth.names)>0)
      sl[oth.names] = NA
    sl = sl[slev]    
    
    out = GRanges(sn, ir, strand = strand(gr), seqlengths = sl)
    values(out) = values(gr);
      
    return(out)
  }

#' gr.stripstrand
#'
#' sets strand to "*"
#' @param gr \code{GRanges} to remove the strand from
#' @name gr.stripstrand
#' @export
gr.stripstrand = function(gr)
  {
    strand(gr) = "*"
    return(gr)
  }

#' gr.flip
#'
#' flips strand on grs
#' optional arg which will determine which ones to flip
#' @param gr \code{GRanges} to flip
#' @param which Default TRUE
#' @name gr.flip
#' @export
gr.flip = function(gr, which = TRUE)
  {    
    if (!is(gr, 'GRanges'))
      stop('GRanges input only')
    
    if (length(gr)==0)
      return(gr)
    
    which = cbind(1:length(gr), which)[,2] == 1

    if (any(which))
      strand(gr)[which] = c('*'='*', '+'='-', '-'='+')[as.character(strand(gr))][which]
    
    return(gr)
  }

#' gr.pairflip
#'
#' "pairs" gr returning a grl with each item consisting
#' of the original gr and its strand flip
#' @param gr \code{GRanges}
#' @name gr.pairflip
#' @export
gr.pairflip = function(gr)
  {
    strand(gr)[strand(gr) =='*'] = '+';
    return(split(c(gr, gr.flip(gr)), rep(c(1:length(gr)), 2)))
  }



#' gr.tostring
#'
#' dumps out a quick text representation of a gr object (ie a character vector)
#' @param gr \code{GRanges}
#' @param places Number of decimal places. Default 2
#' @param interval Default 1e6
#' @param unit Default "MB"
#' @param prefix Default "chr"
#' @return text representation of input
#' @name gr.tostring
#' @export
gr.tostring = function(gr, places = 2, interval = 1e6, unit = 'MB', prefix = 'chr')
{
  p1 = round(start(gr)/interval, places);
  p2 = round(end(gr)/interval, places);
  return(paste(prefix, as.character(seqnames(gr)), ':', p1, '-', p2, ' ', unit, sep = ''));
}

#' gr.tile
#'
#' tiles interval (or whole genome) with segments of <= specified width.  Returns strandless gr
#' "tiles". 
#'
#' input can be seqinfo object (in which case whole genome will be tiled);
#' if inputted grs overlap, will first reduce then tile.
#' @name gr.tile
#' @param gr GRanges to tile, also can be seqlengths or seqinfo
#' @param w Width of the binds. Default 1e3
#' @export
gr.tile = function(gr, w = 1e3)
  {
    if (!is(gr, 'GRanges'))
      gr = seqinfo2gr(gr);

    ix = which(width(gr)>0)
    gr = gr[ix] 

    if (length(gr)==0)
        return(gr[c()][, c()])
    sn = as.character(seqnames(gr))
    str = as.character(strand(gr))
    pos = lapply(1:length(gr), function(x) c(seq(start(gr)[x], end(gr)[x], w), end(gr)[x]+1));
    starts = lapply(pos, function(x) x[1:(length(x)-1)])
    ends = lapply(pos, function(x) x[2:length(x)]-1);
    chr = lapply(1:length(gr), function(x) rep(sn[x], length(starts[[x]])))
    strs = lapply(1:length(gr), function(x) rep(str[x], length(starts[[x]])))
    query.id = lapply(1:length(gr), function(x) rep(x, length(starts[[x]])))
    out = GRanges(unlist(chr), IRanges(unlist(starts), unlist(ends)), strand = unlist(strs), seqlengths = seqlengths(gr))
    out$query.id = ix[unlist(query.id)]
    out$tile.id = unlist(lapply(query.id, function(x) 1:length(x)))
    return(out)
  }


#' gr.flatmap
#'
#' Takes gr (Granges object) and maps onto a flattened coordinate system defined by windows (GRanges object)
#' a provided "gap" (in sequence units).  If squeeze == T then will additionally squeeze ranges into xlim.
#'
#' output is list with two fields corresponding to data frames:
#' $grl.segs = data frame of input gr's "lifted" onto new flattened coordinate space (NOTE: nrow of this not necessarily equal to length(gr))
#' $window.segs = the coordinates of input windows in the new flattened (and squeezed) space
#'
#' @param gr \code{GRanges} to flatten
#' @param windows \code{GRanges} to flatten onto
#' @param gap Default 0
#' @param strand.agnostic Default FALSE
#' @param squeeze Default FALSE. If TRUE, then will additionally squeeze ranges into xlim
#' @param xlim Default c(0,1)
#' @param pintersect Default FALSE
#' @return output is list with two fields corresponding to data frames:
#'   $grl.segs = data frame of input gr's "lifted" onto new flattened coordinate space (NOTE: nrow of this not necessarily equal to length(gr))
#'   $window.segs = the coordinates of input windows in the new flattened (and squeezed) space
#' 
#' FIX: turn this into Chain object
#' @name gr.flatmap
#' @export
gr.flatmap = function(gr, windows, gap = 0, strand.agnostic = T, squeeze = F, xlim = c(0, 1), pintersect=FALSE)
  {
    if (strand.agnostic)
      strand(windows) = "*"

    ## now flatten "window" coordinates, so we first map gr to windows
    ## (replicating some gr if necessary)
#    h = findOverlaps(gr, windows)

    h = gr.findoverlaps(gr, windows, pintersect=pintersect);

    window.segs = gr.flatten(windows, gap = gap)

    grl.segs = as.data.frame(gr);
    grl.segs = grl.segs[values(h)$query.id, ];
    grl.segs$query.id = values(h)$query.id;
    grl.segs$window = values(h)$subject.id
    grl.segs$start = start(h);
    grl.segs$end = end(h);
    grl.segs$pos1 = pmax(window.segs[values(h)$subject.id, ]$start,
      window.segs[values(h)$subject.id, ]$start + grl.segs$start - start(windows)[values(h)$subject.id])
    grl.segs$pos2 = pmin(window.segs[values(h)$subject.id, ]$end,
      window.segs[values(h)$subject.id, ]$start + grl.segs$end - start(windows)[values(h)$subject.id])
    grl.segs$chr = grl.segs$seqnames

    if (squeeze)
      {
        min.win = min(window.segs$start)
        max.win = max(window.segs$end)    
        grl.segs$pos1 = affine.map(grl.segs$pos1, xlim = c(min.win, max.win), ylim = xlim)
        grl.segs$pos2 = affine.map(grl.segs$pos2, xlim = c(min.win, max.win), ylim = xlim)
        window.segs$start = affine.map(window.segs$start, xlim = c(min.win, max.win), ylim = xlim)
        window.segs$end = affine.map(window.segs$end, xlim = c(min.win, max.win), ylim = xlim)
       }
    
    return(list(grl.segs = grl.segs, window.segs = window.segs))
  }

#' Affinely maps 1D points
#'
#' affinely maps 1D points in vector x from interval xlim to interval ylim,
#' ie takes points that lie in 
#' interval xlim and mapping onto interval ylim using linear / affine map defined by:
#' (x0,y0) = c(xlim(1), ylim(1)),
#' (x1,y1) = c(xlim(2), ylim(2))
#' (using two point formula for line)
#' useful for plotting.
#'
#' @param x vector of 1D points
#' @param ylim interval to map onto
#' @param xlim interval to map from. Default c(min(x), max(x))
#' @param cap Default FALSE
#' @param cap.min Default \code{cap}
#' @param cap.max Default \code{cap}
#' @param clip Default TRUE
#' @param clip.min Default \code{clip}
#' @param clip.max Default \code{clip.max}
#' if cap.max or cap.min == T then values outside of the range will be capped at min or max
#' @export
affine.map = function(x, ylim = c(0,1), xlim = c(min(x), max(x)), cap = F, cap.min = cap, cap.max = cap, clip = T, clip.min = clip, clip.max = clip)
  {
  #  xlim[2] = max(xlim);
  #  ylim[2] = max(ylim);
    
    if (xlim[2]==xlim[1])
      y = rep(mean(ylim), length(x))
    else
      y = (ylim[2]-ylim[1]) / (xlim[2]-xlim[1])*(x-xlim[1]) + ylim[1]

    if (cap.min)
      y[x<min(xlim)] = ylim[which.min(xlim)]
    else if (clip.min)
      y[x<min(xlim)] = NA;
    
    if (cap.max)
      y[x>max(xlim)] = ylim[which.max(xlim)]
    else if (clip.max)
      y[x>max(xlim)] = NA;
    
    return(y)
  }


#' Faster version of GRanges::findOverlaps
#'
#' (faster) replacement for GRanges version of findOverlaps
#' returns granges of matches with two additional fields
#' $query.id - index of matching query
#' $subject.id - index of matching subject
#'
#' pintersect employs pintersect to find overlaps, in general this is slower, but can be much faster with much lower
#' memory footprint for large ranges sets with many different seqnames (eg transcriptome)
#' max.chunk controls the maximum number of range pairs that compared at once 
#'
#' Optional "by" field is a character scalar that specifies a metadata column present in both query and subject
#' that will be used to additionally restrict matches, i.e. to pairs of ranges that overlap and also
#' have the same values of their "by" fields
#'
#' @param query Query \code{GRanges}
#' @param subject Subject \code{GRanges}
#' @param ignore.strand [Default TRUE] Ignore the strand when doing overlap queries
#' @param first [Default FALSE]
#' @param qcol any query metadata columns to add to result
#' @param scol any subject metadata columns to add to result
#' @param max.chunk [Default 1e13] If query is bigger than this, chunks into smaller pieces and sends to multi-cores
#' @param foverlaps Should we use data.table::foverlaps? Auto detects this
#' @param pintersect Should we use IRanges::pintersect? Auto determines this if NA
#' @param verbose Default FALSE
#' @param type Default "any"
#' @param by Do overlaps within groups
#' @param mc.cores Default 1. Only active if exceeded max.chunk (ideally should not use)
#' @param return.type Default "same"
#' @param ... = additional args for findOverlaps (IRanges version)
#' @name gr.findoverlaps
#' @export 
gr.findoverlaps = function(query, subject, ignore.strand = TRUE, first = FALSE,
    qcol = NULL, ## any query meta data columns to add to result
    scol = NULL, ## any subject meta data columns to add to resultx
    max.chunk = 1e13,
    foverlaps = ifelse(is.na(as.logical(Sys.getenv('GRFO_FOVERLAPS'))), TRUE, as.logical(Sys.getenv('GRFO_FOVERLAPS'))) & exists('foverlaps'),
    pintersect = NA,
    verbose = F,
    type = 'any', 
    by = NULL, 
    mc.cores = 1,
    return.type = 'same',
    ...)
  {

     if (type != 'any')
         {
             foverlaps = FALSE
             pintersect = FALSE
         }
      
  if (nchar(foverlaps)==0)
      foverlaps = TRUE

  if (is.na(foverlaps))
      foverlaps = TRUE
  
  isdt <- any(class(query) == 'data.table' )
  if (return.type == 'same')
    return.type <- ifelse(isdt, 'data.table', 'GRanges')
  
  if (!((inherits(subject, 'GRanges') | inherits(subject, 'data.table')) & (inherits(query, 'GRanges') | inherits(query, 'data.table'))))
      stop('both subject and query have to be GRanges or data.table')
  
  if (is.na(pintersect))
    if (isdt)
      pintersect <- length(unique(query$seqnames)) > 50 & length(unique(subject$seqnames)) > 50
    else
      pintersect <- seqlevels(query) > 50 && seqlevels(subject) > 50
  if (is.na(pintersect))
    pintersect <- FALSE

  
  if (!is.null(qcol))
      if (!all(qcol %in% names(values(query))))
          stop('Some qcol are not present in meta data of query')

  if (!is.null(scol))
      if (!all(scol %in% names(values(subject))))
          stop('Some scol are not present in meta data of subject')
          
  if (!is.null(by))
    if (!(by %in% names(values(query)) & by %in% names(values(subject))))
      stop('"by" field must be meta data column of both query and subject')
    
    if ((as.numeric(length(query)) * as.numeric(length(subject))) > max.chunk)
      {
        if (verbose) 
          cat('Overflow .. computing overlaps in chunks.  Adjust max.chunk parameter to gr.findoverlaps to avoid chunked computation\n')
        chunk.size = floor(sqrt(max.chunk));
        ix1 = c(seq(1, length(query), chunk.size), length(query)+1)
        ix2 = c(seq(1, length(subject), chunk.size), length(subject)+1)
        ij = cbind(rep(1:(length(ix1)-1), length(ix2)-1), rep(1:(length(ix2)-1), each = length(ix1)-1))
        if (verbose)
          print(paste('Number of chunks:', nrow(ij)))

        out = do.call('c', mclapply(1:nrow(ij),
            function(x)
                        {
                          if (verbose)
                            cat(sprintf('chunk i = %s-%s (%s), j = %s-%s (%s)\n', ix1[ij[x,1]], ix1[ij[x,1]+1]-1, length(query),
                                        ix2[ij[x,2]], (ix2[ij[x,2]+1]-1), length(subject)))
                          i.chunk = ix1[ij[x,1]]:(ix1[ij[x,1]+1]-1)
                          j.chunk = ix2[ij[x,2]]:(ix2[ij[x,2]+1]-1)
                          out = gr.findoverlaps(query[i.chunk], subject[j.chunk],  ignore.strand = ignore.strand, first = first, pintersect=pintersect, by = by, qcol = qcol, verbose = verbose, foverlaps = foverlaps, scol = scol, type = type, ...)
                          out$query.id = i.chunk[out$query.id]
                          out$subject.id = j.chunk[out$subject.id]
                          return(out)
                        }, mc.cores=mc.cores))

        convert = FALSE
        if ((return.type == 'same' & is(query, 'data.table')) | return.type == 'data.table')
            out = grdt(out)
        return(out)            
      }

  if (foverlaps)
      {
          if (verbose)
              print('overlaps by data.table::foverlaps')
          if (ignore.strand)
              by = c(by, 'seqnames',  'start', 'end')
          else
              by = c(by, 'seqnames', 'strand', 'start', 'end')

          if (!is.data.table(query))
              {
                  names(query) = NULL
                  querydt = grdt(query[, setdiff(by, c('seqnames', 'start', 'end', 'strand'))])
              }
          else
              {
                  if (!all(by %in% names(query)))
                      stop(paste('the following columns are missing from query:',
                                 paste(by, collapse = ',')))
                      
                  querydt = query[, by, with = FALSE]
              }
          
          if (!is.data.table(subject))
              {
                  names(subject) = NULL
                  subjectdt = grdt(subject[, setdiff(by, c('seqnames', 'start', 'end', 'strand'))])
              }
          else
              {
                  if (!all(by %in% names(subject)))
                      stop(paste('the following columns are missing from subejct:',
                                 paste(by, collapse = ',')))
                  subjectdt = subject[, by, with = FALSE]
              }
          
          
          ix1 = querydt$query.id = 1:nrow(querydt)
          ix2 = subjectdt$subject.id = 1:nrow(subjectdt)

          querydt = querydt[start<=end, ]
          subjectdt = subjectdt[start<=end, ]
          
          querydt = querydt[, c('query.id', by), with = F]
          subjectdt = subjectdt[, c('subject.id', by), with = F]
          setkeyv(querydt, by)
          setkeyv(subjectdt, by)

         
          h.df = foverlaps(querydt, subjectdt, by.x = by, by.y = by, mult = 'all', type = 'any', verbose = verbose)
          h.df = h.df[!is.na(subject.id) & !is.na(query.id), ]
          h.df[, start := pmax(start, i.start)]
          h.df[, end := pmin(end, i.end)]

          if (verbose)
              cat(sprintf('Generated %s overlaps\n', nrow(h.df)))          
      }
  else
      {

          if (isdt) {
              sn1 <- query$seqnames
              sn2 <- subject$seqnames
          } else {
              sn1 = as.character(seqnames(query))
              sn2 = as.character(seqnames(subject))
          }
          if (is.null(by))
              {
                  ix1 = which(sn1 %in% sn2)
                  ix2 = which(sn2 %in% sn1)
              }
          else
              {
                  by1 = values(query)[, by]
                  by2 = values(subject)[, by]
                  ix1 = which(sn1 %in% sn2 & by1 %in% by2)
                  ix2 = which(sn2 %in% sn1 & by2 %in% by1)
                  by1 = by1[ix1]
                  by2 = by2[ix2]
              }
          
          query.ix = query[ix1]
          subject.ix = subject[ix2]
          sn1 = sn1[ix1]
          sn2 = sn2[ix2]
          
          
          if (pintersect)
              {
                  if (verbose)
                      print('overlaps by pintersect')
                  if (length(sn1)>0 & length(sn2)>0)
                      {

                          if (is.null(by))
                              {
                                  dt1 <- data.table(i=seq_along(sn1), sn=sn1, key="sn")
                                  dt2 <- data.table(j=seq_along(sn2), sn=sn2, key="sn")                
                                  ij <- merge(dt1, dt2, by = 'sn', allow.cartesian=TRUE)
                              }
                          else
                              {
                                  dt1 <- data.table(i=seq_along(sn1), sn=sn1, by = by1, key=c("sn", "by"))
                                  dt2 <- data.table(j=seq_along(sn2), sn=sn2, by = by2, key=c("sn", "by"))
                                  ij <- merge(dt1, dt2, by = c('sn', 'by'), allow.cartesian=TRUE)
                              }

                          if (ignore.strand && isdt)
                              subject$strand <- '*'
                          else if (ignore.strand)
                              strand(subject) = '*'

                          qr <- query.ix[ij$i]
                          sb <- subject.ix[ij$j]
                          if (!isdt) {
                              seqlengths(qr) <- rep(NA, length(seqlengths(qr)))
                              seqlengths(sb) <- rep(NA, length(seqlengths(sb)))
                          }
                          
                          if (!isdt && any(as.character(seqnames(qr)) != as.character(seqnames(sb))))
                              warning('gr.findoverlaps: violated pintersect assumption')

                          ## changed to ranges(qr) etc rather than just GRanges call. Major problem if too many seqlevels
                          if (isdt) {
                              rqr <- IRanges(start=qr$start, end=qr$end)
                              rsb <- IRanges(start=sb$start, end=sb$end)              
                          } else {
                              rqr <- ranges(qr)
                              rsb <- ranges(sb)              
                          }
                          tmp <- pintersect(rqr, rsb, resolve.empty = 'start.x', ...)
                          names(tmp) = NULL
                          non.empty = which(width(tmp)!=0)
                          h.df = as.data.frame(tmp[non.empty])
                          if (isdt)
                              h.df$seqnames <- qr$seqnames[non.empty]
                          else
                              h.df$seqnames <- as.character(seqnames(qr))[non.empty]
                          h.df$query.id = ij$i[non.empty]
                          h.df$subject.id = ij$j[non.empty]
                      }
                  else
                      h.df = data.frame()
              }
          else
              {
                  if (verbose)
                      print('overlaps by findOverlaps')
                  if (isdt) {
                      rqr <- IRanges(start=query.ix$start, end=query.ix$end)
                      rsb <- IRanges(start=subject.ix$start, end=subject.ix$end)              
                  } else {
                      rqr <- ranges(query.ix)
                      rsb <- ranges(subject.ix)              
                  }

                  h <- findOverlaps(rqr, rsb, type = type)
                  r <- ranges(h, rqr, rsb)   
                  h.df <- data.frame(start = start(r), end = end(r), query.id = queryHits(h), subject.id = subjectHits(h), stringsAsFactors = F);
                                        #        sn.query = as.character(seqnames(query))[h.df$query.id]        
                                        #        sn.subject = as.character(seqnames(subject))[h.df$subject.id]
                  sn.query <- sn1[h.df$query.id]
                  sn.subject <- sn2[h.df$subject.id]

                  if (is.null(by))
                      keep.ix <- sn.query == sn.subject
                  else
                      {
                          by.query <- by1[h.df$query.id]
                          by.subject <- by2[h.df$subject.id]
                          keep.ix <- sn.query == sn.subject & by.query == by.subject
                      }
                  
                  h.df <- h.df[keep.ix, ]
                  h.df$seqnames <- sn.query[keep.ix];
              }

          if (!ignore.strand)
              {
                  h.df$strand <- str.query <- as.character(strand(query)[ix1[h.df$query.id]])
                  str.subject <- as.character(strand(subject)[ix2[h.df$subject.id]])
                  h.df <- h.df[which(str.query == str.subject | str.query == '*' | str.subject == '*'),]            
              }
          else if (nrow(h.df)>0)
              h.df$strand = '*'
      }
    
    if (first)
      h.df = h.df[!duplicated(h.df$query.id), ]

     if (return.type=='GRanges') 
       if (nrow(h.df)>0)           
           {
               if (('strand' %in% names(h.df)))
                   out.gr = GRanges(h.df$seqnames, IRanges(h.df$start, h.df$end),
                       query.id = ix1[h.df$query.id], subject.id = ix2[h.df$subject.id], strand = h.df$strand, seqlengths = seqlengths(query))
               else
                   out.gr = GRanges(h.df$seqnames, IRanges(h.df$start, h.df$end),
                       query.id = ix1[h.df$query.id], subject.id = ix2[h.df$subject.id], seqlengths = seqlengths(query))
                   
               if (!is.null(qcol))
                   values(out.gr) = cbind(values(out.gr), values(query)[out.gr$query.id, qcol, drop = FALSE])

               if (!is.null(scol))
                   values(out.gr) = cbind(values(out.gr), values(subject)[out.gr$subject.id, scol, drop = FALSE])               

               return(out.gr)
           }
       else
         return(GRanges(seqlengths = seqlengths(query)))
     else 
         if (nrow(h.df)>0) {
             
             if (!is.data.table(h.df))
                 h.df = as.data.table(h.df)
             h.df$query.id <- ix1[h.df$query.id]
             h.df$subject.id <- ix2[h.df$subject.id]

             if (!is.null(qcol))
                 h.df = cbind(h.df, as.data.table(as.data.frame(values(query))[h.df$query.id, qcol, drop = FALSE]))
             
             if (!is.null(scol))
                 h.df = cbind(h.df, as.data.table(as.data.frame(values(subject))[h.df$subject.id, scol, drop = FALSE]))
             
             if ('i.start' %in% colnames(h.df))
                 h.df[, i.start := NULL]
                          
             if ('i.end' %in% colnames(h.df))
                 h.df[, i.end := NULL]
             
             return(h.df)
       } else {
         return(data.table())
       }
   }


#' Faster implementatin of GenomicRangs::match 
#'
#' Faster implementation of GRanges match (uses gr.findoverlaps)
#' returns indices of query in subject or NA if none found
#' @param query \code{GRanges} object as query
#' @param subject \code{GRanges} object as subject
#' @param max.slice Default Inf. If query is bigger than this, chunk into smaller on different cores
#' @param verbose Default FALSE
#' @param mc.cores Default 1. Only works if exceeded max.slice
#' @param ... arguments to be passed to \link{gr.findoverlaps}
#' @return returns indices of query in subject or NA if none found
#' @name gr.match
#' @export
gr.match = function(query, subject, max.slice = Inf, verbose = FALSE, mc.cores = 1, ...)
  {
      if (length(query)>max.slice)
          {
              verbose = TRUE
              ix.l = split(1:length(query), ceiling(as.numeric((1:length(query)/max.slice))))
              return(do.call('c', mclapply(ix.l, function(ix) {
                  if (verbose)
                      cat(sprintf('Processing %s to %s\n', min(ix), max(ix)))                
                  gr.match(query[ix, ], subject, verbose = TRUE, ...)
              }, mc.cores = mc.cores)))
          }
      
    tmp = gr.findoverlaps(query, subject, ...)
    tmp = tmp[!duplicated(tmp$query.id)]
    out = rep(NA, length(query))
    out[tmp$query.id] = tmp$subject.id
    return(out)    
   }


#' gr.tile.map
#'
#' Given two tilings of the genome (eg at different resolution)
#' query and subject outputs a length(query) list whose items are integer vectors of indices in subject
#' overlapping that overlap that query (strand non-specific)
#'
#' @note Assumes that input query and subject have no gaps (including at end) or overlaps, i.e. ignores end()
#' coordinates and only uses "starts"
#' @param query Query
#' @param subject Subject
#' @param mc.cores number of cores
#' @param verbose Default FALSE
#' @export
gr.tile.map = function(query, subject, mc.cores = 1, verbose = FALSE)
  {
    ix.q = order(query)
    ix.s = order(subject)
        
    q.chr = as.character(seqnames(query))[ix.q]
    s.chr = as.character(seqnames(subject))[ix.s]

    ql = split(ix.q, q.chr)
    sl = split(ix.s, s.chr)

    tmp = mcmapply(
      function(x,y)
      {
        if (length(y)==0)
          return(NULL)
        all.pos = c(start(query)[x], start(subject)[y])
        is.q = c(rep(T, length(x)), rep(F, length(y)))
        all.ix = c(x, y)
        ord.ix = order(all.pos)
        all.pos = all.pos[ord.ix]
        is.q = is.q[ord.ix]
        all.ix = all.ix[ord.ix]        
        out = matrix(NA, nrow = length(all.pos), ncol = 2)
        last.x = last.y = NA
        for (i in 1:length(all.pos))
          {
            if (verbose)
              if (i %% 100000 == 0) cat('Iteration', i, 'of', length(all.pos), '\n')
            if (is.q[i])
              {
                out[i, ] = c(all.ix[i], last.y)

                if (i<length(all.pos)) ## edge case where subject and query intervals share a start point, leading to two consecutive all.pos
                  if (all.pos[i] == all.pos[i+1])
                    out[i, ] = NA
                      
                last.x = all.ix[i]                
              }            
            else
              {
                out[i, ] = c(last.x, all.ix[i])
                
                if (i<length(all.pos)) ## edge case where subject and query intervals share a start point, leading to two consecutive all.pos
                  if (all.pos[i] == all.pos[i+1])
                    out[i, ] = NA                
                
                last.y = all.ix[i]
              }           
          }
        out = out[rowSums(is.na(out))==0, ]
        return(out)          
      }, ql, sl[names(ql)], mc.cores = mc.cores, SIMPLIFY = FALSE)
    
    m = munlist(tmp)[, -c(1:2), drop = FALSE]
    out = split(m[,2], m[,1])[as.character(1:length(query))]
    names(out) = as.character(1:length(query))
    return(out)
  }

#' @name gr.in
#' @title gr.in
#' @description
#'
#' faster implementation of GRanges %over%  (uses gr.findoverlaps)
#'
#' returns T / F vector if query range i is found in any range in subject
#'
#' by = column name in query and subject that we additionally control for a match (passed on to gr.findoverlaps)
#' @param query \code{GRanges} of query
#' @param subject \code{GRanges} of subject
#' @param by column name in query and subject that we additionally control for a match (passed on to gr.findoverlaps)
#' @param ... params to send to \code{\link{gr.findoverlaps}}
#' @export
gr.in = function(query, subject, by = NULL, pintersect=FALSE,...)
  {
    tmp = gr.findoverlaps(query, subject, by = by, pintersect=pintersect, ...)
    out = rep(FALSE, length(query))
    out[tmp$query.id] = TRUE

    return(out)    
   }

# gr.duplicated
#
# more flexible version of gr.duplicated that allows to restrict duplicates
# using "by" columns and allows in exact matching 
gr.duplicated = function(query, by = NULL, type = 'any')
    {        
        return(duplicated(gr.match(query, query, by = by , type = type)))
    }

#' gr.collapse
#'
#' like "reduce" except only collapses <<adjacent>> ranges in the input
#' returning the collapsed ranges
#' @param gr \code{GRanges} to collapse
#' @param pad padding to place around ranges. Default 1
#' @name gr.collapse
#' @export
gr.collapse = function(gr, pad = 1)
  {
    tmp = gr.findoverlaps(gr + pad, gr + pad, ignore.strand = FALSE)
    m = rep(FALSE, length(gr))
    m[tmp$query.id[tmp$query.id == (tmp$subject.id-1)]] = TRUE
    
    ## will not collapse if two intersecting ranges are in the wrong "order" (ie not increasing (decreasing) on pos (neg) strand
    m[which((strand(gr)[-length(gr)] == '+' & (start(gr)[-length(gr)] > start(gr)[-1])) |
            (strand(gr)[-length(gr)] == '-' & (end(gr)[-length(gr)] < end(gr)[-1])))] = FALSE
    
    m = as(m, 'IRanges')

    if (length(m)>0)
      {
        end(m) = end(m)+1
        tmp = cbind(start(gr)[start(m)], end(gr)[start(m)], start(gr)[end(m)], end(gr)[end(m)])
        s = pmin(start(gr)[start(m)], end(gr)[start(m)], start(gr)[end(m)], end(gr)[end(m)])
        e = pmax(start(gr)[start(m)], end(gr)[start(m)], start(gr)[end(m)], end(gr)[end(m)])
        return(GRanges(seqnames(gr)[start(m)], IRanges(s, e), strand = strand(gr)[start(m)], seqlengths = seqlengths(gr)))
      }
    else
      return(gr[c()])                
  }

#' gr.val
#'
#' annotates gr's in "query" with aggregated values of gr's in "target" in field "val"
#'
#' If "val" is numeric: given "target" gr with value column "val" representing ranged data (ie segment intensities) computes the value
#' in each "query" gr as the weighted mean of its intersection with target (ie the target values weighted by the width 
#' of the intersections).  
#'
#' applications include querying the average value of target across a given query interval (eg exon to gene pileup)
#' or recasting a high res tiling in terms of low res intervals.  Usually query intervals are bigger than the target intervals.
#'
#' If "val" is a character field: then aggregation will paste together the (unique), verlapping values, collapsing by comma
#' 
#' returns query with the "val" field populated
#'
#' Optional
#'
#' query and target can be GRangesLists's, in which case val will refer to GRangesList level values fields
#' @name gr.val
#' @param query GRanges of query ranges whose "val" column we will populate with aggregated values of target
#' @param target GRanges of target ranges that already have "val" column populated
#' @param val dummy
#' @param mean scalar logical flag if FALSE then will return sum instead of mean, only applies if target "val" column i snumeric
#' @param weighted Default \code{mean}
#' @param na.rm Default FALSE
#' @param by.prefix Default \code{val}
#' @param merge if merge = FALSE then will cross every range in query with every level of "by" in target (and create data matrix), otherwise will assume query has "by" and merge only ranges that have matching "by" values in both query and target
#' @param verbose Default FALSE
#' @param FUN takes two  arguments (value, na.rm = TRUE) if weighted = FALSE, and three (value, width, na.rm = TRUE) if weighted = TRUE
#' @param ignore.strand Default TRUE
#' @param default.val dummy
#' @param max.slice Default Inf
#' @param mc.cores Number of cores (only if query exceed \code{max.slice})
#' @param ... params to be passed to gr.findoverlaps
#' @param sep scalar character, specifies character to use as separator when aggregating character "vals" from target, only applies if target is numeric
#' @param by scalar character, specifies additional "by" column of query AND target that will be used to match up query and target pairs (i.e. in addition to pure GRanges overlap), default is NULL
#' @export
gr.val = function(query, target, val = NULL,
    mean = TRUE, # if false then will return (weighted) <sum> instead of <mean>, only applies if target is numeric
    weighted = mean, # if false will return unweighted sum / mean
    na.rm = F, # only applies if val column of target is numeric
    by = NULL,
    by.prefix = val,
    merge = FALSE, # if merge = FALSE then will cross every range in query with every level of "by" in target (and create data matrix), otherwise will assume query has "by" and merge only ranges that have matching "by" values in both query and target
    verbose = FALSE,
    FUN = NULL, ## takes two  arguments (value, na.rm = TRUE) if weighted = FALSE, and three (value, width, na.rm = TRUE) if weighted = TRUE
    ignore.strand = TRUE,
    default.val = NA, 
    max.slice = Inf, ## max row slice to process at a time
    mc.cores = 1,  ## is slicing how many cores to split across
  ..., 
  sep = ', ' # only applies if val column of target is character  
  )
  {
      if (is.null(val))
          val = 'value'

      if (!(all(ix <- val %in% names(values(target)))))
          values(target)[, val[!ix]] = 1
            
      if (length(query)>max.slice)
          {
              verbose = TRUE
              ix.l = split(1:length(query), ceiling(as.numeric((1:length(query)/max.slice))))
              return(do.call('grbind', mclapply(ix.l, function(ix) {
                  if (verbose)
                      cat(sprintf('Processing %s to %s of %s\n', min(ix), max(ix), length(query)))               
                  gr.val(query[ix, ], target = target, val= val, mean = mean, weighted = weighted, na.rm = na.rm, verbose = TRUE, by = by, FUN = FUN, merge = merge, ignore.strand = ignore.strand, ...)
              }, mc.cores = mc.cores)))
          }
      
    if (inherits(target, 'GRangesList'))
      {
        target.was.grl = T;
        target.grl.id = as.data.frame(target)$element        
        val.vec = lapply(val, function(x) values(target)[, x])
        target = unlist(target)
        val.vec = lapply(val.vec, function(X) val.vec[target.grl.id])
      }
    else
      {
        if (!is.null(val))
          val.vec = lapply(val, function(x) values(target)[, x])
        else
          val.vec = list(rep(1, length(target)))
          
        target.grl.id = 1:length(target);
        target.was.grl = F;
      }

    if (inherits(query, 'GRangesList'))
      {
        query.was.grl = T;
        query.grl.id = rep(1:length(query), elementLengths(query))
        query = unlist(query)
      }
    else
      {
        query.grl.id = 1:length(query);
        query.was.grl = F;
      }

      if (!is.null(FUN))
          {
              args = names(formals(FUN))[1:3]
         
              if (!is.null(args))
                  {
                      if (any(is.na(args)))
                          args[is.na(args)] = ''
                      
                      if (weighted)
                          {
                              if (any(!(c('x', 'w', 'na.rm') %in% args)))
                                  warning('FUN input must be function with three arguments: "x" = value, "w" = interval width, "na.rm" = na.rm flag')
                          }
                      else
                          {
                              if (any(!(c('x', 'na.rm') %in% args)))
                                  warning('FUN input must be function with two arguments: "x" = value, "na.rm" = na.rm flag')
                          }
                  }
          }
      
    if (!merge)
        hits = gr.findoverlaps(query, target, scol = by, ignore.strand = ignore.strand, verbose = verbose, return.type = 'data.table', ...)
    else
        hits = gr.findoverlaps(query, target, by = by, ignore.strand = ignore.strand, verbose = verbose, return.type = 'data.table', ...)
      
      if (verbose)
          cat(sprintf('aggregating hits\n'))

    vals = val
    val.vecs = val.vec
      
    for (vix in 1:length(vals))
        {
            val = vals[[vix]]
            val.vec = val.vecs[[vix]]
            is.char = is.character(values(target)[, val]);
            if (is.null(by) | merge == TRUE)
                {
                    if (nrow(hits)>0)
                        {
                            values(query)[, val] = NA;
                            hits[, width := as.numeric(end - start)+1]
                            setkey(hits,  query.id)                    
                            if (is.char)
                                {
                                    values(query)[, val] = '';
                                    hits$id = 1:nrow(hits);
                                    tmp = hits[, list(val = paste(setdiff(val.vec[subject.id], NA), collapse = sep)), by = query.id]
                                    if (!is.na(default.val))
                                        tmp[is.na(tmp)] = default.val
                                    values(query)[tmp[,query.id], val] = tmp[,val]
                                }
                            else
                                {           

                                        #            val.vec = as.numeric(values(target)[, val]);
                                    val.vec = as.numeric(val.vec)
                                    if (weighted)
                                        {
                                            if (!is.null(FUN))
                                                tmp = hits[, list(val = do.call(FUN, list(val.vec[subject.id], width, na.rm = na.rm))), by = query.id]
                                            else if (mean)
                                                tmp = hits[, list(val = sum(width * val.vec[subject.id], na.rm = na.rm)/sum(width)), by = query.id]
                                            else
                                                tmp = hits[, list(val = sum(width * val.vec[subject.id], na.rm = na.rm)), by = query.id]
                                        }
                                    else
                                        {
                                            if (!is.null(FUN))
                                                tmp = hits[, list(val = do.call(FUN, list(val.vec[subject.id], na.rm = na.rm))), by = query.id]
                                            else if (mean)
                                                tmp = hits[, list(val =  mean(val.vec[subject.id], na.rm = na.rm)), by = query.id]
                                            else
                                                tmp = hits[, list(val = sum(val.vec[subject.id], na.rm = na.rm)), by = query.id]
                                        }

                                    if (!is.na(default.val))
                                        tmp[is.na(tmp)] = default.val
                                    
                                    values(query)[tmp[,query.id], val] = tmp[,val]
                                }
                        }
                }
            else ## by is not null
                {
                    
                    if (!is.null(by.prefix))                  
                        if (is.na(by.prefix))
                            by.prefix = NULL
                        else if (nchar(by.prefix)==0)
                            by.prefix = NULL
                    
                    if (nrow(hits)>0)
                        {
                            hits[, width := as.numeric(end - start)+1]
                            if (is.char)
                                {
                                    hits$id = 1:nrow(hits);
                                    tmp = hits[, list(val = paste(setdiff(val.vec[subject.id], NA), collapse = sep)), keyby = list(query.id, bykey = eval(parse(text=by)))]
                                    
                                    tmp2 = dcast.data.table(tmp, query.id ~ bykey, value.var = 'val')
                                    setkey(tmp2, query.id)
                                    new.df = as.data.frame(tmp2[list(1:length(query)), ])[ ,-1]
                                    
                                    if (!is.na(default.val))
                                        new.df[is.na(new.df)] = default.val
                                    
                                    if (!is.null(by.prefix))
                                        colnames(new.df) =  paste(by.prefix, names(tmp2)[-1], sep = '.')
                                    else
                                        colnames(new.df) =  names(tmp2)[-1]
                                    new.names = c(colnames(values(query)), colnames(new.df))
                                    values(query) = cbind(values(query), new.df)
                                    colnames(values(query)) = new.names
                                }
                            else
                                {           
                                    
                                        #            val.vec = as.numeric(values(target)[, val]);
                                    val.vec = as.numeric(val.vec)
                                    if (weighted)
                                        {
                                            if (!is.null(FUN))
                                                {
                                                    tmp = hits[, list(val = do.call(FUN, list(val.vec[subject.id], width, na.rm = na.rm))), keyby = list(query.id, bykey = eval(parse(text=by)))]
                                                }
                                            else if (mean)
                                                tmp = hits[, list(val = sum(width * val.vec[subject.id], na.rm = na.rm)/sum(width)), keyby = list(query.id, bykey = eval(parse(text=by)))]
                                            else
                                                tmp = hits[, list(val = sum(width * val.vec[subject.id], na.rm = na.rm)), keyby = list(query.id, bykey = eval(parse(text=by)))]
                                        }
                                    else
                                        {
                                            if (!is.null(FUN))
                                                tmp = hits[, list(val = do.call(FUN, list(val.vec[subject.id], na.rm = na.rm))), keyby = list(query.id, bykey = eval(parse(text=by)))]
                                            else if (mean)
                                                tmp = hits[, list(val =  mean(val.vec[subject.id], na.rm = na.rm)), keyby = list(query.id, bykey = eval(parse(text=by)))]
                                            else
                                                tmp = hits[, list(val = sum(val.vec[subject.id], na.rm = na.rm)), keyby = list(query.id, bykey = eval(parse(text=by)))]
                                        }                            
                                    
                                    tmp2 = dcast.data.table(tmp, query.id ~ bykey, value.var = 'val')
                                    setkey(tmp2, query.id)
                                    new.df = as.data.frame(tmp2[list(1:length(query)), ])[ ,-1, drop = FALSE]

                                    if (!is.na(default.val))
                                        new.df[is.na(new.df)] = default.val
                                    
                                    if (!is.null(by.prefix))
                                        colnames(new.df) =  paste(by.prefix, names(tmp2)[-1], sep = '.')
                                    else
                                        colnames(new.df) =  names(tmp2)[-1]
                                    
                                    new.names = c(colnames(values(query)), colnames(new.df))
                                    values(query) = cbind(values(query), new.df)
                                    colnames(values(query)) = new.names                            
                                }
                        }        
                }
        }

    if (query.was.grl)
      query = split(query, query.grl.id)    

    return(query)
  }

#' gr.dist
#'
#' computes pairwise distance between elements of two gr objects of length n and m
#' returning n by m matrix of distances between item i of gr1 and item j of gr2
#'
#' distances are computed as follows:
#' NA for ranges on different seqnames
#' 0 for overlapping ranges)
#' min(abs(end1-end2), abs(end1-start2), abs(start1-end2), abs(start1-end1),) for all others
#'
#' if only gr1 is provided, then will return n x n matrix of gr's vs themselves
#'
#' if max.dist = T then will replace min with max above
#' @name gr.dist
#' @export
gr.dist = function(gr1, gr2 = NULL, ignore.strand = T, ...)
{
  if (is.null(gr2))
    gr2 = gr1;

  if (ignore.strand)
    {
      strand(gr1) = '*'
      strand(gr2) = '*'
    }
  
  ix1 = rep(1:length(gr1), length(gr2))
  ix2 = rep(1:length(gr2), each = length(gr1))

  out = matrix(suppressWarnings(distance(gr1[ix1], gr2[ix2], ...)), nrow = length(gr1), ncol = length(gr2))

  return(out)
}

#' Provide transparency to colors
#'
#' takes provided colors and gives them the specified alpha (ie transparency) value
#'
#' @param col color string
#' @param alpha alpha value between 0 and 1
#' @return rgb color like the input, but with transparency added
#' @export
alpha = function(col, alpha)
  {    
    col.rgb = col2rgb(col)
    return(rgb(red = col.rgb['red', ]/255, green = col.rgb['green', ]/255, blue = col.rgb['blue', ]/255, alpha = alpha))
  }

#' Filters GRangesList to only include ranges in the specified window
#' 
#' (this is different from %in% which does not remove non matching ranges from the grls)
#'
#' does not return list in necessarily same order
# @param grl \link{GRangesList} to filter
# @param windows \link{GRanges} windows to keep
#' @name grl.filter
#' @export
grl.filter = function(grl, windows)
  {
    tmp = as.data.frame(grl);
    tmp = tmp[seg.on.seg(tmp, windows), ]
    FORBIDDEN = c('seqnames', 'start', 'end', 'strand', 'ranges', 'seqlevels', 'seqlengths', 'isCircular', 'genome', 'width', 'element');
    gr.metadata = tmp[, setdiff(colnames(tmp), FORBIDDEN)];

    if (!is.null(dim(gr.metadata)))
      out.grl = split(GRanges(tmp$seqnames, IRanges(tmp$start, tmp$end), seqlengths = seqlengths(grl), gr.metadata,
        strand = tmp$strand), tmp$element)
    else
      out.grl = split(GRanges(tmp$seqnames, IRanges(tmp$start, tmp$end), seqlengths = seqlengths(grl),
        strand = tmp$strand), tmp$element);

    values(out.grl) = values(grl)[match(names(out.grl), names(grl)), ]
    return(out.grl);
  }

#' grl.allin
#'
#' Like %in% for grl but now will return a logical vector that is true at position if i
#' only if the ranges in grl[i] intersect <<all>>, <<some>>, <<only>>  windows in the subject
#'
#' eg can use to identify read pairs whose ends are contained inside two genes)
#' @param grl \link{GRangesList} object to query
#' @param windows \link{GRanges} windows to check against
#' @param some [Default FALSE]
#' @param only [Default FALSE]
#' @name grl.allin
#' @export
grl.in = grl.allin = function(grl, windows, some = FALSE, only = FALSE)
  {
    if (length(grl)==0)
      return(logical())

    if (length(windows)==0)
      return(rep(T, length(grl)))
                 
    numwin = length(windows);    
    gr = grl.unlist(grl)
    m = as.data.frame(gr.findoverlaps(gr, windows))

    out = rep(FALSE, length(grl))
    if (nrow(m)==0)
      return(out)
    m$grl.id = gr$grl.ix[m$query.id]

    if (some)
      tmp = aggregate(formula = subject.id ~ grl.id, data = m, FUN = function(x) length(intersect(1:numwin, x))>0)
    else if (only)
      return(mapply(function(x, y) length(setdiff(x, y))==0,
                    split(1:length(gr), factor(gr$grl.ix, 1:length(grl))),
                    split(m$query.id, factor(m$grl.id, 1:length(grl)))))   
    else
      tmp = aggregate(formula = subject.id ~ grl.id, data = m, FUN = function(x) length(setdiff(1:numwin, x))==0)
    
    out = rep(FALSE, length(grl))
    out[tmp[,1]] = tmp[,2]
    return(out)    
  }

#' grl.split
#'
#' splits GRL's with respect to their seqnames and strand (default), returning
#' new grl whose items only contain ranges with a single seqname / strand
#'
#' can also split by arbitrary (specified) genomic ranges value fields
#' @param grl \code{GRangesList} to split
#' @param seqname Default TRUE
#' @param strand Default TRUE
#' @param values columns of values field in grl
#' @name grl.split
#' @export
grl.split = function(grl, seqname = TRUE, strand = TRUE,
  values = c() # columns of values field in grl
  )
  {
    ele = tryCatch(as.data.frame(grl)$element, error = function(e) e)    
    if (inherits(ele, 'error'))
      {
        if (is.null(names(grl)))
          nm = 1:length(names(grl))
        else
          nm = names(grl)
        
        ele = unlist(lapply(1:length(grl), function(x) rep(nm[x], length(grl[[x]]))))
      }
    
    gr = unlist(grl)
    names(gr) = NULL;

    by = ele;
    if (seqname)
      by = paste(by, seqnames(gr))

    if (strand)
      by = paste(by, strand(gr))

    values = intersect(names(values(gr)), values);
    if (length(values)>0)
      for (val in values)
        by = paste(by, values(gr)[, val])

    out = split(gr, by);
    names(out) = ele[!duplicated(by)]

    values(out) = values(grl[ele[!duplicated(by)]])

    return(out)
  }

#' grl.stripnames
#'
#' get rid of gr names inside a grl
#' @param grl \link{GRangesList} to string names from
#' @name grl.stripnames
#' @export
grl.stripnames = function(grl)
  {
    ele = tryCatch(as.data.frame(grl)$element, error = function(e) e)    
    if (inherits(ele, 'error'))
      ele = unlist(lapply(1:length(grl), function(x) rep(x, length(grl[[x]]))))
      
    gr = unlist(grl);
    names(gr) = NULL;

    out = split(gr, ele);
    values(out) = values(grl)
    names(out) = names(grl)
    
    return(out)
  }


#' grl.unlist
#'
#' Does a "nice" unlist of a grl object adding a field "grl.ix" denoting which element of the grl
#' each gr corresponds to
#'
#' and a field grl.iix which saves the (local) index that that gr was in its corresponding grl item
#' @param grl GenomicRangesList to unlist
#' @name grl.unlist
#' @export
grl.unlist = function(grl)
  {
    if (length(grl) == 0) ## JEREMIAH
      return(GRanges())
#      return(grl) 
    names(grl) = NULL

    as.df = as.data.frame(grl)
    
    el = as.df$element
    if (is.null(el))
        el = as.df$group
       
    out = unlist(grl)
    out$grl.ix = el
    tmp = rle(el)
    out$grl.iix = unlist(sapply(tmp$lengths, function(x) 1:x))
    values(out) = cbind(values(grl)[out$grl.ix, , drop = FALSE], values(out))
    return(out)
  }

#' grl.span
#'
#' Returns GRanges object representing the left / right extent of each GRL item.  In case of "chimeric" GRL items (ie that map
#' to two chromosomes) there are two options:
#' (1) specify "chr" chromosome as argument to subset GRL's that are on that chromosome, and compute GRL extents from this, any GRL
#'     full outside of that chromosome will get a 0 width GRL 
#' (2) (default) allow chimeric GRL items to get an extent that is with respect to the first chromosome in that GRL 
#'
#' If a grl item contains ranges that lie on different chromosomes, then corresponding grange will have chromosome "NA" and IRange(0, 0)
#' @param grl \link{GRangesList} to query
#' @param chr [Default NULL]
#' @param ir [Default FALSE]
#' @param keep.strand [Default TRUE]
#' @name grl.span
#' @export
grl.span = function(grl, chr = NULL, ir = FALSE, keep.strand = TRUE)
  {
    if (is.null(names(grl)))
      names(grl) = 1:length(grl);

    tmp = tryCatch(as.data.frame(grl), error = function(e) e)
    
    if (inherits(tmp, 'error')) ## gr names are screwy so do some gymnastics    
      {
        if (is.null(names(grl)))
          names.grl = 1:length(grl)
        else
          names.grl = names(grl);
        
        element = as.character(Rle(names.grl, sapply(grl, length)))
        tmp.gr = unlist(grl)
        names(tmp.gr) = NULL;
        tmp = as.data.frame(tmp.gr);
        tmp$element = element;        
      }
    
    if (is.null(chr))
      {
        chrmap = aggregate(formula = seqnames ~ element, data = tmp, FUN = function(x) x[1]);
        chrmap = structure(as.character(chrmap[,2]), names = chrmap[,1])

        if (keep.strand)
          {
            strmap = aggregate(formula = as.character(strand) ~ element, data = tmp, FUN =
              function(x) {y = unique(x); if (length(y)>1) return('*') else y[1]})
            strmap = structure(as.character(strmap[,2]), names = strmap[,1])
            str = strmap[names(grl)]; 
          }
        else
          str = '*'
        
        tmp = tmp[tmp$seqnames == chrmap[tmp$element], ]; ## remove all gr from each GRL item that don't map to the chr of the first gr
        chr = chrmap[names(grl)];
        out.gr = GRanges(chr, IRanges(1,0), seqlengths = seqlengths(grl), strand = str)
      }
    else 
      {
        if (length(chr)>1)
          warning('chr has length greater than 1, only the first element will be used')
        tmp = tmp[tmp$seqnames == chr[1], ]
        out.gr = rep(GRanges(chr, IRanges(1, 0)), length(grl)) # missing values
      }

    if (nrow(tmp)>0)
      {
        tmp = split(GRanges(tmp$seqnames, IRanges(tmp$start, tmp$end)), tmp$element)
        out.gr[match(names(tmp), names(grl))] = GRanges(chr[names(tmp)],
                IRanges(sapply(start(tmp), min), sapply(end(tmp), max)), strand = strand(out.gr)[match(names(tmp), names(grl))]);
        names(out.gr) = names(grl)
      }
    return(out.gr)
  }


#' grl.pivot
#'
#' "pivots" grl object "x" by returning a new grl "y" whose
#' kth item is gr object of ranges x[[i]][k] for all i in 1:length(x)
#'
#' Assumes all grs in "x" are of equal length
#' @param x GenomicRangesList object to pivot
#' @name grl.pivot
#' @export
grl.pivot = function(x)
  {
    if (length(x) == 0)
      return(GRangesList(GRanges(seqlengths = seqlengths(x)), GRanges(seqlengths = seqlengths(x))))
    return(split(unlist(x), rep(1:length(x[[1]]), length(x))))
  }


#' get.var.col
#'
#' simple function storing default
#' variant color scheme
#' @name get.var.col
get.varcol = function()
  {
    VAR.COL = c('XA' = 'green', 'XG' = 'brown', 'XC' = 'blue', 'XT' = 'red', 'D'= alpha('lightblue', 0.4),
    'I'= 'purple', 'N' = alpha('white', 0.8), 'S' = alpha('pink', 0.9))
    return(VAR.COL)
  }

##
##
## $$$$$$$\                                   $$\                        $$\                 $$\   $$\   $$\     $$\ $$\ 
## $$  __$$\                                  $$ |                       $$ |                $$ |  $$ |  $$ |    \__|$$ |
## $$ |  $$ | $$$$$$$\ $$$$$$\  $$$$$$\$$$$\$$$$$$\   $$$$$$\   $$$$$$\  $$ | $$$$$$$\       $$ |  $$ |$$$$$$\   $$\ $$ |
## $$$$$$$  |$$  _____|\____$$\ $$  _$$  _$$\_$$  _| $$  __$$\ $$  __$$\ $$ |$$  _____|      $$ |  $$ |\_$$  _|  $$ |$$ |
## $$  __$$< \$$$$$$\  $$$$$$$ |$$ / $$ / $$ |$$ |   $$ /  $$ |$$ /  $$ |$$ |\$$$$$$\        $$ |  $$ |  $$ |    $$ |$$ |
## $$ |  $$ | \____$$\$$  __$$ |$$ | $$ | $$ |$$ |$$\$$ |  $$ |$$ |  $$ |$$ | \____$$\       $$ |  $$ |  $$ |$$\ $$ |$$ |
## $$ |  $$ |$$$$$$$  \$$$$$$$ |$$ | $$ | $$ |\$$$$  \$$$$$$  |\$$$$$$  |$$ |$$$$$$$  |      \$$$$$$  |  \$$$$  |$$ |$$ |
## \__|  \__|\_______/ \_______|\__| \__| \__| \____/ \______/  \______/ \__|\_______/        \______/    \____/ \__|\__|
##
##
## Rsamtools util
##
## wrapper functions around Rsamtools and rtracklayer
## to help extract read info from bam file, mutations pileups, and coverage from wig / bigwig files
##
##
                                                                                                                      
#' Read BAM file into GRanges or data.table
#'
#' Wrapper around Rsamtools bam scanning functions,
#' by default, returns GRangesList of read pairs for which <at least one> read lies in the supplied interval
#' @param bam Input bam file. Advisable to make "bam" a BamFile instance instead of a plain string, so that the index does not have to be reloaded.
#' @param bami Input bam index file.
#' @param gr GRanges of intervals to retrieve
#' @param intervals GRanges of intervals to retrieve
#' @param stripstrand Flag to ignore strand information on the query intervals. Default TRUE
#' @param what What fields to pull down from BAM. Default \code{scanBamWhat()}
#' @param unpack.flag Add features corresponding to read flags. Default FALSE
#' @param verbose Increase verbosity
#' @param tag Additional tags to pull down from the BAM (e.g. 'R2')
#' @param isPaired See documentation for \code{scanBamFlag}. Default NA
#' @param isProperPair See documentation for \code{scanBamFlag}. Default NA
#' @param isUnmappedQuery See documentation for \code{scanBamFlag}. Default NA
#' @param hasUnmappedMate See documentation for \code{scanBamFlag}. Default NA
#' @param isNotPassingQualityControls See documentation for \code{scanBamFlag}. Default NA
#' @param isDuplicate See documentation for \code{scanBamFlag}. Default FALSE
#' @param isValidVendorRead See documentation for \code{scanBamFlag}. Default TRUE
#' @param as.grl Return reads as GRangesList. Controls whether \code{get.pairs.grl} does split. Default TRUE
#' @param as.data.table Return reads in the form of a data.table rather than GRanges/GRangesList
#' @param ignore.indels messes with cigar to read BAM with indels removed. Useful for breakpoint mapping on contigs
#' @param size.limit Default 1e6
#' @param ... passed to \code{scanBamFlag}
#' @return Reads in one of GRanges, GRangesList or data.table
#' @import Rsamtools
#' @import data.table
#' @export
read.bam = function(bam, intervals = NULL,## GRanges of intervals to retrieve
    gr = intervals,
    all = FALSE, 
    bami = NULL,  
  pairs.grl = TRUE, # if TRUE will return GRangesList of read pairs for whom at least one read falls in the supplied interval
#  paired = F, # if TRUE, will used read bam gapped alignment pairs warning: will throw out pairs outside of supplied window
#  gappedAlignment = T, # if false just read alignments using scanbam
  stripstrand = TRUE, 
  what = scanBamWhat(),
  unpack.flag = FALSE, # will add features corresponding to read flags
  verbose = FALSE,
  tag = NULL,
  isPaired = NA, ## if these features are NA, then reads satisfying both T and F will be returned
  isProperPair = NA, 
  isUnmappedQuery = NA,
  hasUnmappedMate = NA,
  isNotPassingQualityControls = NA,
  isDuplicate = F,
  isValidVendorRead = TRUE,
  as.grl=TRUE, ## return pairs as grl, rather than GRanges .. controls whether get.pairs.grl does split (t/c rename to pairs.grl.split)
  as.data.table=FALSE, ## returns reads in the form of a data table rather than GRanges/GRangesList
  ignore.indels=FALSE, ## messes with cigar to read BAM with indels removed. Useful for breakpoint mapping on contigs
  size.limit = 1e6,
  ... # passed to scanBamFlag (
  )
{
  if (!inherits(bam, 'BamFile'))
    {
      if (is.null(bami))
        {
          if (file.exists(bai <- gsub('.bam$', '.bai', bam)))
            bam = BamFile(bam, bai)
          else if (file.exists(bai <- paste(bam, '.bai', sep = '')))
            bam = BamFile(bam, bai)
          else
            bam = BamFile(bam)
        }
      else
        bam = BamFile(bam, index = bami)
    }
      

  # if intervals unspecified will try to pull down entire bam file (CAREFUL)

  if (length(intervals)==0)
      intervals = NULL

  if (is.null(intervals))
      intervals = gr
  
  if (is.null(intervals))
      {
          if (all)              
              intervals = seqinfo2gr(seqinfo(bam))
          else
              stop('Must provide non empty interval list')
      }
      
  if (class(intervals) == 'data.frame')
    intervals = seg2gr(intervals);

  if (inherits(intervals, 'GRangesList'))
    intervals = unlist(intervals);
  
  if (stripstrand)
    strand(intervals) = '*'

  intervals = reduce(intervals);

  now = Sys.time();

  if (pairs.grl)
    paired = F
  
  flag = scanBamFlag(isPaired = isPaired, isProperPair = isProperPair, isUnmappedQuery = isUnmappedQuery,
    hasUnmappedMate = hasUnmappedMate, isNotPassingQualityControls = isNotPassingQualityControls,
    isDuplicate = isDuplicate, ...)

  tag = unique(c('MD', 'MQ', tag))
  param = ScanBamParam(which = gr.fix(intervals, bam, drop = T), what = what, flag = flag, tag = tag)

  if (verbose)
      cat('Reading bam file\n')
   if (class(bam) == 'BamFile')
     out <- scanBam(bam, param=param)
   else      
     out <- scanBam(bam, index=bami, param=param)
   if (verbose) {
     print(Sys.time() - now)
     print('BAM read. Making into data.frame')
   }

   out <- out[sapply(out, function(x) length(x$qname)>0)]

   if (length(out)>0)
     {
       if (verbose) {
         print(Sys.time() - now)
         print('combining lists')
       }
       out <- as.data.frame(rbindlist(lapply(out, function(x) 
         {
           x <- c(x[-match('tag', names(x))], x$tag)
           
           x <- x[sapply(x, length)>0]
           conv <- which(!(sapply(x, class) %in% c('integer', 'numeric', 'character')))
           x[conv] <- lapply(x[conv], as.character)

           for (t in tag)
             if (!(t %in% names(x)))
               x[[t]] = rep(NA, length(x$qname))

           if (!('R2' %in% names(x)) && 'R2' %in% tag)
             x$R2 <- rep(NA, length(x$qname))
           if (!('Q2' %in% names(x)) && 'Q2' %in% tag)
             x$Q2 <- rep(NA, length(x$qname))
           x
        })))

       ## faster CIGAR string parsing with vectorization and data tables
       if (verbose) {
         print(Sys.time() - now)
         print('filling pos2 from cigar')
       }
       if (ignore.indels) {
         cigar <- gsub('[0-9]+D', '', gsub('([0-9]+)I', '\\1M', out$cigar))  ## Remove deletions, turn insertions to matches
         cig <- splitCigar(cigar)
         torun=sapply(cig, function(y) any(duplicated((y[[1]][y[[1]]==M]))))
         M <- charToRaw('M')
         new.cigar <- sapply(cig[torun], function(y) {
                 lets <- y[[1]][!duplicated(y[[1]])]
                 vals <- y[[2]][!duplicated(y[[1]])]
                 vals[lets==M] <- sum(y[[2]][y[[1]]==M])
                 lets <- strsplit(rawToChar(lets), '')[[1]]
                 paste(as.vector(t(matrix(c(vals, lets), nrow=length(vals), ncol=length(lets)))), collapse='')
                 })
         out$cigar[torun] <- new.cigar
       }
       cigs <- countCigar(out$cigar)
       out$pos2 <- out$pos + cigs[, "M"]

       if (verbose) {
         print(Sys.time() - now)
         print('fixing seqdata')
       }
       out$qwidth = nchar(out$seq)
       unm = is.na(out$pos)          
       if (any(unm))
         {
           out$pos[unm] = 1
           out$pos2[unm] = 0
           out$strand[unm] = '*'
         }                    
       gr.fields = c('rname', 'strand', 'pos', 'pos2');
       vals = out[, setdiff(names(out), gr.fields)]

       if (!as.data.table) {
         out <- GRanges(out$rname, IRanges(out$pos, pmax(0, out$pos2-1)), strand = out$strand, seqlengths = seqlengths(intervals))
         values(out) <- vals;
       } else {
         out <- data.table(seqnames=out$rname, start=out$pos, end= pmax(out$pos2-1, 0), strand=out$strand)
         val <- data.table(vals)
         out <- cbind(out, val)
       }
       #out$uname = paste(out$qname, ifelse(bamflag(out$flag)[, 'isFirstMateRead'], '_r1', '_r2'), sep = '')
     }
   else {
     if (!as.data.table)
       return(GRanges(seqlengths = seqlengths(intervals)))
     else
       return(data.table())
   }

  if (verbose)
    {
      if (as.data.table)
        cat(sprintf('Extracted %s reads\n', nrow(out)))
      else
        cat(sprintf('Extracted %s reads\n', length(out)))        
      print(Sys.time() - now)
    }

  if (pairs.grl)
    {
      if (verbose)
        cat('Pairing reads\n')
      out <- get.pairs.grl(out, as.grl=as.grl)
      if (verbose)
        {
          cat('done\n')
          print(Sys.time() - now)
        }
      if (as.grl && !as.data.table) {
        names(out) = NULL;
        values(out)$col = 'gray';
        values(out)$border = 'gray';
      } 
    }

    return(out)
}

#' Quick way to get tiled coverage via piping to samtools (~10 CPU-hours for 100bp tiles, 5e8 read pairs)
#'
#' Gets coverage for window size "window", pulling "chunksize" records at a time and incrementing bin
#' corresponding to midpoint or overlaps of corresponding (proper pair) fragment (uses TLEN and POS for positive strand reads that are part of a proper pair)
#'
#' @param bam.file character scalar input bam file
#' @param window integer scalar window size (in bp)
#' @param chunksize integer scalar, size of window
#' @param min.mapq integer scalar, minimim map quality reads to consider for counts
#' @param verbose dummy
#' @param max.tlen max paired-read insert size to consider
#' @param st.flag samtools flag to filter reads on [Default: -f 0x02 -F 0x10]
#' @param fragments dummy
#' @param region dummy
#' @param do.gc dummy
#' @param midpoint if TRUE will only use the fragment midpoint, if FALSE will count all bins that overlap the fragment
#' @return GRanges of "window" bp tiles across seqlengths of bam.file with meta data field $counts specifying fragment counts centered
#' in the given bin.
#' @import Rsamtools
#' @import data.table
#' @export
bam.cov = function(bam.file, window = 1e2, chunksize = 1e5, min.mapq = 30, verbose = TRUE,
    max.tlen = 1e4, ## max insert size to consider
    st.flag = "-f 0x02 -F 0x10",
    fragments = TRUE,
    region = NULL, 
    do.gc = FALSE, 
    midpoint = TRUE ## if TRUE will only use the fragment midpoint, if FALSE will count all bins that overlap the fragment
    )
  {
    cmd = 'samtools view %s %s -q %s | cut -f "3,4,9"' ## cmd line to grab the rname, pos, and tlen columns
    
    sl = seqlengths(BamFile(bam.file))
    
    counts = lapply(sl, function(x) rep(0, ceiling(x/window)))
    numwin = sum(sapply(sl, function(x) ceiling(x/window)))

    if (!is.null(region))
        {
            cat(sprintf('Limiting to region %s\n', region))
            cmd = 'samtools view %s %s -q %s %s | cut -f "3,4,9"' ## cmd line to grab the rname, pos, and tlen columns
            if (!file.exists(paste(bam.file, '.bam', sep = '')))
                if (file.exists(bai.file <- gsub('.bam$', '.bai', bam.file)))                    
                    {                
                        .TMP.DIR = '~/temp/.samtools'
                        system(paste('mkdir -p', TMP.DIR))
                        tmp.fn = paste(normalizePath(TMP.DIR), '/tmp', runif(1), sep = '')
                        system(sprintf('ln -s %s %s.bam', bam.file, tmp.fn))
                        system(sprintf('ln -s %s %s.bam.bai', bai.file, tmp.fn))
                    }

            cat('Calling', sprintf(cmd, st.flag, paste(tmp.fn, 'bam', sep = '.'), min.mapq, region), '\n')
            p = pipe(sprintf(cmd, st.flag, paste(tmp.fn, 'bam', sep = '.'), min.mapq, region), open = 'r')
        }
    else
        {
            cat('Calling', sprintf(cmd, st.flag, bam.file, min.mapq), '\n')
            p = pipe(sprintf(cmd, st.flag, bam.file, min.mapq), open = 'r')
        }
    
    i = 0    
    sl.dt = data.table(chr = names(sl), len = sl)
    counts = sl.dt[, list(start = seq(1, len, window)), by = chr]
    counts = counts[, bin := 1:length(start), by = chr]
    counts[, end := pmin(start + window-1, sl[chr])]
    counts[, count := 0]
    counts[, rowid := 1:length(count)]
    setkeyv(counts, c("chr", "bin")) ## now we can quickly populate the right entries
    totreads = 0

    st = Sys.time()
    if (verbose)
        cat('Starting fragment count on', bam.file, 'with bin size', window, 'and min mapQ', min.mapq, 'and insert size limit', max.tlen, 'with midpoint set to', midpoint, '\n')

    while (length(chunk <- readLines(p, n = chunksize))>0)
      {
          i = i+1

          if (fragments)
              {
                  chunk = fread(paste(chunk, collapse = "\n"), header = F)[abs(V3)<=max.tlen, ]
                  if (midpoint) ## only take midpionts              
                      chunk[, bin := 1 + floor((V2 + V3/2)/window)] ## use midpoint of template to index the correct bin
                  else ## enumerate all bins containing fragment i.e. where fragments overlap multiple bins  (slightly slower)
                      {
                          if (verbose)
                              cat('!!!! Counting all overlapping bins !!!\n')                  
                          chunk[, ":="(bin1 = 1 + floor((V2)/window), bin2 = 1 + floor((V2+V3)/window))]
                          chunk = chunk[, list(V1, bin = bin1:bin2), by = list(ix = 1:length(V1))]
                      }
              }
          else ## just count reads
              {
                  cat('counting reads\n')
                  chunk = fread(paste(chunk, collapse = "\n"), header = F)
                  chunk[, bin := 1 + floor((V2)/window)]
              }

          tabs = chunk[, list(newcount = length(V1)), by = list(chr = as.character(V1), bin)] ## tabulate reads to bins data.table style
          counts[tabs, count := count + newcount] ## populate latest bins in master data.table

         ## should be no memory issues here since we preallocate the data table .. but they still appear
          if (do.gc)
              {
                  print('GC!!')
                  print(gc())
              }
#          print(tables())
          
          ## report timing
          if (verbose)
              {
                  cat('bam.cov.tile.st ', bam.file, 'chunk', i, 'num fragments processed', i*chunksize, '\n')
                  timeelapsed = as.numeric(difftime(Sys.time(), st, units = 'hours'))
                  meancov = i * chunksize / counts[tabs[nrow(tabs),], ]$rowid  ## estimate comes from total reads and "latest" bin filled
                  totreads = meancov * numwin
                  tottime = totreads*timeelapsed/(i*chunksize)
                  rate = i*chunksize / timeelapsed / 3600
                  cat('mean cov:', round(meancov,1), 'per bin, estimated tot fragments:', round(totreads/1e6,2), 'million fragments, processing', rate,
                      'fragments/second\ntime elapsed:', round(timeelapsed,2), 'hours, estimated time remaining:', round(tottime - timeelapsed,2), 'hours', ', estimated total time', round(tottime,2), 'hours\n')
              }
      }

    gr = GRanges(counts$chr, IRanges(counts$start, counts$end), count = counts$count, seqinfo = Seqinfo(names(sl), sl))    
    if (verbose)
        cat("Finished computing coverage, and making GRanges\n")
    close(p)

    if (!is.null(region))
        system(sprintf('rm %s.bam %s.bam.bai', tmp.fn, tmp.fn))
   
    return(gr)
}

#' Compute rpkm counts from counts
#'
#' takes countbam (or bam.cov.gr) output "counts" and computes rpkm by aggregating across "by" variable
#' @param counts GRanges, data.table or data.frame with records, width fields
#' @param by Field to group counts by 
#' @note The denominator (ie total reads) is just the sum of counts$records
#' @export
counts2rpkm = function(counts, by)
  {
    out = aggregate(1:nrow(counts), by = list(by), FUN = function(x) sum(counts$records[x])/ sum(counts$width[x]/1000));
    out[,2] = out[,2]/sum(counts$records)*1e6;
    names(out) = c('by', 'rpkm');
    return(out);
  }


#' Create GRanges of read mates from reads
#'
#' @return \code{GRanges} corresponding to mates of reads
#' @name get.mate.gr
#' @export
get.mate.gr = function(reads)
  {

    if (inherits(reads, 'GRanges')) {
      mpos = values(reads)$mpos
      mrnm = as.vector(values(reads)$mrnm)
      mapq = values(reads)$MQ
      bad.chr = !(mrnm %in% seqlevels(reads)); ## these are reads mapping to chromosomes that are not in the current "genome"
      mrnm[bad.chr] = as.character(seqnames(reads)[bad.chr]) # we set mates with "bad" chromosomes to have 0 width and same seqnames (ie as if unmapped)
    } else if (inherits(reads, 'data.table')) {
      mpos <- reads$mpos
      mrnm <- reads$mrnm
      mapq = reads$MQ
      bad.chr <- !(mrnm %in% c(seq(22), 'X', 'Y', 'M'))
      mrnm[bad.chr] <- reads$seqnames[bad.chr]
    }
    
    if (inherits(reads, 'GappedAlignments'))
      mwidth = qwidth(reads)
    else
      {
        mwidth = reads$qwidth
        mwidth[is.na(mwidth)] = 0
      }
        
    mwidth[is.na(mpos)] = 0
    mwidth[bad.chr] = 0;  # we set mates with "bad" chromosomes to have 0 width
    mpos[is.na(mpos)] = 1;
    
    if (inherits(reads, 'GappedAlignments'))
      GRanges(mrnm, IRanges(mpos, width = mwidth), strand = c('+', '-')[1+bamflag(reads)[, 'isMateMinusStrand']], seqlengths = seqlengths(reads), qname = values(reads)$qname, mapq = mapq)
    else if (inherits(reads, 'GRanges'))
      GRanges(mrnm, IRanges(mpos, width = mwidth), strand = c('+', '-')[1+bamflag(reads$flag)[, 'isMateMinusStrand']], seqlengths = seqlengths(reads), qname = values(reads)$qname, mapq = mapq)
    else if (inherits(reads, 'data.table'))
      ab=data.table(seqnames=mrnm, start=mpos, end=mpos + mwidth - 1, strand=c('+','-')[1+bamflag(reads$flag)[,'isMateMinusStrand']], qname=reads$qname, mapq = mapq)
  }

#' Takes reads object and returns grl with each read and its mate (if exists)
#'
#' @param reads \code{GRanges} holding reads
#' @param as.grl Default TRUE. Return as a \code{GRangesList}
#' @param verbose Default FALSE
#' @name get.pairs.grl
#' @export
get.pairs.grl = function(reads, as.grl = TRUE, verbose = F)
  {

    isdt <- inherits(reads, 'data.table')

    bad.col = c("seqnames", "ranges", "strand", "seqlevels",
      "seqlengths", "isCircular", "genome", "start", "end", "width", "element")

    if (verbose)
      cat('deduping\n')

    if (is(reads, 'GappedAlignmentPairs'))
      reads = unlist(reads)

    if (inherits(reads, 'GRanges')) {
      d <- duplicated(values(reads)$qname) ## duplicates are already paired up
      qpair.ix <- values(reads)$qname %in% unique(values(reads)$qname[d])
    } else if (isdt) {
      d <- duplicated(reads$qname)
      qpair.ix <- reads$qname %in% unique(reads$qname[d])
    }
    
    if (!inherits(reads, 'GenomicRanges') && !inherits(reads, 'data.table'))
      {
        if (verbose)
          cat('converting to granges\n')
        r.gr = granges(reads)
      }
    else if (!isdt)
      r.gr = reads[, c()]
    else
      r.gr <- reads

    if (verbose)
      cat('grbinding\n')
    
    m.gr = get.mate.gr(reads[!qpair.ix]);

    if (inherits(reads, 'GRanges')) {
      m.val <- values(m.gr)
      values(m.gr) = NULL;
      r.gr = c(r.gr, m.gr);
      mcols(r.gr) <- rrbind2(mcols(reads)[, setdiff(colnames(values(reads)), bad.col)], m.val)
    } else if (isdt) {
      m.gr <- m.gr[, setdiff(colnames(reads), colnames(m.gr)) := NA, with=FALSE]
      r.gr <- rbind(reads, m.gr, use.names=TRUE)
      setkey(r.gr, qname)
    }

    if (as.grl && !isdt) {
      if (verbose)
        cat('splitting\n')
      return(split(r.gr, as.character(r.gr$qname)))
    } else {
      return(r.gr)
    }
  }

#' chunk
#'
#' takes same input as seq (from, to, by, length.out) and outputs a 2 column matrix of indices 
#' corresponding to "chunks"
#'
#' @param from dummy
#' @param to dummy
#' @param by Default 1
#' @param length.out Default NULL
#' @return 2-column matrix of indices corresponding to "chunks"
#' @export
chunk = function(from, to = NULL, by = 1, length.out = NULL)
  {
    if (is.null(to))
      {
        to = from;
        from = 1;
      }

    if (is.null(length.out))
      tmp = c(seq(from = from, to = to, by = by), to + 1)
    else
      tmp = c(seq(from = from, to = to, length.out = length.out), to + 1)
    
    out = floor(cbind(tmp[-length(tmp)], tmp[-1]-1))
    
    return(out)
  }



#' rrbind = function(df1, df2, [df3 ... etc], )
#'
#' like rbind, but takes the intersecting columns of the dfs
#'
#' if union flag is used then will take union of columns (and put NA's for columns of df1 not in df2 and vice versa)
#' @param ... list of data frames to concatenate
#' @param union if union flag is used then will take union of columns (and put NA's for columns of df1 not in df2 and vice versa). Default TRUE
#' @name rrbind
rrbind = function(..., union = TRUE)
  {     
    dfs = list(...);  # gets list of data frames
    if (any(ix <- sapply(dfs, function(x) class(x)[1])!='data.frame'))
        dfs[ix] = lapply(dfs[ix], as.data.frame)

    dfs = dfs[!sapply(dfs, is.null)]    
    dfs = dfs[sapply(dfs, ncol)>0]

    ## defactorize (need to do to cat without introducing NA's in weird places)
    dfs = lapply(dfs, function(x) { for (y in names(x)) if (is.factor(x[,y])) x[, y] = as.character(x[, y]); return(x)})
    
    names.list = lapply(dfs, names);
    classes = unlist(lapply(dfs, function(x) sapply(names(x), function(y) class(x[, y]))))
    cols = unique(unlist(names.list));
    unshared = lapply(names.list, function(x) setdiff(cols, x));
    unshared.u = unique(unlist(unshared))
    ix = which(sapply(dfs, nrow)>0)
    expanded.dfs = lapply(ix, function(x)
      {
        dfs[[x]][, unshared[[x]]] = as.character(NA);
        return(dfs[[x]][, cols, drop = F])
      })
    
    out = do.call('rbind', expanded.dfs);
    
    if (any(uix <<- which(classes[unshared.u] != 'character')))
      {
          ix = match(unshared.u, names(out))
          for (j in uix) ### HACK to prevent stupid class mismatches leading to NA BS
              out[, ix[j]] = as(out[, ix[j]], classes[unshared.u[j]])
      }
    
    if (!union)
      {
        shared = setdiff(cols, unique(unlist(unshared)))
        out = out[, shared];
      }    
    
   return(out)
}

#' count.clips
#'
#' takes gr or gappedalignment object and uses cigar field (or takes character vector of cigar strings)
#' and returns data frame with fields (for character input)
#' $right.clips number of "right" soft clips (eg cigar 89M12S)
#' #left.clips number of "left" soft clips (eg cigar 12S89M)
#' or appends these fields to the reads object
#'
#' @param reads GenomicRanges holding the reads
#' @param hard [Default TRUE] option counts hard clips
#' @name count.clips
#' @export
count.clips = function(reads, hard = FALSE)
{
  if (length(reads) == 0)
    return(reads)
  if (inherits(reads, 'GRanges') | inherits(reads, 'GappedAlignments'))
    cigar = values(reads)$cigar
  else
    cigar = reads;

  if (!inherits(cigar, 'character') & !inherits(cigar, 'factor'))
    stop('Input must be GRanges, GappedAlignments, or character vector')
  
  out = data.frame(left.clips = rep(0, length(cigar)), right.clips = rep(0, length(cigar)));

  re.left = '^(\\d+)S.*'
  re.right = '.*[A-Z](\\d+)S$'

  lclip.ix = grep(re.left, cigar)
  rclip.ix = grep(re.right, cigar)

  if (length(lclip.ix)>0)
    {
      left.clips = gsub(re.left, '\\1', cigar[lclip.ix])
      out$left.clips[lclip.ix] = as.numeric(left.clips)
    }

  if (length(rclip.ix)>0)
    {
      right.clips = gsub(re.right, '\\1', cigar[rclip.ix])
      out$right.clips[rclip.ix] = as.numeric(right.clips)
    }

  if (inherits(reads, 'GRanges') | inherits(reads, 'GappedAlignments'))
    {
      values(reads)$right.clips = out$right.clips
      values(reads)$left.clips = out$left.clips
      out = reads
    }

  return(out)  
}  

#' varbase
#'
#' takes gr or gappedalignment object "reads" and uses cigar, MD, seq fields
#' to return variant bases and ranges
#'
#' returns grl (of same length as input) of variant base positions with character vector $varbase field populated with variant bases
#' for each gr item in grl[[k]], with the following handling for insertions, deletions, and substitution gr's:
#'
#' substitutions: nchar(gr$varbase) = width(gr) of the corresponding var
#' insertions: nchar(gr$varbase)>=1, width(gr) ==0
#' deletions: gr$varbase = '', width(gr)>=1
#'
#' Each gr also has $type flag which shows the cigar string code for the event ie
#' S = soft clip --> varbase represents clipped bases
#' I = insertion --> varbase represents inserted bases
#' D = deletion --> varbase is empty
#' X = mismatch --> varbase represents mismatched bases
#' @param reads GenomicRanges to extract variants from
#' @param soft [Default TRUE]
#' @param verbose [Default TRUE]
#' @name varbase
#' @export
varbase = function(reads, soft = TRUE, verbose = TRUE)
{
  nreads = length(reads)
  if (inherits(reads, 'GRangesList'))
    {
      was.grl = TRUE
      r.id = as.data.frame(reads)$element
      reads = unlist(reads)      
    }
  else if (inherits(reads, 'data.frame'))
    {      
      r.id = 1:nrow(reads)
      nreads = nrow(reads)
      was.grl = FALSE
    }
  else    
    {      
      r.id = 1:length(reads)
      was.grl = FALSE
    }

  if (!inherits(reads, 'GRanges') & !inherits(reads, 'GappedAlignments') & !inherits(reads, 'data.frame'))
    stop('Reads must be either GRanges, GRangesList, or GappedAlignment object')
  else if (inherits(reads, 'data.frame'))
    {
      if (is.null(reads$cigar) | is.null(reads$seq))
        stop('Reads must have cigar and seq fields specified')
    }
  else if (is.null(values(reads)$cigar) | is.null(values(reads)$seq))
    stop('Reads must have cigar and seq fields specified')

  if (is.data.frame(reads))
    {
      sl = NULL
      sn =  reads$seqnames
      cigar = as.character(reads$cigar)
      seq = as.character(reads$seq)     
      str = reads$strand
            
      if (!is.null(reads$MD))
          md = as.character(reads$MD)
      else
        md = rep(NA, length(cigar))
    }  
  else
    {
      sl = seqlengths(reads)
      sn =  seqnames(reads)
      cigar = as.character(values(reads)$cigar)
      seq = as.character(values(reads)$seq)      
      str = as.character(strand(reads))
      
      if (!is.null(values(reads)$MD))
          md = as.character(values(reads)$MD)
      else
        md = rep(NA, length(cigar))
    }
          
  if (!inherits(cigar, 'character') & !inherits(cigar, 'character') & !inherits(md, 'character'))
    stop('Input must be GRanges with seq, cigar, and MD fields populated or GappedAlignments object')

  ix = which(!is.na(cigar))

  if (length(ix)==0)
    return(rep(GRangesList(GRanges()), nreads))
  
  cigar = cigar[ix]
  seq = seq[ix]
  md = md[ix]
  str = str[ix]

    if (is.data.frame(reads))
    {
      r.start = reads$start[ix]
      r.end = reads$end[ix]
    }
  else
    {
      r.start = start(reads)[ix]
      r.end = end(reads)[ix];
    }

  flip = str == '-'
  
  seq = strsplit(seq, '')
 
  cigar.vals = lapply(strsplit(cigar, "\\d+"), function(x) x[2:length(x)])
  cigar.lens = lapply(strsplit(cigar, "[A-Z]"), as.numeric)

  clip.left = sapply(cigar.vals, function(x) x[1] == 'S')
  clip.right = sapply(cigar.vals, function(x) x[length(x)] == 'S')

  if (any(clip.left))    
    r.start[clip.left] = r.start[clip.left]-sapply(cigar.lens[which(clip.left)], function(x) x[1])

  if (any(clip.right))
    r.end[clip.right] = r.end[clip.right]+sapply(cigar.lens[which(clip.right)], function(x) x[length(x)])

  # split md string into chars after removing "deletion" signatures and also
  # any soft clipped base calls (bases followed by a 0
  md.vals = strsplit(gsub('([ATGC])', '|\\1|', gsub('\\^[ATGC]+', '|', md)), '\\|')

  # ranges of different cigar elements relative to query ie read-centric coordinates
  starts.seq = lapply(1:length(cigar.lens), function(i)
    {
      x = c(0, cigar.lens[[i]])
      x[which(cigar.vals[[i]] %in% c('D', 'H', 'N'))+1] = 0  ## deletions have 0 width on query
      cumsum(x[1:(length(x)-1)])+1
    })

  ends.seq = lapply(1:length(cigar.lens), function(i)
    {
      x = cigar.lens[[i]];
      x[which(cigar.vals[[i]] %in% c('D', 'H', 'N'))] = 0
      cumsum(x)
    })
 
  # ranges of different cigar elements relative to reference coordinatse
  starts.ref = lapply(1:length(cigar.lens), function(i)
    {
      x = c(0, cigar.lens[[i]]);
      x[which(cigar.vals[[i]] %in% c('I'))+1] = 0 ## insertions have 0 width on reference / subject
      cumsum(x[1:(length(x)-1)]) + r.start[i]
    })

  ends.ref = lapply(1:length(cigar.lens), function(i)
    {
      x = cigar.lens[[i]];
      x[which(cigar.vals[[i]] %in% c('I'))] = 0 
      cumsum(x) + r.start[i] - 1
    })

  # now using MD find coordinates of mismatched bases (using starts and ends of M regions)
  
  # find coords of subs on genome
  tmp.pos = lapply(1:length(md.vals), function(i)
    {
      x = md.vals[[i]]
      nix = grepl('[ATGC]', x);
      if (!any(nix))
        return(c())
      p = rep(0, length(x))
      p[!nix] = as.numeric(x[!nix])
      p[nix] = 1
      s.pos.m = cumsum(p)[nix] ## position of subs in read
      mix = cigar.vals[[i]]=='M' 
      m.st = cumsum(c(1, ends.seq[[i]][mix]-starts.seq[[i]][mix]+1))
      m.st.g = starts.ref[[i]][mix]
      m.st.r = starts.seq[[i]][mix]
      s.match = rep(NA, length(s.pos.m)) ## indices of matching "M" cigar element for each sub
      for (ii in 1:length(s.pos.m))
        {
          j = 0;
          done = FALSE
          for (j in 0:(length(m.st)-1))
            if (s.pos.m[ii] < m.st[j+1])
              break
          s.match[ii] = j
        }
      
      s.pos.g = m.st.g[s.match] + s.pos.m-m.st[s.match]
      s.pos.r = m.st.r[s.match] + s.pos.m-m.st[s.match]
      
      return(rbind(s.pos.g, s.pos.r))      
    })
  subs.pos = lapply(tmp.pos, function(x) x[1,])
  subs.rpos = lapply(tmp.pos, function(x) x[2,])  
#  subs.base = lapply(md.vals, grep, pattern = '[ATGC]', value = T)
  subs.base = lapply(1:length(seq), function(x) seq[[x]][subs.rpos[[x]]])
  
  # make sure md and cigar are consistent
  # (for some reason - sometimes soft clipped mismatches are included in MD leading to a longer MD string)
  # also some MD are NA
  mlen.cigar = sapply(1:length(ends.seq), function(x) {mix = cigar.vals[[x]]=='M'; sum(ends.seq[[x]][mix]-starts.seq[[x]][mix]+1)})
  mlen.md = sapply(md.vals, function(x) {ix = grepl('[ATGC]', x); sum(as.numeric(x[!ix])) + sum(nchar(x[ix]))})
  good.md = which(!is.na(md))
#  good.md = which(mlen.md == mlen.cigar & !is.na(md))
 
  if (any(na <- is.na(md)))
    {
      warning('MD field absent from one or more input reads')
      good.md = which(!na)
    }
  
  # now make granges of subs
  if (length(good.md)>0)
    {
      if (any(mlen.md[good.md] != mlen.cigar[good.md]))
        warning('the lengths of some MD strings do not match the number of M positions on the corresponding CIGAR string: some variants may not be correctly mapped to the genome')      
      
      iix.md = unlist(lapply(good.md, function(x) rep(x, length(subs.pos[[x]]))))
      tmp = unlist(subs.pos[good.md])
      if (!is.null(tmp))
        {
          subs.gr = GRanges(sn[ix][iix.md], IRanges(tmp, tmp), strand = '*')
          values(subs.gr)$varbase = unlist(subs.base[good.md])
          values(subs.gr)$type = 'X'
          values(subs.gr)$iix = ix[iix.md]
        }
      else
        subs.gr = GRanges()
    }
  else
    subs.gr = GRanges()
  
  iix = unlist(lapply(1:length(cigar.vals), function(x) rep(x, length(cigar.vals[[x]]))))
  cigar.vals = unlist(cigar.vals)
  cigar.lens = unlist(cigar.lens)
  starts.seq = unlist(starts.seq)
  ends.seq = unlist(ends.seq)
  starts.ref = unlist(starts.ref)
  ends.ref = unlist(ends.ref)
 
  ## pick up other variants (including soft clipped and indel)
  is.var = cigar.vals != 'M'
  iix = iix[is.var]
  cigar.vals = cigar.vals[is.var]
  cigar.lens = cigar.lens[is.var]
  starts.ref = starts.ref[is.var]
  ends.ref = ends.ref[is.var]
  starts.seq = starts.seq[is.var]
  ends.seq = ends.seq[is.var]
  str <- str[iix] # JEREMIAH

  if (length(cigar.vals)>0)
    {
      var.seq = lapply(1:length(cigar.vals),
        function(i)
        {
          if (ends.seq[i]<starts.seq[i])
            return('') # deletion
          else
            seq[[iix[i]]][starts.seq[i]:ends.seq[i]] #insertion
        })
      other.gr = GRanges(sn[ix][iix], IRanges(starts.ref, ends.ref), strand = str, seqlengths = sl)
      values(other.gr)$varbase = sapply(var.seq, paste, collapse = '')
      values(other.gr)$type = cigar.vals
      values(other.gr)$iix = ix[iix];
      
      out.gr = sort(c(subs.gr, other.gr))
    }
  else
    out.gr = subs.gr

  
  # add default colors to out.gr
  VAR.COL = get.varcol()
  
  col.sig = as.character(out.gr$type)
  xix = out.gr$type == 'X'
  col.sig[xix] = paste(col.sig[xix], out.gr$varbase[xix], sep = '')
  out.gr$col = VAR.COL[col.sig]
  out.gr$border = out.gr$col

  if (!soft)
    if (any(soft.ix <<- out.gr$type == 'S'))
      out.gr = out.gr[-which(soft.ix)]
  
  out.grl = rep(GRangesList(GRanges()), nreads)  
  out.iix = r.id[values(out.gr)$iix]
  values(out.gr)$iix = NULL
  
  tmp.grl = split(out.gr, out.iix)
  out.grl[as.numeric(names(tmp.grl))] = tmp.grl
  values(out.grl)$qname[r.id] = reads$qname
  
  return(out.grl)  
}

#' splice.cigar
#'
#' takes gr or gappedalignment object "reads" and parses cigar fields
#' to return grl corresponding to spliced alignments on the genome corresponding to
#' portions of the cigar
#'
#' ie each outputted grl item contains the granges corresponding to all non-N portions of cigar string
#'
#' if grl provided as input (e.g. paired reads) then all of the spliced ranges resulting from each
#' input grl item will be put into the corresponding output grl item 
#'
#' NOTE: does not update MD tag
#'
#' if use.D = TRUE, then will treat "D" (deletion) in addition to "N" flags as indicative of deletion event.
#' @param reads \code{GRanges} reads
#' @param verbose Default TRUE
#' @param fast Default TRUE
#' @param use.D Default TRUE
#' @param rem.soft Default TRUE
#' @param get.seq Default FALSE
#' @param return.grl Default TRUE
#' @name splice.cigar
#' @export
splice.cigar = function(reads, verbose = TRUE, fast = TRUE, use.D = TRUE, rem.soft = TRUE, get.seq = FALSE, return.grl = TRUE)
{
  nreads = length(reads)

  if (nreads==0)
      if (return.grl)
          return(GRangesList())
      else
          return(GRanges)  
  
  if (inherits(reads, 'GRangesList'))
    {
      was.grl = TRUE
      r.id = as.data.frame(reads)$element
      reads = unlist(reads)
    }
  else
    {
      r.id = 1:length(reads)
      was.grl = FALSE
    }


  if (is.data.frame(reads))
      {
          sl = NULL
          sn =  reads$seqnames
          cigar = as.character(reads$cigar)
          seq = as.character(reads$seq)     
          str = reads$strand
          
          if (!is.null(reads$MD))
              md = as.character(reads$MD)
          else
              md = rep(NA, length(cigar))
      }  
  else
      {
          sl = seqlengths(reads)
          sn =  seqnames(reads)
          cigar = as.character(values(reads)$cigar)
          seq = as.character(values(reads)$seq)      
          str = as.character(strand(reads))
          
          if (!is.null(values(reads)$MD))
              md = as.character(values(reads)$MD)
          else
              md = rep(NA, length(cigar))
      }

  
  if (!inherits(reads, 'GRanges') & !inherits(reads, 'GappedAlignments') & !inherits(reads, 'data.frame'))
    stop('Reads must be either GRanges, GRangesList, or GappedAlignment object')
  else if (is.null(values(reads)$cigar) | is.null(values(reads)$seq))
    stop('Reads must have cigar and seq fields specified')

  if (!inherits(cigar, 'character') & !inherits(cigar, 'character') & !inherits(md, 'character'))
    stop('Input must be GRanges with seq, cigar, and MD fields populated or GappedAlignments object')

  
  ix = which(!is.na(reads$cigar))
  
  if (length(ix)==0)
      if (return.grl)
          return(rep(GRangesList(GRanges()), nreads))
      else
          return(GRanges())

  if (fast)
      {
          ir = cigarRangesAlongReferenceSpace(reads[ix]$cigar, N.regions.removed = FALSE, with.ops = TRUE, reduce.ranges = FALSE)
          irul = unlist(ir)
          out.gr = GRanges(rep(seqnames(reads)[ix], elementLengths(ir)), shift(IRanges(irul), rep(start(reads)[ix]-1, elementLengths(ir))),
              strand = rep(strand(reads)[ix], elementLengths(ir)), seqlengths = seqlengths(reads))
          out.gr$type = names(irul)
          out.gr$rid = ix[rep(1:length(ir), elementLengths(ir))]
          out.gr$riid = unlist(lapply(elementLengths(ir), function(x) 1:x))
          out.gr$fid = r.id[out.gr$rid]          
          out.gr$qname = reads$qname[out.gr$rid]
                    
          if (return.grl)
              {
                  out.grl = rep(GRangesList(GRanges()), nreads)
                  tmp.grl = split(out.gr, out.gr$fid)
                  out.grl[as.numeric(names(tmp.grl))] = tmp.grl
                  return(out.grl)
              }
          else
              return(out.gr)          
      }
  else
      {
                
          cigar = cigar[ix]
          str = str[ix]

          if (is.data.frame(reads))
              {
                  r.start = reads$start[ix]
                  r.end = reads$end[ix]
              }
          else
              {
                  r.start = start(reads)[ix]
                  r.end = end(reads)[ix];
              }
          
          flip = str == '-'   
          cigar.vals = lapply(strsplit(cigar, "\\d+"), function(x) x[2:length(x)])
          cigar.lens = lapply(strsplit(cigar, "[A-Z]"), as.numeric)

          clip.left = sapply(cigar.vals, function(x) x[1] == 'S')
          clip.right = sapply(cigar.vals, function(x) x[length(x)] == 'S')

                                        # ranges of different cigar elements relative to query ie read-centric coordinates
          starts.seq = lapply(1:length(cigar.lens), function(i)
              {
                  x = c(0, cigar.lens[[i]])
                  x[which(cigar.vals[[i]] %in% c('D', 'H', 'N'))+1] = 0  ## deletions have 0 width on query
                  cumsum(x[1:(length(x)-1)])+1
              })

          ends.seq = lapply(1:length(cigar.lens), function(i)
              {
                  x = cigar.lens[[i]];
                  x[which(cigar.vals[[i]] %in% c('D', 'H', 'N'))] = 0
                  cumsum(x)
              })
          
                                        # ranges of different cigar elements relative to reference coordinatse
          starts.ref = lapply(1:length(cigar.lens), function(i)
              {
                  x = c(0, cigar.lens[[i]]);
                  x[which(cigar.vals[[i]] %in% c('I'))+1] = 0 ## insertions have 0 width on reference / subject
                  cumsum(x[1:(length(x)-1)]) + r.start[i]
              })

          ends.ref = lapply(1:length(cigar.lens), function(i)
              {
                  x = cigar.lens[[i]];
                  x[which(cigar.vals[[i]] %in% c('I'))] = 0 
                  Cumsum(x) + r.start[i] - 1
              })
          
          iix = unlist(lapply(1:length(cigar.vals), function(x) rep(x, length(cigar.vals[[x]]))))
          cigar.vals = unlist(cigar.vals)
          cigar.lens = unlist(cigar.lens)
          starts.seq = unlist(starts.seq)
          ends.seq = unlist(ends.seq)
          starts.ref = unlist(starts.ref)
          ends.ref = unlist(ends.ref)
          
          ## pick up splice (including soft clipped and indel)
          splice.char = 'N'
          
          if (use.D)
              splice.char = c(splice.char, 'D')

          if (rem.soft)
              splice.char = c(splice.char, 'S')
          
          is.splice = !(cigar.vals %in% splice.char)
          iix = iix[is.splice]
          cigar.vals = cigar.vals[is.splice]
          cigar.lens = cigar.lens[is.splice]
          starts.ref = starts.ref[is.splice]
          ends.ref = ends.ref[is.splice]
         starts.seq = starts.seq[is.splice]
          ends.seq = ends.seq[is.splice]
          str <- str[iix] # 

          out.gr = GRanges()
          if (length(cigar.vals)>0)
              {
                  other.gr = GRanges(sn[ix][iix], IRanges(starts.ref, ends.ref), strand = str, seqlengths = sl)
                  
                  if (get.seq)
                      {
                          var.seq = lapply(1:length(cigar.vals),
                              function(i)
                                  {
                                      if (ends.seq[i]<starts.seq[i])
                                          return('') # deletion
                                      else
                                          seq[[iix[i]]][starts.seq[i]:ends.seq[i]] #insertion
                                  })
                          values(other.gr)$seq = sapply(var.seq, paste, collapse = '')
                      }
                  
                  values(other.gr)$type = cigar.vals
                  values(other.gr)$iix = ix[iix];

                  out.gr = c(other.gr, out.gr)
              }
          
          out.gr$rid = out.gr$iix
          out.iix = r.id[out.gr$rid]
          values(out.gr)$iix = NULL
          out.gr$fid = out.iix
          out.gr$qname = reads$qname[out.gr$rid]
          
          if (return.grl)
              {
                  out.grl = rep(GRangesList(GRanges()), nreads)
                  tmp.grl = split(out.gr, out.iix)
                  out.grl[as.numeric(names(tmp.grl))] = tmp.grl
                  return(out.grl)
              }
          else
              return(out.gr)
      }
}
   
#' Get GC content from reference genome
#'
#' Uses BSgenome package to compute gc content for a collection of segments in seg data frame ($chr, $start, $end or $chr, $pos1, $pos2 or $chr, $begin, $end)
#' Returns vector of gc content of length nrow(segs).
#' @param segs Segment data frame to pull gc from
#' @param bs_genome A \code{\link{BSgenome}} object. Perhaps \code{BSgenome.Hsapiens.UCSC.hg19::Hsapiens}
#' @import BSgenome
#' @export
#' @name gc_content
gc_content = function(segs, bs_genome) ##build = 'hg19')
  {
    segs = standardize_segs(segs, chr = TRUE);

## NEW 
        tmp = getSeq(bs_genome, segs$chr, segs$pos1, segs$pos2, as.character = TRUE)
##     if (build == 'hg19') {
##       if (requireNamespace("BSgenome.Hsapiens.UCSC.hg19", quietly=TRUE)) {
##         tmp = getSeq(BSgenome.Hsapiens.UCSC.hg19::Hsapiens, segs$chr, segs$pos1, segs$pos2, as.character = TRUE)
##       }
##     }
##     else if (build == 'hg18') {
##       if (requireNamespace("BSgenome.Hsapiens.UCSC.hg18", quietly=TRUE)) {
##         tmp = getSeq(BSgenome.Hsapiens.UCSC.hg18::Hsapiens, segs$chr, segs$pos1, segs$pos2, as.character = TRUE)
##       }
##     }
##     else
##       stop('gc_content: hg build not recognized');

    ## OLD
##    if (build == 'hg19')
##      library(BSgenome.Hsapiens.UCSC.hg19)
##    else if (build == 'hg18')
 ##     library(BSgenome.Hsapiens.UCSC.hg18)
  ##  else
   ##   stop('gc_content: hg build not recognized');

   ## tmp = getSeq(Hsapiens, segs$chr, segs$pos1, segs$pos2, as.character = T)

    return(as.numeric(sapply(gregexpr('[GC]', tmp), length)/sapply(tmp, nchar)))
  }


#' bamflag
#'
#' shortcut .. assumes reads are GappedAlignments with flag variable or actual integers representing bam flag
#' @param reads GenomicRanges holding the reads
#' @name bamflag
#' @export
bamflag = function(reads)
  {
    if (inherits(reads, 'GappedAlignments') | inherits(reads, 'data.frame') | inherits(reads, 'GRanges'))
      bf = reads$flag
    else
      bf = reads

    out = matrix(as.numeric(intToBits(bf)), byrow = T, ncol = 32)[, 1:12, drop = FALSE]
    colnames(out) = c('isPaired', 'isProperPair', 'isUnmappedQuery', 'hasUnmappedMate', 'isMinusStrand', 'isMateMinusStrand', 'isFirstMateRead', 'isSecondMateRead', 'isNotPrimaryRead', 'isNotPassingQualityControls', 'isDuplicate', 'isSupplementary')

    return(out)
#    if (inherits(reads, 'GappedAlignments'))      
#      return(bamFlagAsBitMatrix(values(reads)$flag))
#    else
#      return(bamFlagAsBitMatrix(reads))
  }


#' bamtag
#'
#' outputs a tag that cats qname, first vs first second mate +/- secondary alignment +/- gr.string
#' to give an identifier for determine duplicates in a read pile
#' @param reads GenomicRanges holding the reads
#' @name bamflag
#' @export
bamtag = function(reads, secondary = F, gr.string = F)
  {
    grs = sec = NULL
    if (secondary)
      sec = bamflag(read$flag[, 'isNotPrimaryRead'])
    
    if (gr.string)
      grs = gr.string(reads, mb  = F)
      
    return(paste(reads$qname, ifelse(bamflag(reads$flag)[, 'isFirstMateRead'], '1', '2'), grs, sec, sep = '_'))
  }


#' import.ucsc
#'
#' wrapper around rtracklayer import that
#' (1) handles "as" formats
#' (2) has additional flag chrsub to sub in 'chr' in selection, and then sub it out of the output
#' @name import.ucsc
#' @importClassesFrom "rtracklayer" WIGFile BEDFile BigWigFile
#' @export
import.ucsc = function(con, selection = NULL, text, chrsub = TRUE, verbose = FALSE, as = NULL, ...)
  {
    si = NULL;

    if (verbose)
        cat('importing', as.character(con), '\n')
    
    if (grepl('(\\.bw)|(\\.bigwig)', con, ignore.case = TRUE))
      {
          if (is.null(as))
              as = 'RleList'
          
        if (is.character(con))
          f = BigWigFile(normalizePath(con))
        else
          f = con
        
        si = tryCatch(seqinfo(f), error = function(con) NULL)
      }
    else if (grepl('\\.wig', con, ignore.case = TRUE))
      {
          if (is.null(as))
              as = 'RleList'
          
        if (is.character(con))
          f = WIGFile(normalizePath(con))
        else
          f = con
        
        si = tryCatch(seqinfo(f), error = function(con) NULL)
      }
    else if (grepl('\\.bed', con, ignore.case = TRUE))
      {
          if (is.null(as))
              as = 'GRanges'
          
        if (is.character(con))
          f = BEDFile(normalizePath(con))
        else
          f = con
                                        #                            si = tryCatch(seqinfo(f), error = function(con) NULL)
        bed.style = T
      }
    else if (grepl('\\.gff', con, ignore.case = TRUE))
      {
          if (is.null(as))
              as = 'GRangesList'
          
        if (is.character(con))
          f = GFFFile(normalizePath(con))
        else
          f = con
        
        si = tryCatch(seqinfo(f), error = function(con) NULL)
      }
    else if (grepl('\\.2bit', con, ignore.case = T))
      {
          if (is.null(as))
              as = 'RleList'
        if (is.character(con))
          f = TwoBitFile(normalizePath(con))
        else
          f = con
        
        si = tryCatch(seqinfo(f), error = function(con) NULL)
      }
    else if (grepl('\\.bedgraph', con, ignore.case = T))
      {
          if (is.null(as))
              as = 'GRanges'
          
        if (is.character(con))
          f = BEDGraphFile(normalizePath(con))
        else
          f = con
        #                           si = tryCatch(seqinfo(f), error = function(con) NULL)
        bed.style = T
      }
    else
      f = con
    
    if (chrsub & !is.null(si) & !is.null(selection))
      selection = gr.fix(gr.chr(selection), si, drop = T)

    if (class(f) %in% c('BEDFile'))
        {
            if (!is.null(selection))
                out = import(f, selection = selection, asRangedData = FALSE, ... )
            else
                out = import(f, asRangedData = FALSE, ...)
        }
    else
        {            
            if (!is.null(selection))
                out = import(f, selection = selection, as = as, ... )
            else
                out = import(f, as = as, ...)
        }

    if (!is(out, 'GRanges'))
        out = as(out, 'GRanges')

#    if (chrsub & !is.null(si))
    if (chrsub)
      out = gr.sub(out, 'chr', '')

    if (verbose)
        cat('Finished importing', as.character(con), '\n')
    
    return(out)
  }

#####################################################################
# 
#
# $$\      $$\ $$\                                      $$\     $$\ $$\ 
# $$$\    $$$ |\__|                                     $$ |    \__|$$ |
# $$$$\  $$$$ |$$\  $$$$$$$\  $$$$$$$\       $$\   $$\$$$$$$\   $$\ $$ |
# $$\$$\$$ $$ |$$ |$$  _____|$$  _____|      $$ |  $$ \_$$  _|  $$ |$$ |
# $$ \$$$  $$ |$$ |\$$$$$$\  $$ /            $$ |  $$ | $$ |    $$ |$$ |
# $$ |\$  /$$ |$$ | \____$$\ $$ |            $$ |  $$ | $$ |$$\ $$ |$$ |
# $$ | \_/ $$ |$$ |$$$$$$$  |\$$$$$$$\       \$$$$$$  | \$$$$  |$$ |$$ |
# \__|     \__|\__|\_______/  \_______|       \______/   \____/ \__|\__|
#
#
# Misc util
################

#' Improved rbidn for intersecting columns of data.frames or data.tables
#'
#' like rbind, but takes the intersecting columns of the dfs
#' rrbind = function(df1, df2, [df3 ... etc], )
#' @param ... list of data frames to concatenate
#' @param union if union flag is used then will take union of columns (and put NA's for columns of df1 not in df2 and vice versa). Default TRUE
#' @param as.data.table [Default FALSE] return as a \link{data.table}
#' @export
rrbind2 = function(..., union = T, as.data.table = FALSE)
  {
    dfs = list(...);  # gets list of data frames
    dfs = dfs[!sapply(dfs, is.null)]    
    dfs = dfs[sapply(dfs, ncol)>0]
    names.list = lapply(dfs, names);
    cols = unique(unlist(names.list));
    unshared = lapply(names.list, function(x) setdiff(cols, x));
    ix = which(sapply(dfs, nrow)>0)

    ## only call expanded dfs if needed
    if (any(sapply(unshared, length) != 0)) 
      expanded.dts <- lapply(ix, function(x) {
          tmp = dfs[[x]]
          if (is.data.table(dfs[[x]]))
              tmp = as.data.frame(tmp)
          tmp[, unshared[[x]]] = NA;
          return(as.data.table(as.data.frame(tmp[, cols])))
      })
    else
      expanded.dts <- lapply(dfs, as.data.table)

    ## convert data frames (or DataFrame) to data table.
    ## need to convert DataFrame to data.frmae for data.table(...) call. 
    ## structure call is way faster than data.table(as.data.frame(...))
    ## and works on data.frame and DataFrame
#    dts <- lapply(expanded.dfs, function(x) structure(as.list(x), class='data.table'))       
 #   rout <- data.frame(rbindlist(dts))

    rout <- rbindlist(expanded.dts)
    if (!as.data.table)
        rout = as.data.frame(rout)
    
    if (!union)
      {
        shared = setdiff(cols, unique(unlist(unshared)))
        rout = rout[, shared];
      }    
    
   return(rout)
  }

#' vaggregate
#'
#' same as aggregate except returns named vector
#' with names as first column of output and values as second
#' @param ... things to aggregate
#' @export
vaggregate = function(...)
  {
    out = aggregate(...);
    return(structure(out[,ncol(out)], names = do.call(paste, lapply(names(out)[1:(ncol(out)-1)], function(x) out[,x]))))
  }


#' Identify matches between query and dictionary
#'
#' Wrapper around matchPdict to identify matches between a query
#' string query and dictionary dict (both BString objects or subclasses)
#'
#' @param query Query
#' @param dict Dictionary
#' @param midpoint Flag for output the coordinates of the match as the location,
#'   where the midpoint of the dict string matches the given query. Default FALSE
#' @return a vector of indices of length width(query) that contains
#' indices of the (starting) dictionary in the query string
#' @export
match.bs = function(query, dict, midpoint = FALSE)
  {
    names(dict) = as.character(1:length(dict))
    
    tmp = sort(unlist(matchPDict(dict, query)))
    out = rep(NA, length(query))
    
    if (!midpoint)
      out[start(tmp)] = as.numeric(names(tmp))
    else
      out[floor((start(tmp)+end(tmp))/2)] = as.numeric(names(tmp))

    return(out)
  }

#' Output standard human genome seqlengths
#'
#' Outputs a standard seqlengths for human genome +/- "chr"
#' @param hg19 Flag for dealing with hg19 (vs hg18). Default TRUE
#' @param chr Flag for whether to keep "chr". Default FALSE
#' @param include.junk Flag for whether to not trim to only 1-22, X, Y, M. Default FALSE
#' @return Seqlengths
#' @import BSgenome
#' @export
hg_seqlengths = function(hg19 = TRUE, chr = FALSE, include.junk = FALSE)
  {
    hg = read_hg(hg19)

    sl = seqlengths(hg)

    if (!include.junk)
      sl = sl[c(paste('chr', 1:22, sep = ''), 'chrX', 'chrY', 'chrM')]
    
    if (!chr)
      names(sl) = gsub('chr', '', names(sl))

    return(sl)          
  }

#' Wrapper around BSgenome call
#'
#' Retreives either the BSgenome hg18 or hg19 genome by default.  Requires packages
#' BSgenome.Hsapiens.UCSC.hg19 for hg19 and BSgenome.Hsapiens.UCSC.hg19 for hg18.
#'
#' If fft = TRUE, can also also return the hg19 ffTrack (requires that the file exists)
#' Requires the existence of environment variable HG.FFT pointing to ffTrack .rds file.. 
#' 
#' @param hg19 Logical whether to return hg18 or hg19 BSgenome. Default TRUE
#' @param fft Logical whether to return an ffTrack. Default FALSE
#' @return BSgenome or ffTrack of the genome
#' @export
read_hg = function(hg19 = T, fft = F)
  {

    if (file.exists(Sys.getenv('HG.FFT')))
      REFGENE.FILE.HG19.FFT = Sys.getenv('HG.FFT')
    else if (file.exists('/home/unix/marcin/DB/ffTracks/hg19.rds'))
      REFGENE.FILE.HG19.FFT = '/home/unix/marcin/DB/ffTracks/hg19.rds'      
    else
      stop("Need to supply environment variable to FFtracked genome. Env Var: HG.FFT")
      
    if (fft)
      return(readRDS(REFGENE.FILE.HG19.FFT))
##     else
##       {
##         require(BSgenome)
##         if (hg19)
##           library(BSgenome.Hsapiens.UCSC.hg19)
##         else
##           library(BSgenome.Hsapiens.UCSC.hg18)
##       }
    return(Hsapiens)
  }

#' Retrieve genomic sequenes
#'
#' Wrapper around getSeq which does the "chr" and seqnames conversion if necessary
#' also handles GRangesList queries
#'
#' @param hg A BSgenome or and ffTrack object with levels = c('A','T','G','C','N')
#' @param gr GRanges object to define the ranges
#' @param unlist logical whether to unlist the final output into a single DNAStringSet. Default TRUE
#' @param mc.cores Optional multicore call. Default 1
#' @param mc.chunks Optional define how to chunk the multicore call. Default mc.cores
#' @param verbose Increase verbosity
#' @return DNAStringSet of sequences
#' @importClassesFrom "Biostrings" DNAString AAString AAStringSet DNAStringSet
#' @export
get_seq = function(hg, gr, unlist = TRUE, mc.cores = 1, mc.chunks = mc.cores,
     as.data.table = FALSE, verbose = FALSE)
{
  if (inherits(gr, 'GRangesList'))
    {
      grl = gr;
      old.names = names(grl);
      gr = unlist(grl);
      names(gr) = unlist(lapply(1:length(grl), function(x) rep(x, length(grl[[x]]))))
      seq = get_seq(hg, gr, mc.cores = mc.cores, mc.chunks = mc.chunks, verbose = verbose)
      cl = class(seq)
      out = split(seq, names(gr))
      out = out[order(as.numeric(names(out)))]
      if (unlist)
        out = do.call('c', lapply(out, function(x) do.call(cl, list(unlist(x)))))
      names(out) = names(grl)
      return(out)
    }
  else
    {
      if (is(hg, 'ffTrack'))
        {
          if (!all(sort(levels(hg)) == sort(c('A', 'T', 'G', 'C', 'N'))))
            cat("ffTrack not in correct format for get_seq, levels must contain only: 'A', 'T', 'G', 'C', 'N'\n")
        }
      else ## only sub in 'chr' if hg is a BSenome        
        if (!all(grepl('chr', as.character(seqnames(gr)))))
          gr = gr.chr(gr)
      
      gr = gr.fix(gr, hg)
      if (mc.cores>1)
        {
          ix = suppressWarnings(split(1:length(gr), 1:mc.chunks))

          if (is(hg, 'ffTrack'))
            {
              mcout <- mclapply(ix, function(x) 
                {
                  tmp = hg[gr[x]]
                  
                  if (any(is.na(tmp)))
                    stop("ffTrack corrupt: has NA values, can't convert to DNAString")

                  if (!as.data.table) {
                    bst = DNAStringSet(sapply(split(tmp, as.vector(Rle(1:length(x), width(gr)[x]))), function(y) paste(y, collapse = '')))
                    names(bst) = names(gr)[x]
                  } else {
                    bst <- data.table(seq=sapply(split(tmp, as.vector(Rle(1:length(x), width(gr)[x]))), function(y) paste(y, collapse='')))
                    bst[, names:=names(gr)[x]]
                  }
                  
                  if (any(strand(gr)[x]=='-'))
                    {
                      ix.tmp = as.logical(strand(gr)[x]=='-')
                      if (!as.data.table)
                        bst[ix.tmp] = Biostrings::complement(bst[ix.tmp])
                      else
                        bst$seq[ix.tmp] <- as.character(Biostrings::complement(DNAStringSet(bst$seq[ix.tmp])))
                    }
                                                     
                  if (verbose)
                    cat('.')
                  
                  return(bst)
                }
                , mc.cores = mc.cores)

            if (!as.data.table)
              {
                if (length(mcout)>1)
                  tmp = c(mcout[[1]], mcout[[2]])
                out <- do.call('c', mcout)[order(unlist(ix))]
              }
            else
               out <- rbindlist(mcout)
            }
          else
            {
              out = do.call(c, mclapply(ix, function(x)
                {
                  if (verbose)
                    cat('.')
                  return(getSeq(hg, gr[x]))
                }
                ,mc.cores = mc.cores))[order(unlist(ix))]
              if (verbose)
                cat('\n')
            }
        }
      else
        {
          if (is(hg, 'ffTrack'))
            {
              tmp = hg[gr]

              tmp[is.na(tmp)] = 'N'

              if (any(is.na(tmp)))
                stop("ffTrack corrupt: has NA values, can't convert to DNAString")
              
              if (as.data.table) {
                bst <- data.table(seq=sapply(split(tmp, as.vector(Rle(1:length(gr), width(gr)))), function(x) paste(x, collapse='')))
                bst[, names:=names(gr)]
              } else {
                bst = DNAStringSet(sapply(split(tmp, as.numeric(Rle(1:length(gr), width(gr)))), function(x) paste(x, collapse = '')))
                names(bst) = names(gr)
              }

              if (any(as.character(strand(gr))=='-'))
                {
                  ix = as.logical(strand(gr)=='-')
                  if (!as.data.table) {
                    bstc <- as.character(bst)
                    bstc[ix] <- as.character(Biostrings::complement(bst[ix]))
                    bst <- DNAStringSet(bstc)  ## BIZARRE bug with line below
                    #bst[ix] = Biostrings::complement(bst[ix])
                  } else { 
                    bst$seq[ix] <- as.character(Biostrings::complement(DNAStringSet(bst$seq[ix])))
                  }
                }
              
              return(bst)
            }
          else            
            out = getSeq(hg, gr)
        }      
      return(out)
    }
}

#' Deduplicate character vectors
#
#' Relabels duplicates in a character vector with .1, .2, .3
#' (where "." can be replaced by any user specified suffix)
#' @param x Character vector to deduplicate.
#' @param suffix User defined suffix
#' @return dedupclicated character vector
#' @export
dedup = function(x, suffix = '.')
{
  dup = duplicated(x);
  udup = unique(x[dup])
  udup.ix = lapply(udup, function(y) which(x==y));
  udup.suffices = lapply(udup.ix, function(y) c('', paste(suffix, 2:length(y), sep = '')))
  out = x;
  out[unlist(udup.ix)] = paste(out[unlist(udup.ix)], unlist(udup.suffices), sep = '');
  return(out)  
}

#' levapply
#
#' Applies FUN locally to levels of x and returns vector of length()
#' (eg can do a "local" order within levels)
#' @param x factor over which to apply the funciton
#' @param by group the factors?
#' @param FUN function to apply. Default \code{order}
#' @export
levapply = function(x, by, FUN = 'order')
  {
    if (!is.list(by))
      by = list(by)
    
    f = factor(do.call('paste', c(list(sep = '|'), by)))
    ixl = split(1:length(x), f);
    ixv = lapply(ixl, function(y) x[y])    
    res = structure(unlist(lapply(ixv, FUN)), names = unlist(ixl))
    out = rep(NA, length(x))
    out[as.numeric(names(res))] = res;
    return(out)
  }

#' Convert from chrXX to numeric format
#'
#' Convert from chrXX to numeric format
#' @param x factor, Rle or character vector with chromosome names
#' @param xy Flag to convert M to 25, Y to 24 and X to 23. Default FALSE
#' @return character vector with xy=FALSE, or numeric vector with xy=TRUE
#' @export
chr2num = function(x, xy = FALSE)
  {
    if (inherits(x, 'factor') | inherits(x, 'Rle'))
      x = as.character(x)
    
     out = gsub('chr', '', x);

     if (!xy)
       out = as.numeric(gsub('M', '25', gsub('Y', '24', gsub('X', '23', out))))
     
     return(out)
  }

#' Count bases in cigar string
#' 
#' Counts the total number of bases, per cigar, that fall into D, I, M, S categories.
#' countCigar makes no distinction between, for instance 1S2M2S, 2S2M1S, or 3S2M
#' @param cigar character vector of cigar strings
#' @return a 4-column, length(cigar)-row matrix with the total counts for each type
#' @export
countCigar <- function(cigar) {

  cigar.vals <- unlist(strsplit(cigar, "\\d+"))
  cigar.lens <- strsplit(cigar, "[A-Z]")
  lens <- nchar(gsub('\\d+', '', cigar))
  lens[is.na(cigar)] <- 1
  
  cigar.lens <- as.numeric(unlist(cigar.lens))
  cigar.vals <- cigar.vals[cigar.vals != ""]
  repr       <- rep(seq_along(cigar), lens)
  dt         <- data.table(val=cigar.vals, lens=cigar.lens, group=repr, key="val")

  smr.d      <- dt["D"][, sum(lens), by=group]
  smr.i      <- dt["I"][, sum(lens), by=group]  
  smr.m      <- dt["M"][, sum(lens), by=group]
  smr.s      <- dt["S"][, sum(lens), by=group]  
  
  out <- matrix(nrow=length(cigar), ncol=4, 0)
  out[smr.d$group,1] <- smr.d$V1
  out[smr.i$group,2] <- smr.i$V1
  out[smr.m$group,3] <- smr.m$V1
  out[smr.s$group,4] <- smr.s$V1
  colnames(out) <- c('D','I','M','S')

  return(out)
}

#' Filter reads by average PHRED score
#' Defines a cutoff score for the mean PHRED quality of a read
#' in a GRanges.
#' @param gr GRanges or data.table of reads that has a \code{qname} and \code{qual} field
#' @param cutoff cutoff score for mean PHRED quality. Default "+"
#' @return GRanges or data.table where reads have mean quality score >= cutoff
#' @export
gr.readfilter <- function(gr, cutoff = '+') {

  cutoff <- as.numeric(charToRaw(cutoff))
  qual   <- as.character(gr$qual)
  logvec <- sapply(qual, function(x) mean(as.numeric(charToRaw(x))) < cutoff)

  qn <- unique(gr$qname[logvec])

  gr <- gr[!(gr$qname %in% qn)]
  
  return(gr)
}

#' Check if reads are clipped
#' 
#' Returns a logical vector of length of the input GRanges that
#' that classifies a read as clipped or not. The user can specify
#' a cutoff value for how many bases need to be clipped.
#' @param gr Granges OR data.table that has \code{cigar} field and \code{qname} field
#' @param clip.cutoff Minimum number of bases that are clipped to call the reads as clipped
#' @return logical of length of input, denoting whether that read is part of a clipped read pair.
#' @export
gr.isclip <- function(gr, clip.cutoff=10) {
	if (inherits(gr, 'GRanges') && length(gr)==0)
	  return(logical(0))
        if (inherits(gr, 'data.table') && nrow(gr) == 0)
          return(logical(0))
          
        if (inherits(gr, 'GRanges'))
          nm <- names(mcols(gr))
        else
          nm <- colnames(gr)
	if (any(!('cigar' %in% nm)))
	  stop('gr.isclip: reads need flag and cigar')
        cig <- countCigar(gr$cigar)
        return(cig[,"S"] >= clip.cutoff)        
##	logvec <- grepl('[0-9][0-9]S', gr$cigar)
##	logvec[is.na(logvec)] <- FALSE
##	return(logvec)
}

#' Checks if reads are discordant
#'
#' Returns a logical vector denoting if a read is discordant.
#' There is only a minimum absolute isize, and any read below this is
#' not considered discordant. This will return logicals based on read pairs
#' @param gr Granges OR data.table that has \code{isize} field and \code{qname} field
#' @param isize Minimum insert size required to call discordant. Default 1000
#' @param unmap.only Find only pairs with an unmapped read
#' @return logical vector of length of input, denoting each read as discordant or not
#' @export
gr.isdisc <- function(gr, isize=1000, unmap.only=FALSE) {

	if (inherits(gr, 'GRanges') && length(gr)==0)
   	  return(logical(0))
	if (inherits(gr, 'data.table') && nrow(gr)==0)
   	  return(logical(0))

        if (inherits(gr, 'GRanges'))
          nm <- names(mcols(gr))
        else
          nm <- colnames(gr)

	if (any(!(c('isize') %in% nm )))
	  stop('gr.isdigsc: reads need flag and cigar')

        if (inherits(gr, 'GRanges'))
          st <- start(gr) == 1
        else
          st <- gr$start == 1
        if (unmap.only)
          logvec <- bitAnd(gr$flag, 8) != 0 | st  # last two get reads with unmapped mates and reads that are unmapped, respectively
        else
          logvec <- abs(gr$isize) >= isize | gr$isize==0 | bitAnd(gr$flag, 8) != 0 | st  # last two get reads with unmapped mates and reads that are unmapped, respectively        
	logvec[is.na(logvec)] <- FALSE  
        qn <- gr$qname[logvec]
        isdisc <- gr$qname %in% qn
        return(isdisc)
}

#' Minimal overlaps for GRanges/GRangesList
#' 
#' Takes any number of GRanges or GRangesList and reduces them to the minimal
#' set of overlapping windows, ignoring strand.
#' @param ... \code{GRanges} or \code{GRangesList}
#' @return GRanges
#' @export
gr.reduce <- function(...) {

	input <- list(...)
	for (i in seq_along(input)) {
          if (inherits(input[[i]], 'GRanges'))
	    input[[i]] <- reduce(gr.stripstrand(input[[i]]))
	  else if (inherits(input[[i]], 'GRangesList'))
   	    input[[i]] <- reduce(gr.stripstrand(unlist(input[[i]])))
          else 
   	    stop('reduce.window: Need to input GRanges or GRangesList objects')
          seqlengths(input[[i]]) <- seqlengths(input[[i]])*NA
        }
        
	output <- do.call('c', input)

	return(sort(reduce(output)))

}

#' Return windows with minimal coverage
#'
#' Takes a set of GRanges and removes any ranges that
#' don't have a minimal coverage value. If you give it
#' a GRangesList, you will get back an unlisted GRanges.
#'
#' @param gr \code{GRanges} to filter
#' @param min.cov Minimum number of overlaps to keep. Default 2
#' @param buffer Add a buffer to the ranges when computing overlaps. Default 0
#' @param ignore.strand Ignore the strand when comparing overlaps. Default TRUE
#' @param pintersect Force the pintersect option for \link{gr.findoverlaps}
#' @return GRanges
#' @export
gr.mincov <- function(gr, min.cov=2, buffer=0, ignore.strand=TRUE, pintersect=FALSE) {

  if (inherits(gr, 'GRangesList'))
    gr <- unlist(gr)
  if (!inherits(gr, 'GRanges'))
    stop('gr.mincov: Requires a GRanges input')

  gr2 <- gr.pad(gr, buffer)

  tab <- table(gr.findoverlaps(gr2, gr2, ignore.strand=ignore.strand, pintersect=pintersect)$subject.id)
  winkeep <- as.numeric(names(tab[tab >= min.cov]))

  return(gr[winkeep])
  
}


#' Nice padding
#' 
#' @return GRanges
#' @keywords internal
gr.pad = function(gr, pad)
    {
        start(gr) = pmax(1, start(gr)-pad)
        end(gr) = pmin(seqlengths(gr)[as.character(seqnames(gr))], end(gr)+pad)
        return(gr)
    }

#' gr2grl 
#' Quick way to make grl from list of indices into a GRanges gr
#'
#' @param gr \code{GRanges} to split
#' @param ix vector to split on
#' @export
gr2grl = function(gr, ix)
  {    
     out = split(gr[unlist(ix)], sapply(1:length(ix), function(x) rep(x, length(ix[[x]]))))
     if (!is.null(names(ix)))
       names(out) = names(ix)
     return(out)
  }

#' Remove chr prefix from GRanges seqlevels
#'
#' @param gr GRanges with chr seqlevel prefixes
#' @return GRanges without chr seqlevel prefixes
#' @export
gr.nochr = function(gr) {
    if (grepl('^chr', seqlevels(gr)[1]))
	seqlevels(gr) = gsub('^chr','', seqlevels(gr))
    return(gr)
}

#' Wrapper to base \code{system} function to call system (e.g. perl) from R.
#' The only benefit to this wrapper is a more controlled verbose argument.
#' 
#' @author Jeremiah Wala \email{jwala@@broadinstitute.org}
#' @param syscall string containing the system call
#' @param verbose print the syscall to screen, and it's stdout
#' @export
#' @examples
#' # system.call('perl s/[0-9]+//g file1 > file2')
system.call <- function(syscall, verbose=T) {
  if (verbose)
    print(syscall)
  if (verbose)
    system(syscall)
  else {
    system(syscall, intern=TRUE) #, stderr=FALSE, stdin=FALSE)
  }
    #system(syscall, intern=T, ignore.stdout=TRUE, ignore.stderr=TRUE)
}

#' @name seg2gr
#' @title Convert GRange like data.frames into GRanges
#' @description
#'
#' Take data frame of ranges "segs" and converts into granges object porting over additional value columns
#' "segs" data frame can obey any number of conventions to specify chrom, start, and end of ranges
#' (eg $pos1, $pos2, $Start_position, $End_position) --> see "standardize_segs" for more info
#'
#' @param segs data frame of segments with fields denoting chromosome, start, end, and other metadata (see standardized segs for seg data frame input formats)
#' @param seqlengths seqlengths of output GRanges object
#' @param seqinfo seqinfo of output GRanges object
#' @import GenomicRanges
#' @export
seg2gr = function(segs, key = NULL, seqlengths = hg_seqlengths(), seqinfo = Seqinfo())
  {
    if (is(segs, 'data.table'))
        segs = as.data.frame(segs)

    if (!inherits(seqinfo, 'Seqinfo'))
      seqinfo = seqinfo(seqinfo);

    if (is.null(seqlengths))
      seqlengths = seqlengths(seqinfo)

    if (!inherits(segs, 'GRanges'))
        {
            if (nrow(segs)==0)
                return(GRanges(seqlengths = seqlengths(seqinfo)))
        }
    else if (length(segs) == 0)
        return(GRanges(seqlengths = seqlengths(seqinfo)))
    
    segs = standardize_segs(segs);

    if (any(bix <- (is.na(segs$chr) | is.na(segs$pos1) | is.na(segs$pos2))))
        {
            warning('Segments with NA values for chromosome, start, and end position detected .. removing')
            segs = segs[!bix, ]
        }
    
    GR.NONO.FIELDS = c('seqnames', 'ranges', 'strand', 'seqlevels', 'seqlengths', 'isCircular', 'start', 'end', 'width', 'element');

    if (is.null(segs$strand))
      segs$strand = "*"

    if (any(ix <- !(segs$strand %in% c('+', '-', '*'))))
      segs$strand[ix] = "*"
    
    if (length(seqlengths)>0)
      {
        if (length(wtf  <- setdiff(segs$chr, names(seqlengths))))
          {
            warning('some seqnames in seg object were not included in provided seqlengths: ', paste(wtf, collapse = ','))
            seqlengths[as.character(wtf)] = NA
          }
        segs$pos1 <- as.numeric(segs$pos1)
        segs$pos2 <- as.numeric(segs$pos2)        

        out = GRanges(seqnames = segs$chr, ranges = IRanges(segs$pos1, segs$pos2), names = levels(levels), strand = segs$strand, seqlengths = seqlengths)
      }
    else
      out = GRanges(seqnames = as.character(segs$chr), ranges = IRanges(segs$pos1, segs$pos2), names = levels(levels), strand = segs$strand)    
    
    if (length(seqinfo)>0)
      out = gr.fix(out, seqinfo)
    else if (is.null(seqlengths))
      out = gr.fix(out)
    
    values(out) = segs[, setdiff(names(segs), GR.NONO.FIELDS)]

    if (!is.null(key))
       names(out) = key;
    
    return(out);
  }

#' grdt
#'
#' Converts gr to data frame
#'
#' and a field grl.iix which saves the (local) index that that gr was in its corresponding grl item
#' @param x \code{GRanges} to convert
#' @name grl.unlist
#' @export
grdt = function(x)
 {
      ## new approach just directly instantiating data table
      cmd = 'data.table(';
      if (is(x, 'GRanges'))
          {
              was.gr = TRUE
              f = c('seqnames', 'start', 'end', 'strand', 'width')
              f2 = c('as.character(seqnames', 'c(start', 'c(end', 'as.character(strand', 'as.numeric(width')
              cmd = paste(cmd, paste(f, '=', f2, '(x))', sep = '', collapse = ','), sep = '')
              value.f = names(values(x))              
          }
      else          
          {
              was.gr = FALSE
              value.f = names(x)
          }
      
      if (length(value.f)>0)
          {
              if (was.gr)
                  cmd = paste(cmd, ',', sep = '')
              class.f = sapply(value.f, function(f) eval(parse(text=sprintf("class(x$'%s')", f))))

              .StringSetListAsList = function(x) ### why do I need to do this, bioconductor peeps??
                  {
                      tmp1 = as.character(unlist(x))
                      tmp2 = rep(1:length(x), elementLengths(x))
                      return(split(tmp1, tmp2))                          
                  }

              ## take care of annoying S4 / DataFrame / data.frame (wish-they-were-non-)issues
              as.statement = ifelse(grepl('Integer', class.f), 'as.integer',
                  ifelse(grepl('Character', class.f), 'as.character',
                         ifelse(grepl('StringSetList', class.f), '.StringSetListAsList',
                                ifelse(grepl('StringSet$', class.f), 'as.character',
                                       ifelse(grepl('List', class.f), 'as.list',
                                              ifelse(grepl('List', class.f), 'as.list', 'c'))))))
              cmd = paste(cmd, paste(value.f, '=', as.statement, "(x$'", value.f, "')", sep = '', collapse = ','), sep = '')
          }

      cmd = paste(cmd, ')', sep = '')
      return(eval(parse(text =cmd)))
  }

#'
#' identifies events that are in ra1 that do not exist in ra2.
#' Aside from ra1 and ra2, all arguments are sent to ra.overlaps
#'
#' @name ra.overlaps
#' @export
ra.setdiff <- function(ra1, ra2, ...) {

  ro <- ra.overlaps(ra1, ra2, ...)

  in.ra1.only <- setdiff(seq_along(ra1), ro[, 'ra1.ix'])
  return(ra1[unique(in.ra1.only)])

}

#' ra.union
#'
#' returns events in ra1 that are in ra2 also
#'
#' @name ra.union
#' @export
ra.union <- function(ra1, ra2, ...) 
  return(ra1[unique(ra.overlaps(ra1, ra2, ...)[, 'ra1.ix'])])


.ls.objects <- function (pos = 1, pattern, order.by,
                                                 decreasing=FALSE, head=FALSE, n=5) {
      napply <- function(names, fn) sapply(names, function(x)
                                                                                    fn(get(x, pos = pos)))
          names <- ls(pos = pos, pattern = pattern)
          obj.class <- napply(names, function(x) as.character(class(x))[1])
          obj.mode <- napply(names, mode)
          obj.type <- ifelse(is.na(obj.class), obj.mode, obj.class)
          obj.prettysize <- napply(names, function(x) {
                                       capture.output(print(object.size(x), units = "auto")) })
          obj.size <- napply(names, object.size)
          obj.dim <- t(napply(names, function(x)
                                                      as.numeric(dim(x))[1:2]))
          vec <- is.na(obj.dim)[, 1] & (obj.type != "function")
          obj.dim[vec, 1] <- napply(names, length)[vec]
          out <- data.frame(obj.type, obj.size, obj.prettysize, obj.dim)
          names(out) <- c("Type", "Size", "PrettySize", "Rows", "Columns")
          if (!missing(order.by))
                    out <- out[order(out[[order.by]], decreasing=decreasing), ]
          if (head)
                    out <- head(out, n)
          out
    }

#' gr.all
#'
#' Return a GRanges that holds interavals for all of HG19
#'
#' @param unmap [default F] Optinally add a "unmapped" chr
#' @param M [default F] Include mitochondrial chr
#' @param Y [default T] Include Y chr
#' @return \code{GRanges} object with one element per chromosome
gr.all <- function(unmap=FALSE, M=FALSE, Y=TRUE) {
  gr <- si2gr(gr.tfix(GRanges(1, IRanges(1,1))))

  if (!M)
    gr <- gr[seqnames(gr) != 'M']
  if (!Y)
    gr <- gr[seqnames(gr) != 'Y']


  if (!unmap)
    return(gr[!seqnames(gr) %in% "Unmapped"])
  else
    return(gr)
}


setGeneric('%|%', function(gr, ...) standardGeneric('%|%'))
setMethod("%|%", signature(gr = "GRanges"), function(gr, df) {
    if (is.data.table(df))
        df = as.data.frame(df)
    else if (inherits(df, 'GRanges'))
        df = values(df)    
    values(gr) = cbind(values(gr), df)
    return(gr)
})


#' @name %+%
#' @title Nudge GRanges right
#' @description
#' Operator to shift GRanges right "sh" bases
#' 
#' @return shifted granges
#' @import GenomicRanges
#' @rdname gr.nudge
#' @exportMethod %+%
#' @export
#' @author Marcin Imielinski
setGeneric('%+%', function(gr, ...) standardGeneric('%+%'))
setMethod("%+%", signature(gr = "GRanges"), function(gr, sh) {
    end(gr) = end(gr)+sh
    start(gr) = start(gr)+sh
    return(gr)
})

#' @name %-%
#' @title Shift GRanges left
#' @description
#' Operator to shift GRanges left "sh" bases
#'
#' df %!% c('string.*to.*match', 'another.string.to.match')
#' 
#' @return shifted granges
#' @import GenomicRanges
#' @rdname gr.nudge
#' @exportMethod %-%
#' @export
#' @author Marcin Imielinski
setGeneric('%-%', function(gr, ...) standardGeneric('%-%'))
setMethod("%-%", signature(gr = "GRanges"), function(gr, sh) {
    start(gr) = start(gr)-sh
    end(gr) = end(gr)-sh
    return(gr)
})


#' @name %*%
#' @title gr.findoverlaps (strand agnostic)
#' @description
#' Shortcut for gr.findoverlaps (standard arguments)
#'
#' gr1 %*% gr2
#' 
#' @return new granges containing every pairwise intersection of ranges in gr1 and gr2 with a join of the corresponding metadata (strand agnostic)
#' @rdname grfo
#' @exportMethod %*%
#' @export
#' @author Marcin Imielinski
#setGeneric('%*%', function(x, ...) standardGeneric('%*%'))
setMethod("%*%", signature(x = "GRanges"), function(x, y) {
    gr = gr.findoverlaps(x, y, qcol = names(values(x)), scol = names(values(y)))
    return(gr)
})


#' @name %^%
#' @title gr.in shortcut
#' @description
#' Shortcut for gr.in (standard arguments)
#'
#' gr1 %^% gr2
#' 
#' @return logical vector of length gr1 which is TRUE at entry i only if gr1[i] intersects at least one interval in gr2 (strand agnostic)
#' @rdname gr.in
#' @exportMethod %^%
#' @export 
#' @author Marcin Imielinski
setGeneric('%^%', function(x, ...) standardGeneric('%^%'))
setMethod("%^%", signature(x = "GRanges"), function(x, y) {
    return(gr.in(x, y))
})


#' @name %$% 
#' @title gr.val shortcut to get mean values of subject "x" meta data fields in query "y" (strand agnostic)
#' @description
#' Shortcut for gr.val (using val = names(values(y)))
#'
#' gr1 %$% gr2
#' 
#' @return logical vector of length gr1 which is TRUE at entry i only if gr1[i] intersects at least one interval in gr2
#' @rdname gr.val
#' @exportMethod %$%
#' @export 
#' @author Marcin Imielinski
setGeneric('%$%', function(x, ...) standardGeneric('%$%'))
setMethod("%$%", signature(x = "GRanges"), function(x, y) {
    return(gr.val(x, y, val = names(values(y))))
})


#' @name %&%
#' @title gr.intersect shortcut
#' @description
#' Shortcut for gr.intersect
#'
#' gr1 %&% gr2
#' 
#' @return granges representing intersection of input intervals
#' @rdname gr.intersect
#' @exportMethod %&%
#' @export 
#' @author Marcin Imielinski
setGeneric('%&%', function(x, ...) standardGeneric('%&%'))
setMethod("%&%", signature(x = "GRanges"), function(x, y) {
    return(reduce(gr.findoverlaps(x[, c()], y[, c()])))
})


#' @name %_%
#' @title setdiff shortcut (strand agnostic)
#' @description
#' Shortcut for union
#'
#' gr1 %_% gr2
#' 
#' @return granges representing setdiff of input interval
#' @rdname gr.setdiff
#' @exportMethod %_%
#' @export 
#' @author Marcin Imielinski
setGeneric('%_%', function(x, ...) standardGeneric('%_%'))
setMethod("%_%", signature(x = "GRanges"), function(x, y) {
    setdiff(gr.stripstrand(x[, c()]), gr.stripstrand(y[, c()]))
})


#' @name gr.union.unstranded
#' @title union shortcut
#' @description
#' Shortcut for union
#'
#' gr1 %|% gr2
#' 
#' @return granges representing setdiff of input interval
#' @rdname gr.union
#' @exportMethod %|%
#' @export 
#' @author Marcin Imielinski
setGeneric('%|%', function(x, ...) standardGeneric('%|%'))
setMethod("%|%", signature(x = "GRanges"), function(x, y) {
    return(reduce(grbind(gr.stripstrand(x[, c()]), gr.stripstrand(y[, c()]))))
})

#' @name %**%
#' @title gr.findoverlaps (respects strand)
#' @description
#' Shortcut for gr.findoverlaps
#'
#' gr1 %**% gr2
#' 
#' @return new granges containing every pairwise intersection of ranges in gr1 and gr2 with a join of the corresponding  metadata
#' @rdname grfo
#' @exportMethod %**%
#' @export
#' @author Marcin Imielinski
setGeneric('%**%', function(x, ...) standardGeneric('%**%'))
setMethod("%**%", signature(x = "GRanges"), function(x, y) {
    gr = gr.findoverlaps(x, y, qcol = names(values(x)), scol = names(values(y)), ignore.strand = FALSE)
    return(gr)
})


#' @name %^^%
#' @title gr.in shortcut (respects strand)
#' @description
#' Shortcut for gr.in 
#'
#' gr1 %^^% gr2
#' 
#' @return logical vector of length gr1 which is TRUE at entry i only if gr1[i] intersects at least one interval in gr2
#' @rdname gr.in
#' @exportMethod %^^%
#' @export 
#' @author Marcin Imielinski
setGeneric('%^^%', function(x, ...) standardGeneric('%^^%'))
setMethod("%^^%", signature(x = "GRanges"), function(x, y) {
    return(gr.in(x, y, ignore.strand = FALSE))
})


#' @name %$$%
#' @title gr.val shortcut to get mean values of subject "x" meta data fields in query "y" (respects strand)
#' @description
#' Shortcut for gr.val (using val = names(values(y)))
#'
#' gr1 %$$% gr2
#' 
#' @return logical vector of length gr1 which is TRUE at entry i only if gr1[i] intersects at least one interval in gr2
#' @rdname gr.val
#' @exportMethod %$$%
#' @export 
#' @author Marcin Imielinski
setGeneric('%$$%', function(x, ...) standardGeneric('%$$%'))
setMethod("%$$%", signature(x = "GRanges"), function(x, y) {
    return(gr.val(x, y, val = names(values(y)), ignore.strand = FALSE))
})


#' @name %&&%
#' @title gr.intersect shortcut (respects strand)
#' @description
#' Shortcut for gr.intersect
#'
#' gr1 %&&% gr2
#' 
#' @return granges representing intersection of input intervals
#' @rdname gr.intersect
#' @exportMethod %&&%
#' @export 
#' @author Marcin Imielinski
setGeneric('%&&%', function(x, ...) standardGeneric('%&&%'))
setMethod("%&&%", signature(x = "GRanges"), function(x, y) {
    return(reduce(gr.findoverlaps(x[, c()], y[, c()], ignore.strand = FALSE)))
})


#' @name %__% 
#' @title setdiff shortcut (respects strand)
#' @description
#' Shortcut for setdiff
#'
#' gr1 %__% gr2
#' 
#' @return granges representing setdiff of input interval
#' @rdname gr.setdiff
#' @exportMethod %__%
#' @export 
#' @author Marcin Imielinski
setGeneric('%__%', function(x, ...) standardGeneric('%__%'))
setMethod("%__%", signature(x = "GRanges"), function(x, y) {
    setdiff(x[, c()], y[, c()])
})


#' @name gr.union.stranded
#' @title setdiff shortcut (respects strand)
#' @description
#' Shortcut for setdiff
#'
#' gr1 %||% gr2
#' 
#' @return granges representing setdiff of input interval
#' @rdname gr.union
#' @exportMethod %||%
#' @export 
#' @author Marcin Imielinski
setGeneric('%||%', function(x, ...) standardGeneric('%||%'))
setMethod("%||%", signature(x = "GRanges"), function(x, y) {
    return(reduce(grbind(x[, c()], y[, c()])))
})

.toggle_grfo = function()
    {
        old.val = as.logical(Sys.getenv('GRFO_FOVERLAPS'))
        if (is.na(old.val))
            old.val = FALSE            
        Sys.setenv(GRFO_FOVERLAPS = !old.val)
        cat('GRFO_FOVERLAPS is', Sys.getenv('GRFO_FOVERLAPS'), '\n\t...Default gr.findoverlaps behavior will use',
            ifelse(Sys.getenv('GRFO_FOVERLAPS'), 'data.table foverlaps', 'IRanges findOverlaps'), '\n')
    }

