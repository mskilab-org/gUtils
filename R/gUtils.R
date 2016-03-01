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

#' Converts \code{\link[GenomicRanges]{GRanges}} to \code{\link[base]{data.frame}}
#' @import data.table
#' @import GenomicRanges
#' @name grdt
grdt = function(gr)
  {
    return(as.data.table(as.data.frame(gr)))
  }

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
#' gr.start(GRanges(1, IRanges(start=c(1,2,3), width=101), seqinfo=Seqinfo("1", 3)), width=20)
#' gr.start(GRanges(1, IRanges(start=c(1,2,3), width=101), seqinfo=Seqinfo("1", 3)), width=200, clip=T)
#' @export
gr.start = function(x, width = 1, force = F, ignore.strand = T, clip = F)
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

###############
#' Convert data.table to GRanges
#'
#' Takes as input a data.table which must have the fields: start, end, strand, seqnames.
#' All of the remaining fields are added as meta data to the GRanges
#' @param dt data.table to convert to GRanges
#' @return GRanges object of length = nrow(dt)
#' @import data.table
#' @import GenomicRanges
#' @examples
#' gr <- dtgr(data.table(start=c(1,2), seqnames=c("X", "1"), end=c(10,20), strand = c('+', '-')))
#' @export
dtgr <- function(dt) {
  library(data.table)
  
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

#' Get GRanges corresponding to beginning of end
#'
#' Alternative to \code{flank} that will provide end positions *within* intervals
#' 
#' @param x \code{GRanges} object to operate on
#' @param width [default = 1] Specify subranges of greater width including the start of the range.
#' @param force [default = F] Allows returned \code{GRanges} to have ranges outside of its \code{Seqinfo} bounds.
#' @param clip [default = F] Trims returned \code{GRanges} so that it does not extend beyond bounds of the input \code{GRanges}
#' @param ignore.strand [default = T] If set to \code{FALSE}, will extend '-' strands from the other direction.
#' @return GRanges object of width 1 ranges representing end of each genomic range in the input.
#' @examples
#' gr.end(GRanges(1, IRanges(start=c(1,2,3), width=101), seqinfo=Seqinfo("1", 3)), width=20)
#' gr.end(GRanges(1, IRanges(start=c(1,2,3), width=101), seqinfo=Seqinfo("1", 3)), width=200, clip=T)
#' @export
gr.end = function(x, width = 1, force = F, ignore.strand = T, clip = T)
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
#' @export
#' @examples
#' ##gr.mid(GRanges(1, IRanges(1000,2000), seqinfo=Seqinfo("1", 2000)), pad=10)
gr.mid = function(x, pad = 0)
  {
      start(x) = end(x) = rowMeans(cbind(start(x), end(x)))
      if (!is.null(pad))
          if (!is.na(pad))
              if (pad>0)
                  x = x + pad

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
#' @param up [default T] See description.
#' @param parallel [default F] If \code{TRUE}, assumes Q and S are same length and this analysis is only performed between the corresponding Q and S pairs. 
#' @return Rounded \code{GRanges}
#' @examples
#' query   <- GRanges(1, IRanges(c(100,110),width=201), seqinfo=Seqinfo("1", 500))
#' subject <- GRanges(1, IRanges(c(160,170),width=201), seqinfo=Seqinfo("1", 500))
#' gr.round(query, subject)
#' @export

gr.round = function(Q, S, up = T, parallel = F)
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


##################################################
#' Generate random GRanges on genome
#'
#' Randomly generates non-overlapping GRanges with supplied widths on supplied genome
#'
#' @param w Vector of widths (length of w determines length of output)
#' @param genome Genome which can be a \code{GRanges}, \code{GRangesList}, or \code{Seqinfo} object. Default is "hg19" from the \code{\link{BSGenome}} package.
#' @param seed [default NA] Optionally specify a seed for the RNG. Defualt behavior is random seed.
#' @return \code{GRanges} with random intervals on the specifed "chromosomes"
#' @note This function is currently quite slow, needs optimization
#' @examples
#' ## Generate a single random interval of width 10, on "chr" of length 1000
#' gr.rand(10, Seqinfo("1", 1000))
#' ## Generate 5 non-overlapping regions of width 10 on hg19
#' gr.rand(rep(10,5))
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


#' Trim a pile of GRanges
#'
#' trims pile of granges relative to the specified <local> coordinates of each range
#' (ie the first coordinate of every gr is 1 and the last is width(gr))
#' if end is larger than the width of the corresponding gr, then the corresponding output will only have end(gr) as its coordinate.
#'
#' this is a role not currently provided by the standard granges fxns
#' (eg shift, reduce, restrict, shift, resize, flank etc)
#' @param gr GRanges to trime
#' @param starts number. Default 1
#' @param ends number. Default 1
#' @export
gr.trim = function(gr, starts=1, ends=1, fromEnd=FALSE, ignore.strand = T)
    {
    starts = cbind(1:length(gr), starts)[, 2]
    ends = cbind(1:length(gr), ends)[, 2]
    
    if (!ignore.strand)
        {
        ix = as.logical(strand(gr)=='-')
        if (any(ix))
            {
            if (fromEnd)
                {
                tmp = starts[ix]
                starts[ix] = ends[ix]
                ends[ix] = tmp-1
            }
            else
                {
                starts[ix] = width(gr)[ix]-starts[ix]+1
                ends[ix] = width(gr)[ix]-ends[ix]+1
            }
        }
    }
      
    if (fromEnd) {
      en = pmax(starts, end(gr)-ends);
  } else {
      ends = pmax(starts, ends);
      ends = pmin(ends, width(gr));
      en = start(gr) + ends - 1;
  }
    
    st = start(gr)+starts-1;
    st = pmin(st, en);
            
    out = GRanges(seqnames(gr), IRanges(st, en),
                   seqlengths = seqlengths(gr), strand = strand(gr))
    values(out) = values(gr)
    return(out)
}

##########################
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
gr.sample = function(gr, k, len = 100, replace = T)
{
  if (!inherits(gr, 'GRanges'))
    gr = seqinfo2gr(gr)
  
  if (length(k)==1)
    {
      gr.f = gr.flatten(gr.trim(gr, starts = 1, ends = width(gr)-len), gap = 0);
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
#' si <- Seqinfo(names(hg_seqlength(), hg_seqlengths()))
#' si2gr(si)
#' @export
si2gr <- seqinfo2gr <- function(si, strip.empty = F)
#' \dontrun{si <- Seqinfo(names(hg_seqlength(), hg_seqlengths()))
#' seqinfo2gr(si)}
#' @export
seqinfo2gr <- function(si, strip.empty = FALSE)
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


#' Shortcut GRanges from Seqinfo
#'
#' Creates a genomic ranges from seqinfo object
#' ie a pile of ranges spanning the genome
#' @param si Seqinfo object
#' @param strip.empty Don't know. Default FALSE
#' @examples
#' \dontrun{si <- Seqinfo(names(hg_seqlength(), hg_seqlengths()))
#' si2gr(si)}
#' @export
si2gr <- function(...) seqinfo2gr(...)


#' rle.query
#'
#' Queries an \code{\link{RleList}} representing genomic data (ie a list whose names represent
#' seqnames ie chromosomes, and lengths represent seqlengths)
#' via \code{GRanges} object
#'
#' @return Rle representing the (concatenated) vector of data (reversing order in case of negative strand input)
#' @note Throws warning if seqlengths(gr) do not correspond to the lengths of the \code{RleList} components
#' @export
rle.query = function(subject.rle, query.gr, verbose = F, mc.cores = 1, chunksize = 1e9) ## mc.cores only relevant if there are over 1e9 bases to be queried from subject.rle
  {
    was.grl = F
    
    if (is(query.gr, 'GRangesList'))
    {
      was.grl = T        
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
              require(data.table)
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

    force.rrbind = F
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
#' @return Concatenated GRangesList 
#' @examples
#' ## Create some dummy data
#' gr1 <- GRanges(1, IRanges(100,1000), my.gene='X', seqinfo=Seqinfo("1", 1000))
#' gr2 <- GRanges(1, IRanges(200,2000), my.annot='Y', seqinfo=Seqinfo("1",2000))
#' gr3 <- GRanges(2, IRanges(1,3000), my.annot='Z', seqinfo=Seqinfo("2",3000))
#' grl1 <- GRangesList(grbind(gr1, gr2))
#' grl2 <- GRangesList(grbind(gr1, gr3))
#' ## Add unique annotation to just one \code{mcols}.
#' mcols(gr1)$my.new.annot=1
#' ## Concatenate
#' grlbind(grl1, grl2)
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
gr2gatk = function(gr, file, add.chr = F)
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
streduce = function(gr, pad = 0, sort = T)
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
#' @param include.val scalar logical, will include in out gr values field of first matching record in input gr
#' @return Simplified GRanges with "field" populated with uniquely contiguous values
#' @export
gr.simplify = function(gr, field = NULL, val = NULL, include.val = T, split = F, pad = 1)
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
#' @name gr.string
#' @export
gr.string = function(gr, add.chr = FALSE, mb = FALSE, round = 3, other.cols = c(), pretty = FALSE)
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
        if (pretty)
            return(paste(sn, ':', str_trim(prettyNum(start(gr), big.mark = ',')), '-', str_trim(prettyNum(end(gr), big.mark = ',')), str, other.str, sep = ''))
        else
            return(paste(sn, ':', start(gr), '-', end(gr), str, other.str, sep = ''))
  }


#' grl.string
#'
#' return ucsc style interval string corresponding to each gr in each grl, one line per grl item, with gr's in 
#' grl separated by sep
#'
#' if mb will return as MB and round to "round"
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
    if (any(ix <- sapply(spl, length)==1))
        spl[ix] = strsplit(gr.string(seqinfo2gr(seqlengths)[sapply(spl[ix], function(x) x[[1]])], mb = F), '[\\:\\-\\+]', perl = T)
    else if (any(ix <- sapply(spl, length)==2))
        spl[ix] = lapply(spl, function(x) x[c(1,2,2)])
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
#' @name parse.gr
#' @export
parse.gr = function(...)
  {
    return(unlist(parse.grl(...)))
  }


#' gstring
#'
#' quick function to parse gr from character vector IGV / UCSC style strings of format gr1;gr2;gr3 where each gr is of format chr:start-end[+/-]
#'
#' @name gstring
#' @export
gstring = function(...)
  {
    return(unlist(parse.grl(...)))
  }

#################
#' ra_breaks
#'
#' takes in either file or data frame from dranger or snowman or path to BND / SV type vcf file
#' and returns junctions in VCF format.
#' 
#' The default output is GRangesList each with a length two GRanges whose strands point AWAY from the break.  If get.loose = TRUE (only relevant for VCF)
#'
#' @name ra_breaks
#' @export
ra_breaks = function(rafile, keep.features = T, seqlengths = hg_seqlengths(), chr.convert = T, snowman = FALSE,  breakpointer = FALSE, seqlevels = NULL,
    get.loose = FALSE ## if TRUE will return a list with fields $junctions and $loose.ends
  )
  {
      if (is.character(rafile))
          {
              if (grepl('(vcf$)|(vcf.gz$)', rafile))
                  {
                      library(VariantAnnotation)
                      vcf = readVcf(rafile, Seqinfo(seqnames = names(seqlengths), seqlengths = seqlengths))
                      if (!('SVTYPE' %in% names(info(vcf)))) {
                        warning('Vcf not in proper format.  Is this a rearrangement vcf?')
                          return(GRangesList());
                        }
                      
                      vgr = rowData(vcf) ## parse BND format

                      ## no events
                      if (length(vgr) == 0)
                        return (GRangesList())

                      ## fix mateids if not included
                      if (!"MATEID"%in%colnames(mcols(vgr))) {
                        nm <- vgr$MATEID <- names(vgr)
                        ix <- grepl("1$",nm)
                        vgr$MATEID[ix] = gsub("(.*?)(1)$", "\\12", nm[ix])
                        vgr$MATEID[!ix] = gsub("(.*?)(2)$", "\\11", nm[!ix])
                        vgr$SVTYPE="BND"
                      }
                      
                      if (!any(c("MATEID", "SVTYPE") %in% colnames(mcols(vgr))))
                        stop("MATEID or SVTYPE not included. Required")
                      vgr$mateid = info(vcf)$MATEID
                      vgr$svtype = info(vcf)$SVTYPE

                      if (!is.null(info(vcf)$SCTG))
                          vgr$SCTG = info(vcf)$SCTG
                      
                      if (sum(vgr$svtype == 'BND')==0)
                          stop('Vcf not in proper format.  Will only process rearrangements in BND format')

                      if (!all(vgr$svtype == 'BND'))
                          warning(sprintf('%s rows of vcf do not have svtype BND, ignoring these'), sum(vgr$svtype != 'BND'))

                      bix = which(vgr$svtype == "BND")
                      vgr = vgr[bix]
                      vgr$first = !grepl('^(\\]|\\[)', vgr$ALT) ## ? is this row the "first breakend" in the ALT string (i.e. does the ALT string not begin with a bracket)
                      vgr$right = grepl('\\[', vgr$ALT) ## ? are the (sharp ends) of the brackets facing right or left
                      vgr$coord = as.character(paste(seqnames(vgr), ':', start(vgr), sep = ''))
                      vgr$mcoord = as.character(gsub('.*(\\[|\\])(.*\\:.*)(\\[|\\]).*', '\\2', vgr$ALT))
                      vgr$mcoord = gsub('chr', '', vgr$mcoord)

                      if (all(is.na(vgr$mateid)))
                          if (!is.null(names(vgr)) & !any(duplicated(names(vgr))))
                              {
                                  warning('MATEID tag missing, guessing BND partner by parsing names of vgr')
                                  vgr$mateid = paste(gsub('::\\d$', '', names(vgr)), (sapply(strsplit(names(vgr), '\\:\\:'), function(x) as.numeric(x[length(x)])))%%2 + 1, sep = '::')
                              }
                          else if (!is.null(vgr$SCTG))
                              {
                                  warning('MATEID tag missing, guessing BND partner from coordinates and SCTG')
                                  require(igraph)
                                  ucoord = unique(c(vgr$coord, vgr$mcoord))
                                  vgr$mateid = paste(vgr$SCTG, vgr$mcoord, sep = '_')

                                  if (any(duplicated(vgr$mateid)))
                                      {
                                          warning('DOUBLE WARNING! inferred mateids not unique, check VCF')
                                          bix = bix[!duplicated(vgr$mateid)]
                                          vgr = vgr[!duplicated(vgr$mateid)]
                                      }
                              }
                          else
                              stop('MATEID tag missing')

                      vgr$mix = as.numeric(match(vgr$mateid, names(vgr)))

                      pix = which(!is.na(vgr$mix))

                      vgr.pair = vgr[pix]

                      if (length(vgr.pair)==0)
                          stop('No mates found despite nonzero number of BND rows in VCF')
                      vgr.pair$mix = match(vgr.pair$mix, pix)
                      vix = which(1:length(vgr.pair)<vgr.pair$mix )
                      vgr.pair1 = vgr.pair[vix]
                      vgr.pair2 = vgr.pair[vgr.pair1$mix]

                      ## now need to reorient pairs so that the breakend strands are pointing away from the breakpoint
                      
                      ## if "first" and "right" then we set this entry "-" and the second entry "+"
                      tmpix = vgr.pair1$first & vgr.pair1$right
                      if (any(tmpix))
                          {
                              strand(vgr.pair1)[tmpix] = '-'
                              strand(vgr.pair2)[tmpix] = '+'
                          }

                      ## if "first" and "left" then "-", "-"
                      tmpix = vgr.pair1$first & !vgr.pair1$right
                      if (any(tmpix))
                          {
                              strand(vgr.pair1)[tmpix] = '-'
                              strand(vgr.pair2)[tmpix] = '-'
                          }

                      ## if "second" and "left" then "+", "-"
                      tmpix = !vgr.pair1$first & !vgr.pair1$right
                      if (any(tmpix))
                          {
                              strand(vgr.pair1)[tmpix] = '+'
                              strand(vgr.pair2)[tmpix] = '-'
                          }

                      ## if "second" and "right" then "+", "+"
                      tmpix = !vgr.pair1$first & vgr.pair1$right
                      if (any(tmpix))
                          {
                              strand(vgr.pair1)[tmpix] = '+'
                              strand(vgr.pair2)[tmpix] = '+'
                          }

                      pos1 = as.logical(strand(vgr.pair1)=='+') ## positive strand junctions shift left by one (ie so that they refer to the base preceding the break for these junctions
                      if (any(pos1))
                          {
                              start(vgr.pair1)[pos1] = start(vgr.pair1)[pos1]-1
                              end(vgr.pair1)[pos1] = end(vgr.pair1)[pos1]-1
                          }

                      pos2 = as.logical(strand(vgr.pair2)=='+') ## positive strand junctions shift left by one (ie so that they refer to the base preceding the break for these junctions
                      if (any(pos2))
                          {
                              start(vgr.pair2)[pos2] = start(vgr.pair2)[pos2]-1
                              end(vgr.pair2)[pos2] = end(vgr.pair2)[pos2]-1
                          }
                      ra = grl.pivot(GRangesList(vgr.pair1[, c()], vgr.pair2[, c()]))

                      this.inf = info(vcf)[bix[pix[vix]], ]

                      if (is.null(this.inf$POS))
                          this.inf = cbind(data.frame(POS = ''), this.inf)
                      if (is.null(this.inf$CHROM))
                          this.inf = cbind(data.frame(CHROM = ''), this.inf)

                      if (is.null(this.inf$MATL))
                          this.inf = cbind(data.frame(MALT = ''), this.inf)
                      
                      this.inf$CHROM = seqnames(vgr.pair1)
                      this.inf$POS = start(vgr.pair1)
                      this.inf$MATECHROM = seqnames(vgr.pair2)
                      this.inf$MATEPOS = start(vgr.pair2)
                      this.inf$MALT = vgr.pair2$ALT
                      
                      values(ra) = cbind(fixed(vcf)[bix[pix[vix]],], this.inf)
                      
                      if (is.null(values(ra)$TIER))
                          values(ra)$tier = ifelse(values(ra)$FILTER == "PASS", 2, 3) ## baseline tiering of PASS vs non PASS variants
                      else
                          values(ra)$tier = values(ra)$TIER

                      if (!get.loose)
                          return(ra)
                      else
                          {
                              npix = is.na(vgr$mix)
                              vgr.loose = vgr[npix, c()] ## these are possible "loose ends" that we will add to the segmentation
                              values(vgr.loose) = cbind(fixed(vcf)[bix[npix], ], info(vcf)[bix[npix], ])

                              return(list(junctions = ra, loose.ends = vgr.loose))
                          }
                  }
              else
                  rafile = read.delim(rafile)
          }
            
     if (is.data.table(rafile))
         {
             require(data.table)
             rafile = as.data.frame(rafile)
         }

    if (nrow(rafile)==0)
        {
            out = GRangesList()
            values(out) = rafile
            return(out)
        }
  
    if (snowman) ## flip breaks so that they are pointing away from junction
      {
        rafile$str1 = ifelse(rafile$strand1 == '+', '-', '+')
        rafile$str2 = ifelse(rafile$strand2 == '+', '-', '+')
      }
      
    if (!is.null(seqlevels)) ## convert seqlevels from notation in tab delim file to actual
      {
        rafile$chr1 = seqlevels[rafile$chr1]
        rafile$chr2 = seqlevels[rafile$chr2]        
      }

     
    if (is.null(rafile$str1))
      rafile$str1 = rafile$strand1

    if (is.null(rafile$str2))
      rafile$str2 = rafile$strand2
     if (!is.null(rafile$pos1) & !is.null(rafile$pos2))
         {
             if (breakpointer)
                 {
                     rafile$pos1 = rafile$T_BPpos1
                     rafile$pos2 = rafile$T_BPpos2
                 }
             
             if (!is.numeric(rafile$pos1))
                 rafile$pos1 = as.numeric(rafile$pos1)

             if (!is.numeric(rafile$pos2))
                 rafile$pos2 = as.numeric(rafile$pos2)

             ## clean the parenthesis from the string

             rafile$str1 <- gsub('[()]', '', rafile$str1)
             rafile$str2 <- gsub('[()]', '', rafile$str2)

             if (is.character(rafile$str1) | is.factor(rafile$str1))
                 rafile$str1 = gsub('0', '-', gsub('1', '+', rafile$str1))
             
             if (is.character(rafile$str2) | is.factor(rafile$str2))
                 rafile$str2 = gsub('0', '-', gsub('1', '+', rafile$str2))
             
             if (is.numeric(rafile$str1))
                 rafile$str1 = ifelse(rafile$str1>0, '+', '-')

             if (is.numeric(rafile$str2))
                 rafile$str2 = ifelse(rafile$str2>0, '+', '-')
             
             rafile$rowid = 1:nrow(rafile)

             bad.ix = is.na(rafile$chr1) | is.na(rafile$chr2) | is.na(rafile$pos1) | is.na(rafile$pos2) | is.na(rafile$str1) | is.na(rafile$str2) | rafile$str1 == '*'| rafile$str2 == '*' | rafile$pos1<0 | rafile$pos2<0
             
             rafile = rafile[which(!bad.ix), ]
             
             if (nrow(rafile)==0)
                 return(GRanges())
             
             seg = rbind(data.frame(chr = rafile$chr1, pos1 = rafile$pos1, pos2 = rafile$pos1, strand = rafile$str1, ra.index = rafile$rowid, ra.which = 1, stringsAsFactors = F),
                 data.frame(chr = rafile$chr2, pos1 = rafile$pos2, pos2 = rafile$pos2, strand = rafile$str2, ra.index = rafile$rowid, ra.which = 2, stringsAsFactors = F))

             if (chr.convert)
                 seg$chr = gsub('25', 'M', gsub('24', 'Y', gsub('23', 'X', seg$chr)))
             
             out = seg2gr(seg, seqlengths = seqlengths)[, c('ra.index', 'ra.which')];
             out = split(out, out$ra.index)
         }
     else if (!is.null(rafile$start1) & !is.null(rafile$start2) & !is.null(rafile$end1) & !is.null(rafile$end2))
         {
             ra1 = gr.flip(GRanges(rafile$chr1, IRanges(rafile$start1, rafile$end1), strand = rafile$str1))
             ra2 = gr.flip(GRanges(rafile$chr2, IRanges(rafile$start2, rafile$end2), strand = rafile$str2))
             out = grl.pivot(GRangesList(ra1, ra2))             
         }
     
     
     if (keep.features)
         values(out) = rafile[, ]

     return(out)
 }

#' Trim a pile of GRanges
#'
#' trims pile of granges relative to the specified <local> coordinates of each range
#' (ie the first coordinate of every gr is 1 and the last is width(gr))
#' if end is larger than the width of the corresponding gr, then the corresponding output will only have end(gr) as its coordinate.
#'
#' this is a role not currently provided by the standard granges fxns
#' (eg shift, reduce, restrict, shift, resize, flank etc)
#' @param gr GRanges to trime
#' @param starts number. Default 1
#' @param ends number. Default 1
#' @export
gr.trim = function(gr, starts=1, ends=1, fromEnd=FALSE, ignore.strand = T)
    {
    starts = cbind(1:length(gr), starts)[, 2]
    ends = cbind(1:length(gr), ends)[, 2]
    
    if (!ignore.strand)
        {
        ix = as.logical(strand(gr)=='-')
        if (any(ix))
            {
            if (fromEnd)
                {
                tmp = starts[ix]
                starts[ix] = ends[ix]
                ends[ix] = tmp-1
            }
            else
                {
                starts[ix] = width(gr)[ix]-starts[ix]+1
                ends[ix] = width(gr)[ix]-ends[ix]+1
            }
        }
    }
      
    if (fromEnd) {
      en = pmax(starts, end(gr)-ends);
  } else {
      ends = pmax(starts, ends);
      ends = pmin(ends, width(gr));
      en = start(gr) + ends - 1;
  }
    
    st = start(gr)+starts-1;
    st = pmin(st, en);
            
    out = GRanges(seqnames(gr), IRanges(st, en),
                   seqlengths = seqlengths(gr), strand = strand(gr))
    values(out) = values(gr)
    return(out)
}


#' ra.merge
#'
#' merges rearrangements input as GRangesList
#' 
ra.merge = function(..., pad = 0, ignore.strand = FALSE)
    {
        ra = list(...)
        nm = names(ra)
        if (is.null(nm))
            nm = paste('ra', 1:length(ra), sep = '')
        nm = paste('seen.by', nm, sep = '.')        
        if (length(nm)==0)
            return(NULL)        
        out = ra[[1]]
        values(out) = cbind(as.data.frame(matrix(FALSE, nrow = length(out), ncol = length(nm), dimnames = list(NULL, nm))), values(out))
        values(out)[, nm[1]] = TRUE
        if (length(ra)>1)
            {
                for (i in 2:length(ra))
                    {                
                        this.ra = ra[[i]]                        
                        if (length(this.ra)>0)
                            {
                                values(this.ra) = cbind(as.data.frame(matrix(FALSE, nrow = length(this.ra), ncol = length(nm), dimnames = list(NULL, nm))), values(this.ra))
                                ovix = ra.overlaps(out, this.ra, pad = pad, ignore.strand = ignore.strand)
                                values(this.ra)[[nm[i]]] = TRUE
                                if (!all(is.na(ovix)))
                                    values(out)[, nm[i]][ovix[,1]] = TRUE 
                                if (!all(is.na(ovix)))
                                    nix = setdiff(1:length(this.ra), ovix[,2])
                                else
                                    nix = 1:length(this.ra)                
                                if (length(nix)>0)    
                                    {    
                                        val1 = values(out)
                                        val2 = values(this.ra)
                                        values(out) = NULL
                                        values(this.ra) = NULL
                                        out = grlbind(out, this.ra[nix])    
                                        values(out) = rrbind2(val1, val2[nix, ])                                                                            
                                    }                            
                            }
                    }
            }
        return(out)
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
#' @name ra.overlaps
#' @export
ra.overlaps = function(ra1, ra2, pad = 0, arr.ind = T, ignore.strand=FALSE, ...)
  {    
    bp1 = grl.unlist(ra1) + pad
    bp2 = grl.unlist(ra2) + pad 
    ix = gr.findoverlaps(bp1, bp2, ignore.strand = ignore.strand, ...)

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
#' @name gr.fix
#' @export
gr.fix = function(gr, genome = NULL, gname = NULL,  drop = F)
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
                require(data.table)
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
#' @name gr.refactor
#' @export
gr.refactor = function(gr, sn, gap = 0, rev = F)
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
#' @name gr.flip
#' @export
gr.flip = function(gr, which = T)
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
#' @export
gr.tile = function(gr, w = 1e3)
  {
    if (!is(gr, 'GRanges'))
      gr = seqinfo2gr(gr);

    ix = which(width(gr)>0)
    gr = gr[ix] 

    if (length(gr)==0)
        return(gr[c()][, c()])   

    ws = as.numeric(ceiling((width(gr))/w))

    st = rep(start(gr), ws)
    en = rep(end(gr), ws)
    strand = rep(as.character(strand(gr)), ws)
    dt = data.table(query.id = rep(1:length(gr), ws), tile.id = 1:sum(ws))
    dt[, numw := 0:(length(tile.id)-1), by = query.id]   
    start = pmin(st+w*dt$numw, en)
    end = pmin(st+w*(dt$numw+1)-1, en)
    
    out = GRanges(rep(as.character(seqnames(gr)), ws), IRanges(start, end), strand = strand, query.id = ix[dt$query.id], tile.id = 1:length(start), seqinfo = seqinfo(gr))    
    
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


#' gr.findoverlaps
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
#' ... = additional args for findOverlaps (IRanges version)
#' @name gr.findoverlaps
#' @export 
gr.findoverlaps = function(query, subject, ignore.strand = T, first = F,
    qcol = NULL, ## any query meta data columns to add to result
    scol = NULL, ## any subject meta data columns to add to resultx
    max.chunk = 1e13,
    foverlaps = ifelse(is.na(as.logical(Sys.getenv('GRFO_FOVERLAPS'))), FALSE, as.logical(Sys.getenv('GRFO_FOVERLAPS'))) & exists('foverlaps'),
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
                  require(data.table)
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


#' gr.match
#' More robust and faster implementation of GenomicRangs::setdiff
#'
#' Robust to common edge cases of setdiff(gr1, gr2)  where gr2 ranges are contained inside gr1's (yieldings
#' setdiffs yield two output ranges for some of the input gr1 intervals.  
#' 
#' @param query \code{GRanges} object as query
#' @param subject \code{GRanges} object as subject
#' @param max.slice Default Inf. If query is bigger than this, chunk into smaller on different cores
#' @param verbose Default FALSE
#' @param mc.cores Default 1. Only works if exceeded max.slice
#' @param ... arguments to be passed to \link{gr.findoverlaps}
#' @return returns indices of query in subject or NA if none found
#' @name gr.match
#' @export
gr.setdiff = function(query, subject, ignore.strand = TRUE, by = NULL,  ...)
    {
        if (!is.null(by)) ## in this case need to be careful about setdiffing only within the "by" level
            {
                tmp = grdt(subject)
                tmp$strand = factor(tmp$strand, c('+', '-', '*'))
                sl = seqlengths(subject)
                gp = seg2gr(tmp[, as.data.frame(gaps(IRanges(start, end), 1, sl[seqnames][1])), by = c('seqnames', 'strand', by)], seqinfo = seqinfo(subject))                
            }
        else ## otherwise easier
            {
                if (ignore.strand)
                    gp = gaps(gr.stripstrand(subject)) %Q% (strand == '*')
                else
                    gp = gaps(subject)
            }
        
        out = gr.findoverlaps(query, gp, qcol = names(values(query)), ignore.strand = ignore.strand, by = by, ...)        
        return(out)        
    }

#' Faster implementatin of GenomicRangs::match 
#'
#' Faster implementation of GRanges match (uses gr.findoverlaps)
#' returns indices of query in subject or NA if none found
#' ... = additional args for findOverlaps (IRanges version)
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
#' @export
gr.tile.map = function(query, subject, mc.cores = 1, verbose = F)
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


#' gr.in
#'
#' faster implementation of GRanges %over%  (uses gr.findoverlaps)
#'
#' returns T / F vector if query range i is found in any range in subject
#'
#' by = column name in query and subject that we additionally control for a match (passed on to gr.findoverlaps)
#' @name gr.in
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
#
gr.duplicated = function(query, by = NULL, type = 'any')
    {        
        return(duplicated(gr.match(query, query, by = by , type = type)))
    }

#' gr.collapse
#'
#' like "reduce" except only collapses <<adjacent>> ranges in the input
#' returning the collapsed ranges
#' @name gr.collapse
#' @export
gr.collapse = function(gr, pad = 1)
  {
    tmp = gr.findoverlaps(gr + pad, gr + pad, ignore.strand = F)
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
#' @param mean scalar logical flag if FALSE then will return sum instead of mean, only applies if target "val" column i snumeric
#' @param sep scalar character, specifies character to use as separator when aggregating character "vals" from target, only applies if target is numeric
#' @param by scalar character, specifies additional "by" column of query AND target that will be used to match up query and target pairs (i.e. in addition to pure GRanges overlap), default is NULL
#' @export
gr.val = function(query, target, val = NULL,
    mean = T, # if false then will return (weighted) <sum> instead of <mean>, only applies if target is numeric
    weighted = mean, # if false will return unweighted sum / mean
    na.rm = F, # only applies if val column of target is numeric
    by = NULL,
    by.prefix = val,
    merge = FALSE, # if merge = FALSE then will cross every range in query with every level of "by" in target (and create data matrix), otherwise will assume query has "by" and merge only ranges that have matching "by" values in both query and target
    verbose = FALSE,
    FUN = NULL, ## takes two  arguments (value, na.rm = TRUE) if weighted = FALSE, and three (value, width, na.rm = TRUE) if weighted = TRUE
    ignore.strand = T,
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

#' gr.dice
#'
#' "dice" up gr into width 1 ranges spanning the input (warning can produce a very large object)
#'
#' if grl = T, then input gr structure will be reflected in the output grl structure
#' @name gr.dice
#' @export
gr.dice = function(gr, grl = F)
  {
    sn = Rle(as.character(seqnames(gr)), width(gr))
    tmp.st = start(gr)
    tmp.en = end(gr)    
    st = unlist(lapply(1:length(gr), function(x) tmp.st[x]:tmp.en[x]))
    str = Rle(as.character(strand(gr)), width(gr))

    ## reverse for negative strand ranges
    if (any(ix <- as.logical(strand(gr) == '-')))
      {
        w = width(gr)
        q.id = unlist(lapply(1:length(gr), function(x) rep(x, w[x])))
        q.l = split(1:length(st), q.id)
        
        for (j in q.l[ix])
          st[j] = rev(st[j])
      }
                
    ix = Rle(1:length(gr), width(gr))        
    out = GRanges(sn, IRanges(st, st), strand = str, seqlengths = seqlengths(gr))
    values(out) = values(gr)[ix, ]

    if (grl)
      out = split(out, ix)
    
    return(out)      
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


#' gr2circos
#'
#' dumps rearrangemnts and cn data stored in GRanges into a simple circos format
#' 
#' ra = GRangesList output of ra_breaks
#' cn = outer ring of segment data
#' cn.raw = outer ring of scatter data
#' labels = GRanges of labels to draw in outer, outer ring
#'
#' the plotted field in cn and cn.raw is specified by data.field
#' $col field of cn and cn.raw will specify color of segment / scatter point
#' $col field of text 
#' 
#' if cn does not contain $col field then colors are automatically populated (blue for intervals with value below cn.low
#' and red for intervals with values above cn.high)
#'
#' ra.col vector or $col field of ra (GRangesList value field) can specify colors (with ra.col taking precedence)
#' @name gr2circos
#' @export
gr2circos = function(outdir, ra = NULL, cn = NULL, cn.raw = NULL, data.field = 'cn',
  text = NULL, ## GRanges of text
  text.field = 'label', ## the column of text from which to pull the label
  text.size = 24,
  text.r0 = 1.3, text.r1 = 1.5,    
  ra.col = NULL,
  karyotype.fn = paste(CIRCOS.DIR, 'data', 'karyotype.hg19.txt',  sep = '/'), png = T, svg = T,
  cn.low = 1, cn.high = 1.6,
  cn.line = FALSE, ## if true will draw line plot instead of bar
  cn.cex = 1, ## thickness of line for CN
  cn.ylim = NULL, cn.tick = 1, cn.r0 = 1.0, cn.r1 = 1.3,
  chroms = NULL)
  {        
    # dump out seg data
    # interval data is chr start end numval

    if (!file.exists(Sys.getenv('CIRCOS_DIR')))
      stop('Need to set CIRCOS_DIR environment variable with path to circos installation')
    
    cn.dat = c()

    if (!is.null(chroms))
      if (length(chroms)>1)
        chroms = paste(chroms, collapse = ';')
                   
    if (!is.null(cn))
      cn.dat = c(cn.dat, values(cn)[, data.field])

    if (!is.null(cn.raw))
      cn.dat = c(cn.dat, values(cn.raw)[, data.field])

    if (is.null(cn.ylim))
      {
        cn.dat = cn.dat[!is.infinite(cn.dat)]
        if (length(cn.dat)>1)
          {
            cn.d = diff(range(cn.dat, na.rm = T))
            cn.ylim = c(min(cn.dat, na.rm = T) - 0.1*cn.d, max(cn.dat, na.rm = T) + 0.1*cn.d)
            if (any(is.infinite(cn.ylim)))
              browser()
          }        
        else
          cn.ylim = c(0, 4)
      }
    else if (length(cn.ylim) != 2)
      stop('cn.ylim must be length 2 vector specifying y limits of cn plot')
    
    data.dir = paste(outdir, 'data', sep = '/');
    dir.create(data.dir, recursive = TRUE, showWarnings = FALSE)

    circos.conf = readLines(paste(CIRCOS.DIR, 'circos.conf', sep = '/'))
    circos.out.fn = paste(outdir, 'circos.conf', sep = '/')
    
    f = file(circos.out.fn, 'w')

    if (!svg)
      circos.conf = gsub('svg. = yes', 'svg* = no', circos.conf)
    
    if (!png)
      circos.conf = gsub('png. = yes', 'png* = no', circos.conf)

    if (!is.null(chroms))
      {
        circos.conf = c(paste('chromosomes = ', chroms), circos.conf)

        if (any(grepl("\\-", chroms)))
          circos.conf = gsub('\\$\\$CHROM.DISPLAY', 'yes', circos.conf)
        else
          circos.conf = gsub('\\$\\$CHROM.DISPLAY', 'no', circos.conf)
      }
    else
      circos.conf = gsub('\\$\\$CHROM.DISPLAY', 'yes', circos.conf)
        
    circos.conf = gsub('\\$\\$CN.R0', cn.r0, circos.conf)
    circos.conf = gsub('\\$\\$CN.R1', cn.r1, circos.conf)
    circos.conf = gsub('\\$\\$TEXT.R0', text.r0, circos.conf)
    circos.conf = gsub('\\$\\$TEXT.R1', text.r1, circos.conf)
    circos.conf = gsub('\\$\\$TEXT.SIZE', text.size, circos.conf)
    circos.conf = gsub('\\$\\$CN.MIN', cn.ylim[1], circos.conf)
    circos.conf = gsub('\\$\\$CN.MAX', cn.ylim[2], circos.conf)
    circos.conf = gsub('\\$\\$AXIS.SPACING', cn.tick, circos.conf)
    circos.conf = gsub('\\$\\$CN.THICKNESS', cn.cex*3, circos.conf)
    circos.conf = gsub('\\$\\$CN.TYPE', c('line', 'line')[as.numeric(cn.line)+1], circos.conf)
    
    writeLines(circos.conf, f)
    close(f)

    circos.conf = readLines(karyotype.fn)
    writeLines(circos.conf, paste(outdir, 'data', 'karyotype.hg19.txt', sep = '/'))

    tmp = read.delim(karyotype.fn, sep = ' ', header = F)
    chroms = as.character(unique(tmp[,3][tmp[,1] == 'chr']))
    
    
    if (!is.null(cn))
      {
        if (any(ix <- as.logical(!(seqnames(cn) %in% chroms))))
          {
            warning(sprintf('Removing segments in  these nonexistent chroms: %s.  \nTip: alter karyotype file to include these links',
                            paste(setdiff(as.character(seqnames(cn)[ix], chroms)), collapse = ',')))
            cn = cn[!ix]
          }                          
        
        cn = cn[!is.na(values(cn)[, data.field])]
        cn = cn[!is.infinite(values(cn)[, data.field])]
                
        if (is.null(cn$col))
          cn$col = c('grey', 'blue', 'red')[1+ 1*as.numeric(values(cn)[, data.field]<cn.low) + 2*as.numeric(values(cn)[, data.field]>cn.high)]
        
        tmp.col = t(col2rgb(cn$col, alpha = TRUE))
        cn$col = paste('(', paste(tmp.col[,1], tmp.col[,2], tmp.col[,3], tmp.col[,4]/255, sep = ','), ')', sep = '')
        
        if (cn.line)
          write.table(data.frame(as.character(seqnames(cn)), start(cn), end(cn), values(cn)[, data.field], paste('color=', cn$col, ',stroke_color=', cn$col, sep = '')), paste(data.dir, 'cn.seg.txt', sep = '/'), col.names = F, sep = ' ', quote = F, row.names = F)
        else
          write.table(data.frame(as.character(seqnames(cn)), start(cn), end(cn), values(cn)[, data.field], paste('color=', cn$col, ',stroke_color=', cn$col, ',fill_color=', cn$col, sep = '')), paste(data.dir, 'cn.seg.txt', sep = '/'), col.names = F, sep = ' ', quote = F, row.names = F)
      }
    else      
      writeLines(c(''), paste(data.dir, 'cn.seg.txt', sep = '/'))

    if (!is.null(text))
      {
        if (any(ix <- as.logical(!(seqnames(text) %in% chroms))))
          {
            warning(sprintf('Removing text at these nonexistent chroms: %s.  \nTip: alter karyotype file to include these links',
                            paste(setdiff(as.character(seqnames(text)[ix], chroms)), collapse = ',')))
            text = text[!ix]
          }                          
        
        if (is.null(text$col))
          text$col = 'gray10'

        if (is.null(text$font))
          text$font = 'italic'

        LEGAL.FONTS = c('light', 'normal', 'default', 'semibold', 'bold', 'italic', 'bolditalic', 'italicbold')
        text$font[!(text$font %in% LEGAL.FONTS)] = 'default'
                
        tmp.col = t(col2rgb(text$col, alpha = TRUE))
        text$col = paste('(', paste(tmp.col[,1], tmp.col[,2], tmp.col[,3], tmp.col[,4]/255, sep = ','), ')', sep = '')
                       
        write.table(data.frame(as.character(seqnames(text)), start(text), end(text), values(text)[, text.field], paste('color=', text$col, ',stroke_color=', text$col, ',fill_color=', text$col, ',label_font=', text$font, sep = '')), paste(data.dir, 'text.txt', sep = '/'), col.names = F, sep = ' ', quote = F, row.names = F)
      }
    else
      writeLines('', paste(data.dir, 'text.txt', sep = '/'))

    if (!is.null(cn.raw))
      {
        if (any(ix <- as.logical(!(seqnames(cn.raw) %in% chroms))))
          {
            warning(sprintf('Removing segments at these nonexistent chroms: %s.  \nTip: alter karyotype file to include these segments',
                            paste(setdiff(as.character(seqnames(c.rawn)[ix], chroms)), collapse = ',')))
            cn.raw = cn.raw[!ix]
          }
        
        cn.raw = cn.raw[!is.na(values(cn.raw)[, data.field])]
        cn.raw = cn.raw[!is.infinite(values(cn.raw)[, data.field])]

        if (!is.null(cn.raw$col))
          {
            tmp.col = t(col2rgb(cn.raw$col, alpha = TRUE))
            cn.raw$col = paste('(', paste(tmp.col[,1], tmp.col[,2], tmp.col[,3], tmp.col[,4]/255, sep = ','), ')', sep = '')
          }
        else
          cn.raw$col = c('grey', 'blue', 'red')[1+ 1*as.numeric(values(cn.raw)[, data.field]<cn.low) + 2*as.numeric(values(cn.raw)[, data.field]>cn.high)]
        
        write.table(data.frame(as.character(seqnames(cn.raw)), start(cn.raw), end(cn.raw), values(cn.raw)[, data.field], paste('color=', cn.raw$col, ',stroke_color=', cn.raw$col, sep = '')), paste(data.dir, 'cn.raw.txt', sep = '/'), col.names = F, sep = ' ', quote = F, row.names = F)
      }
    else      
      writeLines(c(''), paste(data.dir, 'cn.raw.txt', sep = '/'))
    
   # links are linkid chr start end color=col
    if (!is.null(ra))
      {
        ra.og = ra;
        ra = grl.pivot(ra)
        
        if (any(ix <- as.logical(!(seqnames(ra[[1]]) %in% chroms) | !(seqnames(ra[[2]]) %in% chroms))))
          {
            warning(sprintf('Removing links to / from these nonexistent chroms: %s.  \nTip: alter karyotype file to include these links',
                            paste(setdiff(as.character(c(seqnames(ra[[1]][ix]), seqnames(ra[[2]][ix]))), chroms), collapse = ',')))
            ra = grl.pivot(ra.og[!ix])
            ra.og = ra.og[ix]
          }                          
        
        if (is.null(ra.col))
          if (is.null(values(ra.og)$ra.col))
            ra.col = alpha(c('lightblue', 'purple'), 0.5)[1 + as.numeric(as.character(seqnames(ra[[1]])) != as.character(seqnames(ra[[2]])))]
          else
            ra.col = values(ra.og)$ra.col
        

        tmp.col = t(col2rgb(ra.col, alpha = T))        
        ra.col = paste('(', paste(tmp.col[,1], tmp.col[,2], tmp.col[,3], tmp.col[,4]/255, sep = ','), ')', sep =  '')

        if (length(unlist(ra))>0)
          {
            links = rbind(
              data.frame(id  = paste('link', 1:length(ra[[1]]), sep = ''), chr = as.character(seqnames(ra[[1]])), pos1 = start(ra[[1]]), pos2 = end(ra[[1]])), 
              data.frame(id  = paste('link', 1:length(ra[[2]]), sep = ''), chr = as.character(seqnames(ra[[2]])), pos1 = start(ra[[2]]), pos2 = end(ra[[2]]))
              );
            
            links$col = paste('color=', ra.col, sep = '')
            links$col[is.na(links$col)] = 'gray'
            links = links[order(-as.numeric(links$id)), ]
#            links = links[order(links$chr, links$pos1), ]

            write.table(links, paste(data.dir, 'links.txt', sep = '/'), col.names = F, sep = ' ', quote = F, row.names = F)
          }
        else
          writeLines('', paste(data.dir, 'links.txt', sep = '/'))
      }
    else
      {
        writeLines('', paste(data.dir, 'links.txt', sep = '/'))
      }

    cat(sprintf('Outputted CIRCOS plot data to directory %s.  To make plot, go into shell / cmd line, cd into directory %s, and type "circos" at the command line to produce svg and png versions of the plot\n', normalizePath(data.dir), normalizePath(outdir)))
  }
   

#' grl.filter
#'
#' filters grl to only include ranges in the specified window
#' (this is different from %in% which does not remove non matching ranges from the grls)
#'
#' does not return list in necessarily same order
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
#' @name grl.allin
#' @export
grl.in = grl.allin = function(grl, windows, some = FALSE, only = FALSE, logical = TRUE, ...)
  {
    if (length(grl)==0)
      return(logical())

    if (length(windows)==0)
      return(rep(T, length(grl)))
                 
    numwin = length(windows);    
    gr = grl.unlist(grl)
    m = grdt(gr.findoverlaps(gr, windows, ...))

    out = rep(FALSE, length(grl))
    if (nrow(m)==0)
      return(out)

    m$grl.id = gr$grl.ix[m$query.id]
    m$grl.iid = gr$grl.iix[m$query.id]

    if (some)
        tmp = as.data.frame(m[, length(unique(grl.iid)), by = grl.id])
    else if (only)
      return(mapply(function(x, y) length(setdiff(x, y))==0,
                    split(1:length(gr), factor(gr$grl.ix, 1:length(grl))),
                    split(m$query.id, factor(m$grl.id, 1:length(grl)))))   
    else
      tmp = aggregate(formula = subject.id ~ grl.id, data = m, FUN = function(x) length(setdiff(1:numwin, x))==0)
    
    out = rep(FALSE, length(grl))
    out[tmp[,1]] = tmp[,2]
    if (logical)
        out = out!=0
    return(out)    
  }

###################
#' grl.split
#'
#' splits GRL's with respect to their seqnames and strand (default), returning
#' new grl whose items only contain ranges with a single seqname / strand
#'
#' can also split by arbitrary (specified) genomic ranges value fields
#' @name grl.split
###################
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

#######################
#' grl.stripnames
#'
#' get rid of gr names inside a grl
#' @name grl.stripnames
#######################
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

#########################
#' grl.span
#'
#' Returns GRanges object representing the left / right extent of each GRL item.  In case of "chimeric" GRL items (ie that map
#' to two chromosomes) there are two options:
#' (1) specify "chr" chromosome as argument to subset GRL's that are on that chromosome, and compute GRL extents from this, any GRL
#'     full outside of that chromosome will get a 0 width GRL 
#' (2) (default) allow chimeric GRL items to get an extent that is with respect to the first chromosome in that GRL 
#'
#' If a grl item contains ranges that lie on different chromosomes, then corresponding grange will have chromosome "NA" and IRange(0, 0)
#' @name grl.span
########################
grl.span = function(grl, chr = NULL, ir = F, keep.strand = T)
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


################################
#' grl.pivot
#'
#' "pivots" grl object "x" by returning a new grl "y" whose
#' kth item is gr object of ranges x[[i]][k] for all i in 1:length(x)
#'
#' Assumes all grs in "x" are of equal length
#' @name grl.pivot
#' @export
################################
grl.pivot = function(x)
  {
    if (length(x) == 0)
      return(GRangesList(GRanges(seqlengths = seqlengths(x)), GRanges(seqlengths = seqlengths(x))))
    return(split(unlist(x), rep(1:length(x[[1]]), length(x))))
  }


########################
#' get.var.col
#'
#' simple function storing default
#' variant color scheme
#' @name get.var.col
########################
get.varcol = function()
  {
    VAR.COL = c('XA' = 'green', 'XG' = 'brown', 'XC' = 'blue', 'XT' = 'red', 'D'= alpha('lightblue', 0.4),
    'I'= 'purple', 'N' = alpha('white', 0.8), 'S' = alpha('pink', 0.9))
    return(VAR.COL)
  }

