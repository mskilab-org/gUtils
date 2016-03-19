#' DNAaseI hypersensitivity sites for hg19
#'
#' DNAaseI hypersensitivity sites from UCSC Table Browser hg19,
#' subsampled to 10,000 sites
#' @name gr.DNAase
#' @docType data
#' @keywords data
#' @format \code{GRanges}
NULL

#' RefSeq genes from UCSC Table Browser hg19, subsampled to 10,000 genes
#'
#' @name gr.genes
#' @docType data
#' @keywords data
#' @format \code{GRanges}
NULL

#' Fake rearrangement data (set 1)
#'
#' @name grl1
#' @docType data
#' @keywords data
#' @format \code{GRangesList}
NULL

#' Fake rearrangement data (set 2)
#'
#' @name grl2
#' @docType data
#' @keywords data
#' @format \code{GRangesList}
NULL

#' \code{Seqinfo} object for hg19
#'
#' @name si
#' @docType data
#' @keywords data
#' @format \code{Seqinfo}
NULL

#' HiC data for chr14 from Lieberman-Aiden 2009 (in hg19), subsampled
#' to 10,000 interactions
#'
#' @name grl.hiC
#' @docType data
#' @keywords data
#' @format \code{GRangesList}
NULL

#' Converts \code{GRanges} to \code{data.table}
#'
#' @importFrom data.table
#'    as.data.table
#' @importFrom GenomicRanges as.data.frame
#' @name gr2dt
#' @param gr \code{GRanges} pile to convert to \code{data.table}
#' @return \code{data.table} with seqnames, start, end, width, strand and all of the meta data. Width is end-inclusive (e.g. [6,7] width = 2)
#' @examples
#' gr2dt(gr.genes)
#' @export
gr2dt <- function(gr)
  {
    ## as.data.frame gives error if duplicated rownames
    if (any(duplicated(names(gr))))
        names(gr) <- NULL
    out <- as.data.table(GenomicRanges::as.data.frame(gr))
    return(out)
  }

#' Get GRanges corresponding to beginning of range
#'
#' @param x \code{GRanges} object to operate on
#' @param width [default = 1] Specify subranges of greater width including the start of the range.
#' @param force [default = F] Allows returned \code{GRanges} to have ranges outside of its \code{Seqinfo} bounds.
#' @param clip [default = F] Trims returned \code{GRanges} so that it does not extend beyond bounds of the input \code{GRanges}
#' @param ignore.strand If set to \code{FALSE}, will extend '-' strands from the other direction [TRUE].
#' @return \code{GRanges} object of width 1 ranges representing start of each genomic range in the input.
#' @importFrom GenomicRanges GRanges
#' @examples
#' gr.start(gr.DNAase, width=200)
#' gr.start(gr.DNAase, width=200, clip=TRUE)
#' @export
gr.start <- function(x, width = 1, force = FALSE, ignore.strand = TRUE, clip = FALSE)
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
            en = pmin(as.vector(start(x))+width-1, seqlengths(x)[as.character(seqnames(x))], na.rm = TRUE)
          }
        else
          {
            st = ifelse(as.logical(strand(x)=='+'),
              as.vector(start(x)),
              pmax(as.vector(end(x))-width+1, 1)
              )

            en = ifelse(as.logical(strand(x)=='+'),
              pmin(as.vector(start(x))+width-1, seqlengths(x)[as.character(seqnames(x))], na.rm = TRUE),
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
#' Takes as input a data.table which must have the following fields: \code{start}, \code{end}, \code{strand}, \code{seqnames}. Will throw
#' an error if any one of these is not present.
#' All of the remaining fields are added as metadata to the \code{GRanges}.
#' @param dt \code{data.table} to convert to \code{GRanges}
#' @return \code{GRanges} object of \code{length = nrow(dt)}
#' @importFrom data.table data.table
#' @importFrom GenomicRanges GRanges mcols<-
#' @importFrom IRanges IRanges
#' @name dt2gr
#' @examples
#' gr <- dt2gr(data.table(start=c(1,2), seqnames=c("X", "1"), end=c(10,20), strand = c('+', '-')))
#' @export
dt2gr <- function(dt) {

  if (any(!c("seqnames","start","end") %in% colnames(dt)))
    stop("gUtils::dt2gr data.table must have seqnames, start, and end")

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

#' Get the right ends of a \code{GRanges}
#'
#' Alternative to \code{GenomicRanges::flank} that will provide end positions *within* intervals
#'
#' @param x \code{GRanges} object to operate on
#' @param width Specify subranges of greater width including the start of the range. \code{[1]}
#' @param force Allows returned \code{GRanges} to have ranges outside of its \code{Seqinfo} bounds. \code{[FALSE]}
#' @param clip Trims returned \code{GRanges} so that it does not extend beyond bounds of the input \code{GRanges}. \code{[TRUE]}
#' @param ignore.strand If set to \code{FALSE}, will extend '-' strands from the other direction. \code{[TRUE]}
#' @return \code{GRanges} object of width = \code{width} ranges representing end of each genomic range in the input.
#' @examples
#' gr.end(gr.DNAase, width=200, clip=TRUE)
#' @importFrom GenomeInfoDb seqlengths
#' @importFrom GenomicRanges strand seqnames values<- values
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
              pmin(as.vector(start(x))+width-1, seqlengths(x)[as.character(seqnames(x))], na.rm = TRUE)
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

#' Get the midpoints of \code{GRanges} ranges
#'
#' @param x \code{GRanges} object to operate on
#' @return \code{GRanges} of the midpoint, calculated from \code{floor(width(x)/2)}
#' @export
#' @importFrom GenomicRanges start<- end<- start end
#' @examples
#' gr.mid(GRanges(1, IRanges(1000,2000), seqinfo=Seqinfo("1", 2000)))
gr.mid = function(x)
  {
      start(x) = end(x) = rowMeans(cbind(start(x), end(x)))
      return(x)
  }

#' Generate random \code{GRanges} on genome
#'
#' Randomly generates non-overlapping \code{GRanges} with supplied widths on supplied genome.
#' Seed can be supplied with \code{set.seed}
#'
#' @param w Vector of widths (length of \code{w} determines length of output)
#' @param genome Genome which can be a \code{GRanges}, \code{GRangesList}, or \code{Seqinfo} object. Default is "hg19" from the \code{BSGenome} package.
#' @return \code{GRanges} with random intervals on the specifed "chromosomes"
#' @note This function is currently quite slow, needs optimization
#' @importFrom GenomeInfoDb seqinfo seqnames<-
#' @importFrom GenomicRanges gaps ranges ranges<-
#' @examples
#' ## Generate a single random interval of width 10, on "chr" of length 1000
#' gr.rand(10, Seqinfo("1", 1000))
#' ## Generate 5 non-overlapping regions of width 10 on hg19
#' library(BSgenome.Hsapiens.UCSC.hg19)
#' gr.rand(rep(10,5), Hsapiens)
#' @export
gr.rand = function(w, genome)
  {

    if (!is(genome, 'Seqinfo'))
      genome = seqinfo(genome)

    sl = seqlengths(genome);
    available = si2gr(genome);

    out = GRanges(rep(names(sl)[1], length(w)), IRanges(rep(1, length(w)), width = 1), seqlengths = seqlengths(genome));
    for (i in 1:length(w))
      {
        if (i == 1)
          available = si2gr(genome)
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
            rstart = ceiling(stats::runif(1)*starts[length(starts)])-starts
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


#' Trims pile of \code{GRanges} relative to the specified <local> coordinates of each range
#'
#' Example: \code{GRanges} with genomic coordinates 1:1,000,000-1,001,000 can get the first 20 and last 50 bases trimmed off with
#' \code{start = 20, end = 950}.
#' if end is larger than the width of the corresponding gr, then the corresponding output will only have \code{end(gr)} as its coordinate.
#'
#' This is a role not currently provided by the standard \code{GenomicRanges} functions
#' (e.g. \code{shift}, \code{reduce}, \code{restrict}, \code{shift}, \code{resize}, \code{flank})
#' @param gr \code{GRanges} to trim
#' @param starts Number of bases to trim off of the front\code{[1]}
#' @param ends Number of bases to trim off of the back\code{[1]}
#' @examples
#' ## trim the first 20 and last 50 bases
#' gr.trim(GRanges(1, IRanges(1e6, width=1000)), starts=20, ends=950)
#' ## return value: GRanges on 1:1,000,019-1,000,949
#' @export
gr.trim = function(gr, starts=1, ends=1)

  ## note that I deleted ignore.strand and fromEnd options
    {
    starts = cbind(1:length(gr), starts)[, 2]
    ends = cbind(1:length(gr), ends)[, 2]

    ##if (!ignore.strand)
    ##    {
    ##    ix = as.logical(strand(gr)=='-')
    ##    if (any(ix))
          ##  {
            ##if (fromEnd)
            ##    {
            ##    tmp = starts[ix]
            ##    starts[ix] = ends[ix]
            ##    ends[ix] = tmp-1
            ##}
            ##else
            ##    {
      ##          starts[ix] = width(gr)[ix]-starts[ix]+1
      ##          ends[ix] = width(gr)[ix]-ends[ix]+1
            ##}
    ##    }
    ##}

  ##  if (fromEnd) {
  ##    en = pmax(starts, end(gr)-ends);
  ##} else {
      ends = pmax(starts, ends);
      ends = pmin(ends, width(gr));
      en = start(gr) + ends - 1;
  ##}

    st = start(gr)+starts-1;
    st = pmin(st, en);

    out = GRanges(seqnames(gr), IRanges(st, en),
                   seqlengths = seqlengths(gr), strand = strand(gr))
    values(out) = values(gr)
    return(out)
}

#' Randomly sample \code{GRanges} intervals within territory
#'
#' Samples \code{k} intervals of length "len" from a pile of \code{GRanges}.
#' \itemize{
#' \item If k is a scalar then will (uniformly) select k intervals from the summed territory of \code{GRanges}
#' \item If k is a vector of length(gr) then will uniformly select k intervals from each.
#' }
#' @param gr \code{GRanges} object defining the territory to sample from
#' @param k Number of ranges to sample
#' @param len Length of the \code{GRanges} element to produce [100]
#' @param replace If TRUE, will bootstrap, otherwise will sample without replacement. [TRUE]
#' @return GRanges of max length sum(k) [if k is vector) or k*length(gr) (if k is scalar) with labels indicating the originating range.
#'
#' @examples
#' ## sample 5 \code{GRanges} of length 10 each from territory of RefSeq genes
#' gr.sample(reduce(gr.genes), k=5, len=10)
#' @note This is different from \code{GenomicRanges::sample} function, which just samples from a pile of \code{GRanges}
#' @export
gr.sample = function(gr, k, len = 100, replace = TRUE)
{
  if (!inherits(gr, 'GRanges'))
    gr = si2gr(gr)

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
        s = terr*stats::runif(k)

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
                        s = (gr.df$end[i]-gr.df$start[i]-gr.df$len[i])*stats::runif(k[i])+gr.df$start[i]

                      return(data.frame(chr = gr.df$chr[i], start=s, end =s+len-1, strand = as.character(strand(gr)[i]), query.id = i))
                    })
      return(gr.fix(seg2gr(do.call('rbind', tmp)), gr))
    }
}

#' Create \code{GRanges} from \code{Seqinfo} or \code{BSgenome}
#'
#' Creates a genomic ranges from seqinfo object
#' ie a pile of ranges spanning the genome
#' @param si \code{Seqinfo} object or a \code{BSgenome} genome
#' @param strip.empty Don't know. \code{[FALSE]}
#' @return \code{GRanges} representing the range of the input genome
#' @examples
#' \dontrun{libary(BSgenome.Hsapiens.UCSC.hg19); si2gr(Hsapiens)}
#' @export
si2gr <- function(si, strip.empty = FALSE)
  {
    if (is(si, "BSgenome"))
        si <- Seqinfo(names(seqlengths(si)), seqlengths(si))
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

#' Concatenate \code{GRanges}, robust to different \code{mcols}
#'
#' Concatenates \code{GRanges} objects, taking the union of their features if they have non-overlapping features
#' @param x First \code{GRanges}
#' @param ... additional \code{GRanges}
#' @note Does not fill in the \code{Seqinfo} for the output \code{GRanges}
#' @return Concatenated \code{GRanges}
#' grbind(gr.genes, gr.DNAase)
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
    # if ('force.rrbind' %in% names(grs))
    #   {
    #     if (is.logical(grs[['force.rrbind']]))
    #       force.rrbind = grs[['force.rrbind']]
    #     grs = grs[-match('force.rrbind', names(grs))]
    #   }

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
        vals[isDataFrame] = lapply(vals[isDataFrame], gr2dt)

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
        sl.new[names(sl)] = pmax(sl.new[names(sl)], sl, na.rm = TRUE)

    bare.grs = lapply(grs, function(x) gr.fix(x[,c()], sl.new))
    out = tryCatch(do.call('c', bare.grs), error = function(e) NULL)

    # if (is.null(out) | is.list(out)) ## something failed with concatenation, likely some weird ghost with concatenating GRanges with 'c', below is a temp fix
    #     {
    #         getsridofghostsomehow =  c(bare.grs[[1]], bare.grs[[2]])
    #         out = tryCatch(do.call('c', bare.grs), error = function(e) NULL)
    #
    #         if (is.null(out) | is.list(out)) ## now we are really reaching
    #             out = seg2gr(do.call('rrbind', lapply(bare.grs, as.data.frame)), seqlengths = sl.new)[, c()]
    #     }

    # hack to deal with empty value columns
    ix <- (sapply(vals, ncol)==0)
    if (any(ix))
        vals[ix] = lapply(which(ix), function(x) data.frame(col.4214124124124 = rep(NA, length(grs[[x]]))))

    if (!force.rrbind)
      tmp = tryCatch(do.call('rrbind', vals), error = function(e) NULL)
    else
        tmp = NULL

    #if (is.null(tmp) | force.rrbind) ## sometimes rrbind2 gets picky because of type checking (rbindlist) .. so just run rrbind then
    #   values(out) = do.call('rrbind', vals)
    # else
       values(out) = tmp

    if (any(ix))
        out$col.4214124124124 = NULL
    return(out)
  }

#' Concatenate \code{GRangesList} objects
#'
#' Concatenates \code{GRangesList} objects taking the union of their \code{mcols} features if they have non-overlapping features
#' @param ... Any number of \code{GRangesList} to concatenate together
#' @return Concatenated \code{GRangesList} with NA filled in for \code{mcols} fields that are non-overlapping
#' @examples
#' ## Concatenate
#' grl.hiC2 <- grl.hiC[1:20]
#' mcols(grl.hiC2)$test = 1
#' grlbind(grl.hiC2, grl.hiC[1:30])
#' @export
#' @importFrom GenomicRanges mcols<- mcols split
grlbind = function(...)
  {
    ## TODO: make this work for when underlying grs do not have matching features
    ## currently will loose gr level features
    grls = list(...)

    ## annoying acrobatics to reconcile gr and grl level features for heterogenous input gr / grls
    grls.ul = lapply(grls, grl.unlist)
    grls.ul.rb = do.call('grbind', grls.ul)
    sp = base::unlist(lapply(1:length(grls), function(x) rep(x, length(grls.ul[[x]]))))
    gix = base::split(grls.ul.rb$grl.ix, sp)
    gjx = base::split(1:length(grls.ul.rb), sp)
    grls.ul.rb$grl.iix = grls.ul.rb$grl.ix = NULL

    grls.vals = lapply(grls, function(x)
      { if (ncol(mcols(x))>0)  return(as.data.frame(mcols(x))) else return(data.frame(dummy241421 = rep(NA, length(x))))})

    grls.new = mapply(function(x,y) GenomicRanges::split(grls.ul.rb[x],y), gjx, gix)

    ## do.call('c', grls.new) is not working for some reason (gives back list again, not GRangesList)
    ## have to do this instead, not ideal
    if (length(grls.new) > 1) {
      out = grls.new[[1]]
        for (i in 2:length(grls.new))
          out = c(out, grls.new[[i]])
    } else {
      out = grls.new[[1]]
    }

    # out = do.call('c', grls.new)
    #
    # if (is.list(out))
    #   {
    #     if (length(grls.new)>1)
    #       {
    #         bla = c(grls.new[[1]], grls.new[[2]]) ## fix R ghost
    #         out = do.call('c', grls.new)
    #         if (is.list(out)) ## if still is list then do manual 'c'
    #              {
    #                out = grls.new[[1]]
    #                for (i in 2:length(grls.new))
    #                    out = c(out, grls.new[[i]])
    #             }
    #       }
    #     else
    #       out = grls.new[[1]]
    #   }

    out.val = do.call('rrbind', grls.vals)
    out.val$dummy241421 = NULL
    GenomicRanges::mcols(out) <- out.val

    return(out)
  }

#' Prepend "chr" to \code{GRanges seqlevels}
#'
#' @param gr \code{GRanges} object to append 'chr' to
#' @return Identical \code{GRanges}, but with 'chr' prepended to each seqlevel
#' @examples
#' gr <-  gr.chr(GRanges(c(1,"chrX"), IRanges(c(1,2), 1)))
#' seqnames(gr)
#' @importFrom GenomeInfoDb seqlevels seqlevels<-
#' @export
gr.chr = function(gr)
  {
    if (any(ix <- !grepl('^chr', seqlevels(gr))))
      seqlevels(gr)[ix] = paste('chr', seqlevels(gr)[ix], sep = "")
    return(gr)
  }

#' Reduce \code{GRanges} and \code{GRangesList} to miminal footprint
#'
#' Shortcut for \code{reduce(sort(gr.stripstrand(unlist(x))))}
#' @param gr \code{GRanges} or \code{GRangesList}
#' @param pad Expand the input data before reducing. \code{[0]}
#' @param sort Flag to sort the output. \code{[TRUE]}
#' @return \code{GRanges} object with no strand information, representing a minimal footprint
#' @importFrom GenomicRanges reduce
#' @examples
#' streduce(grl.hiC, pad=10)
#' streduce(gr.genes, pad=1000)
#' @export
streduce = function(gr, pad = 0, sort = TRUE)
  {

    if (inherits(gr, 'GRangesList'))
      gr = unlist(gr)

    if (any(is.na(seqlengths(gr))))
      gr = gr.fix(gr)

    #out = suppressWarnings(sort(reduce(gr.stripstrand(gr+pad))))
    out = suppressWarnings(sort(reduce(gr.stripstrand(gr + pad))))
    suppressWarnings(start(out) <-pmax(1, start(out)))
#    out <- gr.tfix(out)
    end(out) = pmin(end(out), seqlengths(out)[as.character(seqnames(out))])


    return(out)
  }

#' Return UCSC style interval string corresponding to \code{GRanges} pile (ie chr:start-end)
#'
#' @param gr \code{GRanges} pile to get intervals from
#' @param add.chr Prepend seqnames with "chr" \code{[FALSE]}
#' @param mb Round to the nearest megabase \code{[FALSE]}
#' @param round If \code{mb} supplied, how many digits to round to. \code{[3]}
#' @param other.cols Names of additional \code{mcols} fields to add to the string (seperated by ";")
#' @name gr.string
#' @examples
#' gr.string(gr.genes, other.cols = c("name", "name2"))
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
      other.str = paste(';', do.call('paste', c(lapply(other.cols, function(x) values(gr)[, x]), list(sep = ';'))))
    else
      other.str = ''

    str = ifelse(as.logical(strand(gr)!='*'), as.character(strand(gr)), '')

    if (mb)
      return(paste(sn, ':', round(start(gr)/1e6, round), '-', round(end(gr)/1e6, round), str, other.str, sep = ''))
    else
      return(paste(sn, ':', start(gr), '-', end(gr), str, other.str, sep = ''))
  }


#' Create string representation of \code{GRangesList}
#'
#' Return ucsc style interval string corresponding to each \code{GRanges} in the \code{GRangesList}.
#' One line per per \code{GRangesList} item. \code{GRanges} elements themselves are separated by \code{sep}
#'
#' @param grl \code{GRangesList} to convert to string vector
#' @param mb Will return as MB and round to "round" [FALSE]
#' @param sep Character to separate single \code{GRanges} ranges [,]
#' @param ... Additional arguments to be passed to \code{gr.string}
#' @return Character vector where each element is a \code{GRanges} pile corresponding to a single \code{GRangesList} element
#' @name grl.string
#' @examples
#' grl.string(grl.hiC, mb=TRUE)
#' @export
grl.string = function(grl, mb= FALSE, sep = ',', ...)
  {

    if (class(grl) == "GRanges")
      return(gr.string(grl, mb=mb, ...))

    if (class(grl) != "GRangesList")
      stop("Input must be GRangesList (or GRanges, which is sent to gr.string)")

    gr = grl.unlist(grl)
    if (!is.null(names(grl)))
      nm = names(grl)
    else
      nm = 1:length(grl)

    grs = gr.string(gr, mb=mb, ...)
    out = sapply(split(grs, gr$grl.ix), paste, collapse = sep)
    names(out) = nm[as.numeric(names(out))]
    return(out)
  }

#' "Fixes" \code{seqlengths} / \code{seqlevels}
#'
#' If "genome" not specified will replace \code{NA} \code{seqlengths} in \code{GRanges} to reflect largest coordinate per \code{seqlevel}
#' and removes all \code{NA seqlevels} after this fix.
#'
#' if "genome" defined (i.e. as \code{Seqinfo} object, or a \code{BSgenome}, \code{GRanges}, \code{GRangesList} object with populated \code{seqlengths}),
#' then will replace \code{seqlengths} in \code{gr} with those for that genome
#' @name gr.fix
#' @param gr \code{GRanges} object to fix
#' @param genome Genome to fix to: \code{Seqinfo}, \code{BSgenome}, \code{GRanges} (w/seqlengths), \code{GRangesList} (w/seqlengths)
#' @param gname Name of the genome (optional, just appends to \code{Seqinfo} of the output) [NULL]
#' @param drop Remove ranges that are not present in the supplied genome [FALSE]
#' @return \code{GRanges} pile with the fixed \code{Seqinfo}
#' @importFrom GenomeInfoDb Seqinfo seqinfo keepSeqlevels seqlevels seqlengths seqlevels<- seqlengths<- genome<- seqnames
#' @export
gr.fix = function(gr, genome = NULL, gname = NULL,  drop = FALSE)
  {
    #### marcin: now it is
    ## if (inherits(gr, "GRangesList"))
    ##   stop("gr.fix not setup to take GRangesList")

    sn = V1 = NULL ## NOTE fix

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
              lens[seqlevels(gr)] = pmax(seqlengths(gr), lens[seqlevels(gr)], na.rm = TRUE)

            }
          else
            lens = structure(seqlengths(genome), names = seqlevels(genome))
        # else
        #   {
        #     if (is.character(genome))
        #       genome = structure(rep(NA, length(genome)), names = genome)
        #
        #     lens = structure(NA, names = union(names(genome), seqlevels(gr)));
        #
        #     lens[seqlevels(gr)] = seqlengths(gr);
        #     lens[names(genome)[!is.na(genome)]] = pmax(lens[names(genome)[!is.na(genome)]], genome[!is.na(genome)], na.rm = TRUE)
        #
        #     if (drop)
        #       lens = lens[names(genome)]
        #   }

        seqlevels(gr, force = TRUE) = names(lens)
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

                tmp.sl = data.table(sn = as.character(seqnames(tmp.gr)), end = end(tmp.gr))[, max(end, na.rm = TRUE), by = sn][,  structure(V1, names = sn)][seqlevels(tmp.gr)]
                names(tmp.sl) = seqlevels(tmp.gr)
                seqlengths(tmp.gr)[!is.na(tmp.sl)] = suppressWarnings(pmax(tmp.sl[!is.na(tmp.sl)], seqlengths(tmp.gr)[!is.na(tmp.sl)], na.rm = TRUE))
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


#' Lay ranges end-to-end onto a derivate "chromosome"
#'
#' Takes pile of \code{GRanges} and returns into a \code{data.frame} with \code{nrow = length(gr)} with each
#' representing the corresponding input range superimposed onto a single "flattened"
#' chromosome, with ranges laid end-to-end
#' @param gr \code{GRanges} to flatten
#' @param gap Number of bases between ranges on the new chromosome \code{[0]}
#' @return \code{data.frame} with start and end coordinates, and all of the original metadata
#' @name gr.flatten
#' @importFrom GenomicRanges mcols
#' @export
gr.flatten = function(gr, gap = 0)
  {
    if (length(gr) == 0)
      return(data.frame())
    else if (length(gr) == 1)
      return(data.frame(start = 1, end = width(gr)))
    else
      {
        starts = as.numeric(cumsum(c(1, width(gr[1:(length(gr)-1)])+gap)))
        ends = as.numeric(starts+width(gr)-1)
		return(cbind(data.frame(start=starts, end=ends), as.data.frame(mcols(gr))))
      }
}

#' gr.stripstrand
#'
#' sets strand to "*"
#' @name gr.stripstrand
#' @keywords internal
gr.stripstrand = function(gr)
  {
    strand(gr) = "*"
    return(gr)
  }

#' Flip strand on \code{GRanges}
#'
#' @name gr.flipstrand
#' @param gr \code{GRanges} pile with strands to be flipped
#' @return \code{GRanges} with flipped strands (+ to -, * to *, - to *)
#' @examples
#' gr.flipstrand(GRanges(1, IRanges(c(10,10,10),20), strand=c("+","*","-")))
#' @export
gr.flipstrand <- function(gr)
  {
    if (!is(gr, 'GRanges'))
      stop('GRanges input only')

    if (length(gr)==0)
      return(gr)

    which = cbind(1:length(gr), TRUE)[,2] == 1

    if (any(which))
      strand(gr)[which] = c('*'='*', '+'='-', '-'='+')[as.character(strand(gr))][which]

    return(gr)
  }

#' Tile ranges across \code{GRanges}
#'
#' Tiles interval (or whole genome) with segments of \code{<=} specified width.
#' @param gr \code{GRanges}, \code{seqlengths} or \code{Seqinfo} range to tile. If has \code{GRanges} has overlaps, will reduce first.
#' @param w Width of each tile
#' @name gr.tile
#' @examples
#' ## 10 tiles of width 10
#' gr1 <- gr.tile(GRanges(1, IRanges(1,100)), w=10)
#' ## make them overlap each other by 5
#' gr1 + 5
#' @export
gr.tile = function(gr, w = 1e3)
  {
    numw = tile.id = query.id = NULL ## getting past NOTE
    if (!is(gr, 'GRanges'))
      gr = si2gr(gr);

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

#' gr.tile.map
#'
#' Given two tilings of the genome (e.g. at different resolution)
#' query and subject outputs a length(query) list whose items are integer vectors of indices in subject
#' overlapping that overlap that query (strand non-specific)
#'
#' @param query Query \code{GRanges} pile, perhaps created from some tile (e.g. \code{gr.tile}), and assumed to have no gaps
#' @param subject Subject \code{GRanges} pile, perhaps created from some tile (e.g. \code{gr.tile}), and assumed to have no gaps
#' @param verbose Increase the verbosity of the output
#' @return \code{list} of length = \code{length(query)}, where each element \code{i} is a vector of indicies in \code{subject} that overlaps element \code{i} of \code{query}
#' @note Assumes that input query and subject have no gaps (including at end) or overlaps, i.e. ignores end()
#' coordinates and only uses "starts"
#' @export
gr.tile.map = function(query, subject, verbose = FALSE)
  {

    mc.cores =1 ## multicore not supported now
  if (length(GenomicRanges::gaps(query)) > 0)
    warning("Query GRanges has gaps. Unexpected behavior may follow")
  if (length(GenomicRanges::gaps(subject)) > 0)
    warning("Subject GRanges has gaps. Unexpected behavior may follow")

    ix.q = order(query)
    ix.s = order(subject)

    q.chr = as.character(seqnames(query))[ix.q]
    s.chr = as.character(seqnames(subject))[ix.s]

    ql = base::split(ix.q, q.chr)
    sl = base::split(ix.s, s.chr)

    tmp <- parallel::mcmapply(
    #tmp <- lapply(
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
      #}, ql, sl[names(ql)], SIMPLIFY = FALSE)

    m = munlist(tmp)[, -c(1:2), drop = FALSE]
    out = split(m[,2], m[,1])[as.character(1:length(query))]
    names(out) = as.character(1:length(query))
    return(out)
  }

#' Dice up \code{GRanges} into \code{width = 1} \code{GRanges} spanning the input (warning can produce a very large object)
#'
#' @param gr \code{GRanges} object to dice
#' @name gr.dice
#' @importFrom S4Vectors Rle
#' @importFrom GenomeInfoDb seqlengths
#' @importFrom GenomicRanges seqnames width strand values<- values strand<- distance
#' @return \code{GRangesList} where kth element is a diced pile of \code{GRanges} from kth input \code{GRanges}
#' @examples
#' gr.dice(GRanges(c(1,4), IRanges(c(10,10),20)))
#' @export
gr.dice = function(gr)
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

    out = S4Vectors::split(out, ix)

    return(out)
  }


#' Pairwise distance between two \code{GRanges}
#'
#' Computes matrix of pairwise distance between elements of two \code{GRanges} objects of length \code{n} and \code{m}.
#'
#' Distances are computed as follows:
#' \itemize{
#' \item NA for ranges on different seqnames
#' \item 0 for overlapping ranges
#' \item min(abs(end1-end2), abs(end1-start2), abs(start1-end2), abs(start1-end1),) for all others
#' }
#' If only \code{gr1} is provided, then will return n x n matrix of \code{gr1} vs itself \cr
#' If \code{max.dist = TRUE}, then will replace \code{min} with \code{max} above
#' @param gr1 First \code{GRanges}
#' @param gr2 Second \code{GRanges}
#' @param ignore.strand Don't required elements be on same strand to avoid \code{NA [FALSE]}
#' @param ... Additional arguments to be supplied to \code{GenomicRanges::distance}
#' @return \code{N} by \code{M} matrix with the pairwise distances, with \code{gr1} on rows and \code{gr2} on cols
#' @name gr.dist
#' @export
gr.dist = function(gr1, gr2 = NULL, ignore.strand = FALSE, ...)
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

#' Check intersection of \code{GRangesList} with windows on genome
#'
#' Like %in% for grl but now will return a logical vector that is true at position if i
#' only if the ranges in grl[i] intersect <<all>>, <<some>>, <<only>>  windows in the subject
#'
#' eg can use to identify read pairs whose ends are contained inside two genes)
#' @param grl \code{GRangesList} object to query for membership in \code{windows}
#' @param windows \code{GRanges} pile of windows
#' @param some Will return \code{TRUE} for \code{GRangesList} elements that intersect at least on window range [FALSE]
#' @param only Will return \code{TRUE} for \code{GRangesList} elements only if there are no elements of query that fail to interesect with windows [FALSE]
#' @param ... Additional parameters to be passed on to \code{GenomicRanges::findOverlaps}
#' @name grl.in
#' @export
grl.in <- function(grl, windows, some = FALSE, only = FALSE, ...)
  {
    grl.iid = grl.id = NULL ## for getting past NOTE

    if (length(grl)==0)
      return(logical())

    if (length(windows)==0)
      return(rep(TRUE, length(grl)))

    numwin = length(windows);
    gr = grl.unlist(grl)
    h = GenomicRanges::findOverlaps(gr, windows, ...)
    m = data.table(query.id = queryHits(h), subject.id = subjectHits(h))
    #m = gr2dt(gr.findoverlaps(gr, windows, ...))

    out = rep(FALSE, length(grl))
    if (nrow(m)==0)
      return(out)

    m$grl.id = gr$grl.ix[m$query.id]
    m$grl.iid = gr$grl.iix[m$query.id]

    if (some)
        tmp = as.data.frame(m[, length(unique(grl.iid)), by = grl.id])
    else if (only)
      return(base::mapply(function(x, y) length(setdiff(x, y))==0,
                    split(1:length(gr), factor(gr$grl.ix, 1:length(grl))),
                    split(m$query.id, factor(m$grl.id, 1:length(grl)))))
    else
      tmp = stats::aggregate(formula = subject.id ~ grl.id, data = m, FUN = function(x) length(setdiff(1:numwin, x))==0)

    out = rep(FALSE, length(grl))
    out[tmp[,1]] = tmp[,2]
    ##if (logical)
        out = out!=0
    return(out)
  }

#' Robust unlisting of \code{GRangesList} that keeps track of origin
#'
#' Does a "nice" unlist of a \code{GRangesList} object adding a field \code{grl.ix} denoting which element of the \code{GRangesList}
#' each \code{GRanges} corresponds to and a field \code{grl.iix} which saves the (local) index that that gr was in its corresponding \code{GRangesList} item
#'
#' In this way, \code{grl.unlist} is reversible, while \code{BiocGenerics::unlist} is not.
#' @name grl.unlist
#' @importFrom BiocGenerics unlist
#' @param grl \code{GRangeList} object to unlist
#' @return \code{GRanges} with added metadata fields \code{grl.ix} and \code{grl.iix}.
#' @examples
#' grl.unlist(grl.hiC)
#' @export
grl.unlist = function(grl)
{
    if (length(grl) == 0) ## JEREMIAH
      return(GRanges())

    names(grl) = NULL
    as.df = as.data.frame(grl)

    el = as.df$element
    if (is.null(el))
        el = as.df$group

    out = BiocGenerics::unlist(grl)
    mcols(out)$grl.ix = el
    tmp = rle(el)

    out$grl.iix = unlist(sapply(tmp$lengths, function(x) 1:x))
    values(out) = cbind(values(grl)[out$grl.ix, , drop = FALSE], values(out))
    return(out)
  }

#' Pivot a \code{GRangesList}, inverting "x" and "y"
#'
#' "Pivots" grl object "x" by returning a new grl "y" whose
#' kth item is gr object of ranges x[[i]][k] for all i in 1:length(x)
#' e.g. If \code{length(grl)} is 50 and length of each \code{GRanges} element inside is 2, then \code{grl.pivot}
#' will produce a length 3 \code{GRangesList} with 50 elements per \code{GRanges}
#'
#' Assumes all grs in "x" are of equal length
#' @name grl.pivot
#' @param x \code{GRangesList} object to pivot
#' @importFrom GenomicRanges GRanges GRangesList split unlist
#' @examples
#' grl.pivot(grl.hiC)
#' @export
grl.pivot = function(x)
  {
    if (length(x) == 0)
      return(GRangesList(GRanges(seqlengths = seqlengths(x)), GRanges(seqlengths = seqlengths(x))))
    #gg <- GenomicRanges::split(GenomicRanges::unlist(x), rep(1:length(x[[1]]), length(x)))
    return(GenomicRanges::split(GenomicRanges::unlist(x), rep(1:length(x[[1]]), length(x))))
}


#' Improved \code{rbind} for intersecting/union columns of \code{data.frames} or \code{data.tables}
#'
#' Like \code{rbind}, but takes the intersecting columns of the data.
#' @param ... Any number of \code{data.frame} or \code{data.table} objects
#' @param union Take union of columns (and put NA's for columns of df1 not in df2 and vice versa). \code{[TRUE]}
#' @param as.data.table Return the binded data as a \code{data.table}. \code{[FALSE]}
#' @return \code{data.frame} or \code{data.table} of the \code{rbind} operation
#' @export
#' @importFrom data.table data.table rbindlist
rrbind = function(..., union = TRUE, as.data.table = FALSE)
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

#' @name seg2gr
#' @title Convert GRange like data.frames into GRanges
#' @description
#'
#' Take data frame of ranges "segs" and converts into granges object porting over additional value columns
#' "segs" data frame can obey any number of conventions to specify chrom, start, and end of ranges
#' (eg $pos1, $pos2, $Start_position, $End_position) --> see "standardize_segs" for more info
#'
#' @importFrom GenomicRanges
#'    GRanges
#' @param segs data frame of segments with fields denoting chromosome, start, end, and other metadata (see standardized segs for seg data frame input formats)
#' @param seqlengths seqlengths of output GRanges object
#' @param seqinfo seqinfo of output GRanges object
#' @keywords internal
seg2gr = function(segs, key = NULL, seqlengths = NULL, seqinfo = Seqinfo())
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

  # if (any(bix <- (is.na(segs$chr) | is.na(segs$pos1) | is.na(segs$pos2))))
  # {
  #   warning('Segments with NA values for chromosome, start, and end position detected .. removing')
  #   segs = segs[!bix, ]
  # }

  GR.NONO.FIELDS = c('seqnames', 'ranges', 'strand', 'seqlevels', 'seqlengths', 'isCircular', 'start', 'end', 'width', 'element');

  if (is.null(segs$strand))
    segs$strand = "*"

  if (any(ix <- !(segs$strand %in% c('+', '-', '*'))))
    segs$strand[ix] = "*"

  if (length(seqlengths)>0)
  {
    # if (length(wtf  <- setdiff(segs$chr, names(seqlengths))))
    # {
    #   warning('some seqnames in seg object were not included in provided seqlengths: ', paste(wtf, collapse = ','))
    #   seqlengths[as.character(wtf)] = NA
    # }
    # segs$pos1 <- as.numeric(segs$pos1)
    # segs$pos2 <- as.numeric(segs$pos2)
    #
    # out = GRanges(seqnames = segs$chr, ranges = IRanges(segs$pos1, segs$pos2), names = levels(levels), strand = segs$strand, seqlengths = seqlengths)
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

#' standardize_segs
#'
#' (data frame seg function)
#'
#' Takes and returns segs data frame standardized to a single format (ie $chr, $pos1, $pos2)
#'
#' if chr = TRUE will ensure "chr" prefix is added to chromossome(if does not exist)
#' @keywords internal
standardize_segs = function(seg, chr = FALSE)
{
  #if (inherits(seg, 'IRangesList'))
  #  seg = irl2gr(seg);

  if (is(seg, 'matrix'))
    seg = as.data.frame(seg, stringsAsFactors = FALSE)

  # if (inherits(seg, 'RangedData') | inherits(seg, 'GRanges') | inherits(seg, 'IRanges'))
  # {
  #   val = as.data.frame(values(seg));
  #   values(seg) = NULL;
  #   seg = as.data.frame(seg, row.names = NULL);  ## returns compressed iranges list
  #   seg$seqnames = as.character(seg$seqnames)
  # }
  # else
    val = NULL;

  field.aliases = list(
    ID = c('id', 'patient', 'Sample'),
    chr = c('seqnames', 'chrom', 'Chromosome', "contig", "seqnames", "seqname", "space", 'chr', 'Seqnames'),
    pos1 = c('start', 'loc.start', 'begin', 'Start', 'start', 'Start.bp', 'Start_position', 'pos', 'pos1', 'left', 's1'),
    pos2 =  c('end', 'loc.end', 'End', 'end', "stop", 'End.bp', 'End_position', 'pos2', 'right', 'e1'),
    strand = c('strand', 'str', 'strand', 'Strand', 'Str')
  );

  if (is.null(val))
    val = seg[, setdiff(names(seg), unlist(field.aliases))]

  seg = seg[, intersect(names(seg), unlist(field.aliases))]

  for (field in setdiff(names(field.aliases), names(seg)))
    if (!(field %in% names(seg)))
      names(seg)[names(seg) %in% field.aliases[[field]]] = field;

  if (chr)
    if (!is.null(seg$chr))
      if (!grepl('chr', seg$chr[1]))
        seg$chr = paste('chr', seg$chr, sep = "");

  if (is.null(seg$pos2))
    seg$pos2 = seg$pos1;

  missing.fields = setdiff(names(field.aliases), c(names(seg), c('chr', 'ID', 'strand')));

  if (length(missing.fields)>0)
    warning(sprintf('seg file format problem, missing an alias for the following fields:\n\t%s',
                    paste(sapply(missing.fields, function(x) paste(x, '(can also be', paste(field.aliases[[x]], collapse = ', '), ')')), collapse = "\n\t")));

  if (!is.null(val))
    seg = cbind(seg, val)

  return(seg)
}

#' Remove chr prefix from GRanges seqlevels
#'
#' @param gr \code{GRanges} with chr seqlevel prefixes
#' @return GRanges without chr seqlevel prefixes
#' @export
gr.nochr = function(gr) {
  if (grepl('^chr', seqlevels(gr)[1]))
    seqlevels(gr) = gsub('^chr','', seqlevels(gr))
  return(gr)
}

# munlist
#
# unlists a list of vectors, matrices, data frames into a n x k matrix
# whose first column specifies the list item index of the entry
# and second column specifies the sublist item index of the entry
# and the remaining columns specifies the value(s) of the vector
# or matrices.
#
# force.cbind = T will force concatenation via 'cbind'
# force.rbind = T will force concatenation via 'rxsbind'
munlist = function(x, force.rbind = FALSE, force.cbind = FALSE, force.list = FALSE)
{
  if (!any(c(force.list, force.cbind, force.rbind)))
  {
    if (any(sapply(x, function(y) is.null(dim(y)))))
      force.list = TRUE
    if (length(unique(sapply(x, function(y) dim(y)[2]))) == 1)
      force.rbind = TRUE
    if ((length(unique(sapply(x, function(y) dim(y)[1]))) == 1))
      force.cbind = TRUE
  }
  else
    force.list = TRUE

  if (force.list)
    return(cbind(ix = unlist(lapply(1:length(x), function(y) rep(y, length(x[[y]])))),
                 iix = unlist(lapply(1:length(x), function(y) if (length(x[[y]])>0) 1:length(x[[y]]) else NULL)),
                 unlist(x)))
  else if (force.rbind)
    return(cbind(ix = unlist(lapply(1:length(x), function(y) rep(y, nrow(x[[y]])))),
                 iix = unlist(lapply(1:length(x), function(y) if (nrow(x[[y]])>0) 1:nrow(x[[y]]) else NULL)),
                 do.call('rbind', x)))
  else if (force.cbind)
    return(t(rbind(ix = unlist(lapply(1:length(x), function(y) rep(y, ncol(x[[y]])))),
                   iix = unlist(lapply(1:length(x), function(y) if (ncol(x[[y]])>0) 1:ncol(x[[y]]) else NULL)),
                   do.call('cbind', x))))
}

#' Wrapper to \code{GenomicRanges::findOverlaps} with added functionality
#'
#' Returns \code{GRanges} of matches with two additional fields:
#' \itemize{
#' \code{$query.id} - index of matching query
#' \code{$subject.id} - index of matching subject
#' }
#' Optional \code{"by"} field is a character scalar that specifies a metadata column present in both query and subject
#' that will be used to additionally restrict matches, i.e. to pairs of ranges that overlap and also
#' have the same values of their \code{"by"} fields
#'
#' @name gr.findoverlaps
#' @importFrom GenomeInfoDb seqlengths seqlengths<-
#' @importFrom GenomicRanges values ranges width strand values<- strand<- seqnames
#' @importFrom data.table is.data.table := setkeyv
#' @importFrom IRanges findOverlaps
#' @importFrom S4Vectors queryHits subjectHits
#' @param query Query \code{GRanges} pile
#' @param subject Subject \code{GRanges} pile
#' @param ignore.strand Don't consider strand information during overlaps. \code{[TRUE]}
#' @param first Restrict to only the first match of the subject (default is to return all matches). \code{[FALSE]}
#' @param qcol \code{character} vector of query meta-data columns to add to results
#' @param scol \code{character} vector of subject meta-data columns to add to results
#' @param type \code{type} argument as defined by \code{IRanges::findOverlaps} (\code{"any"}, \code{"start"}, \code{"end"}, \code{"within"}, \code{"equal"}). \code{["any"]}
#' @param by Meta-data column to consider when performing overlaps [NULL]
#' @param return.type Select data format to return (supplied as character): \code{"same"}, \code{"data.table"}, \code{"GRanges"}. \code{["same"]}
#' @param ... Additional arguments sent to \code{IRanges::findOverlaps}.
#' @return \code{GRanges} pile of the intersection regions, with \code{query.id} and \code{subject.id} marking sources
#' @export
gr.findoverlaps = function(query, subject, ignore.strand = TRUE, first = FALSE,
                           qcol = NULL, ## any query meta data columns to add to result
                           scol = NULL, ## any subject meta data columns to add to resultx
                           type = 'any',
                           by = NULL,
                           return.type = 'same',
                           ...)
{

  subject.id = query.id = i.start = i.end = NULL ## for NOTE

  isdt <- any(class(query) == 'data.table' )
  if (return.type == 'same')
    return.type <- ifelse(isdt, 'data.table', 'GRanges')

  if (!return.type %in% c("data.table", "GRanges"))
    stop("return.type must be one of: same, data.table, GRanges")

  if (!((inherits(subject, 'GRanges') | inherits(subject, 'data.table')) & (inherits(query, 'GRanges') | inherits(query, 'data.table'))))
    stop('both subject and query have to be GRanges or data.table')

  if (!is.null(qcol))
    if (!all(qcol %in% names(values(query))))
      stop('Some qcol are not present in meta data of query')

  if (!is.null(scol))
    if (!all(scol %in% names(values(subject))))
      stop('Some scol are not present in meta data of subject')

  if (!is.null(by))
    if (!(by %in% names(values(query)) & by %in% names(values(subject))))
      stop('"by" field must be meta data column of both query and subject')

  if (is.data.table(query))
    query   <- dt2gr(query)
  if (is.data.table(subject))
    subject <- dt2gr(subject)

  ## perform the actual overlaps
  h <- GenomicRanges::findOverlaps(query, subject, type = type, ignore.strand = ignore.strand, ...)
  r <- ranges(h, ranges(query), ranges(subject))
  h.df <- data.table(start = start(r), end = end(r), query.id = queryHits(h),
                     subject.id = subjectHits(h), seqnames = as.character(seqnames(query)[queryHits(h)]))

  ## add the seqnames, and subset if have a "by"
  if (nrow(h.df) > 0) {
    if (!is.null(by)) {
      by.query <- values(query)[h.df$query.id, by]
      by.subject <- values(subject)[h.df$subject.id, by]
      keep.ix <- by.query == by.subject
      h.df <- h.df[keep.ix, ]
    }
  }

  ## if empty, return now
  if (nrow(h.df) == 0) {
    if (return.type == "GRanges")
      return(GRanges(seqlengths = seqlengths(query)))
    else
      return(data.table())
  }

  ## write the strand
  if (!ignore.strand)
    h.df$strand <- as.character(strand(query)[h.df$query.id])
  else
    h.df$strand = '*'

  ## limit to first hits?
  if (first)
    h.df = h.df[!duplicated(h.df$query.id), ]

  ## format into correct output format
  if (return.type=='GRanges') {
    ## this takes a while, I think because of validity check. Way to hack it with direct slot access via @?
    # out.gr = GRanges()
    # out.gr@seqnames = S4Vectors::Rle(factor(h.df$seqnames))
    # out.gr@ranges = IRanges::IRanges(h.df$start, h.df$end)
    # out.gr@elementMetadata = S4Vectors::DataFrame(query.id = h.df$query.id, subject.id = h.df$subject.id)
    # out.gr@strand = S4Vectors::Rle(h.df$strand)
    # out.gr@seqinfo = seqinfo(query)
    out.gr = GRanges(h.df$seqnames, IRanges(h.df$start, h.df$end), query.id = h.df$query.id, subject.id = h.df$subject.id, seqlengths = seqlengths(query),
                     strand = h.df$strand)

    if (!is.null(qcol))
      values(out.gr) = cbind(values(out.gr), values(query)[out.gr$query.id, qcol, drop = FALSE])

    if (!is.null(scol))
      values(out.gr) = cbind(values(out.gr), values(subject)[out.gr$subject.id, scol, drop = FALSE])

    return(sort(out.gr))
  } else { ## return data.table

    if (!is.null(qcol))
      h.df = cbind(h.df, as.data.table(as.data.frame(values(query))[h.df$query.id, qcol, drop = FALSE]))

    if (!is.null(scol))
      h.df = cbind(h.df, as.data.table(as.data.frame(values(subject))[h.df$subject.id, scol, drop = FALSE]))

    if ('i.start' %in% colnames(h.df))
      h.df[, i.start := NULL]

    if ('i.end' %in% colnames(h.df))
      h.df[, i.end := NULL]

    return(h.df)
  }
}


#' Alternative \code{GenomicRanges::match} that accepts \code{"by"} argument and \code{data.table} inputs
#'
#' Wrapper to \code{GenomicRanges::match} (uses \code{\link{gr.findoverlaps}})
#' @return Vector of length = \code{length(query)} with subject indices of *first* subject in query, or NA if none found.
#' This behavior is different from \code{\link{gr.findoverlaps}}, which will
#' return *all* indicies of subject in query (in the case of one query overlaps with multiple subject)
#' ... = additional args for findOverlaps (IRanges version)
#' @name gr.match
#' @param query Query \code{GRanges} pile
#' @param subject Subject \code{GRanges} pile
#' @param ... Additional arguments to be passed along to \code{\link{gr.findoverlaps}}.
#' @export
gr.match = function(query, subject, ...)
{

  ## if no need to call gr.findoverlaps, go to findOverlaps directly
  # if (!is.data.table(query) && !is.data.table(subject) && length(list(...)) == 0) {
  #   h = findOverlaps(query, subject)
  #   tmp = data.table(subject.id = subjectHits(h), query.id = queryHits(h))
  # } else {
    tmp = gr.findoverlaps(query, subject, ...)
  # }
  tmp = tmp[!duplicated(tmp$query.id)]
  out = rep(NA, length(query))
  out[tmp$query.id] = tmp$subject.id
  return(out)
}

#' @name %*%
#' @title Metadata join with coordinates as keys (wrapper to \code{\link{gr.findoverlaps}})
#' @description
#' Shortcut for gr.findoverlaps with \code{qcol} and \code{scol} filled in with all the query and subject metadata names.
#' This function is useful for piping \code{GRanges} operations together. Another way to think of %*% is as a
#' join of the metadata, with genomic coordinates as the keys. \cr
#' Example usage: \cr
#' x %*% y
#' @param x \code{GRanges}
#' @param y \code{GRanges}
#' @return \code{GRanges} containing every pairwise intersection of ranges in \code{x} and \code{y} with a join of the corresponding  metadata
#' @rdname grfo
#' @exportMethod %*%
#' @export
#' @importFrom methods setMethod
#' @author Marcin Imielinski
#' @docType methods
#' @aliases %*%,GRanges-method
#' @examples
#' gr.genes %*% gr.DNAase
#setGeneric('%*%', function(gr, ...) standardGeneric('%*%'))
setMethod("%*%", signature(x = "GRanges"), function(x, y) {
    gr = gr.findoverlaps(x, y, qcol = names(values(x)), scol = names(values(y)))
   return(gr)
})

#' Find overlaps between rearrangements represented by \code{GRangesList} objects
#'
#' Determines overlaps between two piles of rearrangement junctions, each \code{GRangesLists} of signed locus pairs,
#' against each other. returning a sparseMatrix that is T at entry ij if junction i overlaps junction j.
#'
#' @param ra1 \code{GRangesList} of pairs of signed ranges representing a rearrangement set
#' @param ra2 \code{GRangesList} of pairs of signed ranges representing a rearrangement set
#' @param pad Pad each breakpoint when considering overlaps. \code{[0]}
#' @param ... Additional arguments to be sent to \code{\link{findOverlaps}} (e.g. \code{ignore.strand})
#' @return \code{matrix} with the indices of \code{ra1} that overlap with \code{ra2} and vice-versa
#' @name ra.overlaps
#' @export
ra.overlaps = function(ra1, ra2, pad = 0, ...)
{
  bp1 = grl.unlist(ra1) + pad
  bp2 = grl.unlist(ra2) + pad
  h = GenomicRanges::findOverlaps(bp1, bp2, ...) ## was gr.findoverlaps
  ix = data.table(query.id = queryHits(h), subject.id = subjectHits(h))
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
    # if (arr.ind)
      return(matrix())
    # else
    #   return(sparseMatrix(length(ra1), length(ra2), x = 0))
  }

  rownames(tmp) = NULL

  colnames(tmp) = c('ra1.ix', 'ra2.ix')

  # if (arr.ind) {
    ro <- tmp[order(tmp[,1], tmp[,2]), ]
    if (class(ro)=='integer')
      ro <- matrix(ro, ncol=2, nrow=1, dimnames=list(c(), c('ra1.ix', 'ra2.ix')))
    return(ro)
  # } else {
  #   ro <- sparseMatrix(tmp[,1], tmp[,2], x = 1, dims = c(length(ra1), length(ra2)))
  #   return(ro)
  # }
}

#' Simplify granges by collapsing all non-overlapping adjacent ranges that share a given "field" value
#' (adjacent == adjacent in the input GRanges object)
#'
#' @param gr takes in gr or grl
#' @param field character scalar, corresponding to value field of gr. \code{[NULL]}
#' @param val \code{[NULL]}
#' @param include.val scalar logical, will include in out gr values field of first matching record in input gr. \code{[TRUE]}
#' @param split Split the output into \code{GRangesList} split by \code{"field"}. \code{[FALSE]}
#' @param pad Pad ranges by this amount before doing merge. [1], which merges contiguous but non-overlapping ranges.
#' @return Simplified GRanges with "field" populated with uniquely contiguous values
#' @export
gr.simplify = function(gr, field = NULL, include.val = TRUE, split = FALSE, pad = 1)
{
  tmp = as.logical(suppressWarnings(width(IRanges::pintersect(ranges(gr[-length(gr)]), ranges(gr[-1]+pad), resolve.empty = 'max.start'))>0) &
                     seqnames(gr[-length(gr)]) == seqnames(gr[-1]) & strand(gr[-length(gr)]) == strand(gr[-1]))

  tmp = as.vector(c(0, cumsum(!tmp)))

  if (!is.null(field))
    tmp = paste(tmp, values(gr)[, field])

  # if (!is.null(val))
  #   tmp = paste(tmp, val)

  r = base::rle(tmp)

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
