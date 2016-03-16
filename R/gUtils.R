#' DNAaseI hypersensitivity sites from UCSC Table Browser hg19,
#' subsampled to 10,000 sites
#'
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
#' library(gUtils)
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
#' library(gUtils)
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
#' Takes as input a data.table which must have the fields: start, end, strand, seqnames.
#' All of the remaining fields are added as meta data to the \code{\link[GenomicRanges]{GRanges}}.
#' Will throw an error if \code{data.table} does not contain seqnames, start and end
#' @param dt data.table to convert to GRanges
#' @return \code{\link[GenomicRanges]{GRanges}} object of length = nrow(dt)
#' @importFrom data.table
#'   data.table
#' @importFrom GenomicRanges
#'   GRanges
#'   mcols<-
#' @importFrom IRanges
#'    IRanges
#' @name dt2gr
#' @examples
#' library(gUtils)
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

#' Get GRanges corresponding to beginning of end
#'
#' Alternative to \code{flank} that will provide end positions *within* intervals
#'
#' @param x \code{GRanges} object to operate on
#' @param width [default = 1] Specify subranges of greater width including the start of the range.
#' @param force [default = F] Allows returned \code{GRanges} to have ranges outside of its \code{Seqinfo} bounds.
#' @param clip [default = F] Trims returned \code{GRanges} so that it does not extend beyond bounds of the input \code{GRanges}
#' @param ignore.strand [default = TRUE] If set to \code{FALSE}, will extend '-' strands from the other direction.
#' @return GRanges object of width 1 ranges representing end of each genomic range in the input.
#' @examples
#' library(gUtils)
#' gr.end(gr.DNAase, width=200, clip=TRUE)
#' @importFrom GenomeInfoDb
#'   seqlengths
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

#' Get the midpoint of range
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

# Round a set of GRanges to another set
#
# "rounds" a set of query ranges Q to a given interval set S using the following rule:
# 1) If q in Q is partially / fully within S then return intersection of q with S.
# 2) If q intersects multiple ranges in S and \code{up = F} then return the "first" range, otherwise the last range
# 3) If q in Q is fully outside of S (ie fully inside not S) then return the \code{start-1} (if \code{up = TRUE}) or \code{end+1} (if \code{up = F})
# of the matching range in not S
#
# @param Q Query \code{GRanges} (strand is ignored)
# @param S Subject \code{GRanges} (strand is ignored)
# @param up [default T] See description.
# @param parallel [default F] If \code{TRUE}, assumes Q and S are same length and this analysis is only performed between the corresponding Q and S pairs.
# @return Rounded \code{GRanges}
# @examples
# query   <- GRanges(1, IRanges(c(100,110),width=201), seqinfo=Seqinfo("1", 500))
# subject <- GRanges(1, IRanges(c(160,170),width=201), seqinfo=Seqinfo("1", 500))
# gr.round(query, subject)
# @export
# gr.round = function(Q, S, up = TRUE, parallel = FALSE)
#   {
#     str = strand(Q)
#     Q = gr.stripstrand(Q)
#     S = gr.stripstrand(S)
#     nS = gaps(S);
#     QS = gr.findoverlaps(Q, S)
#     tmp = gr.findoverlaps(Q, nS)
#     QnotS = nS[tmp$subject.id]
#     QnotS$query.id = tmp$query.id
#
#     if (parallel)
#       {
#         QS = QS[QS$query.id==QS$subject.id]
#         QnotS = QnotS[QnotS$query.id==QnotS$subject.id]
#       }
#
#     if (up)
#       suppressWarnings(end(QnotS) <- start(QnotS) <- end(QnotS)+1)
#     else
#       suppressWarnings(start(QnotS) <- end(QnotS) <- start(QnotS)-1)
#
#     suppressWarnings(out <- sort(grbind(QS, QnotS)))
#
#     if (up)
#       {
#         out = rev(out)
#         out = out[!duplicated(out$query.id)]
#       }
#     else
#       out = out[!duplicated(out$query.id)]
#
#     out = out[order(out$query.id)]
#     values(out) = values(Q)
#     names(out) = names(Q)
#     strand(out) = str;
#     return(out)
#   }
#

#' Generate random GRanges on genome
#'
#' Randomly generates non-overlapping GRanges with supplied widths on supplied genome.
#' Seed can be supplied with \code{set.seed}
#'
#' @param w Vector of widths (length of w determines length of output)
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
#' e.g. A \code{GRanges} with genomic coordinates 1:1,000,000-1,001,000 can get the first 20 and last 50 bases trimmed off with
#' \code{start = 20, end = 950}.
#' if end is larger than the width of the corresponding gr, then the corresponding output will only have end(gr) as its coordinate.
#'
#' This is a role not currently provided by the standard \code{GRanges} functions
#' (eg shift, reduce, restrict, shift, resize, flank etc)
#' @param gr \code{GRanges} to trim
#' @param starts number [1]
#' @param ends number [1]
#' @examples
#' ## trim the first 20 and last 50 bases
#' gr.trim(GRanges(1, IRanges(1e6, width=1000)), starts=20, ends=950)
#' ## return GRanges on 1:1,000,019-1,000,949
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
#' Samples k intervals of length "len" from a pile of \code{GRanges}.
#' If k is a scalar then will (uniformly) select k intervals from the summed territory of \code{GRanges}
#' If k is a vector of length(gr) then will uniformly select k intervals from each.
#' from a tiling of the set (and may do fewer than k samples if width(gr[i])<= k[i] *len)
#' If k[i] = NA, will return tiling of that interval, if k = NA will return tiling of the entire
#' gr's (with length len tiles).
#'
#' @param gr \code{GRanges} object defining the territory to sample from
#' @param k Number of ranges to sample
#' @param len Length of the \code{GRanges} element to produce [100]
#' @param replace If TRUE, will bootstrap, otherwise will sample without replacement. [TRUE]
#' @return GRanges of max length sum(k) [if k is vector) or k*length(gr) (if k is scalar) with labels indicating the originating range.
#'
#' @examples
#' ## sample 5 \code{GRanges} of length 10 each from territory of RefSeq genes
#' gr.sample(reduce(gr.genes), k=5, len=10)
#' @note This is different from overloaded sample() function implemented in GenomicRanges class, which just samples from a pile of GRanges
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
#' @param strip.empty Don't know. Default FALSE
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

#' Concatenate GRanges, robust to different \code{mcols}
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
#' @return Concatenated GRangesList
#' @examples
#' ## Concatenate
#' #grl.hiC2 <- grl.hiC[1:20]
#' #mcols(grl.hiC2)$test = 1
#' #grlbind(grl.hiC2, grl.hiC[1:30])
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

#' Prepend "chr" to GRanges seqlevels
#'
#' @param gr GRanges object to append 'chr' to
#' @return Identical \code{GRanges}, but with 'chr' prepended to each seqlevel
#' @examples
#' library(gUtils)
#' gr.chr(GRanges(c(1,"chrX"), IRanges(c(1,2), 1)))
#' @importFrom GenomeInfoDb
#'    seqlevels
#'    seqlevels<-
#' @export
gr.chr = function(gr)
  {
    if (any(ix <- !grepl('^chr', seqlevels(gr))))
      seqlevels(gr)[ix] = paste('chr', seqlevels(gr)[ix], sep = "")
    return(gr)
  }

#' Shortcut for \code{reduce(sort(gr.stripstrand(unlist(x))))}
#'
#' @param gr takes in gr or grl
#' @param pad asdf. Default 0
#' @param sort Flag to sort the output. Default TR#' @return GRanges
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
#' If mb will return as MB and round to "round"
#' @param gr \code{GRanges} pile to get intervals from
#' @param add.chr Prepend seqnames with "chr" [FALSE]
#' @param mb Round to the nearest megabase [FALSE]
#' @param round If \code{mb} supplied, how many digits to round to
#' @param other.cols TODO
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

#' gr.sub
#'
#' will apply gsub to seqlevels of gr, by default removing 'chr', and "0.1" suffixes, and replacing "MT" with "M"
#' @name gr.sub
#' @keywords internal
gr.sub = function(gr, a = c('(^chr)(\\.1$)', 'MT'), b= c('', 'M'))
  {
    tmp = mapply(function(x, y) seqlevels(gr) <<- gsub(x, y, seqlevels(gr)), a, b)
    return(gr)
  }

#' "Fixes" seqlengths / seqlevels
#'
#' If "genome" not specified will replace NA seq lengths in GR to reflect largest coordinate per seqlevel
#' and removes all NA seqlevels after this fix.
#'
#' if "genome" defined (ie as Seqinfo object, or a BSgenome, GRanges, GRnagesList object with populated seqlengths) then will replace
#' seqlengths in gr with those for that genome (and if drop = T, drop all ranges without
#' seqlevels in that genome)
#' @name gr.fix
#' @param gr \code{GRanges} object to fix
#' @param genome Genome to fix to: \code{Seqinfo}, \code{BSgenome}, \code{GRanges} (w/seqlengths), \code{GRangesList} (w/seqlengths)
#' @param gname Name of the genome (optional, just appends to \code{Seqinfo} of the output) [NULL]
#' @param drop Remove ranges that are not present in the supplied genome [FALSE]
#' @return \code{GRanges} pile with the fixed \code{Seqinfo}
#' @importFrom GenomeInfoDb
#'    Seqinfo
#'    seqinfo
#'    keepSeqlevels
#'    seqlevels
#'    seqlengths
#'    seqlevels<-
#'    seqlengths<-
#'    genome<-
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
        else
          {
            if (is.character(genome))
              genome = structure(rep(NA, length(genome)), names = genome)

            lens = structure(NA, names = union(names(genome), seqlevels(gr)));

            lens[seqlevels(gr)] = seqlengths(gr);
            lens[names(genome)[!is.na(genome)]] = pmax(lens[names(genome)[!is.na(genome)]], genome[!is.na(genome)], na.rm = TRUE)

            if (drop)
              lens = lens[names(genome)]
          }

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
#' Takes pile of GRanges and returns into a data.frame with nrow = length(gr) with each
#' representing the corresponding input range superimposed onto a single "flattened"
#' chromosome, with ranges laid end-to-end
#' @param gr \code{GRanges} to flatten
#' @param gap Number of bases between ranges on the new chromosome [0]
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
#' library(gUtils)
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

#' Tile ranges across a genomic range
#'
#' tiles interval (or whole genome) with segments of <= specified width.  Returns strandless gr
#' "tiles".
#'
#' input can be seqinfo object (in which case whole genome will be tiled);
#' if inputted grs overlap, will first reduce then tile.
#' @param gr \code{GRanges}, \code{seqlengths} or \code{seqinfo} range to tile. If has \code{GRanges} has overlaps, will reduce first.
#' @param w Width of each tile
#' @name gr.tile
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

#' Faster replacement for \code{GRanges} version of \code{findOverlaps}
#'
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
#' ... = additional args for \code{findOverlaps} (IRanges version)
#' @name gr.findoverlaps
#' @importFrom GenomeInfoDb
#'   seqlengths
#'   seqlengths<-
#' @importFrom GenomicRanges values ranges width strand values<- strand<- seqnames
#' @importFrom data.table is.data.table := setkeyv
#' @importFrom IRanges findOverlaps
#' @importFrom S4Vectors queryHits subjectHits
#' @param query Query \code{GRanges} pile
#' @param subject Subject \code{GRanges} pile
#' @param ignore.strand Don't consider strand information during overlaps [TRUE]
#' @param first TODO
#' @param qcol \code{character} vector of query meta-data columns to add to results
#' @param scol \code{character} vector of subject meta-data columns to add to results
#' @param foverlaps Use \code{data.table::foverlaps} instead of \code{IRanges::findOverlaps}. Overrules \code{pintersect}
#' @param pintersect Use \code{IRanges::pintersect} function. Useful for overlaps with many, many chromosomes. Default is TRUE if \code{length(unique(seqnames)) > 50}
#' @param verbose Increase the verbosity during runtime [FALSE]
#' @param type TODO
#' @param by Meta-data column to consider when performing overlaps [NULL]
#' @param return.type Select data format to return: \code{same}, \code{data.table}, \code{GRanges}
#' @param ... TODO
#' @return \code{GRanges} pile of the intersection regions, with \code{query.id} and \code{subject.id} marking sources
#' @export
gr.findoverlaps = function(query, subject, ignore.strand = TRUE, first = FALSE,
    qcol = NULL, ## any query meta data columns to add to result
    scol = NULL, ## any subject meta data columns to add to resultx
    ###max.chunk = 1e13, ## broken, see note below
    foverlaps = ifelse(is.na(as.logical(Sys.getenv('GRFO_FOVERLAPS'))), FALSE, as.logical(Sys.getenv('GRFO_FOVERLAPS'))) & exists('foverlaps'),
    pintersect = NA,
    verbose = FALSE,
    type = 'any',
    by = NULL,
    ###mc.cores = 1,
    return.type = 'same',
    ...)
  {

    subject.id = query.id = i.start = i.end = NULL ## for NOTE
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

    ### @param mc.cores Number of cores to use, if ranges exceed \code{max.chunk} [1]
    ### BROKEN @param max.chunk Maximum number of ranges to consider in one chunk (to keep down memory) [1e13]
    ## Jeremiah 3/8/16 - This is broken, as the output of gr.findoverlaps is different
    ## depending on max.chunk. Only support one chunk for now
    # if ((as.numeric(length(query)) * as.numeric(length(subject))) > max.chunk)
    #   {
    #     if (verbose)
    #       cat('Overflow .. computing overlaps in chunks.  Adjust max.chunk parameter to gr.findoverlaps to avoid chunked computation\n')
    #     chunk.size = floor(sqrt(max.chunk));
    #     ix1 = c(seq(1, length(query), chunk.size), length(query)+1)
    #     ix2 = c(seq(1, length(subject), chunk.size), length(subject)+1)
    #     ij = cbind(rep(1:(length(ix1)-1), length(ix2)-1), rep(1:(length(ix2)-1), each = length(ix1)-1))
    #     if (verbose)
    #       print(paste('Number of chunks:', nrow(ij)))
    #
    #     out = do.call('c', mclapply(1:nrow(ij),
    #         function(x)
    #                     {
    #                       if (verbose)
    #                         cat(sprintf('chunk i = %s-%s (%s), j = %s-%s (%s)\n', ix1[ij[x,1]], ix1[ij[x,1]+1]-1, length(query),
    #                                     ix2[ij[x,2]], (ix2[ij[x,2]+1]-1), length(subject)))
    #                       i.chunk = ix1[ij[x,1]]:(ix1[ij[x,1]+1]-1)
    #                       j.chunk = ix2[ij[x,2]]:(ix2[ij[x,2]+1]-1)
    #                       out = gr.findoverlaps(query[i.chunk], subject[j.chunk],  ignore.strand = ignore.strand, first = first, pintersect=pintersect, by = by, qcol = qcol, verbose = verbose, foverlaps = foverlaps, scol = scol, type = type, ...)
    #                       out$query.id = i.chunk[out$query.id]
    #                       out$subject.id = j.chunk[out$subject.id]
    #                       return(out)
    #                     }, mc.cores=mc.cores))
    #
    #     convert = FALSE
    #     if ((return.type == 'same' & is(query, 'data.table')) | return.type == 'data.table')
    #         out = gr2dt(out)
    #     return(out)
    #   }

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
                  querydt = gr2dt(query[, setdiff(by, c('seqnames', 'start', 'end', 'strand'))])
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
                  subjectdt = gr2dt(subject[, setdiff(by, c('seqnames', 'start', 'end', 'strand'))])
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

          querydt = querydt[, c('query.id', by), with = FALSE]
          subjectdt = subjectdt[, c('subject.id', by), with = FALSE]
          setkeyv(querydt, by)
          setkeyv(subjectdt, by)

          h.df = data.table::foverlaps(querydt, subjectdt, by.x = by, by.y = by, mult = 'all', type = 'any', verbose = verbose)
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
                  h.df <- data.frame(start = start(r), end = end(r), query.id = queryHits(h), subject.id = subjectHits(h), stringsAsFactors = FALSE);
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

               return(sort(out.gr))
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

#' Faster \code{GenomicRanges::match}
#'
#' Faster implementation of GRanges match (uses gr.findoverlaps)
#' returns indices of query in subject or NA if none found
#' ... = additional args for findOverlaps (IRanges version)
#' @name gr.match
#' @param query Query \code{GRanges} pile
#' @param subject Subject \code{GRanges} pile
#' @param max.slice Maximum number of ranges to consider at once [Inf]
#' @param verbose Increase the verbosity during runtime
#' @param mc.cores Number of cores to use, if ranges exceed \code{max.slice}
#' @param ... Additional arguments to be passed along to \code{gr.findoverlaps}
#' @importFrom parallel mclapply
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
#' Given two tilings of the genome (e.g. at different resolution)
#' query and subject outputs a length(query) list whose items are integer vectors of indices in subject
#' overlapping that overlap that query (strand non-specific)
#'
#' @param query Query \code{GRanges} pile
#' @param subject Subject \code{GRanges} pile
#' @param verbose Increase the verbosity of the output
#' @return TODO
#' @note Assumes that input query and subject have no gaps (including at end) or overlaps, i.e. ignores end()
#' coordinates and only uses "starts"
#' @export
gr.tile.map = function(query, subject, verbose = FALSE)
  {
    ix.q = order(query)
    ix.s = order(subject)

    q.chr = as.character(seqnames(query))[ix.q]
    s.chr = as.character(seqnames(subject))[ix.s]

    ql = split(ix.q, q.chr)
    sl = split(ix.s, s.chr)

    #tmp = mcmapply(
    tmp <- lapply(
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
      #}, ql, sl[names(ql)], mc.cores = mc.cores, SIMPLIFY = FALSE)
      }, ql, sl[names(ql)], SIMPLIFY = FALSE)


    m = munlist(tmp)[, -c(1:2), drop = FALSE]
    out = split(m[,2], m[,1])[as.character(1:length(query))]
    names(out) = as.character(1:length(query))
    return(out)
  }

#' Faster version of \code{GRanges} \code{over}
#'
#' Uses \code{\link{gr.findoverlaps}} for a faster \code{over}
#' by = column name in query and subject that we additionally control for a match (passed on to gr.findoverlaps)
#' @name gr.in
#' @param query Query \code{GRanges} pile
#' @param subject Subject \code{GRanges} pile
#' @param ... Additional arguments to pass to \code{\link{gr.findoverlaps}}
#' @return logical vector if query range i is found in any range in subject
#' @export
gr.in = function(query, subject, ...)
  {
    tmp = gr.findoverlaps(query, subject, ...)
    out = rep(FALSE, length(query))
    out[tmp$query.id] = TRUE

    return(out)
   }

#' Dice up \code{GRanges} into width 1 \code{GRanges} spanning the input (warning can produce a very large object)
#'
#' @param gr \code{GRanges} object to dice
#' @name gr.dice
#' @importFrom S4Vectors Rle
#' @importFrom GenomeInfoDb seqlengths
#' @importFrom GenomicRanges seqnames width strand values<- values strand<- distance
#' @return \code{GRangesList} where kth element is a diced pile of \code{GRanges} from kth input \code{GRanges}
#' @examples
#' library(gUtils)
#' library(S4Vectors)
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
#' if max.dist = TRUE then will replace min with max above
#' @param gr1 First \code{GRanges}
#' @param gr2 Second \code{GRanges}
#' @param ignore.strand Don't required elements be on same strand to avoid NA [FALSE]
#' @param ... Additional arguments to be supplied to \code{GenomicRanges::distance}
#' @return Matrix with the pairwise distances, with \code{gr1} on rows and \code{gr2} on cols
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
#' @param ... Additional parameters to be passed on to \code{gr.findoverlaps}
#' @name grl.in
#' @export
grl.in <- function(grl, windows, some = FALSE, only = FALSE, ...)
  {
    grl.iid = grl.id = NULL ## for getting past NOTE

    if (length(grl)==0)
      return(logical())

    if (length(windows)==0)
      return(rep(T, length(grl)))

    numwin = length(windows);
    gr = grl.unlist(grl)
    m = gr2dt(gr.findoverlaps(gr, windows, ...))

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

#' grl.split
#'
#' splits GRL's with respect to their seqnames and strand (default), returning
#' new grl whose items only contain ranges with a single seqname / strand
#'
#' can also split by arbitrary (specified) genomic ranges value fields
#' @name grl.split
#' @keywords internal
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

#' Robust unlisting of \code{GRangesList} that keeps track of origin
#'
#' Does a "nice" unlist of a \code{GRangesList} object adding a field "grl.ix" denoting which element of the \code{GRangesList}
#' each \code{GRanges} corresponds to and a field grl.iix which saves the (local) index that that gr was in its corresponding grl item
#'
#' In this way, \code{grl.unlist} is reversible, while \code{unlist} is not.
#' @name grl.unlist
#' @importFrom BiocGenerics unlist
#' @param grl \code{GRangeList} object to unlist
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

# grl.span
#
# Returns GRanges object representing the left / right extent of each GRL item.  In case of "chimeric" GRL items (ie that map
# to two chromosomes) there are two options:
# (1) specify "chr" chromosome as argument to subset GRL's that are on that chromosome, and compute GRL extents from this, any GRL
#     full outside of that chromosome will get a 0 width GRL
# (2) (default) allow chimeric GRL items to get an extent that is with respect to the first chromosome in that GRL
#
# If a grl item contains ranges that lie on different chromosomes, then corresponding grange will have chromosome "NA" and IRange(0, 0)
# @name grl.span
# @keywords internal
# grl.span = function(grl, chr = NULL, ir = FALSE, keep.strand = TRUE)
#   {
#     if (is.null(names(grl)))
#       names(grl) = 1:length(grl);
#
#     tmp = tryCatch(as.data.frame(grl), error = function(e) e)
#
#     if (inherits(tmp, 'error')) ## gr names are screwy so do some gymnastics
#       {
#         if (is.null(names(grl)))
#           names.grl = 1:length(grl)
#         else
#           names.grl = names(grl);
#
#         element = as.character(Rle(names.grl, sapply(grl, length)))
#         tmp.gr = unlist(grl)
#         names(tmp.gr) = NULL;
#         tmp = as.data.frame(tmp.gr);
#         tmp$element = element;
#       }
#
#     if (is.null(chr))
#       {
#         chrmap = stats::aggregate(formula = seqnames ~ element, data = tmp, FUN = function(x) x[1]);
#         chrmap = structure(as.character(chrmap[,2]), names = chrmap[,1])
#
#         if (keep.strand)
#           {
#             strmap = stats::aggregate(formula = as.character(strand) ~ element, data = tmp, FUN =
#               function(x) {y = unique(x); if (length(y)>1) return('*') else y[1]})
#             strmap = structure(as.character(strmap[,2]), names = strmap[,1])
#             str = strmap[names(grl)];
#           }
#         else
#           str = '*'
#
#         tmp = tmp[tmp$seqnames == chrmap[tmp$element], ]; ## remove all gr from each GRL item that don't map to the chr of the first gr
#         chr = chrmap[names(grl)];
#         out.gr = GRanges(chr, IRanges(1,0), seqlengths = seqlengths(grl), strand = str)
#       }
#     else
#       {
#         if (length(chr)>1)
#           warning('chr has length greater than 1, only the first element will be used')
#         tmp = tmp[tmp$seqnames == chr[1], ]
#         out.gr = rep(GRanges(chr, IRanges(1, 0)), length(grl)) # missing values
#       }
#
#     if (nrow(tmp)>0)
#       {
#         tmp = split(GRanges(tmp$seqnames, IRanges(tmp$start, tmp$end)), tmp$element)
#         out.gr[match(names(tmp), names(grl))] = GRanges(chr[names(tmp)],
#                 IRanges(sapply(start(tmp), min), sapply(end(tmp), max)), strand = strand(out.gr)[match(names(tmp), names(grl))]);
#         names(out.gr) = names(grl)
#       }
#     return(out.gr)
#   }


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


#' Improved \code{rbind} for intersecting/union columns of data.frames or data.tables
#'
#' Like \code{rbind}, but takes the intersecting columns of the data
#' rrbind = function(df1, df2, [df3 ... etc], )
#' @param ... Any number of \code{data.frame} or \code{data.table} objects
#' @param union Take union of columns (and put NA's for columns of df1 not in df2 and vice versa). [TRUE]
#' @param as.data.table Return the binded data as a \code{data.table}
#' @return A\code{data.frame} or \code{data.table} of the rbind operation
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


#' munlist
#'
#' unlists a list of vectors, matrices, data frames into a n x k matrix
#' whose first column specifies the list item index of the entry
#' and second column specifies the sublist item index of the entry
#' and the remaining columns specifies the value(s) of the vector
#' or matrices.
#'
#' force.cbind = TRUE will force concatenation via 'cbind'
#' force.rbind = TRUE will force concatenation via 'rxsbind'
#' @keywords internal
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

  if (inherits(seg, 'RangedData') | inherits(seg, 'GRanges') | inherits(seg, 'IRanges'))
  {
    val = as.data.frame(values(seg));
    values(seg) = NULL;
    seg = as.data.frame(seg, row.names = NULL);  ## returns compressed iranges list
    seg$seqnames = as.character(seg$seqnames)
  }
  else
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


