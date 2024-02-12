rowid=.N=colid=V1=inner=startFrom=breakAt=upTo=qid2=tmp.start=subject.id=query.id=transition=ix=nacluster=clusterlength=ix.new=slen=slev=clilp=val=query.id=subject.id=lev=len=cs=type=group=run=chr=ra1.ix=ra2.ix=pos2=chr=V1=NULL

#' DNAaseI hypersensitivity sites for hg19A
#'
#' DNAaseI hypersensitivity sites from UCSC Table Browser hg19,
#' subsampled to 10,000 sites
#' @name example_dnase
#' @docType data
#' @keywords data
#' @format \code{GRanges}
NULL




#' RefSeq genes for hg19
#'
#' RefSeq genes with exon count and name
#' @name example_genes
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




#' @name hg_seqlengths
#' @title Output standard human genome seqlengths
#' @description
#'
#' Outputs a standard seqlengths for human genome +/- "chr".
#'
#' @note A default genome can be set with the environment variable DEFAULT_GENOME. This
#' can be the full namespace of the genome  e.g.: \code{DEFAULT_GENOME=BSgenome.Hsapiens.UCSC.hg19::Hsapiens} OR  a URL / file path pointing to a chrom.sizes text file (e.g. http://genome.ucsc.edu/goldenpath/help/hg19.chrom.sizes) specifying a genome definition
#' @param genome A \code{BSgenome} or object with a \code{seqlengths} accessor. Default is hg19, but loads with warning unless explicitly provided
#' @param chr boolean Flag for whether to keep "chr". (default = FALSE)
#' @param include.junk boolean Flag for whether to not trim to only 1-22, X, Y, M. (default = FALSE)
#' @return Named integer vector with elements corresponding to the genome seqlengths
#' @importFrom utils read.delim
#' @importFrom stats setNames
#' @importFrom utils relist
#' @importFrom methods as is
#' @importFrom S4Vectors elementNROWS
#' @author Marcin Imielinski
#' @export
hg_seqlengths = function(genome = NULL, chr = TRUE, include.junk = FALSE)
{
  sl = NULL
  dbs = ''

  if (!is.null(genome))
    dbs = genome
  else
  {
    if (nchar(Sys.getenv("DEFAULT_GENOME"))>0)
       dbs = Sys.getenv("DEFAULT_GENOME")
    else if (nchar(Sys.getenv("DEFAULT_BSGENOME"))>0)
      dbs = Sys.getenv("DEFAULT_BSGENOME")   
  }

  if (nchar(dbs) == 0)
  {
    warning('hg_seqlengths: supply genome seqlengths or set default with env variable DEFAULT_GENOME (e.g. Sys.setenv(DEFAULT_GENOME = "BSgenome.Hsapiens.UCSC.hg19::Hsapiens").  DEFAULT_BSGENOME can also be set to a path or URL of a tab delimited text *.chrom.sizes file')
    return(NULL)
  } else{
    tmp = suppressWarnings(tryCatch(read.delim(dbs, header = FALSE), error= function(e) NULL))
    if (is.null(tmp))
    {
      genome = tryCatch(eval(parse(text=dbs)), error = function(e) NULL)
      if (is.null(genome)){
        stop(sprintf("Error loading %s as chrom.sizes or BSGenome library ...\nPlease check DEFAULT_GENOME setting and set to either an R library BSGenome object or a valid http URL or filepath pointing to a chrom.sizes tab delimited text file.", dbs))
                }
    }else{
      sl = structure(tmp[,2], names = as.character(tmp[,1]))
    }
  }
 
    if (is.null(sl)){
        sl = seqlengths(genome)
    }

    if (!chr){
        names(sl) = gsub('chr', '', names(sl))
    }

    if (!include.junk){
        sl = sl[nchar(names(sl))<=8]
    }

    return(sl)
}





#' @name gr2dt
#' @title Converts \code{GRanges} to \code{data.table}
#' @description
#'
#' Converts \code{GRanges} to \code{data.table}
#' and a field grl.iix which saves the (local) index that that gr was in its corresponding grl item
#'
#' @param x \code{GRanges} to convert
#' @return data.table of GRanges columns ('seqnames', 'start', 'end', 'strand', 'width') and metadata columns
#' @export
gr2dt = function(x)
{
    ## new approach just directly instantiating data table
    cmd = 'data.frame(';

    if (is(x, 'GRanges'))
    {
        ## as.data.table complains if duplicated row names
        if (any(duplicated(names(x)))){
            names(x) <- NULL
        }

        was.gr = TRUE
        f = c('seqnames', 'start', 'end', 'strand', 'width')
        f2 = c('as.character(seqnames', 'c(start', 'c(end', 'as.character(strand', 'as.numeric(width')
        cmd = paste(cmd, paste(f, '=', f2, '(x))', sep = '', collapse = ','), sep = '')
        value.f = names(values(x))
    } else {
        was.gr = FALSE
        value.f = names(x)
    }

    if (length(value.f)>0)
    {
        if (was.gr){
            cmd = paste(cmd, ',', sep = '')
        }
        class.f = sapply(value.f, function(f) eval(parse(text=sprintf("class(x$'%s')", f))))

        .StringSetListAsList = function(x){
            tmp1 = as.character(unlist(x))
            tmp2 = rep(1:length(x), S4Vectors::elementNROWS(x))
            return(split(tmp1, tmp2))
        }

        ## take care of annoying S4 / DataFrame / data.frame (wish-they-were-non-)issues
        as.statement = ifelse(grepl('Integer', class.f), 'as.integer',
                       ifelse(grepl('Character', class.f), 'as.character',
                       ifelse(grepl('((StringSet)|(Compressed.*))List', class.f), '.StringSetListAsList',
                       ifelse(grepl('StringSet$', class.f), 'as.character',
                       ifelse(grepl('factor$', class.f), 'as.character',
                       ifelse(grepl('List', class.f), 'as.list',
                       ifelse(grepl('factor', class.f), 'as.character',
                       ifelse(grepl('List', class.f), 'as.list', 'c'))))))))
        cmd = paste(cmd, paste(value.f, '=', as.statement, "(x$'", value.f, "')", sep = '', collapse = ','), sep = '')
    }

    cmd = paste(cmd, ')', sep = '')


    out = tryCatch(data.table::as.data.table(eval(parse(text =cmd))), error = function(e) NULL)    

    nr = 0
    if (!is.null(out)) { nr = nrow(out) }

    if (nr != length(x))
    {
        out = as.data.table(x)
    }

    return(out)
}




#' @name gr.start
#' @title Get GRanges corresponding to beginning of range
#' @description
#'
#' Get GRanges corresponding to beginning of range
#'
#' @param x \code{GRanges} object to operate on
#' @param width integer Specify subranges of greater width including the start of the range. (default = 1)
#' @param force boolean Allows returned \code{GRanges} to have ranges outside of its \code{Seqinfo} bounds. (default = FALSE)
#' @param ignore.strand boolean If set to \code{FALSE}, will extend '-' strands from the other direction (default = TRUE)
#' @param clip boolean Trims returned \code{GRanges} so that it does not extend beyond bounds of the input \code{GRanges} (default = TRUE)
#' @return \code{GRanges} object of width 1 ranges representing start of each genomic range in the input.
#' @importFrom GenomicRanges GRanges
#' @examples
#'
#' gr.start(example_dnase, width=200)
#' gr.start(example_dnase, width=200, clip=TRUE)
#'
#' @export
gr.start = function(x, width = 1, force = FALSE, ignore.strand = TRUE, clip = TRUE)
{
    if (length(x)==0){
        return(x)
    }

    width = pmax(width, 1)

    if (any(seqlengths(x)==0) | any(is.na(seqlengths(x)))){
        warning('Warning: Check or fix seqlengths, some are equal 0 or NA, may lead to negative widths')
    }

    .grstart = function(x)
    {
        if (force)
        {
            if (ignore.strand)
            {
                st = as.vector(start(x))
                en = as.vector(start(x))+width-1
            } else {
                st = ifelse(as.character(strand(x)) %in% c('*', '+'),
                            as.vector(start(x)),
                            as.vector(end(x))-width+1)

                en = ifelse(as.character(strand(x)) %in% c('*', '+'),
                            as.vector(start(x))+width-1,
                            as.vector(end(x))
                            )
            }
        } else{
            if (ignore.strand){
                st = start(x)
                en = pmin(as.vector(start(x))+width-1, seqlengths(x)[as.character(seqnames(x))], na.rm = TRUE)
            } else{
                st = ifelse(as.character(strand(x)) %in% c('*', '+'),
                            as.vector(start(x)),
                            pmax(as.vector(end(x))-width+1, 1)
                            )

                en = ifelse(as.character(strand(x)) %in% c('*', '+'),
                            pmin(as.vector(start(x))+width-1, seqlengths(x)[as.character(seqnames(x))], na.rm = TRUE),
                            as.vector(end(x)))
            }
        }

        if (clip){
            en = pmin(en, end(x))
            st = pmax(st, start(x))
        }

        ir = IRanges(st, en)
    }

    ir = tryCatch(.grstart(x), error = function(e) NULL)

    if (is.null(ir)){
        warning("Warning: One or more ranges are out of bounds on seqlengths, fixing and rerunning")
        x = gr.fix(x)
        ir = .grstart(x)
    }

    out = GRanges(seqnames(x), ir, seqlengths = seqlengths(x), strand = strand(x))


    values(out) = values(x)
    return(out)
}




#' @name dt2gr
#' @title Convert data.table to GRanges
#' @description
#'
#' Takes as input a data.table which must have the following fields: \code{start}, \code{end}, \code{strand}, \code{seqnames}. Will throw
#' an error if any one of these is not present.
#' All of the remaining fields are added as metadata to the \code{GRanges}.
#'
#' @param dt data.table or data.frame to convert to \code{GRanges}
#' @param seqlengths named integer vector representing genome (default = hg_seqlengths())
#' @param seqinfo seqinfo of output GRanges object
#' @return \code{GRanges} object of \code{length = nrow(dt)}
#' @importFrom data.table data.table
#' @importFrom GenomicRanges GRanges mcols
#' @importFrom IRanges IRanges
#' @examples
#' converted_gr = dt2gr(data.table(start=c(1,2), seqnames=c("X", "1"), end=c(10,20), strand = c('+', '-')))
#' @export
dt2gr = function(dt, key = NULL, seqlengths = NULL, seqinfo = Seqinfo()) {

  if (!inherits(dt, 'data.frame') & !inherits(dt, 'data.table')){
    stop("Error: Input needs to be data.table or data.frame")
  }

  if (!is(dt, 'data.table'))
    dt = as.data.table(dt)

  if (is.null(seqlengths))
    seqlengths = seqlengths(Seqinfo())

  if ((!is.integer(seqlengths) && !is.numeric(seqlengths)) || is.null(names(seqlengths)))
    stop('seqlengths must be a named integer vector')

  out = NULL;
  tryCatch({
    dt$seqnames = as.character(dt$seqnames)
    sl = dt[, max(end), keyby = seqnames]
    nms = setdiff(union(sl$seqnames, names(seqlengths)), NA)
    seqlengths = structure(pmax(sl[nms, V1], seqlengths[nms], na.rm = TRUE), names = nms)

    rr <- IRanges(dt$start, dt$end)
    if (!'strand' %in% colnames(dt)){
      dt$strand <- '*'
    }
    sf <- factor(dt$strand, levels=c('+', '-', '*'))
    ff <- factor(dt$seqnames, levels=unique(dt$seqnames))

    out <- GRanges(seqnames=ff, ranges=rr, strand=sf, seqlengths = seqlengths)
    if (inherits(dt, 'data.table')){
      mc <- as.data.frame(dt[, setdiff(colnames(dt),
                                       c('start', 'end', 'seqnames', 'strand',
                                         'width', 'seqlevels', 'seqlengths', 'element')),
                             with=FALSE])
    } else if (inherits(dt, 'data.frame')){
      mc <- as.data.frame(dt[, setdiff(colnames(dt),
                                       c('start', 'end', 'seqnames', 'strand',
                                         'width', 'seqlevels', 'seqlengths', 'element')),
                             drop = FALSE])
    }

    if (nrow(mc)){
      mcols(out) <- mc
    }
  }, error = function(e) NULL)

  if (is.null(out)){
    warning('Warning: Coercing to GRanges via non-standard columns')
    out = seg2gr(dt, seqlengths)
  }
  if ("width" %in% names(values(out))){
    out$width = NULL
  }

  return(out)
}




#' @name gr.end
#' @title Get the right ends of a \code{GRanges}
#' @description
#'
#' Alternative to \code{GenomicRanges::flank} that will provide end positions *within* intervals
#'
#' @param x \code{GRanges} object to operate on
#' @param width integer Specify subranges of greater width including the start of the range. (default = 1)
#' @param force boolean Allows returned \code{GRanges} to have ranges outside of its \code{Seqinfo} bounds. (default = FALSE)
#' @param ignore.strand boolean If set to \code{FALSE}, will extend '-' strands from the other direction. (default = TRUE)
#' @param clip boolean Trims returned \code{GRanges} so that it does not extend beyond bounds of the input (default = TRUE)
#' @return GRanges object of width = \code{width} ranges representing end of each genomic range in the input.
#' @examples
#' gr.end(example_dnase, width=200, clip=TRUE)
#' @importFrom GenomeInfoDb seqlengths
#' @importFrom GenomicRanges strand seqnames values<- values
#' @author Marcin Imielinski
#' @export
gr.end = function(x, width = 1, force = FALSE, ignore.strand = TRUE, clip = TRUE)
{
    if (length(x)==0){
        return(x)
    }

    if (any(seqlengths(x)==0) | any(is.na(seqlengths(x)))){
        warning('Warning: Check or fix seqlengths, some are equal 0 or NA, may lead to negative widths')
    }

    width = pmax(width, 1)

    .grend = function(x)
    {
        if (force)
        {
            if (ignore.strand)
            {
                st = as.vector(end(x))-width+1
                en = as.vector(end(x))
            } else{
                st = ifelse(as.character(strand(x)) %in% c('*', '+'),
                            as.vector(end(x))-width+1,
                            as.vector(start(x)))

                en = ifelse(as.character(strand(x)) %in% c('*', '+'),
                            as.vector(end(x)),
                            as.vector(start(x))+width-1)
            }
            out = GRanges(seqnames(x), IRanges(st, en), seqlengths = seqlengths(x), strand = strand(x))
        } else{
            if (ignore.strand)
            {
                st = pmax(as.vector(end(x))-width+1, 1)
                en = as.vector(end(x))
            } else{
                st = ifelse(as.character(strand(x)) %in% c('*', '+'),
                            pmax(as.vector(end(x))-width+1, 1),
                            as.vector(start(x)))

                en = ifelse(as.character(strand(x)) %in% c('*', '+'),
                            as.vector(end(x)),
                            pmin(as.vector(start(x))+width-1, seqlengths(x)[as.character(seqnames(x))], na.rm = TRUE))
            }

            out = GRanges(seqnames(x), IRanges(st, en), seqlengths = seqlengths(x), strand = strand(x))
        }

        if (clip){
            en = pmin(en, end(x))
            st = pmax(st, start(x))
        }

        return(IRanges(st, en))
    }

    ir = tryCatch(.grend(x), error = function(e) NULL)

    if (is.null(ir)){
        warning("Warning: One or more ranges are out of bounds on seqlengths, fixing and rerunning")
        x = gr.fix(x)
        ir = .grend(x)
    }

    out = GRanges(seqnames(x), ir, seqlengths = seqlengths(x), strand = strand(x))

    values(out) = values(x)
    return(out)
}




#' @name gr.mid
#' @title Get the midpoints of \code{GRanges} ranges
#' @description
#'
#' Get the midpoints of \code{GRanges} ranges
#'
#' @param x \code{GRanges} object to operate on
#' @return \code{GRanges} of the midpoint, calculated from \code{floor(width(x)/2)}
#' @importFrom GenomicRanges start<- end<- start end
#' @examples
#' gr.mid(GRanges(1, IRanges(1000,2000), seqinfo=Seqinfo("1", 2000)))
#' @export
gr.mid = function(x)
{
    start(x) = end(x) = rowMeans(cbind(start(x), end(x)))
    return(x)
}




#' @name gr.rand
#' @title Generate random \code{GRanges} on genome
#' @description
#'
#' Randomly generates non-overlapping \code{GRanges} with supplied widths on supplied genome.
#' Seed can be supplied with \code{set.seed}
#'
#' @param w vector of widths (length of \code{w} determines length of output)
#' @param genome GRanges, GRangesList, or Seqinfo genome. Default is "hg19" from the \code{BSGenome} package.
#' @return \code{GRanges} with random intervals on the specifed "chromosomes"
#' @note This function is currently quite slow, needs optimization
#' @importFrom GenomeInfoDb seqinfo seqnames<-
#' @importFrom GenomicRanges gaps ranges ranges<-
#' @examples
#'
#' Sys.setenv(DEFAULT_GENOME = "http://mskilab.com/gUtils/hg19/hg19.chrom.sizes")
#'
#' ## Generate 5 non-overlapping regions of width 10 on hg19
#' gr.rand(rep(10,5), si2gr(hg_seqlengths()))
#'
#' @author Marcin Imielinski
#' @export
gr.rand = function(w, genome)
{
    if (!is(genome, 'Seqinfo')){
        genome = seqinfo(genome)
    }

    sl = seqlengths(genome);
    available = si2gr(genome);

    out = GRanges(rep(names(sl)[1], length(w)), IRanges(rep(1, length(w)), width = 1), seqlengths = seqlengths(genome));
    for (i in 1:length(w))
    {
        if (i == 1){
            available = si2gr(genome)
        } else{
            available = gaps(out[1:(i-1)])
            available = available[strand(available)=='*']
        }

        available = available[width(available)>w[i]]

        if (length(available)>0){
            end(available) = end(available)-w[i]
            starts = c(1, cumsum(as.numeric(width(available))+1))
            rstart = ceiling(stats::runif(1)*starts[length(starts)])-starts
            rind = max(which(rstart>0))
            new.chr = seqnames(available[rind])
            new.ir = IRanges(rstart[rind]+start(available[rind])-1, width = w[i])

            ## FIX: this is the slowest part
            seqnames(out)[i] = new.chr;
            ranges(out)[i] = new.ir;
        } else{
            stop('Error: Allocation failed. Supplied widths are likely too large')
        }
    }

    return(out)
}




#' @name gr.trim
#' @title Trims pile of \code{GRanges} relative to the specified <local> coordinates of each range
#' @description
#'
#' Example: \code{GRanges} with genomic coordinates 1:1,000,000-1,001,000 can get the first 20 and last 50 bases trimmed off with
#' \code{start = 20, end = 950}.
#' if end is larger than the width of the corresponding gr, then the corresponding output will only have \code{end(gr)} as its coordinate.
#'
#' This is a role not currently provided by the standard \code{GenomicRanges} functions
#' (e.g. \code{shift}, \code{reduce}, \code{restrict}, \code{shift}, \code{resize}, \code{flank})
#'
#' @param gr \code{GRanges} to trim
#' @param starts integer Beginning of interval trimmed; Number of bases to trim off of the front\code{[1]}
#' @param ends integer End of interval trimmed
#' @examples
#'
#' ## trim the first 20 and last 50 bases
#' gr.trim(GRanges(1, IRanges(1e6, width=1000)), starts=20, ends=950)
#' ## return value: GRanges on 1:1,000,019-1,000,949
#'
#' @return GRanges with trimmed intervals relative to the specified <local> coordinates of each range
#' @export
gr.trim = function(gr, starts=1, ends=1)
{
    starts = cbind(1:length(gr), starts)[, 2]
    ends = cbind(1:length(gr), ends)[, 2]

    ends = pmax(starts, ends);
    ends = pmin(ends, width(gr));
    en = start(gr) + ends - 1;

    st = start(gr)+starts-1;
    st = pmin(st, en);

    out = GRanges(seqnames(gr), IRanges(st, en), seqlengths = seqlengths(gr), strand = strand(gr))

    values(out) = values(gr)

    return(out)
}




#' @name gr.sample
#' @title Randomly sample \code{GRanges} intervals within territory
#' @description
#'
#' Samples \code{k} intervals of length "len" from a pile of \code{GRanges}.
#' \itemize{
#' \item If k is a scalar then will (uniformly) select k intervals from the summed territory of \code{GRanges}
#' \item If k is a vector of length(gr) then will uniformly select k intervals from each.
#' }
#'
#' @param gr Granges defining the territory to sample from
#' @param k integer Number of ranges to sample
#' @param wid integer Length of the \code{GRanges} element to produce (default = 100)
#' @param replace boolean If TRUE, will bootstrap. If FALSE, otherwise will sample without replacement. (default = TRUE)
#' @return GRanges of max length sum(k) [if k is vector) or k*length(gr) (if k is scalar) with labels indicating the originating range.
#' @examples
#'
#' ## sample 5 \code{GRanges} of length 10 each from territory of RefSeq genes
#' gr.sample(reduce(example_genes), k=5, wid=10)
#'
#' @note This is different from \code{GenomicRanges::sample} function, which just samples from a pile of \code{GRanges}
#' @author Marcin Imielinski
#' @export
gr.sample = function(gr, k, wid = 100, replace = TRUE)
{
    if (!inherits(gr, 'GRanges')){
        gr = si2gr(gr)
    }

    if (length(k)==1){
        gr$ix.og = 1:length(gr)
        gr = gr[width(gr)>=wid]
        if (length(gr)==0){
            stop('Error: Input territory has zero regions of sufficient width')
        }
        gr.f = as.data.table(gr.flatten(gr.trim(gr, starts = 1, ends = width(gr)-wid), gap = 0))
        terr = sum(gr.f$end-gr.f$start + 1)
        st = gr.f$start;

        gr.f$ix = 1:nrow(gr.f)


        if (!replace){
            if (!is.na(k)){
                s = sort(wid*sample(floor(terr/wid), k, replace = FALSE))
            } else{
                s = sort(seq(1, terr, wid))
            }
        } else{
            s = sort(terr*stats::runif(k))
        }

        si = rep(NA, length(s))
        ## concatenate input and random locations
        gr.r = data.table(start = s, end = s+wid-1, ix = as.numeric(NA))
        tmp = rbind(gr.f,  gr.r, fill = TRUE)[order(start), ]

        ## match random events to their last "ix"
        ## find all non NA to NA transitions
        ## and tag all NA runs with the same number
        tmp[, transition := is.na(c(NA, ix[-length(ix)])) != is.na(ix) & is.na(ix)]
        tmp[, nacluster := ifelse(is.na(ix), cumsum(transition), 0)]
        tmp[, clusterlength := length(start), by = nacluster]
        tmp[, ix.new := ifelse(transition, c(NA, ix[-length(ix)]), ix)]
        tmp[, ix.new := ix.new[1], by = nacluster]

        ## now lift up random ranges to their original coordinates
        tmp[is.na(ix), ":="(seqnames = as.character(seqnames(gr)[ix.new]),
                            pos1 = start + start(gr)[ix.new]-gr.f$start[ix.new],
                            pos2 = end + start(gr)[ix.new]-gr.f$start[ix.new])]

        tmp = tmp[is.na(ix), ]

        out = GRanges(tmp$seqnames, IRanges(tmp$pos1, tmp$pos2), strand = '*', seqlengths = seqlengths(gr))
        out$query.id = gr.f$ix.og[tmp$ix.new]
        return(out)
    } else{
        gr.df = data.frame(chr = as.character(seqnames(gr)), start = start(gr), end = end(gr))
        ##gr.df$k = k;
        gr.df$length = wid
        gr.df$replace = replace
        tmp = lapply(1:length(gr), function(i){
            if (!gr.df$replace[i]){
                if (!is.na(k[i])){
                    w = floor(width(gr)[i]/wid)
                    k[i] = min(k[i], w)
                    if (k[i]>0) {
                        s = wid*sample(w, k[i], replace = FALSE) + gr.df$start[i]
                    } else{
                        warning("Warning: trying to sample range of length > width of supplied GRanges element. Returning NULL for this element.")
                        return(NULL)
                    }
                } else{
                    s = seq(gr.df$start[i], gr.df$end[i], wid)
                }
            } else{
                s = (gr.df$end[i]-gr.df$start[i]-gr.df$length[i])*stats::runif(k[i])+gr.df$start[i]
            }

            return(data.frame(chr = gr.df$chr[i], start=s, end =s+wid-1, strand = as.character(strand(gr)[i]), query.id = i))
        })

        ## add check that not all widths are zero
        if (all(sapply(tmp, is.null))){
            stop("Error: Could not sample any ranges. Check that width of input is greater than request width of output.")
        }

        return(gr.fix(seg2gr(do.call('rbind', tmp)), gr))
    }
}




#' @name si2gr
#' @title Create \code{GRanges} from \code{Seqinfo} or \code{BSgenome}
#' @description
#'
#' Creates a genomic ranges from seqinfo object
#' i.e. a pile of ranges spanning the genome
#'
#' @param si \code{Seqinfo} object or a \code{BSgenome} genome
#' @param strip.empty boolean Flag to output non-zero GRanges only (default = FALSE)
#' @return \code{GRanges} representing the range of the input genome
#' @examples
#'
#' si2gr(hg_seqlengths())
#'
#' @export
si2gr = function(si, strip.empty = FALSE)
{
    if (is(si, 'BSgenome')){
        si = Seqinfo(names(seqlengths(si)), seqlengths(si))
    }

    ## treat si as seqlengths if vector
    ## examples
    ## 'si2gr(seqlengths(si))'
    ## 
    ## gr=GRanges('2:2000-3000')
    ## si2gr(seqlengths(gr))
    ## 
    if (is(si, 'vector')){
        si = Seqinfo(seqlengths = si, seqnames = names(si))
    } else if (!is(si, 'Seqinfo')){
        si = seqinfo(si)
    }

    sl = seqlengths(si)
    sn = seqnames(si);
    sl[is.na(sl)] = 0;

    if (strip.empty){
        sn = sn[sl!=0];
        sl = sl[sl!=0];
    }

    sigr = GRanges(sn, IRanges(rep(1, length(sl)), width = sl), seqlengths = seqlengths(si), strand = rep('+', length(sl)))
    names(sigr) = sn;

    return(sigr)
}




#' @name grbind
#' @title Concatenate \code{GRanges}, robust to different \code{mcols}
#' @description
#'
#' Concatenates \code{GRanges} objects, taking the union of their features if they have non-overlapping features
#'
#' @param x GRanges input GRanges
#' @param ... additional input GRanges
#' @note Does not fill in the \code{Seqinfo} for the output \code{GRanges}
#' @return Concatenated \code{GRanges}
#' grbind(example_genes, example_dnase)
#' @export
grbind = function(x, ...)
{
    if (missing('x')){
        grs = list(...)
    } else if (class(x) != 'list'){
        grs <- c(list(x), list(...))
    } else{
        grs <- c(x, list(...))
    }

    force.rrbind = FALSE

    keep = sapply(grs, length)>0 & sapply(grs, function(x) inherits(x, 'GRanges'))
    grs = grs[keep]

    if (length(grs)==0){
        return(NULL)
    }

    if (length(grs)==1){
        return(grs[[1]])
    }

    vals = lapply(grs, function(x) values(x))

    ## DataFrame from IRanges package can hold XStringSets. Convert first

    isDataFrame <- sapply(vals, class) == 'DataFrame'
    if (any(isDataFrame))
    {
        tmp = tryCatch(lapply(vals[isDataFrame], gr2dt), error = function(e) NULL) ## sometimes works, sometimes doesn't
        if (is.null(tmp)){
            tmp = lapply(vals[isDataFrame, as.data.frame])
        }

        vals[isDataFrame] = tmp
    }


    ### FIXING seqlengths when not exactly matching
    sls = lapply(grs, seqlengths)
    names(sls) = NULL
    levs = unique(names(unlist(sls)))
    sl.new = structure(rep(0, length(levs)), names = levs)
    for (sl in sls){
        sl.new[names(sl)] = pmax(sl.new[names(sl)], sl, na.rm = TRUE)
    }

    bare.grs = lapply(grs, function(x) gr.fix(x[,c()], sl.new))

    ##out = tryCatch(do.call('c', bare.grs), error = function(e) NULL) ## this is annoyingly not working
    out <- dt2gr(rbindlist(lapply(bare.grs, gr2dt)), seqlengths = sl.new)

    ix <- (sapply(vals, ncol)==0)
    if (any(ix)){
        vals[ix] = lapply(which(ix), function(x) data.frame(col.4214124124124 = rep(NA, length(grs[[x]]))))
    }

    if (!force.rrbind){
        tmp = tryCatch(do.call('rrbind', vals), error = function(e) NULL)
    } else{
        tmp = NULL
    }

    values(out) = tmp

    if (any(ix)){
        out$col.4214124124124 = NULL
    }

    return(out)
}




#' @name grl.bind
#' @title Concatenate \code{GRangesList} objects.
#' @description
#'
#' Concatenates \code{GRangesList} objects taking the union of their \code{mcols} features if they have non-overlapping features
#'
#' @param ... GRangesList Any number of \code{GRangesList} to concatenate together
#' @return Concatenated \code{GRangesList} with NA filled in for \code{mcols} fields that are non-overlapping. Note that the
#' elements are re-named with sequential numbers
#' @examples
#'
#' ## Concatenate
#' grl.hiC2 <- grl.hiC[1:20]
#' mcols(grl.hiC2)$test = 1
#' grl.bind(grl.hiC2, grl.hiC[1:30])
#'
#' @export
#' @author Marcin Imielinski
#' @importFrom GenomicRanges mcols<- mcols split
grl.bind = function(...)
{
  ## TODO: make this work for when underlying grs do not have matching features
  ## currently will loose gr level features
  grls = list(...)

  ## check the input

  if(any(sapply(grls, function(x) !inherits(x, "GRangesList")))){
    stop("Error: All inputs must inherit GRangesList")
  }

  ## annoying acrobatics to reconcile gr and grl level features for heterogenous input gr / grls
  grls.ul = lapply(grls, grl.unlist)
  grls.ul.rb = do.call('grbind', grls.ul)

  if (length(grls.ul.rb)==0)
    {
      out = GRangesList()
      seqinfo(out) = seqinfo(grls[[1]])
      return(out)
    }

  sp = base::unlist(lapply(1:length(grls), function(x) rep(x, length(grls.ul[[x]]))))
  gix = base::split(grls.ul.rb$grl.ix, sp)
  gjx = base::split(1:length(grls.ul.rb), sp)
  grls.ul.rb$grl.iix = grls.ul.rb$grl.ix = NULL

  grls.vals = lapply(grls, function(x){
    if (ncol(mcols(x))>0){
      return(as.data.frame(mcols(x)))
    } else{
      return(data.frame(dummy241421 = rep(NA, length(x))))
    }
  })

  grls.new = mapply(function(x,y) GenomicRanges::split(grls.ul.rb[x],y), gjx, gix)

  ## do.call('c', grls.new) is not working for some reason (gives back list again, not GRangesList)
  ## have to do this instead, not ideal
  if (length(grls.new) > 1) {
    out = grls.new[[1]]
    for (i in 2:length(grls.new)){
      out = c(out, grls.new[[i]])
    }
  } else {
    out = grls.new[[1]]
  }

  out.val = do.call('rrbind', grls.vals)
  out.val$dummy241421 = NULL
  GenomicRanges::mcols(out) <- out.val

  return(out)
}




#' @name gr.chr
#' @title Prepend "chr" to \code{GRanges seqlevels}
#' @description
#'
#' Prepend "chr" to \code{GRanges seqlevels}
#'
#' @param gr \code{GRanges} object to append 'chr' to
#' @return Identical \code{GRanges}, but with 'chr' prepended to each seqlevel
#' @examples
#'
#' gr <- gr.chr(GRanges(c(1,"chrX"), IRanges(c(1,2), 1)))
#' seqnames(gr)
#'
#' @details implementation looks for strings that begin with 0-9, X,Y,M, or Un for
#' the unmapped chromosome regions/decoys.  
#' regex pattern: "^([0-9,X, Y, M,x,y,m, Un, {EBV}]+)"
#' replacement pattern: "chr\\1"
#' 
#' @importFrom GenomeInfoDb seqlevels seqlevels<-
#' @author Max Chao
#' @export
gr.chr  = function(gr){
  if (grepl('^chr', seqlevels(gr)[1]) == FALSE){
    seqlevels(gr) = gsub(seqlevels(gr), 
                         pattern = "^([0-9,X, Y, M,x,y,m, Un, {EBV}]+)", 
                         replacement = "chr\\1")
  }
  return(gr)
}


#' @name streduce
#' @title Reduce \code{GRanges} and \code{GRangesList} to miminal footprint
#' @description
#'
#' Reduce \code{GRanges} and \code{GRangesList} to miminal footprint
#'
#' Shortcut for \code{reduce(sort(gr.stripstrand(unlist(x))))}
#'
#' @param gr \code{GRanges} or \code{GRangesList}
#' @param pad integer Expand the input data before reducing. (default = 0)
#' @param sort boolean Flag to sort the output. (default = TRUE)
#' @return \code{GRanges} object with no strand information, representing a minimal footprint
#' @importFrom GenomicRanges reduce
#' @examples
#'
#' streduce(grl.hiC, pad=10)
#' streduce(example_genes, pad=1000)
#'
#' @export
streduce = function(gr, pad = 0, sort = TRUE)
{

    if (inherits(gr, 'GRangesList')){
        gr = unlist(gr)
    }

    if (any(is.na(seqlengths(gr)))){
        gr = gr.fix(gr)
    }

    out = suppressWarnings(sort(reduce(gr.stripstrand(gr + pad))))
    suppressWarnings(start(out) <-pmax(1, start(out)))
                                        #    out <- gr.tfix(out)
    end(out) = pmin(end(out), seqlengths(out)[as.character(seqnames(out))])

    return(out)
}




#' @name gr.string
#' @title Return UCSC style interval string corresponding to \code{GRanges} pile (ie chr:start-end)
#' @description
#'
#' Return UCSC style interval string corresponding to \code{GRanges} pile (ie chr:start-end)
#'
#' @param gr \code{GRanges} pile to get intervals from
#' @param add.chr boolean Flag to prepend seqnames with "chr" (default = FALSE)
#' @param mb boolean Flag to round to the nearest megabase (default = FALSE)
#' @param round integer If \code{mb} supplied, the number of digits to round to. (default = 3)
#' @param other.cols character vector Names of additional \code{mcols} fields to add to the string (seperated by ";")
#' @param pretty boolean Flag to output interval string in more readable format
#' @examples
#' gr.string(example_genes, other.cols = c("name", "name2"))
#' @return UCSC style interval string corresponding to \code{GRanges} pile
#' @export
#' @author Marcin Imielinski
gr.string = function(gr, add.chr = FALSE, mb = FALSE, round = 3, other.cols = c(), pretty = FALSE)
{
    if (length(gr)==0){
        return(as.character(NULL))
    }

    sn = as.character(seqnames(gr))

    if (add.chr){
        sn = paste('chr', sn, sep = '');
    }

    other.cols = intersect(names(values(gr)), other.cols)
    if (length(other.cols)>0){
        other.str = paste(' ', do.call('paste', c(lapply(other.cols, function(x) values(gr)[, x]), list(sep = ' '))))
    } else{
        other.str = ''
    }

    str = ifelse(as.logical(strand(gr)!='*'), as.character(strand(gr)), '')

    
    if (mb){
        return(paste(sn, ':', round(start(gr)/1e6, round), '-', round(end(gr)/1e6, round), str, other.str, sep = ''))
    }
    else{
        if (pretty){
            return(paste(sn, ':', stringr::str_trim(format(start(gr), big.mark = ',')), '-', stringr::str_trim(format(end(gr), big.mark = ',')), str, other.str, sep = ''))
        }
        else{
            return(paste(sn, ':', start(gr), '-', end(gr), str, other.str, sep = ''))
        }
    }
}




#' @name grl.reduce
#' @title grl.reduce
#' @description
#'
#' Quickly computes GRanges +/- padding inside a GRangesList
#' Can use with split / unlist
#'
#' @param grl \code{GRangesList} input
#' @param pad integer Padding to add to ranges inside GRangesList before reducing (default = 0)
#' @param clip boolean Flag to add to ranges inside GRangesList before reducing (default = FALSE)
#' @return \code{GRangesList} with GRanges of intervals "original GRanges +/- padding"
#' @export
#' @author Marcin Imielinski
grl.reduce = function(grl, pad = 0, clip = FALSE)
{
    ## function also in skitools
    label.runs = function(x){
        as.integer(ifelse(x, cumsum(diff(as.numeric(c(FALSE, x))) > 0), NA))
    }

    sl = data.table(lev = seqlevels(grl), len = seqlengths(grl))
    setkey(sl, lev)
    gr = grl.unlist(grl)

    ## clip to seqlengths
    if (clip){
        start(gr) = pmax(1, start(gr)-pad)
        end(gr) = pmin(sl[as.character(seqnames(gr)), len], end(gr)+clip, na.rm = TRUE)
    } else{
        gr = gr + pad
    }

    gr$group = gr$grl.ix

    grd = data.table(chr = as.character(seqnames(gr)), pos = c(start(gr), end(gr)), type = rep(as.numeric(c(1,-1)), each = length(gr)), group = rep(gr$group, 2))
    setkeyv(grd, c('group', 'chr', 'pos'))

    grd[, cs := cumsum(type), by = group]
    grd[, run := label.runs(cs!=0), by = group][, run := ifelse(is.na(run), c(NA, run[-length(run)]), run)][, run := paste(group, run)]
    out = grd[, .(chr = chr[1], start = pos[1], end = pos[length(pos)], group = group[1]), by = run]

    out.gr = dt2gr(out)  ## convert to GRanges
    out.grl = split(out.gr, factor(out.gr$group, 1:length(grl)))
    values(out.grl) = values(grl)

    if (!is.null(names(grl))){
        names(out.grl) = names(grl)
    }

    return(out.grl)
}




#' @name grl.string
#' @title Create string representation of \code{GRangesList}
#' @description
#'
#' Return UCSC style interval string corresponding to each \code{GRanges} in the \code{GRangesList}.
#' One line per per \code{GRangesList} item. \code{GRanges} elements themselves are separated by \code{sep}
#'
#' @param grl \code{GRangesList} to convert to string vector
#' @param mb boolean Flag to return as MB and round to "round" (default = FALSE)
#' @param sep Character to separate single \code{GRanges} ranges (default = ',')
#' @param ... Additional arguments to be passed to \code{gr.string}
#' @return Character vector where each element is a \code{GRanges} pile corresponding to a single \code{GRangesList} element
#' @examples
#' grl.string(grl.hiC, mb=TRUE)
#' @author Marcin Imielinski
#' @export
grl.string = function(grl, mb= FALSE, sep = ',', ...)
{

    if (class(grl) == "GRanges"){
        return(gr.string(grl, mb=mb, ...))
    }


    if (!inherits(grl, "GRangesList")){
        stop("Error: Input must be GRangesList (or GRanges, which is sent to gr.string)")
    }

    gr = grl.unlist(grl)
    if (!is.null(names(grl))){
        nm = names(grl)
    } else{
        nm = 1:length(grl)
    }

    grs = gr.string(gr, mb=mb, ...)
    out = sapply(split(grs, gr$grl.ix), paste, collapse = sep)
    names(out) = nm[as.numeric(names(out))]
    return(out)
}





#' @name gr.flipstrand
#' @title Flip strand on \code{GRanges}
#' @description
#'
#' Flip strand on \code{GRanges}
#'
#' @param gr \code{GRanges} pile with strands to be flipped
#' @return \code{GRanges} with flipped strands (+ to -, * to *, - to *)
#' @examples
#'
#' gr.flipstrand(GRanges(1, IRanges(c(10,10,10),20), strand=c("+","*","-")))
#'
#' @export
gr.flipstrand =function(gr){

    if (!is(gr, 'GRanges')){
        stop('Error: GRanges input only. Please check the documentation.')
    }

    if (length(gr)==0){
        return(gr)
    }

    which = cbind(1:length(gr), TRUE)[,2] == 1

    if (any(which)){
        strand(gr)[which] = c('*'='*', '+'='-', '-'='+')[as.character(strand(gr))][which]
    }

    return(gr)
}


#' @name gr.fix
#' @title "Fixes" \code{seqlengths} / \code{seqlevels}
#' @description
#'
#' If "genome" not specified will replace \code{NA} \code{seqlengths} in GRanges to reflect largest coordinate per \code{seqlevel}
#' and removes all \code{NA seqlevels} after this fix.
#'
#' If "genome" defined (i.e. as \code{Seqinfo} object, or a \code{BSgenome}, \code{GRanges}, \code{GRangesList} object with populated \code{seqlengths}),
#' then will replace \code{seqlengths} in \code{gr} with those for that genome
#'
#' @name gr.fix
#' @param gr \code{GRanges} object to fix
#' @param genome Genome to fix to: \code{Seqinfo}, \code{BSgenome}, \code{GRanges} (w/seqlengths), \code{GRangesList} (w/seqlengths) (default = NULL)
#' @param gname string Name of the genome (optional, just appends to \code{Seqinfo} of the output) (default = NULL)
#' @param drop boolean Remove ranges that are not present in the supplied genome (default = FALSE)
#' @param pruning.mode string Controls how to prune x. Four pruning modes are currently defined: "error", "coarse", "fine", and "tidy". See GenomeInfoDb documentation.
#' @return \code{GRanges} pile with the fixed \code{Seqinfo}
#' @importFrom GenomeInfoDb Seqinfo seqinfo keepSeqlevels seqlevels seqlengths seqlevels<- seqlengths<- genome<- seqnames
#' @export
gr.fix = function(gr, genome = NULL, gname = NULL, drop = FALSE, pruning.mode = "coarse")
{
    sn = V1 = NULL ## NOTE fix

    if (!is.null(genome))
    {
        ## assume seqlengths was provided (ie named vector of lengths)
        if (is.vector(genome)){
            genome = Seqinfo(names(genome), seqlengths = genome)
        } else if (!(is(genome, 'character') | inherits(genome, 'GRanges') | inherits(genome, 'BSgenome') | inherits(genome, 'GRangesList') | inherits(genome, 'Seqinfo')) & !is.vector(genome)){
            genome = seqinfo(genome)
        }

        if (!is.vector(genome)){
            if (!drop){
                levs = union(seqlevels(genome), as.character(seqnames(seqinfo(gr))))
                lens = structure(rep(NA, length(levs)), names = levs)
                lens[seqlevels(genome)] = seqlengths(genome);
                lens[seqlevels(gr)] = pmax(seqlengths(gr), lens[seqlevels(gr)], na.rm = TRUE)
            } else{
                lens = structure(seqlengths(genome), names = seqlevels(genome))
            }
        }

        seqlevels(gr, pruning.mode = "coarse") = names(lens)
        seqlengths(gr) = lens;
    } else{
        if (length(gr)>0){
            if (is(gr, 'GRangesList')){
                tmp.gr = unlist(gr)
            } else{
                tmp.gr = gr
            }

            tmp.sl = data.table(sn = as.character(seqnames(tmp.gr)), end = end(tmp.gr))[, max(end, na.rm = TRUE), by = sn][,  structure(V1, names = sn)][seqlevels(tmp.gr)]
            names(tmp.sl) = seqlevels(tmp.gr)
            seqlengths(tmp.gr)[!is.na(tmp.sl)] = suppressWarnings(pmax(tmp.sl[!is.na(tmp.sl)], seqlengths(tmp.gr)[!is.na(tmp.sl)], na.rm = TRUE))
            gr = tmp.gr
        }
        
        if (drop){
            gr = keepSeqlevels(gr, seqlevels(gr)[!is.na(seqlengths(gr))], pruning.mode = "coarse")
        }

    }

    ## hack to get rid of annoying "genome"
    if (!is.null(gname)){
        si = seqinfo(gr)
        genome(si) = gname
        gr@seqinfo = si
    }

    return(gr)
}




#' @name gr.flatten
#' @title Lay ranges end-to-end onto a derivate "chromosome"
#' @description
#'
#' Takes pile of \code{GRanges} and returns into a \code{data.frame} with \code{nrow = length(gr)} with each
#' representing the corresponding input range superimposed onto a single "flattened"
#' chromosome, with ranges laid end-to-end
#'
#' @param gr \code{GRanges} to flatten
#' @param gap integer Number of bases between ranges on the new chromosome (default = 0)
#' @return \code{data.frame} with start and end coordinates, and all of the original metadata
#' @name gr.flatten
#' @importFrom GenomicRanges mcols
#' @export
gr.flatten = function(gr, gap = 0)
{
    if (length(gr) == 0){
        return(data.frame())
    } else if (length(gr) == 1){
        return(data.frame(start = 1, end = width(gr)))
    } else{
        starts = as.numeric(cumsum(c(1, width(gr[1:(length(gr)-1)])+gap)))
        ends = as.numeric(starts+width(gr)-1)
        return(cbind(data.frame(start=starts, end=ends), as.data.frame(mcols(gr))))
    }
}




#' @name gr.stripstrand
#' @title gr.stripstrand
#' @description
#'
#' Sets strand to "*"
#'
#' @name gr.stripstrand
#' @param gr \code{GRanges} to remove strand information from
#' @return \code{GRanges} with strand set to \code{*}
#' @export
gr.stripstrand = function(gr)
{
    strand(gr) = "*"
    return(gr)
}



#' @name gr.flipstrand
#' @title Flip strand on \code{GRanges}
#' @description
#'
#' Flip strand on \code{GRanges}
#'
#' @param gr \code{GRanges} pile with strands to be flipped
#' @return \code{GRanges} with flipped strands (+ to -, * to *, - to *)
#' @examples
#' gr.flipstrand(GRanges(1, IRanges(c(10,10,10),20), strand=c("+","*","-")))
#' @export
gr.flipstrand= function(gr)
{

    if (!is(gr, 'GRanges')){
        stop('Error: GRanges input only')
    }

    if (length(gr)==0){
        return(gr)
    }

    which = cbind(1:length(gr), TRUE)[,2] == 1

    if (any(which)){
        strand(gr)[which] = c('*'='*', '+'='-', '-'='+')[as.character(strand(gr))][which]
    }

    return(gr)
}



#' @name gr.pairflip
#' @title Create pairs of ranges and their strand-inverse
#' @description
#'
#' From a \code{GRanges} returns a \code{GRangesList} with each item consisting
#' of the original \code{GRanges} and its strand flip
#'
#' @name gr.pairflip
#' @param gr \code{GRanges}
#' @return \code{GRangesList} with each element of length 2
#' @export
gr.pairflip = function(gr)
{
    strand(gr)[strand(gr) =='*'] = '+';
    return(split(c(gr, gr.flipstrand(gr)), rep(c(1:length(gr)), 2)))
}




#' @name gr.tile
#' @title Tile ranges across \code{GRanges}
#' @description
#'
#' Tiles interval (or whole genome) with segments of \code{<=} specified width.
#'
#' @param gr \code{GRanges}, \code{seqlengths} or \code{Seqinfo} range to tile. If has \code{GRanges} has overlaps, will reduce first.
#' @param width integer Width of each tile (default = 1e3)
#' @examples
#'
#' ## 10 tiles of width 10
#' gr1 <- gr.tile(GRanges(1, IRanges(1,100)), width=10)
#'
#' ## make them overlap each other by 5
#' gr1 + 5
#'
#' @return GRanges with tiled intervals
#' @export
gr.tile = function(gr, width = 1e3)
{
    numw = tile.id = query.id = NULL ## getting past NOTE
    if (is(gr, 'data.table')){
#        message('Executing dt2gr() on input \n');
        gr = dt2gr(gr);
    } else if (!is(gr, 'GRanges')){
 #       message('Executing si2gr() on input \n');
        gr = si2gr(gr);
    }

    ix = which(width(gr) > 0)
    gr = gr[ix]

    if (length(gr)==0){
        return(gr[c()][, c()])
    }

    ws = as.numeric(ceiling((width(gr))/width))

    st = rep(start(gr), ws)
    en = rep(end(gr), ws)
    strand = rep(as.character(strand(gr)), ws)
    dt = data.table(query.id = rep(1:length(gr), ws), tile.id = 1:sum(ws))
    dt[, numw := 0:(length(tile.id)-1), by = query.id]
    start = pmin(st+width*dt$numw, en)
    end = pmin(st+width*(dt$numw+1)-1, en)

    out = GRanges(rep(as.character(seqnames(gr)), ws), IRanges(start, end), strand = strand, query.id = ix[dt$query.id], tile.id = 1:length(start), seqinfo = seqinfo(gr))

    return(out)
}




#' @name gr.tile.map
#' @title gr.tile.map
#' @description
#'
#' Given two tilings of the genome (e.g. at different resolution)
#' query and subject outputs a length(query) list whose items are integer vectors of indices in subject
#' overlapping that overlap that query (strand non-specific)
#'
#' @param query Query \code{GRanges} pile, perhaps created from some tile (e.g. \code{gr.tile}), and assumed to have no gaps
#' @param subject Subject \code{GRanges} pile, perhaps created from some tile (e.g. \code{gr.tile}), and assumed to have no gaps
#' @param verbose Increase the verbosity of the output (default = FALSE)
#' @return \code{list} of length = \code{length(query)}, where each element \code{i} is a vector of indicies in \code{subject} that overlaps element \code{i} of \code{query}
#' @note Assumes that input query and subject have no gaps (including at end) or overlaps, i.e. ignores end()
#' coordinates and only uses "starts"
#' @author Marcin Imielinski
#' @export
gr.tile.map = function(query, subject, verbose = FALSE)
{

    ## munlist (also in skitools, but added here to uncouple from this dependency)
    ##
    ## unlists a list of vectors, matrices, data frames into a n x k matrix
    ## whose first column specifies the list item index of the entry
    ## and second column specifies the sublist item index of the entry
    ## and the remaining columns specifies the value(s) of the vector
    ## or matrices.
    ##
    ## force.cbind = T will force concatenation via 'cbind'
    ## force.rbind = T will force concatenation via 'rxsbind'

    munlist = function(x, force.rbind = FALSE, force.cbind = FALSE, force.list = FALSE){
        if (!any(c(force.list, force.cbind, force.rbind))){
            if (any(sapply(x, function(y) is.null(dim(y))))){
                force.list = TRUE
            }
            if (length(unique(sapply(x, function(y) dim(y)[2]))) == 1){
                force.rbind = TRUE
            }
            if ((length(unique(sapply(x, function(y) dim(y)[1]))) == 1)){
                force.cbind = TRUE
            }
        } else{
            force.list = TRUE
        }

        if (force.list){
            return(cbind(ix = unlist(lapply(1:length(x), function(y) rep(y, length(x[[y]])))),
                iix = unlist(lapply(1:length(x), function(y) if (length(x[[y]])>0) 1:length(x[[y]]) else NULL)),
                unlist(x)))
        } else if (force.rbind){
            return(cbind(ix = unlist(lapply(1:length(x), function(y) rep(y, nrow(x[[y]])))),
                iix = unlist(lapply(1:length(x), function(y) if (nrow(x[[y]])>0) 1:nrow(x[[y]]) else NULL)),
                do.call('rbind', x)))
        } else if (force.cbind){
            return(t(rbind(ix = unlist(lapply(1:length(x), function(y) rep(y, ncol(x[[y]])))),
                iix = unlist(lapply(1:length(x), function(y) if (ncol(x[[y]])>0) 1:ncol(x[[y]]) else NULL)),
                do.call('cbind', x))))
        }
    }

    mc.cores =1 ## multicore not supported now

    if (length(GenomicRanges::gaps(query)) > 0){
        warning("Warning: Query GRanges has gaps. Unexpected behavior may follow")
    }

    if (length(GenomicRanges::gaps(subject)) > 0){
        warning("Warning: Subject GRanges has gaps. Unexpected behavior may follow")
    }

    ix.q = order(query)
    ix.s = order(subject)

    q.chr = as.character(seqnames(query))[ix.q]
    s.chr = as.character(seqnames(subject))[ix.s]

    ql = base::split(ix.q, q.chr)
    sl = base::split(ix.s, s.chr)

    tmp <- parallel::mcmapply(
        function(x,y){
            if (length(y)==0){
                return(NULL)
            }
            all.pos = c(start(query)[x], start(subject)[y])
            is.q = c(rep(T, length(x)), rep(F, length(y)))
            all.ix = c(x, y)
            ord.ix = order(all.pos)
            all.pos = all.pos[ord.ix]
            is.q = is.q[ord.ix]
            all.ix = all.ix[ord.ix]
            out = matrix(NA, nrow = length(all.pos), ncol = 2)
            last.x = last.y = NA
            for (i in 1:length(all.pos)){
                if (verbose){
                    if (i %% 100000 == 0){
                        cat('Iteration', i, 'of', length(all.pos), '\n')
                    }
                }
                if (is.q[i]){
                    out[i, ] = c(all.ix[i], last.y)

                    ## edge case where subject and query intervals share a start point, leading to two consecutive all.pos
                    if (i<length(all.pos)){
                        if (all.pos[i] == all.pos[i+1]){
                            out[i, ] = NA
                        }
                    }

                    last.x = all.ix[i]
                } else{
                    out[i, ] = c(last.x, all.ix[i])

                    ## edge case where subject and query intervals share a start point, leading to two consecutive all.pos
                    if (i<length(all.pos)){
                        if (all.pos[i] == all.pos[i+1]){
                            out[i, ] = NA
                        }
                    }
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




#' @name gr.val
#' @title Annotate \code{GRanges} with values from another \code{GRanges}
#' @description
#'
#' Annotates \code{GRanges} in \code{query} with aggregated values of \code{GRanges} in \code{target} in field \code{val}.
#' If \code{val} is numeric: given \code{target} with value column \code{target} representing ranged data
#' (i.e. segment intensities), thn computes the value
#' in each \code{query} \code{GRanges} as the weighted mean of its intersection with target
#' (ie the target values weighted by the width of the intersections).
#'
#' Applications include (among others):
#' \itemize{
#' \item Querying the average value of target across a given query interval (e.g. exon to gene pileup)
#' \item recasting a high res tiling in terms of low res intervals.
#' }
#' Usually query intervals are bigger than the target intervals.
#'
#' @note \code{query} and \code{target} can be \code{GRangesList} object, in which case val will refer to \code{GRangesList} level values fields
#' @param query \code{GRanges} of query ranges whose \code{val} column we will populate with aggregated values of \code{target}
#' @param target \code{GRanges} of target ranges that already have "val" column populated
#' @param val If a character field: then aggregation will paste together the (unique), overlapping values, collapsing by comma. (default = NULL)
#' @param mean boolean If \code{FALSE} then will return sum instead of mean, only applies if target \code{val} column is numeric. (default = TRUE)
#' @param weighted Calculate a weighted mean. If \code{FALSE}, calculates unweighted mean. (default = 'mean')
#' @param na.rm boolean Remove NA values when calulating means. only applies if val column of target is numeric (default = FALSE)
#' @param by scalar character, specifies additional "by" column of query AND target that will be used to match up query and target pairs (i.e. in addition to pure GRanges overlap). (default = NULL)
#' @param by.prefix Choose a set of \code{val} fields by a shared prefix. (default = 'val')
#' @param merge boolean If merge = FALSE then will cross every range in query with every level of "by" in target (and create data matrix), otherwise will assume query has "by" and merge only ranges that have matching "by" values in both query and target (default = FALSE)
#' @param FUN Optional different function to call than mean. Takes two arguments (value, na.rm = TRUE) if weighted = FALSE, and three (value, width, na.rm = TRUE) if weighted = TRUE. (default = NULL)
#' @param default.val If no hit in \code{target} found in \code{query}, fill output \code{val} field with this value. (default = NA)
#' @param max.slice integer Maximum number of query ranges to consider in one memory chunk. (default = Inf)
#' @param mc.cores integer Number of cores to use when running in chunked mode (default = 1)
#' @param sep string Specifies character to use as separator when aggregating character "vals" from target, only applies if target is character (default = ', ')
#' @param verbose boolean Increase the verbosity of the output (default = FALSE)
#' @param ... Additional arguments to be sent to \code{\link{gr.findoverlaps}}.
#' @return \code{query} with the \code{val} field populated
#' @author Marcin Imielinski
#' @export
gr.val = function(query, target, val = NULL, mean = TRUE, weighted = mean, na.rm = FALSE, by = NULL, by.prefix = val, merge = FALSE,
    FUN = NULL, default.val = NA, max.slice = Inf, mc.cores = 1,  sep = ', ', verbose = FALSE, ...)
{
    query.id = subject.id = NULL ## fix NOTE

    if (is.null(val)){
        val = 'value'
    }

    if (!(all(ix <- val %in% names(values(target))))){
        values(target)[, val[!ix]] = 1
    }

    if (length(query) > max.slice){
        ix.l = split(1:length(query), ceiling(as.numeric((1:length(query)/max.slice))))
        return(do.call('grbind', parallel::mclapply(ix.l, function(ix) {
            if (verbose){
                cat(sprintf('Processing intervals %s to %s of %s\n', min(ix), max(ix), length(query)))
            }
            gr.val(query[ix, ], target = target, val= val, mean = mean, weighted = weighted, na.rm = na.rm, verbose = verbose, by = by, FUN = FUN, merge = merge, ...)
        }, mc.cores = mc.cores)))
    }

    if (inherits(target, 'GRangesList')){
        target.was.grl = TRUE;
        target.grl.id = as.data.frame(target)$element
        val.vec = lapply(val, function(x) values(target)[, x])
        target = unlist(target)
        val.vec = lapply(val.vec, function(X) val.vec[target.grl.id])
    } else{
        if (!is.null(val)){
            val.vec = lapply(val, function(x) values(target)[, x])
        } else{
            val.vec = list(rep(1, length(target)))
        }

        target.grl.id = 1:length(target);
        target.was.grl = FALSE;
    }

    if (inherits(query, 'GRangesList')){
        query.was.grl = TRUE;
        query.grl.id = rep(1:length(query), S4Vectors::elementNROWS(query))
        query = unlist(query)
    } else{
        query.grl.id = 1:length(query);
        query.was.grl = FALSE;
    }

    if (!is.null(FUN)){
        args = names(formals(FUN))[1:3]

        if (!is.null(args)){
            if (any(is.na(args))){
                args[is.na(args)] = ''
            }

            if (weighted){
                if (any(!(c('x', 'w', 'na.rm') %in% args))){
                    warning('Warning: FUN input must be function with three arguments: "x" = value, "w" = interval width, "na.rm" = na.rm flag')
                }
            } else{
                if (any(!(c('x', 'na.rm') %in% args))){
                    warning('Warning: FUN input must be function with two arguments: "x" = value, "na.rm" = na.rm flag')
                }
            }
        }
    }

    if (!merge){
        hits = gr.findoverlaps(query, target, scol = by, verbose = verbose, return.type = 'data.table', ...)
    } else{
        hits = gr.findoverlaps(query, target, by = by, verbose = verbose, return.type = 'data.table', ...)
    }

    if (verbose){
        cat(sprintf('aggregating hits\n'))
    }

    vals = val
    val.vecs = val.vec

    for (vix in 1:length(vals)){
        val = vals[[vix]]
        val.vec = val.vecs[[vix]]
        is.char = is.character(values(target)[, val]);
        if (is.null(by) | merge == TRUE)
        {
            if (nrow(hits)>0)
            {
                values(query)[, val] = NA;
                hits[, width := as.numeric(end - start)+1]
                data.table::setkey(hits,  query.id)
                if (is.char)
                {
                    values(query)[, val] = '';
                    hits$id = 1:nrow(hits);
                    if (!is.null(FUN)){
                        tmp = hits[, list(val = do.call(FUN, list(val.vec[subject.id], width, na.rm = na.rm))), by = query.id]
                    } else{
                        tmp = hits[, list(val = paste(setdiff(val.vec[subject.id], NA), collapse = sep)), by = query.id]
                    }
                    if (!is.na(default.val)){
                        tmp[is.na(tmp)] = default.val
                    }
                    values(query)[tmp[,query.id], val] = tmp[,val]
                } else{

                    val.vec = as.numeric(val.vec)
                    if (weighted) {
                        if (!is.null(FUN)){
                            tmp = hits[, list(val = do.call(FUN, list(val.vec[subject.id], width, na.rm = na.rm))), by = query.id]
                        } else if (mean){
                            tmp = hits[, list(val = sum(width * val.vec[subject.id], na.rm = na.rm)/sum(width)), by = query.id]
                        } else{
                            tmp = hits[, list(val = sum(width * val.vec[subject.id], na.rm = na.rm)), by = query.id]
                        }
                    } else{
                        if (!is.null(FUN)){
                            tmp = hits[, list(val = do.call(FUN, list(val.vec[subject.id], na.rm = na.rm))), by = query.id]
                        } else if (mean){
                            tmp = hits[, list(val =  mean(val.vec[subject.id], na.rm = na.rm)), by = query.id]
                        } else{
                            tmp = hits[, list(val = sum(val.vec[subject.id], na.rm = na.rm)), by = query.id]
                        }
                    }

                    if (!is.na(default.val)){
                        tmp[is.na(tmp)] = default.val
                    }

                    values(query)[tmp[,query.id], val] = tmp[,val]
                }
            } else{
                values(query)[, val] = NA
            }
        } else {
            ## by is not null
            if (!is.null(by.prefix)){

                if (is.na(by.prefix)){
                    by.prefix = NULL
                } else if (nchar(by.prefix)==0){
                    by.prefix = NULL
                }
            }

            if (nrow(hits)>0) {
                hits[, width := as.numeric(end - start)+1]
                if (is.char) {
                    hits$id = 1:nrow(hits);
                    tmp = hits[, list(val = paste(setdiff(val.vec[subject.id], NA), collapse = sep)), keyby = list(query.id, bykey = eval(parse(text=by)))]

                    tmp2 = data.table::dcast.data.table(tmp, query.id ~ bykey, value.var = 'val')
                    data.table::setkey(tmp2, query.id)
                    new.df = as.data.frame(tmp2[list(1:length(query)), ])[ ,-1]

                    if (!is.na(default.val)){
                        new.df[is.na(new.df)] = default.val
                    }

                    if (!is.null(by.prefix)){
                        colnames(new.df) =  paste(by.prefix, names(tmp2)[-1], sep = '.')
                    } else{
                        colnames(new.df) =  names(tmp2)[-1]
                    }
                    new.names = c(colnames(values(query)), colnames(new.df))
                    values(query) = cbind(values(query), new.df)
                    colnames(values(query)) = new.names
                } else{
                    val.vec = as.numeric(val.vec)
                    if (weighted)
                    {
                        if (!is.null(FUN))
                        {
                            tmp = hits[, list(val = do.call(FUN, list(val.vec[subject.id], width, na.rm = na.rm))), keyby = list(query.id, bykey = eval(parse(text=by)))]
                        } else if (mean){
                            tmp = hits[, list(val = sum(width * val.vec[subject.id], na.rm = na.rm)/sum(width)), keyby = list(query.id, bykey = eval(parse(text=by)))]
                        } else{
                            tmp = hits[, list(val = sum(width * val.vec[subject.id], na.rm = na.rm)), keyby = list(query.id, bykey = eval(parse(text=by)))]
                        }
                    } else{
                        if (!is.null(FUN)){
                            tmp = hits[, list(val = do.call(FUN, list(val.vec[subject.id], na.rm = na.rm))), keyby = list(query.id, bykey = eval(parse(text=by)))]
                        } else if (mean){
                            tmp = hits[, list(val =  mean(val.vec[subject.id], na.rm = na.rm)), keyby = list(query.id, bykey = eval(parse(text=by)))]
                        } else{
                            tmp = hits[, list(val = sum(val.vec[subject.id], na.rm = na.rm)), keyby = list(query.id, bykey = eval(parse(text=by)))]
                        }
                    }

                    tmp2 = data.table::dcast.data.table(tmp, query.id ~ bykey, value.var = 'val')
                    data.table::setkey(tmp2, query.id)
                    new.df = as.data.frame(tmp2[list(1:length(query)), ])[ ,-1, drop = FALSE]

                    if (!is.na(default.val)){
                        new.df[is.na(new.df)] = default.val
                    }

                    if (!is.null(by.prefix)){
                        colnames(new.df) =  paste(by.prefix, names(tmp2)[-1], sep = '.')
                    } else{
                        colnames(new.df) =  names(tmp2)[-1]
                    }

                    new.names = c(colnames(values(query)), colnames(new.df))
                    values(query) = cbind(values(query), new.df)
                    colnames(values(query)) = new.names
                }
            } else{
                if(is.null(by.prefix)){
                    for (val in levels(factor(values(target)[, by]))){
                        values(query)[, val] = NA
                    }
                } else {
                    for (val in levels(factor(values(target)[, by]))){
                        values(query)[, paste(by.prefix, val, sep='.')] = NA
                    }
                }
            }
        }
    }

    if (query.was.grl){
        query = split(query, query.grl.id)
    }

    return(query)
}




#' @name gr.duplicated
#' @title Allows to restrict duplicates using "by" columns and allows in exact matching
#' @description
#'
#' Allows to restrict duplicates using "by" columns and allows in exact matching
#'
#' @param query GRanges to query
#' @param by string Column 'by' used to restrict duplicates. See the 'by' argument for function gr.match()
#' @param type string 'type' used. See the 'type' argument for function gr.match()
#' @return boolean vector of match status
#' @export
gr.duplicated = function(query, by = NULL, type = 'any')
{
    return(duplicated(gr.match(query, query, by = by, type = type)))
}




#' @name gr.dice
#' @title  Dice up \code{GRanges} into \code{width = 1} \code{GRanges} spanning the input (warning can produce a very large object)
#' @description
#'
#' Dice up \code{GRanges} into \code{width = 1} \code{GRanges} spanning the input (warning can produce a very large object)
#'
#' @param gr \code{GRanges} object to dice
#' @importFrom S4Vectors Rle
#' @importFrom GenomeInfoDb seqlengths
#' @importFrom GenomicRanges seqnames width strand values<- values strand<- distance
#' @return \code{GRangesList} where kth element is a diced pile of \code{GRanges} from kth input \code{GRanges}
#' @examples
#' gr.dice(GRanges(c(1,4), IRanges(c(10,10),20)))
#' @author Marcin Imielinski
#' @export
gr.dice = function(gr)
{
    sn = Rle(as.character(seqnames(gr)), width(gr))
    tmp.st = start(gr)
    tmp.en = end(gr)
    st = unlist(lapply(1:length(gr), function(x) tmp.st[x]:tmp.en[x]))
    str = Rle(as.character(strand(gr)), width(gr))

    ## reverse for negative strand ranges
    if (any(ix <- as.logical(strand(gr) == '-'))) {
        w = width(gr)
        q.id = unlist(lapply(1:length(gr), function(x) rep(x, w[x])))
        q.l = split(1:length(st), q.id)

        for (j in q.l[ix]){
            st[j] = rev(st[j])
        }
    }

    ix = Rle(1:length(gr), width(gr))
    out = GRanges(sn, IRanges(st, st), strand = str, seqlengths = seqlengths(gr))
    values(out) = values(gr)[ix, ]

    out = GenomicRanges::split(out, ix)

    return(out)
}




#' @name gr.dist
#' @title Pairwise distance between two \code{GRanges}
#' @description
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
#'
#' @param gr1 First \code{GRanges}
#' @param gr2 Second \code{GRanges}
#' @param ignore.strand Don't required elements be on same strand to avoid \code{NA [FALSE]}
#' @param ... Additional arguments to be supplied to \code{GenomicRanges::distance}
#' @return \code{N} by \code{M} matrix with the pairwise distances, with \code{gr1} on rows and \code{gr2} on cols
#' @author Marcin Imielinski
#' @export
gr.dist = function(gr1, gr2 = NULL, ignore.strand = FALSE, ...)
{
    if (is.null(gr2)){
        gr2 = gr1
    }

    if (ignore.strand){
        strand(gr1) = '*'
        strand(gr2) = '*'
    }

    ix1 = rep(1:length(gr1), length(gr2))
    ix2 = rep(1:length(gr2), each = length(gr1))

    out = matrix(suppressWarnings(distance(gr1[ix1], gr2[ix2], ...)), nrow = length(gr1), ncol = length(gr2))

    return(out)
}


#' @name rle.query
#' @title Queries an \code{\link{RleList}} representing genomic data
#' @description
#'
#' Queries an \code{\link{RleList}} representing genomic data
#'
#' (i.e. a list whose names represent seqnames ie chromosomes, and lengths represent seqlengths)
#' via \code{GRanges} object
#'
#' @param subject.rle \code{Rle}
#' @param query.gr integer GRangeslist or GRanges
#' @param chunksize integer Number of \code{query.gr} ranges to consider in one memory chunk. (default = 1e9)
#' @param mc.cores Number of cores to apply when doing chunked operation (default = 1)
#' @param verbose boolean Set the verbosity of the output (default = FALSE)
#' @return Rle representing the (concatenated) vector of data (reversing order in case of negative strand input)
#' @note Throws warning if seqlengths(gr) do not correspond to the lengths of the \code{RleList} components
#' @export
rle.query = function(subject.rle, query.gr, chunksize = 1e9, mc.cores = 1, verbose = FALSE)
{
    was.grl = FALSE

    if (is(query.gr, 'GRangesList'))
    {
        was.grl = TRUE
        query.gr = grl.unlist(query.gr)
    }

    chunksize = pmin(1e9, chunksize)

    if ((sum(as.numeric(width(query.gr)))) > chunksize) ## otherwise integer overflow
    {
        tmp = rle(ceiling(cumsum(as.numeric(width(query.gr)))/chunksize))
        chunks = cbind(cumsum(c(1, tmp$lengths[-length(tmp$lengths)])), cumsum(c(tmp$lengths)))

        if (verbose){
            cat(sprintf('chunking up into %s chunks \n', nrow(chunks)))
        }

        out = do.call('c', parallel::mclapply(1:nrow(chunks), function(x) rle.query(subject.rle, query.gr[chunks[x,1]:chunks[x,2]]), mc.cores = mc.cores))
    } else {
        out = Rle(NA, sum(as.numeric(width(query.gr))));

        if (length(query.gr)>1){
            st.ix = cumsum(c(1, width(query.gr)[1:(length(query.gr)-1)]))
        } else{
            st.ix = 1
        }
        out.ix = IRanges(st.ix, st.ix + width(query.gr)-1) ## ranges in out corresponding to query

        for (chr in intersect(names(subject.rle), unique(as.character(seqnames(query.gr)))))
        {
            ix = which(as.character(seqnames(query.gr)) == chr)
            rix = ranges(query.gr)[ix]
            m = max(end(rix))

            ## pad subject rle if not long enough
            if (length(subject.rle[[chr]]) < m)
            {
                subject.rle[[chr]] = c(subject.rle[[chr]], Rle(NA, m - length(subject.rle[[chr]])))
            }
            ## out[unlist(as.integer(out.ix[ix]))] = subject.rle[[chr]][rix] ### previously (before as.integer of IRanges ceased to be supported)
            out[out.ix[ix]] = subject.rle[[chr]][rix]
        }

        if('-' %in% as.character(strand(query.gr)))
        {
            id = NULL; V1 = NULL ## NOTE fix
            tmp = data.table(ix = 1:sum(width(out.ix)), id = rep(1:length(out.ix), width(out.ix)), strand = rep(as.character(strand(query.gr)), width(out.ix)), key = 'ix')
            out = out[tmp[, rev(ix), by = id][, V1]]
        }
    }

    if (was.grl){
        out = split(out, Rle(query.gr$grl.ix, width(query.gr)))
    }

    return(out)
}




#' @name grl.in
#' @title Check intersection of \code{GRangesList} with windows on genome
#' @description
#'
#' Check intersection of \code{GRangesList} with windows on genome
#'
#' Like %in% for grl but now will return a logical vector that is true at position if i
#' only if the ranges in grl[i] intersect <<all>>, <<some>>, <<only>>  windows in the subject
#'
#' e.g. can use to identify read pairs whose ends are contained inside two genes)
#'
#' @param grl \code{GRangesList} object to query for membership in \code{windows}
#' @param windows \code{GRanges} pile of windows
#' @param some boolean Will return \code{TRUE} for \code{GRangesList} elements that intersect at least on window range (default = FALSE)
#' @param only boolean Will return \code{TRUE} for \code{GRangesList} elements only if there are no elements of query that fail to intersect with windows (default = FALSE)
#' @param logical boolean Will return logical otherwise will return numeric vector of number of windows overlapping each grl (default = TRUE)
#' @param exact boolean Will return exact intersection (default = FALSE)
#' @param ignore.strand boolean Will return exact intersection (default = TRUE)
#' @param maxgap integer maxgap is interpreted as the max gap allowed between two ranges to be considered as overlapping. See documentation for IRanges::findOverlaps() (default = -1L)
#' @param ... Additional parameters to be passed on to \code{GenomicRanges::findOverlaps}
#' @return boolean vector of match status
#' @export
grl.in = function(grl, windows, some = FALSE, only = FALSE, logical = TRUE, exact = FALSE, ignore.strand = TRUE, ...)
{

    grl.iid = grl.id = NULL ## for getting past NOTE

    if (length(grl)==0){

        if (logical){
            return(logical())
        } else{
            return(numeric())
        }
    }

    if (length(windows)==0){
        if (logical){
            return(rep(FALSE, length(grl)))
        } else{
            return(rep(0, length(grl)))
        }
    }


    numwin = length(windows);

    gr = grl.unlist(grl)
    if (logical){
        h = tryCatch(GenomicRanges::findOverlaps(gr, windows, ignore.strand = ignore.strand, ...), error = function(e) NULL)
        if (!is.null(h)){
            m = data.table(query.id = queryHits(h), subject.id = subjectHits(h))
        } else{
            m = gr2dt(gr.findoverlaps(gr, windows, ignore.strand = ignore.strand, ...))
        }
    } else{
        some = TRUE
        m = gr2dt(gr.findoverlaps(gr, windows, ignore.strand = ignore.strand, ...))
    }
    
    if (exact){
        m = m[start == start(gr)[query.id] & start == start(windows)[subject.id] & end == end(gr)[query.id] & end == end(windows)[subject.id], ]
    }

    out = rep(FALSE, length(grl))
    if (nrow(m)==0){
        return(out)
    }

    m$grl.id = gr$grl.ix[m$query.id]
    m$grl.iid = gr$grl.iix[m$query.id]

    if (some){
        tmp = as.data.frame(m[, length(unique(grl.iid)), by = grl.id])
    }
    else if (only){
        return(base::mapply(function(x, y) length(setdiff(x, y))==0,
                            split(1:length(gr), factor(gr$grl.ix, 1:length(grl))),
                            split(m$query.id, factor(m$grl.id, 1:length(grl)))))
    }
    else{
        tmp = stats::aggregate(subject.id ~ grl.id, data = m, FUN = function(x) numwin-length(setdiff(1:numwin, x)))
    }

    out = rep(FALSE, length(grl))
    out[tmp[,1]] = tmp[,2]

    if (logical){
        out = out!=0
    }

    return(out)
}






#' @name grl.unlist
#' @title Robust unlisting of \code{GRangesList} that keeps track of origin
#' @description
#'
#' Robust unlisting of \code{GRangesList} that keeps track of origin
#'
#' Does a "nice" unlist of a \code{GRangesList} object adding a field \code{grl.ix} denoting which element of the \code{GRangesList}
#' each \code{GRanges} corresponds to and a field \code{grl.iix} which saves the (local) index that that gr was in its corresponding \code{GRangesList} item
#'
#' In this way, \code{grl.unlist} is reversible, while \code{BiocGenerics::unlist} is not.
#'
#' @importFrom BiocGenerics unlist
#' @param grl \code{GRangeList} object to unlist
#' @param use.group \code{logical} force grl.unlist to use the "group" field for setting grl.ix. By default, grl.unlist would first try to use the "element", but if the input GRangesList contains a GRanges with a field that starts with the substring "element" then grl.unlist would fail.
#' @return \code{GRanges} with added metadata fields \code{grl.ix} and \code{grl.iix}.
#' @examples
#'
#' grl.unlist(grl.hiC)
#'
#' @export
grl.unlist = function(grl, use.group = FALSE)
{
    if (length(grl) == 0){
        return(GRanges())
    }

    if (is(grl, 'GRanges')){
        grl$grl.ix = 1
        grl$grl.iix = 1:length(grl)
        return(grl)
    }

    names(grl) = NULL
    as.df = as.data.frame(grl)

    el = NULL
    if (!use.group == TRUE){
        el = as.df$element
    }

    if (is.null(el)){
        el = as.df$group
    }

    out = BiocGenerics::unlist(grl)
    mcols(out)$grl.ix = el
    tmp = rle(el)

    nm = setdiff(names(values(grl)), c('grl.ix', 'grl.iix'))
    out$grl.iix = as.integer(unlist(sapply(tmp$lengths, function(x) 1:x)))

    tryCatch({
        values(out) = BiocGenerics::cbind(values(grl)[out$grl.ix, nm, drop = FALSE], values(out))
        },
        error = function(e){
            message('An error occured. This could occur if you have a field name in your input GRangesList that starts with "element". If this is indeed the case then please set the flag use.group to TRUE.')
            message('Here are more details about the error:')
            message(e)
        }
    )
    return(out)
}




#' @name grl.pivot
#' @title Pivot a \code{GRangesList}, inverting "x" and "y"
#' @description
#'
#' "Pivots" GRangesList object "x" by returning a new GRangesList "y" whose
#' kth item is GRanges object of ranges x[[i]][k] for all i in 1:length(x)
#'
#' e.g. If the length of a GRangesList `grl` is 50, `length(grl)=50  and length of each \code{GRanges} element inside is 2,
#' then the function \code{grl.pivot} will produce a length 3 \code{GRangesList} with 50 elements per \code{GRanges}
#'
#' Note: Assumes all GRanges in "x" are of equal length
#'
#' @param x \code{GRangesList} object to pivot
#' @importFrom GenomicRanges GRanges GRangesList split
#' @return GRangesList with inverted 'x' and 'y'
#' @examples
#' grl.pivot(grl.hiC)
#' @export
grl.pivot = function(x)
{
    if (length(x) == 0){
        return(GRangesList(GRanges(seqlengths = seqlengths(x)), GRanges(seqlengths = seqlengths(x))))
    }
    return(GenomicRanges::split(BiocGenerics::unlist(x), rep(1:length(x[[1]]), length(x))))
}




#' @name rrbind
#' @title Improved \code{rbind} for intersecting/union columns of \code{data.frames} or \code{data.tables}
#' @description
#'
#' Improved \code{rbind} for intersecting/union columns of \code{data.frames} or \code{data.tables}
#'
#' Like \code{rbind}, but takes the intersecting columns of the data.
#'
#' @param ... Any number of \code{data.frame} or \code{data.table} objects
#' @param union Take union of columns (and put NA's for columns of df1 not in df2 and vice versa). (default = TRUE)
#' @param as.data.table Return the binded data as a \code{data.table}. (default = FALSE)
#' @return \code{data.frame} or \code{data.table} of the \code{rbind} operation
#' @export
#' @importFrom data.table data.table rbindlist
#' @author Marcin Imielinski
rrbind = function (..., union = TRUE, as.data.table = FALSE)
{
    dfs = list(...)
    dfs = dfs[!sapply(dfs, is.null)]
    dfs = dfs[sapply(dfs, ncol) > 0]

    if (any(mix <- sapply(dfs, class) == 'matrix')){
        dfs[mix] = lapply(dfs, as.data.frame)
    }

    names.list = lapply(dfs, names)
    cols = unique(unlist(names.list))
    unshared = lapply(names.list, function(x) setdiff(cols, x))
    ix = which(sapply(dfs, nrow) > 0)
   
    if (any(sapply(unshared, length) != 0)){
        expanded.dts <- lapply(ix, function(x){
            tmp = dfs[[x]]
            if (is.data.table(dfs[[x]])){
                tmp = as.data.frame(tmp)
            }
            tmp[, unshared[[x]]] = NA
            return(data.table::as.data.table(as.data.frame(tmp[,
                                                               cols, drop = FALSE])))
        })
    } else{
        expanded.dts <- lapply(dfs, function(x) as.data.table(as.data.frame(x)[, cols, drop = FALSE]))
    }

    rout = tryCatch(rbindlist(expanded.dts, fill = TRUE), error = function(e) NULL)

    if (is.null(rout)){
        rout = data.table::as.data.table(do.call("rbind", lapply(expanded.dts,
                                                                 as.data.frame)))
    }
    
    if (!as.data.table){
        rout = as.data.frame(rout)
    }
    
    if (!union) {
        shared = setdiff(cols, unique(unlist(unshared)))
        rout = rout[, shared]
    }

    return(rout)
}




#' @name gr.sub
#' @title Apply \code{gsub} to seqlevels of a \code{GRanges}
#' @description
#'
#' Apply gsub to seqlevels of gr, by default removing 'chr', and "0.1" suffixes, and replacing "MT" with "M"
#'
#' @param gr \code{GRanges} to switch out seqlevels for
#' @param a vector of regular expressions of things to substitute out (default = c("^chr", "MT"))
#' @param b vector of values to substitute in (default = c("", "M"))
#' @return GRanges with substitutions
#' @export
gr.sub = function (gr, a = c("^chr", "MT"), b = c("", "M"))
{
    subs = cbind(a, b)
    tmp.gr = tryCatch(
    {
        tmp.gr = gr
        for (i in 1:nrow(subs)){
            seqlevels(tmp.gr) = unique(gsub(subs[i,1], subs[i,2], seqlevels(tmp.gr)))
        }
        tmp.gr
        }, error = function(e) NULL)

    if (is.null(tmp.gr)){
        warning('Warning: gr.sub had to convert GRanges to data.table before replacing seqlevels: check input seqlevels e.g. for mixed chr and non-chr seqlevels')

        is.list = FALSE
        gr.len = length(gr)
        if (is(gr, 'GRangesList')){
            is.list = TRUE
            if (!is.null(names(gr))){
                nm = structure(names(gr), names = as.character(1:length(gr)))
            } else{
                nm = NULL
            }
            gr = grl.unlist(gr)
        }

        tmp = gr2dt(gr)
        tmp$width = NULL
        sl = seqlevels(gr)
        for (i in 1:nrow(subs)){
            sl = gsub(subs[i,1], subs[i,2], sl)
            tmp[, seqnames := gsub(subs[i,1], subs[i,2], seqnames)]
        }
        sln = data.table(slev = sl, slen = seqlengths(gr))[, .(slen = max(slen, na.rm = TRUE)), by = slev][, structure(slen, names = slev)]
        tmp.gr = dt2gr(tmp, seqlengths = sln)

        if (is.list){
            newnm = setdiff(names(values(tmp.gr)), c('grl.ix', 'grl.iix'))
            tmp.gr = split(tmp.gr[, newnm], factor(tmp.gr$grl.ix, as.character(1:gr.len)))
            if (is.null(nm)){
                names(tmp.gr) = NULL
            } else{
                names(tmp.gr) = nm[names(tmp.gr)]
            }
        }
    }
    return(tmp.gr)
}




#' @name seg2gr
#' @title Convert GRanges-like data.frame into GRanges
#' @description
#'
#' Input data.frame of segments "segs" and converts into GRanges object porting over additional value columns
#'
#' "segs" data.frame/data.table can obey any number of conventions to specify chrom, start, and end of ranges
#' (e.g. $pos1, $pos2, $Start_position, $End_position)
#'
#' Please see documentation for function 'standardize_segs()' for more details.
#'
#' @param segs data.frame (or data.table) of segments with fields denoting chromosome, start, end, and other metadata. (See function 'standardize_segs()' for 'seg' data.frame/data.table input formats)
#' @param seqlengths seqlengths of output GRanges object (default = NULL)
#' @param seqinfo seqinfo of output GRanges object (default = Seqinfo())
#' @return GRanges from converted "segs" data.frame/data.table
#' @export
seg2gr = function(segs, seqlengths = NULL, seqinfo = Seqinfo())
{
    if (is(segs, 'data.table')){
        segs = as.data.frame(segs)
    }

    if (!is(segs, 'data.frame')){
        segs = as.data.frame(segs)
    }

    if (!inherits(seqinfo, 'Seqinfo')){
        seqinfo = seqinfo(seqinfo)
    }

    if (is.null(seqlengths)){
        seqlengths = seqlengths(seqinfo)
    }

    if (!inherits(segs, 'GRanges')){
        if (nrow(segs)==0){
            return(GRanges(seqlengths = seqlengths(seqinfo)))
        }
    } else if (length(segs) == 0){
        return(GRanges(seqlengths = seqlengths(seqinfo)))
    }

    segs = standardize_segs(segs)

    GR.NONO.FIELDS = c('seqnames', 'ranges', 'strand', 'seqlevels', 'seqlengths', 'isCircular', 'start', 'end', 'width', 'element', 'pos1', 'pos2', 'chr');

    if (is.null(segs$strand)){
        segs$strand = "*"
    }

    if (any(ix <- !(segs$strand %in% c('+', '-', '*')))){
        segs$strand[ix] = "*"
    }

    segs$chr = as.character(segs$chr)
    sl = data.table::as.data.table(segs)[, max(pos2), keyby = chr ]
    nms = setdiff(union(sl$chr, names(seqlengths)), NA)
    seqlengths = structure(pmax(sl[nms, V1], seqlengths[nms], na.rm = TRUE), names = nms)

    if (length(seqlengths)>0){
        if (length(wtf  <- setdiff(segs$chr, names(seqlengths)))){
            warning('Warning: Some seqnames in seg object were not included in provided seqlengths: ', paste(wtf, collapse = ','))
            seqlengths[as.character(wtf)] = NA
        }
        segs$pos1 = as.numeric(segs$pos1)
        segs$pos2 = as.numeric(segs$pos2)

        out = GRanges(seqnames = segs$chr, ranges = IRanges(segs$pos1, segs$pos2),strand = segs$strand, seqlengths = seqlengths)

    } else{
        out = GRanges(seqnames = as.character(segs$chr), ranges = IRanges(segs$pos1, segs$pos2), strand = segs$strand)
    }

    if (length(seqinfo)>0){
        out = gr.fix(out, seqinfo)
    } else if (is.null(seqlengths)){
        out = gr.fix(out)
    }

    values(out) = segs[, setdiff(names(segs), GR.NONO.FIELDS), drop = FALSE]

    return(out)
}




#' @name standardize_segs
#' @title Takes and returns segs data frame standardized to a single format
#' @description
#'
#' Takes and returns segs data frame standardized to a single format (i.e. $chr, $pos1, $pos2)
#'
#' If chr = TRUE will ensure "chr" prefix is added to chromossome(if does not exist)
#'
#' "segs" data.frame can obey any number of conventions to specify chrom, start, and end of ranges
#' (e.g. $pos1, $pos2, $Start_position, $End_position)
#'
#' Conventions:
#' \itemize{
#' \item \code{ID} - 'id', 'patient', 'sample'
#' \item \code{chr} - 'seqnames', 'chrom', 'chromosome', 'rname', 'space', 'contig'
#' \item \code{pos1} - 'start', 'loc.start', 'start.bp', 'start_position', 'begin', 'pos', 'pos1', 'left', 's1'
#' \item \code{pos2} - 'end', 'loc.end', 'stop', 'end.bp', 'end_posiiton', 'pos2', 'right', 'e1'
#' \item s\code{strand} - 'strand', 'str'
#' }
#'
#' @import GenomicRanges
#' @param seg data.frame or data.table or GRanges of segments with fields denoting chromosome, start, end, and other metadata.
#' @param chr boolean Flag to force add chromosomes (default = FALSE)
#' @return data.frame or data.table with standardized segments
#' @export
standardize_segs = function(seg, chr = FALSE)
{

    if (is(seg, 'matrix')){
        seg = as.data.frame(seg, stringsAsFactors = FALSE)
    }

    if (is(seg, 'GRanges')){
        seg = as.data.frame(gr2dt(seg))
    }

    if (is(seg, 'data.table')){
        seg = as.data.frame(seg)
    }

    if (!is(seg, 'data.frame') & !is(seg, 'data.table') & !is(seg, 'GRanges')){
        stop("Error: input 'segs' must be type data.frame, data.table or GRanges. Please see the documentation for information.")
    }

    val = NULL;

    field.aliases = list(
        ID = c('id', 'patient', 'Sample', 'ID', 'Patient', 'sample'),
        chr = c('seqnames', "seqname", 'chrom', 'Chromosome', 'chr', 'chromosome', 'Seqnames', 'Seqname', "contig", "rname", "space"),
        pos1 = c('start', 'loc.start', 'begin', 'Start', 'start', 'Start.bp', 'start_position', 'Start_position', 'Start_Position', 'start_Position', 'pos', 'pos1', 'left', 's1'),
        pos2 =  c('end', 'loc.end', 'End', 'end', "stop", 'End.bp', 'end_position', 'End_position', 'end_Position', 'End_Position', 'pos2', 'right', 'e1'),
        strand = c('strand', 'str', 'strand', 'Strand', 'Str')
    )

    
    if (is.null(val)){
        val = seg[, setdiff(names(seg), unlist(field.aliases)), drop = FALSE] 
    }

    seg = seg[, intersect(names(seg), unlist(field.aliases)), drop = FALSE]

    for (field in setdiff(names(field.aliases), names(seg))){
        if (!(field %in% names(seg))){
            names(seg)[names(seg) %in% field.aliases[[field]]] = field;
        }
    }
    

    if (chr){
        if (!is.null(seg$chr)){
            if (!grepl('chr', seg$chr[1])){
                seg$chr = paste('chr', seg$chr, sep = "");
            }
        } else{
            message('seg$chr field already exists. As seg$chr != NULL, flag "chr=TRUE" had no effect.')
        }
    }

    if (is.null(seg$pos2)){
        seg$pos2 = seg$pos1
    }

    missing.fields = setdiff(names(field.aliases), c(names(seg), c('chr', 'ID', 'strand')));

    if (length(missing.fields)>0){
        warning(sprintf('Warning: seg file format problem, missing an alias for the following fields:\n\t%s',
                        paste(sapply(missing.fields, function(x) paste(x, '(can also be', paste(field.aliases[[x]], collapse = ', '), ')')), collapse = "\n\t")));
    }

    if (!is.null(val)){
        seg = cbind(seg, val)
    }

    return(seg)
}




#' @name gr.nochr
#' @title Remove chr prefix from GRanges seqlevels
#' @description
#'
#' Remove chr prefix from GRanges seqlevels
#'
#' @param gr \code{GRanges} with chr seqlevel prefixes
#' @return GRanges without chr seqlevel prefixes
#' @export
gr.nochr = function(gr) {
    if (grepl('^chr', seqlevels(gr)[1])){
        seqlevels(gr) = gsub('^chr','', seqlevels(gr))
    }
    return(gr)
}




#' @name gr.findoverlaps
#' @title Wrapper to \code{GenomicRanges::findOverlaps} with added functionality
#' @description
#'
#' Wrapper to \code{GenomicRanges::findOverlaps} with added functionality
#'
#' Returns \code{GRanges} of matches with two additional fields:
#' \itemize{
#' \item \code{$query.id} - index of matching query
#' \item \code{$subject.id} - index of matching subject
#' }
#' Optional \code{"by"} field is a character scalar that specifies a metadata column present in both query and subject
#' that will be used to additionally restrict matches, i.e. to pairs of ranges that overlap and also
#' have the same values of their \code{"by"} fields
#'
#' @importFrom GenomeInfoDb seqlengths seqlengths<-
#' @importFrom GenomicRanges values ranges width strand values<- strand<- seqnames
#' @importFrom data.table is.data.table := setkeyv
#' @importFrom IRanges findOverlaps
#' @importFrom S4Vectors queryHits subjectHits
#' @param query Query \code{GRanges} pile
#' @param subject Subject \code{GRanges} pile
#' @param ignore.strand boolean Strand ignored (default = TRUE)
#' @param first boolean Flag if TRUE restricts to only the first match of the subject. If FALSE will return all matches. (default = FALSE)
#' @param qcol \code{character} vector of query meta-data columns to add to results (default = NULL)
#' @param scol \code{character} vector of subject meta-data columns to add to results (default = NULL)
#' @param type \code{type} argument as defined by \code{IRanges::findOverlaps} (\code{"any"}, \code{"start"}, \code{"end"}, \code{"within"}, \code{"equal"}). (default = 'any')
#' @param by vector Meta-data column to consider when performing overlaps (default = NULL)
#' @param return.type character Select data format to return (supplied as character): \code{"same"}, \code{"data.table"}, \code{"GRanges"}. (default = 'same')
#' @param max.chunk integer Maximum number of \code{query*subject} ranges to consider at once. Lower number increases runtime but decreased memory. If \code{length(query)*length(subject)} is less than \code{max.chunk}, overlaps will run in one batch. (default = 1e3)
#' @param verbose boolean Flag to increase the verbosity. (default = FALSE)
#' @param mc.cores Number of cores to use when running in chunked mode (default = 1)
#' @param ... Additional arguments sent to \code{IRanges::findOverlaps}.
#' @return \code{GRanges} pile of the intersection regions, with \code{query.id} and \code{subject.id} marking sources
#' @export
gr.findoverlaps = function(query, subject, ignore.strand = TRUE, first = FALSE, qcol = NULL, scol = NULL, type = 'any', by = NULL, return.type = 'same', max.chunk = 1e13, verbose = FALSE, mc.cores = 1, maxgap = -1L, ...)
{

    subject.id = query.id = i.start = i.end = NULL ## for NOTE

    isdt <- any(class(query) == 'data.table' )

    if (return.type == 'same'){
        return.type <- ifelse(isdt, 'data.table', 'GRanges')
    }

    if (!return.type %in% c("data.table", "GRanges")){
        stop('Error: argument "return.type" must be one of: same, data.table, GRanges')
    }

    if (!((inherits(subject, 'GRanges') | inherits(subject, 'data.table')) & (inherits(query, 'GRanges') | inherits(query, 'data.table'))))
        stop('Error: Both subject and query have to be GRanges or data.table')

    if (!is.null(qcol)){
        if (!all(qcol %in% names(values(query)))){
            stop('Error: Some qcol are not present in meta data of query')
        }
    }

    if (!is.null(scol)){
        if (!all(scol %in% names(values(subject)))){
          stop('Error: Some scol are not present in meta data of subject')
        }
    }

    if (!is.null(by)){
        if (!(all(by %in% names(values(query))) & all(by %in% names(values(subject))))){
           stop('Error: "by" field must be meta data column of both query and subject')
        }
    }

    if (is.data.table(query)){
       query   = dt2gr(query)
    }
    if (is.data.table(subject)){
       subject = dt2gr(subject)
    }

    ss <- seqinfo(query)
    ## chunked operation


    if ((as.numeric(length(query)) * as.numeric(length(subject))) > max.chunk){
        if (verbose){
            cat('Overflow .. computing overlaps in chunks.  Adjust max.chunk parameter to gr.findoverlaps to avoid chunked computation\n')
        }


        chunk.size = floor(sqrt(max.chunk));
        ix1 = c(seq(1, length(query), chunk.size), length(query)+1)
        ix2 = c(seq(1, length(subject), chunk.size), length(subject)+1)
        ij = cbind(rep(1:(length(ix1)-1), length(ix2)-1), rep(1:(length(ix2)-1), each = length(ix1)-1))
        if (verbose){
            print(paste('Number of chunks:', nrow(ij)))
        }

        out = rbindlist(parallel::mclapply(1:nrow(ij),
            function(x){
                if (verbose){
                    cat(sprintf('chunk i = %s-%s (%s), j = %s-%s (%s)\n', ix1[ij[x,1]], ix1[ij[x,1]+1]-1, length(query),
                        ix2[ij[x,2]], (ix2[ij[x,2]+1]-1), length(subject)))
                }
                i.chunk = ix1[ij[x,1]]:(ix1[ij[x,1]+1]-1)
                j.chunk = ix2[ij[x,2]]:(ix2[ij[x,2]+1]-1)
                out = gr.findoverlaps(query[i.chunk], subject[j.chunk],  ignore.strand = ignore.strand, first = first, by = by, qcol = qcol, verbose = verbose, scol = scol, type = type, max.chunk = Inf, ...)
                out$query.id = i.chunk[out$query.id]
                out$subject.id = j.chunk[out$subject.id]
                return(as.data.table(out))
            }, mc.cores = mc.cores), fill = TRUE)

        out = dt2gr(out, seqinfo = ss)

        ## sort by position, then sort by query, then subject id
        out <- sort(out[order(out$query.id, out$subject.id)])

        convert = FALSE

        if ((return.type == 'same' & is(query, 'data.table')) | return.type == 'data.table') {
            out = gr2dt(out)
        } else{
            out = gr.fix(out, ss)
        }

        return(out)

    }

    ## perform the actual overlaps
    h <- tryCatch(GenomicRanges::findOverlaps(query, subject, type = type, ignore.strand = ignore.strand, maxgap = maxgap, ...), error = function(e) NULL)
    ## if any seqlengths badness happens overrride
    if (is.null(h)){
        warning('findOverlaps applied to ranges with non-identical seqlengths')
        query = gr.fix(query, subject)
        subject = gr.fix(subject, query)
        h <- GenomicRanges::findOverlaps(query, subject, type = type, ignore.strand = ignore.strand, maxgap = maxgap, ...)
    }



    ## r <- ranges(h, ranges(query), ranges(subject))
    r <- overlapsRanges(GenomicRanges::ranges(query), GenomicRanges::ranges(subject), h)
    h.df <- data.table(start = start(r), end = end(r), query.id = queryHits(h),
                     subject.id = subjectHits(h), seqnames = as.character(seqnames(query)[queryHits(h)]))

    ## add the seqnames, and subset if have a "by"
    if (nrow(h.df) > 0) {
        if (!is.null(by)) {
            by.query <- values(query)[h.df$query.id, by]
            by.subject <- values(subject)[h.df$subject.id, by]

            if (length(by)==1){
                keep.ix = by.query == by.subject
            }
            else{
                keep.ix = apply(as.matrix(by.query) == as.matrix(by.subject), 1, all)
            }

        h.df = h.df[keep.ix, ]
        }
    }

    ## if empty, return now
    if (nrow(h.df) == 0) {
        if (return.type == "GRanges"){
            return(GRanges(seqlengths = seqlengths(query)))
        } else{
            return(data.table())
        }
    }

    ## write the strand
    if (!ignore.strand){
        h.df$strand <- as.character(strand(query)[h.df$query.id])
    } else{
        h.df$strand = '*'
    }

    ## limit to first hits?
    if (first){
        h.df = h.df[!duplicated(h.df$query.id), ]
    }

    ## format into correct output format
  if (return.type=='GRanges') {
        out.gr = dt2gr(h.df, seqlengths = seqlengths(query))

        if (!is.null(qcol)){
          new.cols = as.data.frame(values(query)[out.gr$query.id, qcol, drop = FALSE])
          colnames(new.cols) = qcol
          values(out.gr) = cbind(as.data.frame(values(out.gr)), new.cols)
        }

        if (!is.null(scol)){
          new.cols = as.data.frame(values(subject)[out.gr$subject.id, scol, drop = FALSE])
          colnames(new.cols) = scol
          values(out.gr) = cbind(as.data.frame(values(out.gr)), new.cols)
        }

        out.gr <- gr.fix(out.gr, ss)
        ## sort by position, then sort by query, then subject id
        return(sort(out.gr[order(out.gr$query.id, out.gr$subject.id)]))
    } else {
        ## return data.table
        if (!is.null(qcol)){
            h.df = cbind(h.df, data.table::as.data.table(as.data.frame(values(query))[h.df$query.id, qcol, drop = FALSE]))
        }

        if (!is.null(scol)){
            h.df = cbind(h.df, data.table::as.data.table(as.data.frame(values(subject))[h.df$subject.id, scol, drop = FALSE]))
        }

        if ('i.start' %in% colnames(h.df)){
            h.df[, i.start := NULL]
        }

        if ('i.end' %in% colnames(h.df)){
            h.df[, i.end := NULL]
        }

        return(h.df)
    }
}








#' @name grl.eval
#' @title Evaluate and aggregate expression on GRanges column in GRangesList
#' @description
#'
#' Evaluate expression 'expr' on indivdual GRanges inside GRangesList.
#' Expression should result in a single i.e. scalar value per GRangesList item.
#'
#' @param grl GRangesList to evaluate over
#' @param expr Any syntactically valid R expression, on columns of GRanges or GRangesList
#' @param condition Optional: any syntactically valid R expression, on columns of GRanges or GRangesList, on which to subset before evaluating main 'expr' (default = NULL)
#' @return GRangesList evaluated by R expressions
#' @export
grl.eval = function(grl, expr, condition = NULL)
{
    expr = paste(deparse(substitute(expr)), collapse = ' ')

    dt = gr2dt(grl.unlist(grl))

    by = 'grl.ix'

    by = paste(by, collapse = ',')

    condition = substitute(condition)

    if (is.null(condition)){
        condition = ''
    } else{
        condition = paste(deparse(condition), collapse = ' ')
    }

    cmd = sprintf('dt[%s, .(val = %s), keyby = .(%s)]', condition, expr, by)

    dtg = tryCatch(
        eval(parse(text = cmd)), error = function(e) NULL)

    if (is.null(dtg)){
        ## hack to instantiate variables parent environment to make them accessible
        ## (since I don't have a clue about R environemnts)
        pf = as.list(parent.frame())

        tryCatch({
            for (nm in names(pf))
                eval(parse(text =  paste("tryCatch(", paste(nm, ' = pf[[nm]]'), ", error = function(e) NULL")))
        }, error = function(e) NULL)

        ## try again
        dtg = eval(parse(text = cmd))
    }

    return(dtg[.(1:length(grl)), val])

}



#' Minimal overlaps for GRanges/GRangesList
#'
#' Takes any number of GRanges or GRangesList and reduces them to the minimal
#' set of overlapping windows, ignoring strand (optional).  Can also
#' collapse only within levels of a meta data field "by"
#'
#' Will populate output with metadata of first row of input contributing the reduced output range.
#'
#' @param ... \code{GRanges} or \code{GRangesList}
#' @return GRanges
#' @export
gr.reduce = function(..., by = NULL, ignore.strand = TRUE, span = FALSE) {
    input = do.call(grbind, list(...))
    if (length(input)==0)
        return(input)
    input.meta = values(input)
    if (is.null(by))
        values(input) = data.frame(i = 1:length(input), bykey = 1)
    else
        values(input) = data.frame(i = 1:length(input), bykey = values(input)[, by])

    if (span)
    {
        if (ignore.strand)
            out = seg2gr(gr2dt(input)[, data.frame(i = i[1], start = min(start), end = max(end)), keyby = list(seqnames, bykey)])
        else
            out = seg2gr(gr2dt(input)[, data.frame(i = i[1], start = min(start), end = max(end)), keyby = list(seqnames, strand, bykey)])
    }
    else
    {
        if (ignore.strand)
            out = seg2gr(gr2dt(input)[, cbind(i = i[1], as.data.frame(reduce(IRanges(start, end)))), keyby = list(seqnames, bykey)])
        else
            out = seg2gr(gr2dt(input)[, cbind(i = i[1], as.data.frame(reduce(IRanges(start, end)))), keyby = list(seqnames, strand, bykey)])
    }

    values(out) = input.meta[out$i, , drop = FALSE]

    return(out)
}


#' @name gr.merge
#' @title merge GRanges by using coordinates as primary key
#' @description
#'
#' Uses gr.findoverlaps() to enable internal and external joins of GRanges using
#' syntax similar to "merge" where merging is done using coordinates +/- "by" fields
#'
#' Uses gr.findoverlaps() / GRanges::findOverlaps for heavy lifting, but returns outputs with
#' metadata populated as well as query and subject ids.  For external joins,
#' overlaps x with gaps(y) and gaps(x) with y.
#'
#' @param query GRanges Set of GRanges to query. Refer to gr.findoverlaps() and GenomicRanges::findOverlaps()
#' @param subject GRanges Set of GRanges as 'subject' in query. Refer to  gr.findoverlaps() and GenomicRanges::findOverlaps()
#' @param by vector Additional metadata fields to join on
#' @param all boolean Flag whether to include left and right joins
#' @param all.query boolean Flag whether to do a left join (default = all)
#' @param all.subject boolean Flag whether to do a right join (default = all)
#' @param ignore.strand boolean Strand information ignored (default = TRUE)
#' @param verbose boolean 
#' @return GRanges merged on 'by' vector
#' @export
gr.merge = function(query, subject, by = NULL, all = FALSE, all.query = all, all.subject = all, ignore.strand = TRUE, verbose = FALSE, ... )
{
    qcol = names(values(query))
    scol = names(values(subject))

    if (!is.null(by)){
        qcol = setdiff(qcol, by)
    }

    ov = gr.findoverlaps(query, subject, qcol = qcol, scol = scol, by = by, ...)

    if (verbose){
        message('inner join yields ', length(ov), ' overlaps')
    }


    if (all.query){
        if (is.null(by))
        {
            if (ignore.strand){
                sgaps = gaps(gr.stripstrand(subject)) %Q% (strand == '*')
            }
            else{
                sgaps = gaps(subject)
            }
        }
        else ## if by not null then we want to define gaps in a group wise manner ..
        {
            sgaps = unlist(do.call('GRangesList', lapply(split(subject, apply(as.matrix(values(subject)[, by]), 1, paste, collapse = ' ')), function(group)
            {
                if (ignore.strand){
                    this.gap = gaps(gr.stripstrand(group)) %Q% (strand == '*')
                }
                else{
                    this.gap = gaps(group)
                }

                values(this.gap)[, by] = values(group)[1, by]
                return(this.gap)
            })))
        }

        ov.left = gr.findoverlaps(query, sgaps, qcol = qcol, scol = NULL, by = by, ignore.strand = ignore.strand, ...)
        if (length(ov.left)>0){
            ov.left$subject.id = NA
        }

        if (verbose){
            message('left join yields ', length(ov.left), ' overlaps')
        }

        ov = grbind(ov, ov.left)
    }

    if (all.subject)
    {
        if (is.null(by))
        {
            if (ignore.strand){
                qgaps = gaps(gr.stripstrand(query)) %Q% (strand == '*')
            }
            else{
                qgaps = gaps(query)
            }
        }
        ## if by not null then we want to define gaps in a group wise manner ..
        else
        {
            qgaps = unlist(do.call('GRangesList', lapply(split(query, apply(as.matrix(values(query)[, by]), 1, paste, collapse = ' ')), function(group){

                if (ignore.strand){
                    this.gap = gaps(gr.stripstrand(group)) %Q% (strand == '*')
                }
                else{
                    this.gap = gaps(group)
                }

                values(this.gap)[, by] = values(group)[1, by]
                return(this.gap)
            })))
        }

        ov.right = gr.findoverlaps(qgaps, subject, qcol = NULL, scol = scol, by = by, ignore.strand = ignore.strand, ...)

        if (length(ov.right)>0){
            ov.right$query.id = NA
        }

        strand(ov.right) = strand(subject)[ov.right$subject.id]

        if (verbose){
            message('right join yields ', length(ov.right), ' overlaps')
        }
        ov = grbind(ov, ov.right)
    }

    return(ov)
}




#' @name gr.disjoin
#' @title GenomicRanges disjoin with some additional functionality
#' @description
#'
#' Identical to GRanges disjoin, except outputs inherit metadata from first overlapping parent instance on input
#'
#' @param x GRanges to disjoin
#' @param ... arguments to disjoin (e.g. with.revmap=FALSE). Please see documentation for GenomicRanges::disjoin()
#' @param ignore.strand logical scalar (default = TRUE)
#' @return GRanges of non-overlapping ranges with metadata
#' @export
gr.disjoin = function(x, ..., ignore.strand = TRUE)
{
    y = disjoin(x, ...)
    ix = gr.match(y, x, ignore.strand = ignore.strand)
    values(y) = values(x)[ix, , drop = FALSE]
    return(y)
}




#' @name gr.in
#' @title Versatile implementation of \code{GenomicRanges::findOverlaps}
#' @description
#'
#' Versatile implementation of \code{GenomicRanges::findOverlaps}
#'
#' Returns boolean vector. TRUE if query range i is found in any range in subject.
#'
#' @param query GRanges Set of GRanges to query. Refer to gr.findoverlaps() and GenomicRanges::findOverlaps()
#' @param subject GRanges Set of GRanges as 'subject' in query. Refer to gr.findoverlaps() and GenomicRanges::findOverlaps()
#' @param ... Arguments to be passed to \code{\link{gr.findoverlaps}}
#' @return boolean vector whereby TRUE is if query range i is found in any range in subject
#' @export
gr.in = function(query, subject, ...)
{

  if (is(query, 'character'))
    query = parse.gr(query)

  if (is(subject, 'character'))
    subject = parse.gr(subject)
    
  tmp = gr.findoverlaps(query, subject, ...)
  out = rep(FALSE, length(query))
  out[tmp$query.id] = TRUE
  
  return(out)
}




#' @name gr.sum
#' @title gr.sum
#' @description
#'
#' Sums GRanges either by doing coverage and either weighting them equally
#' or using a field "weight".  Will return either sum or average.
#'
#' @param gr \code{GRanges} to sum
#' @param field metadata field from gr to use as a weight
#' @param mean logical scalar specifying whether to divide the output at each interval but the total number of intervals overlapping it (only applies if field == NULL) (default FALSE)
#' @return non-overlapping GRanges spanning the seqlengths of gr with $score (if field is NULL) or $field specifying the sum / mean at that position
#' @export
gr.sum = function(gr, field = NULL, mean = FALSE)
{
  SHIFT = 0
  if (length(gr)>0)
    {
      SHIFT = pmax(SHIFT, -min(start(gr))-1)
    }
      
    if (is.null(field)){
        weight = rep(1, length(gr))
    }
    else{
        weight = values(gr)[, field]
    }

    ## in case we go "off the seqlengths" we recreate a GR with adjusted seqlengths so coverage will process correctly
    tmp.gr = GenomicRanges::shift(gr[, c()], SHIFT)
    tmp.gr = gr.fix(tmp.gr)
    out = GenomicRanges::shift(as(coverage(tmp.gr, weight = weight), 'GRanges'), -SHIFT)

    if (!is.null(field))
    {
        if (mean)
        {
            count = GenomicRanges::shift(as(coverage(tmp.gr, weight = 1), 'GRanges'), -SHIFT)
            out$score = out$score/count$score[gr.match(out, count)] ## divide by total count at each location
        }
        names(values(out))[length(names(values(out)))] = field
    }

    out = gr.fix(out, gr)
    return(out)
}




#' @name gr.collapse
#' @title Collapse adjacent ranges
#' @description
#'
#' Like \code{GenomicRanges::reduce} except only collapses <<adjacent>> ranges in the input
#'
#' @param gr \code{GRanges} to collapse
#' @param pad Padding that allows for not quite adjacent elements to be considered overlapping. (default = 1)
#' @return GRanges with collapsed adjacent GRanges
#' @author Marcin Imielinski
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
    else{
        return(gr[c()])
    }
}




#' @name gr.match
#' @title Alternative \code{GenomicRanges::match} that accepts additional \code{\link{gr.findoverlaps}} options
#' @description
#'
#' Alternative \code{GenomicRanges::match} that accepts additional \code{\link{gr.findoverlaps}} options
#'
#' Wrapper to \code{GenomicRanges::match} (uses \code{\link{gr.findoverlaps}}). This allows users to
#' match on additional \code{by} fields, or chunk into smaller pieces for lower memory.
#'
#' @param query Query \code{GRanges} pile
#' @param subject Subject \code{GRanges} pile
#' @param max.slice max slice of query to match at a time
#' @param verbose whether to give verbose outputx
#' @importFrom parallel mclapply
#' @param ... Additional arguments to be passed along to \code{\link{gr.findoverlaps}}.
#' @return Vector of length = \code{length(query)} with subject indices of *first* subject in query, or NA if none found.
#' This behavior is different from \code{\link{gr.findoverlaps}}, which will
#' return *all* indicies of subject in query (in the case of one query overlaps with multiple subject)
#' ... = additional args for findOverlaps (IRanges version)
#' @author Marcin Imielinski
#' @export
gr.match = function(query, subject, max.slice = Inf, verbose = FALSE, ...)
{

    if (length(query)>max.slice)
    {
        ix.l = split(1:length(query), ceiling(as.numeric((1:length(query)/max.slice))))
        return(do.call('c', parallel::mclapply(ix.l, function(ix) {
            if (verbose){
                cat(sprintf('Processing %s to %s\n', min(ix), max(ix)))
            }
 
            gr.match(query[ix, ], subject, verbose = verbose, ...)
        })))
    }

  tmp = gr.findoverlaps(query, subject, ...)

  if (length(tmp)>0)
    {
      tmp = gr2dt(tmp)[order(subject.id), ][!duplicated(query.id), ]
    }
  out = rep(NA, length(query))
  out[tmp$query.id] = tmp$subject.id
  return(out)
}



setGeneric('%+%', function(gr, x) standardGeneric('%+%'))

#' @name %+%
#' @title Nudge GRanges right
#' @description
#'
#' Operator to shift GRanges right "sh" bases
#'
#' @return shifted granges
#' @rdname gr.nudge-shortcut
#' @exportMethod %+%
#' @aliases %+%,GRanges-method
#' @author Marcin Imielinski
#' @export
setMethod("%+%", signature(gr = 'GRanges'), function(gr, x) {
  end(gr) = end(gr)+x
  start(gr) = start(gr)+x
  return(gr)
})


setGeneric('%-%', function(gr, x) standardGeneric('%-%'))

#' @name %-%
#' @title Shift GRanges left
#' @description
#'
#' Operator to shift GRanges left "sh" bases
#'
#' df %!% c('string.*to.*match', 'another.string.to.match')
#'
#' @return shifted GRanges
#' @rdname gr.nudge
#' @export
#' @author Marcin Imielinski
setMethod("%-%", signature(gr = 'GRanges'), function(gr, x) {
    start(gr) = start(gr)-x
    end(gr) = end(gr)-x
    return(gr)
})


setGeneric('%&%', function(x, y) standardGeneric('%&%'))

#' @name %&%
#' @title subset x on y ranges while ignoring strand (strand-agnostic)
#' @description
#'
#' shortcut for x[gr.in(x,y)]
#'
#' gr1 %&% gr2 returns the subsets of gr1 that overlaps gr2
#'
#' @return subset of gr1 that overlaps gr2
#' @rdname gr.in-shortcut
#' @exportMethod %&%
#' @aliases %&%, GRanges-method
#' @author Marcin Imielinski
setMethod("%&%", signature(x = 'GRanges'), function(x, y) {
    if (is.character(y)){
        y = parse.gr(y)
    }
    return(x[gr.in(x, y)])
})

.grlandfun = function(x, y) {
  if (!length(x))
    return(x)

    if (is.character(y)){
        y = parse.gr(y)
    }
    gr = unlist(x, use.names = FALSE)
    if (inherits(y, "GRanges")) {
        w_in = gr.in(gr, y)
    } else if (inherits(y, "GRangesList")) {
        w_in = gr.in(gr, unlist(y))
    }
    split_by = rep(1:length(x), times = elementNROWS(x))[w_in]
    grl_ix = unique(split_by)
    ret_grl = split(gr[w_in], split_by)
    names(ret_grl) = names(x)[grl_ix]
    return(ret_grl)
}


#' @name %&%.GRangesList
#' @title subset x on y ranges while ignoring strand (strand-agnostic)
#' @description
#'
#' shortcut for x[gr.in(x,y)]
#'
#' gr1 %&% gr2 returns the subsets of gr1 that overlaps gr2
#'
#' @return subset of gr1 that overlaps gr2
#' @rdname gr.in-shortcut
#' @exportMethod %&%
#' @aliases %&%, GRangesList-method
#' @author Marcin Imielinski
setMethod("%&%", signature(x = 'GRangesList'), .grlandfun)

#' @name %&%.CompressedGRangesList
#' @title subset x on y ranges while ignoring strand (strand-agnostic)
#' @description
#'
#' shortcut for x[gr.in(x,y)]
#'
#' gr1 %&% gr2 returns the subsets of gr1 that overlaps gr2
#'
#' @return subset of gr1 that overlaps gr2
#' @rdname gr.in-shortcut
#' @exportMethod %&%
#' @aliases %&%, CompressedGRangesList-method
#' @author Marcin Imielinski
setMethod("%&%", signature(x = 'CompressedGRangesList'), .grlandfun)


setGeneric('%&&%', function(x, y) standardGeneric('%&&%'))

#' @name %&&%
#' @title Subset x on y ranges, strand-specific
#' @description
#'
#' shortcut for x[gr.in(x,y)]
#'
#' gr1 %&&% gr2 returns the subsets of gr1 that overlaps gr2
#'
#' @return subset of gr1 that overlaps gr2
#' @rdname gr.in-strand-shortcut
#' @exportMethod %&&%
#' @aliases %&&%,GRanges-method
#' @author Marcin Imielinski
setMethod("%&&%", signature(x = "GRanges"), function(x, y) {
    if (is.character(y)){
        y = parse.gr(y)
    }
    return(x[gr.in(x, y, ignore.strand = FALSE)])
})

.grlandfunsign = function(x, y) {
    if (is.character(y)){
        y = parse.gr(y)
    }
    gr = unlist(x, use.names = FALSE)
    if (inherits(y, "GRanges")) {
        w_in = gr.in(gr, y, ignore.strand = FALSE)
    } else if (inherits(y, "GRangesList")) {
        w_in = gr.in(gr, unlist(y), ignore.strand = FALSE)
    }
    split_by = rep(1:length(x), times = elementNROWS(x))[w_in]
    grl_ix = unique(split_by)
    ret_grl = split(gr[w_in], split_by)
    names(ret_grl) = names(x)[grl_ix]
    return(ret_grl)
}

#' @name %&&%.GRangesList
#' @title Subset x on y ranges, strand-specific
#' @description
#'
#' shortcut for x[gr.in(x,y)]
#'
#' gr1 %&&% gr2 returns the subsets of gr1 that overlaps gr2
#'
#' @return subset of gr1 that overlaps gr2
#' @rdname gr.in-strand-shortcut
#' @exportMethod %&&%
#' @aliases %&&%,GRangesList-method
#' @author Marcin Imielinski
setMethod("%&&%", signature(x = 'GRangesList'), .grlandfunsign)

#' @name %&&%.CompressedGRangesList
#' @title Subset x on y ranges, strand-specific
#' @description
#'
#' shortcut for x[gr.in(x,y)]
#'
#' gr1 %&&% gr2 returns the subsets of gr1 that overlaps gr2
#'
#' @return subset of gr1 that overlaps gr2
#' @rdname gr.in-strand-shortcut
#' @exportMethod %&&%
#' @aliases %&&%,CompressedGRangesList-method
#' @author Marcin Imielinski
setMethod("%&&%", signature(x = 'CompressedGRangesList'), .grlandfunsign)


setGeneric('%O%', function(x, ...) standardGeneric('%O%'))

#' @name %O%
#' @title gr.val shortcut to get fractional overlap of gr1 by gr2, strand-agnostic
#' @description
#'
#' Shortcut for gr.val (using val = names(values(y)))
#'
#' gr1 %O% gr2
#'
#' @return fractional overlap of gr1 with gr2
#' @rdname gr.val-O
#' @exportMethod %O%
#' @aliases %O%,GRanges-method
#' @author Marcin Imielinski
#' @param x See \link{gr.val}
#' @param ... See \link{gr.val}
setMethod("%O%", signature(x = "GRanges"), function(x, y) {

    query.id = NULL; ## NOTE fix
    ov = gr2dt(gr.findoverlaps(x, reduce(gr.stripstrand(y))))
    if (nrow(ov) > 0) {
        ov = ov[ , sum(width), keyby = query.id]
        x$width.ov = 0
        x$width.ov[ov$query.id] = ov$V1
        return(x$width.ov/width(x))
    } else{
        return(rep(0, length(x)))
    }
})




#' @name %OO%
#' @title gr.val shortcut to get fractional overlap of gr1 by gr2, strand-specific
#' @description
#'
#' Shortcut for gr.val (using val = names(values(y)))
#'
#' gr1 %OO% gr2
#'
#' @return fractional overlap  of gr1 with gr2
#' @docType methods
#' @aliases %OO%,GRanges-method
#' @rdname gr.val-fractional
#' @exportMethod %OO%
#' @author Marcin Imielinski
#' @param x See \link{gr.val}
#' @param ... See \link{gr.val}
setGeneric('%OO%', function(x, ...) standardGeneric('%OO%'))
setMethod("%OO%", signature(x = "GRanges"), function(x, y) {
    query.id = NULL; ## NOTE fix
    ov = gr2dt(gr.findoverlaps(x, reduce(y), ignore.strand = FALSE))
    if (nrow(ov)>0){
        ov = ov[ , sum(width), keyby = query.id]
        x$width.ov = 0
        x$width.ov[ov$query.id] = ov$V1
        return(x$width.ov/width(x))
    } else{
        return(rep(0, length(x)))
    }
})




#' @name %o%
#' @title gr.val shortcut to total per interval width of overlap of gr1 with gr2, ignoring strand
#' @description
#'
#' Shortcut for gr.val (using val = names(values(y)))
#'
#' gr1 %o% gr2
#'
#' @return bases overlap of gr1 with gr2
#' @rdname gr.val-total
#' @exportMethod %o%
#' @aliases %o%,GRanges-method
#' @author Marcin Imielinski
#' @param x See \link{gr.val}
#' @param ... See \link{gr.val}
setGeneric('%o%', function(x, ...) standardGeneric('%o%'))
setMethod("%o%", signature(x = "GRanges"), function(x, y) {
    query.id = NULL; ## NOTE fix
    ov = gr2dt(gr.findoverlaps(x, reduce(gr.stripstrand(y))))
    if (nrow(ov)>0){
        ov = ov[ , sum(width), keyby = query.id]
        x$width.ov = 0
        x$width.ov[ov$query.id] = ov$V1
        return(x$width.ov)
    } else{
        return(rep(0, length(x)))
    }
})




#' @name %oo%
#' @title gr.val shortcut to total per interval width of overlap of gr1 with gr2, strand-specific
#' @description
#'
#' gr1 %oo% gr2
#'
#' @param x See \link{gr.val}
#' @param ... See \link{gr.val}
#' @return bases overlap  of gr1 with gr2
#' @rdname gr.val-strand
#' @aliases %oo%,GRanges-method
#' @exportMethod %oo%
#' @author Marcin Imielinski
setGeneric('%oo%', function(x, ...) standardGeneric('%oo%'))
setMethod("%oo%", signature(x = "GRanges"), function(x, y) {
    query.id = NULL; ## NOTE fix
    ov = gr2dt(gr.findoverlaps(x, y, ignore.strand = FALSE))
    if (nrow(ov)>0){
        ov = ov[ , sum(width), keyby = query.id]
        x$width.ov = 0
        x$width.ov[ov$query.id] = ov$V1
        return(x$width.ov)
    } else{
        return(rep(0, length(x)))
    }
})




#' @name %N%
#' @title gr.val shortcut to get total numbers of intervals in gr2 overlapping with each interval in  gr1, ignoring strand
#' @description
#'
#' gr1 %N% gr2
#'
#' @return bases overlap of gr1 with gr2
#' @rdname gr.val-shortcut
#' @docType methods
#' @aliases %N%,GRanges-method
#' @exportMethod %N%
#' @author Marcin Imielinski
#' @param x See \link{gr.val}
#' @param ... See \link{gr.val}
setGeneric('%N%', function(x, ...) standardGeneric('%N%'))
setMethod("%N%", signature(x = "GRanges"), function(x, y) {
    query.id = NULL; ## NOTE fix
    ov = gr2dt(gr.findoverlaps(x, y, ignore.strand = TRUE))
    if (nrow(ov)>0){
        ov = ov[ , length(width), keyby = query.id]
        x$width.ov = 0
        x$width.ov[ov$query.id] = ov$V1
        return(x$width.ov)
    } else{
        return(rep(0, length(x)))
    }
})




#' @name %NN%
#' @title gr.val shortcut to get total numbers of intervals in gr2 overlapping with each interval in  gr1, respecting strand
#' @description
#'
#' gr1 %NN% gr2
#'
#' @return bases overlap  of gr1 with gr2
#' @rdname gr.val-numbers
#' @aliases %NN%,GRanges-method
#' @exportMethod %NN%
#' @author Marcin Imielinski
#' @param x See \link{gr.val}
#' @param ... See \link{gr.val}
setGeneric('%NN%', function(x, ...) standardGeneric('%NN%'))
setMethod("%NN%", signature(x = "GRanges"), function(x, y) {
    query.id = NULL; ## NOTE fix
    ov = gr2dt(gr.findoverlaps(x, y, ignore.strand = FALSE))
    if (nrow(ov)>0){
        ov = ov[ , length(width), keyby = query.id]
        x$width.ov = 0
        x$width.ov[ov$query.id] = ov$V1
        return(x$width.ov)
    } else{
        return(rep(0, length(x)))
    }
    return(x$width.ov)
})




#' @name %_%
#' @title \code{BiocGenerics::setdiff} shortcut (strand agnostic)
#' @description
#'
#' Shortcut for \code{BiocGenerics::setdiff}
#'
#' gr1 <- GRanges(1, IRanges(10,20), strand="+")
#' gr2 <- GRanges(1, IRanges(15,25), strand="-")
#' gr3 <- "1:1-15"
#' gr1 %_% gr2
#' gr1 %_% gr3
#'
#' @param x \code{GRanges} object to to
#' @param ... A \code{GRanges} or a character to be parsed into a \code{GRanges}
#' @return \code{GRanges} representing setdiff of input interval
#' @rdname gr.setdiff
#' @aliases %_%,GRanges-method
#' @exportMethod %_%
#' @author Marcin Imielinski
setGeneric('%_%', function(x, ...) standardGeneric('%_%'))
setMethod("%_%", signature(x = "GRanges"), function(x, y) {
    if (is.character(y)){
        y = parse.gr(y)
    }
    setdiff(gr.stripstrand(x[, c()]), gr.stripstrand(y[, c()]))
})



#' @name %)%
#' @title \code{gr.end} shortcut (strand agnostic)
#' @description
#'
#' Shortcut for \code{gr.end}
#'
#' gr1 <- GRanges(1, IRanges(10,20), strand="+")
#' gr1 %)% 5 # right 5 bases of gr1
#'
#' @param x \code{GRanges} object 
#' @param y num bases to get from right side of interval
#' @return \code{GRanges} right sides of input intervals
#' @rdname gr.end
#' @aliases %)%,GRanges-method
#' @exportMethod %)%
#' @author Marcin Imielinski
setGeneric('%)%', function(x, ...) standardGeneric('%)%'))
setMethod("%)%", signature(x = "GRanges"), function(x, y) {
  gr.end(x, y, ignore.strand = TRUE)
})

#' @name %))%
#' @title \code{gr.end} shortcut (strand agnostic)
#' @description
#'
#' Shortcut for \code{gr.end}
#'
#' gr1 <- GRanges(1, IRanges(10,20), strand="+")
#' gr1 %))% 5 # get three prime 5 bases of gr1
#'
#' @param x \code{GRanges} object 
#' @param y num bases to get from right side of interval
#' @return \code{GRanges} 3' sides of input intervals
#' @rdname gr.end
#' @aliases %))%,GRanges-method
#' @exportMethod %)%
#' @author Marcin Imielinski
setGeneric('%))%', function(x, ...) standardGeneric('%))%'))
setMethod("%))%", signature(x = "GRanges"), function(x, y) {
  gr.end(x, y, ignore.strand = FALSE)
})


#' @name %(%
#' @title \code{gr.end} shortcut (strand agnostic)
#' @description
#'
#' Shortcut for \code{gr.end}
#'
#' gr1 <- GRanges(1, IRanges(10,20), strand="+")
#' gr1 %(% 5 # left 5 bases of gr1
#'
#' @param x \code{GRanges} object 
#' @param y num bases to get from left side of interval
#' @return \code{GRanges} left sides of input intervals
#' @rdname gr.end
#' @aliases %(%,GRanges-method
#' @exportMethod %(%
#' @author Marcin Imielinski
setGeneric('%(%', function(x, ...) standardGeneric('%(%'))
setMethod("%(%", signature(x = "GRanges"), function(x, y) {
  gr.start(x, y, ignore.strand = TRUE)
})

#' @name %((%
#' @title \code{gr.end} shortcut (strand agnostic)
#' @description
#'
#' Shortcut for \code{gr.end}
#'
#' gr1 <- GRanges(1, IRanges(10,20), strand="+")
#' gr1 %((% 5 # get three prime 5 bases of gr1
#'
#' @param x \code{GRanges} object 
#' @param y num bases to get from left side of interval
#' @return \code{GRanges} 5' sides of input intervals
#' @rdname gr.end
#' @aliases %((%,GRanges-method
#' @exportMethod %)%
#' @author Marcin Imielinski
setGeneric('%((%', function(x, ...) standardGeneric('%((%'))
setMethod("%((%", signature(x = "GRanges"), function(x, y) {
  gr.start(x, y, ignore.strand = FALSE)
})

#' @name %Q%
#' @title query ranges by applying an expression to ranges metadata
#' @description
#'
#' gr %Q% query returns the subsets of gr1 that matches meta data statement in query
#'
#' @return subset of gr that matches query
#' @rdname gr.query
#' @docType methods
#' @aliases %Q%,GRanges-method
#' @param x \code{GRanges} to match against a query \code{GRanges}
#' @param y \code{GRanges} with metadata to be queried
#' @export
#' @author Marcin Imielinski
setGeneric('%Q%', function(x, ...) standardGeneric('%Q%'))
setMethod("%Q%", signature(x = "GRanges"), function(x, y) {
    condition_call  = substitute(y)
    ## serious R voodoo gymnastics .. but I think finally hacked it to remove ghosts
    ## create environment that combines the calling env with the granges env
    env = as(c(as.list(parent.frame(2)), as.list(as.data.frame(x))), 'environment')
    parent.env(env) = parent.frame()
    ix = tryCatch(eval(condition_call, env), error = function(e) NULL)
    if (is.null(ix))
    {
        condition_call  = substitute(y)
        ix = eval(condition_call, GenomicRanges::as.data.frame(x))
    }
    return(x[ix])
})

.qgrlfun =  function(x, y) {
    condition_call  = substitute(y)
    ## serious R voodoo gymnastics .. but I think finally hacked it to remove ghosts
    ## create environment that combines the calling env with the granges env
    env = as(c(as.list(parent.frame(2)), as.list(mcols(x))), 'environment')
    parent.env(env) = parent.frame()
    ix = tryCatch(eval(condition_call, env), error = function(e) NULL)
    if (is.null(ix))
    {
        condition_call  = substitute(y)
        ix = eval(condition_call, GenomicRanges::as.data.frame(x))
    }
    return(x[ix])
}

setMethod("%Q%", signature(x = "GRangesList"),.qgrlfun)
setMethod("%Q%", signature(x = "CompressedGRangesList"), .qgrlfun)

.qqgrlfun = function(x, y) {
    tmp_vals = mcols(x)
    tmp_gr = unlist(x, use.names = FALSE)
    tmp_nm = names(tmp_gr)
    tmp_gr = setNames(tmp_gr, NULL)
    mcols(tmp_gr)[["split_by"]] = rep(1:length(x), times = elementNROWS(x))
    condition_call  = substitute(y)
    env = as(c(as.list(parent.frame(2)), as.list(as.data.frame(tmp_gr))), 'environment')
    parent.env(env) = parent.frame()
    ix = tryCatch(eval(condition_call, env), error = function(e) NULL)
    if (is.null(ix))
    {
        condition_call  = substitute(y)
        ix = eval(condition_call, GenomicRanges::as.data.frame(tmp_gr))
    }
    tmp_gr = tmp_gr[ix]
    split_by = rep(1:length(x), times = elementNROWS(x))[ix]
    names(tmp_gr) = tmp_nm[ix]
    ret_grl = split(tmp_gr, split_by)
    grl_id = unique(split_by)
    names(ret_grl) = names(x)[grl_id]
    mcols(ret_grl) = tmp_vals[grl_id,]
    return(ret_grl)
}

#' @name %QQ%
#' @title query ranges by applying an expression to GRanges metadata within a GRangesList
#' @description
#'
#' grl %QQ% query returns the subsets of GRanges elements within grl that matches meta data statement in query
#'
#' @return subset of grl that matches query
#' @rdname grl.query
#' @docType methods
#' @aliases %QQ%,GRanges-method
#' @param x \code{GRangesList}
#' @param y expression querying metadata columns evaluating to a logical or integer matching to a subset of GRanges elements of x 
#' @export
#' @author Kevin Hadi
setGeneric('%QQ%', function(x, ...) standardGeneric('%QQ%'))
setMethod(`%QQ%`, signature(x = "GRangesList"), .qqgrlfun)
setMethod(`%QQ%`, signature(x = "CompressedGRangesList"), .qqgrlfun)



#' @name %^%
#' @title gr.in shortcut
#' @description
#'
#' Shortcut for gr.in (standard arguments)
#'
#' gr1 %^% gr2
#'
#' @return logical vector of length gr1 which is TRUE at entry i only if gr1[i] intersects at least one interval in gr2 (strand agnostic)
#' @rdname gr.in-shortcut
#' @param x \code{GRanges} object
#' @param ... additional arguments to gr.in
#' @export
#' @docType methods
#' @aliases %^%,GRanges-method
#' @param x See \link{gr.in}
#' @param ... See \link{gr.in}
setGeneric('%^%', function(x, ...) standardGeneric('%^%'))
setMethod("%^%", signature(x = "GRanges"), function(x, y) {
    if (is.character(y)){
        y = parse.gr(y)
    }
    return(gr.in(x, y))
})



#' @name %^^%
#' @title grl.in shortcut
#' @description
#'
#' Shortcut for grl.in (standard arguments)
#'
#' grl %^% gr
#'
#' @return logical vector of length length(grl) which is TRUE at entry i only if grl[i] intersects at least one interval in gr (strand agnostic)
#' @rdname grl.in-shortcut
#' @param x \code{GRanges} object
#' @param ... additional arguments to gr.in
#' @export
#' @docType methods
#' @aliases %^%,GRanges-method
#' @param x See \link{gr.in}
#' @param ... See \link{gr.in}
setMethod("%^%", signature(x = "GRangesList"), function(x, y) {
    if (is.character(y)){
        y = parse.gr(y)
    }
    return(grl.in(x, y))
})

setMethod("%^%", signature(x = "CompressedGRangesList"), function(x, y) {
    if (is.character(y)){
        y = parse.gr(y)
    }
    return(grl.in(x, y))
})





#' @name %$%
#' @title gr.val shortcut to get mean values of subject "x" meta data fields in query "y" (strand-agnostic)
#' @description
#'
#' Shortcut for gr.val (using val = names(values(y)))
#'
#' gr1 %$% gr2
#'
#' @return gr1 with extra meta data fields populated from gr2
#' @rdname gr.val-shortcut
#' @docType methods
#' @param x \code{GRanges} object
#' @aliases %$%,GRanges-method
#' @exportMethod %$%
#' @author Marcin Imielinski

setGeneric('%$%', function(x, ...) standardGeneric('%$%'))
setMethod("%$%", signature(x = "GRanges"), function(x, y) {
    return(gr.val(x, y, val = names(values(y))))
})




#' @name %$$%
#' @title gr.val shortcut to get mean values of subject "x" meta data fields in query "y" (strand-specific)
#' @description
#'
#' Shortcut for gr.val (using val = names(values(y)))
#'
#' gr1 %$$% gr2
#'
#' @return gr1 with extra meta data fields populated from gr2
#' @rdname gr.val-mean
#' @exportMethod %$$%
#' @aliases %$$%,GRanges-method
#' @author Marcin Imielinski
#' @param x See \link{gr.val}
#' @param ... See \link{gr.val}
setGeneric('%$$%', function(x, ...) standardGeneric('%$$%'))
setMethod("%$$%", signature(x = "GRanges"), function(x, y) {
    if (is.character(y)){
        y = parse.gr(y)
    }
    return(gr.val(x, y, val = names(values(y)), ignore.strand = FALSE))
})




#' @name %*%
#' @title Metadata join with coordinates as keys (wrapper to \code{\link{gr.findoverlaps}})
#' @description
#'
#' Shortcut for gr.findoverlaps with \code{qcol} and \code{scol} filled in with all the query and subject metadata names.
#' This function is useful for piping \code{GRanges} operations together. Another way to think of %*% is as a
#' join of the metadata, with genomic coordinates as the keys. \cr
#' Example usage: \cr
#' x %*% y
#'
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
#' example_genes %*% example_dnase
## setGeneric('%*%', function(...) standardGeneric('%*%'))
## setGeneric('%*%')
setMethod("%*%", signature(x = "GRanges"), function(x, y) {
    gr = gr.findoverlaps(x, y, qcol = names(values(x)), scol = names(values(y)))
    return(gr)
})




#' @name %**%
#' @title shortcut for gr.findoverlaps (strand-specific)
#' @description
#'
#' Shortcut for gr.findoverlaps
#'
#' gr1 %**% gr2
#'
#' @return new granges containing every pairwise intersection of ranges in gr1 and gr2 with a join of the corresponding metadata
#' @rdname grfo-shortcut
#' @exportMethod %**%
#' @aliases %**%,GRanges-method
#' @author Marcin Imielinski
#' @param x See \link{gr.findoverlaps}
#' @param ... See \link{gr.findoverlaps}
setGeneric('%**%', function(x, ...) standardGeneric('%**%'))
setMethod("%**%", signature(x = "GRanges"), function(x, y) {
    if (is.character(y)){
        y = parse.gr(y)
    }
    gr = gr.findoverlaps(x, y, qcol = names(values(x)), scol = names(values(y)), ignore.strand = FALSE)
    return(gr)
})




#' @name %^^%
#' @title gr.in shortcut (strand-specific)
#' @description
#'
#' Shortcut for gr.in
#'
#' gr1 %^^% gr2
#'
#' @return logical vector of length gr1 which is TRUE at entry i only if gr1[i] intersects at least one interval in gr2
#' @rdname gr.in-shortcut
#' @aliases %^^%,GRanges-method
#' @exportMethod %^^%
#' @author Marcin Imielinski
#' @param x See \link{gr.in}
#' @param ... See \link{gr.in}
setGeneric('%^^%', function(x, ...) standardGeneric('%^^%'))
setMethod("%^^%", signature(x = "GRanges"), function(x, y) {
    if (is.character(y)){
        y = parse.gr(y)
    }
    return(gr.in(x, y, ignore.strand = FALSE))
})


#' @name gr.setdiff
#' @title gr.setdiff
#' @description
#'
#' More robust and faster implementation of GenomicRanges::setdiff()
#'
#' Robust to common edge cases of setdiff(gr1, gr2)  where gr2 ranges are contained inside gr1's (yieldings
#' setdiffs yield two output ranges for some of the input gr1 intervals.
#'
#' @param query \code{GRanges} object as query
#' @param subject \code{GRanges} object as subject
#' @param ignore.strand boolean Flag to ignore strands. Refer to 'gr.findoverlaps()'. (default = TRUE)
#' @param by vector Meta-data column to consider when performing overlaps. Refer to 'gr.findoverlaps()' documentation (default = NULL)
#' @param ... arguments to be passed to \link{gr.findoverlaps}
#' @return Returns indices of query in subject. If none found, NA
#' @export
gr.setdiff = function(query, subject, ignore.strand = TRUE, by = NULL, ...)
{
    ## in this case need to be careful about setdiffing only within the "by" level
  if (!is.null(by)){
    if (ignore.strand)
    {
      query = gr.stripstrand(query)
      subject = gr.stripstrand(subject)
    }
    tmp = gr2dt(subject)
    tmp$strand = factor(tmp$strand, c('+', '-', '*'))
    sl = seqlengths(query)
    gp = dt2gr( ## gaps using by
      tmp[,
          as.data.frame(
            gaps(GRanges(seqnames,             
                         IRanges(start, end),
                         seqlengths = sl,
                         strand = strand))),
          , by = by],      
      seqinfo = seqinfo(query))
    
    ## if any other bydiff need these gaps
    if (length(bydiff <- setdiff(values(query)[[by]],
                values(subject)[[by]]))>0)
    {
      gp = grbind(gp, query[values(query)[[by]] %in% bydiff])
    }

    if (ignore.strand)
      gp = gp %Q% (strand == '*')

  } else {
    ## otherwise easier
    if (ignore.strand){
      gp = gaps(gr.stripstrand(subject)) %Q% (strand == '*')
    } else{
      gp = gaps(subject)
    }
  }
  
  out = gr.findoverlaps(query, gp, qcol = names(values(query)), ignore.strand = ignore.strand, by = by, ...)
    return(out)
}



#' @name gr.simplify
#' @title Calculates pairwise distance for rearrangements represented by \code{GRangesList} objects
#' @description
#'
#' Calculates pairwise distance for rearrangements represented by \code{GRangesList} objects
#'
#' @param gr GRanges or GRangesList input
#' @param field character scalar, corresponding to value field of gr. (default = NULL)
#' @param val character scalar (default = NULL)
#' @param include.val boolean Flag will include in out gr values field of first matching record in input gr. (default = TRUE)
#' @param split boolean Flag to split the output into \code{GRangesList} split by \code{"field"}. (default = FALSE)
#' @param pad integer Pad ranges by this amount before doing merge. [1], which merges contiguous but non-overlapping ranges (default = 1)
#' @return Simplified GRanges with "field" populated with uniquely contiguous values
#' @export
gr.simplify = function(gr, field = NULL, val = NULL, include.val = TRUE, split = FALSE, pad = 1)
{
    tmp = as.logical(suppressWarnings(width(GenomicRanges::pintersect(ranges(gr[-length(gr)]), ranges(gr[-1]+pad), resolve.empty = 'max.start')) > 0) &
                     seqnames(gr[-length(gr)]) == seqnames(gr[-1]) & strand(gr[-length(gr)]) == strand(gr[-1]))

    tmp = as.vector(c(0, cumsum(!tmp)))

    if (!is.null(field)){
        tmp = paste(tmp, values(gr)[, field])
    }

    if (!is.null(val)){
        tmp = paste(tmp, val)
    }

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

    if (include.val){
        values(out) = values(gr)[ix, ]
    }

    if (split){
        if (!is.null(field)){
            out = GenomicRanges::split(out, values(gr)[ix, field])
        } else{
            out = GRangesList(out)
        }
    }

    return(out)
}




#' @name parse.gr
#' @title parse.gr
#' @description
#'
#' Quick function to parse GRanges from character vector IGV-/UCSC-style strings of format gr1;gr2;gr3
#' where each input GRanges is of the format chr:start-end[+/-]
#'
#' @param ... arguments to parse.grl i.e. character strings in UCSC style chr:start-end[+-]
#' @author Marcin Imielinski
#' @return GRanges parsed from IGV-/UCSC-style strings
#' @export
parse.gr = function(...)
{
    grl.unlist(parse.grl(...))
}




#' @name parse.grl
#' @title parse.grl
#' @description
#'
#' Quick function to parse \code{GRangesList} from character vector IGV / UCSC style strings of format gr1;gr2;gr3 where each gr is of format chr:start-end[+/-]
#'
#' @param x  character vector representing a GRangesList with UCSC style coordinates (chr:start-end[+-]) representing a [signed] Granges and  ";" separators within each item of x separating individaul each GRAnges
#' @param seqlengths named integer vector representing genome (default = hg_seqlengths())
#' @author Marcin Imielinski
#' @return GRangesList parsed from IGV-/UCSC-style strings
#' @export
parse.grl = function(x, seqlengths = hg_seqlengths(), meta = NULL)
{
    nm = names(x)
    tmp = strsplit(x, '\\s*[;\\,\\|]\\s*')
    tmp.u = unlist(tmp)
    tmp.u = gsub('\\,', '', tmp.u)
    tmp.id = rep(1:length(tmp), sapply(tmp, length))
    str = gsub('.*([\\+\\-])$','\\1', tmp.u)
    spl = strsplit(tmp.u, "[\\-\\+\\:]", perl = T)

    if (any(ix <- sapply(spl, length)==2)){
        spl[ix] = lapply(which(ix), function(x) spl[[x]][c(1:2,2)])
    }

    if (any(ix <- sapply(spl, length)!=3))
    {
        if (is.null(seqlengths)){
            stop('Error: Need to define genome boundaries to use chromosome only coordinate strings')
        }

        tmp.sigr = gUtils::si2gr(seqlengths)
        tmp.chrnames = sapply(spl[ix], function(x) x[[1]])
        ix.present = tmp.chrnames %in% names(tmp.sigr)

        if (any(ix.present))
        {
          key = tmp.chrnames[ix.present]
          tmp.grs = gUtils::gr.string(tmp.sigr[key], mb = F)
          spl[ix[ix.present]] = strsplit(tmp.grs, '[\\:\\-\\+]', perl = T)
        }
        
        if (any(!ix.present))
          spl = spl[!ix[!ix.present]]
    }

    if (length(spl)==0)
      return(c())

    if (any(ix <- !str %in% c('+', '-'))){
        str[ix] = '*'
    }
    df = cbind(as.data.frame(matrix(unlist(spl), ncol = 3, byrow = T), stringsAsFactors = F), str)
    names(df) = c('chr', 'start', 'end', 'strand')
    df$start = as.numeric(df$start)
    df$end = as.numeric(df$end)
    rownames(df) = NULL

    gr = gUtils::seg2gr(df, seqlengths = seqlengths)[, c()]
    grl = GenomicRanges::split(gr, tmp.id)
    names(grl) = nm
    if (!is.null(meta) && nrow(meta) == length(grl))
      values(grl) = meta
    return(grl)
}




#' @name anchorlift
#' @title anchorlift
#' @description
#'
#' "lifts" all queries with respect to subject in coordinates that are within "pad"
#' i.e. puts the queries into subject-centric coordinates, which is a new genome with label "Anchor" (default)
#'
#' Respects strand of subject (i.e. if subject strand gr is "-" then will lift all queries to the left of it
#' into positive subject-centric coordinates). Keeps
#' track of subject and query id for later deconvolution if need be.
#'
#' @param query  GRanges that will be lifted around the subject
#' @param subject GRanges around which the queries will be lifted
#' @param window integer specifying how far around each subject to gather query intervals to lift (default = 1e9)
#' @param by character vector specifying additional columms (e.g. sample id) around which to restrict overlaps (via gr.findoverlaps()). Refer to `gr.findoverlaps()` documentation. (default = NULL)
#' @param seqname String specifying the name of the output sequence around which to anchor (default = "Anchor")
#' @param include.values Boolean Flag whether to include values from query and subject (default = TRUE)
#' @return anchorlifted GRanges
#' @author Marcin Imielinski
#' @export
anchorlift = function(query, subject, window = 1e9, by = NULL, seqname = "Anchor", include.values = TRUE, verbose = TRUE, ...)
{               
    if (window > 1e9)
    {
        stop("Window can't be bigger than 1e9")
    }
    
    if (as.numeric(length(query))*as.numeric(length(subject))==0){
        return(NULL)
    }

    if (verbose)
    {
        message('Computing overlaps')
    }
    ov = gr.findoverlaps(query, subject+window, by = by, verbose = verbose, ...)

    if (length(ov) == 0){
        return(NULL)
    }

    if (verbose)
    {
        message('Lifting ranges')
    }
    
    ## new coordinate is query relative to <midpoint> of subject   
    pad = start(subject)[ov$subject.id] + round(width(subject)[ov$subject.id]/2)
    nov = GRanges(seqnames(query)[ov$query.id], IRanges(start(query)[ov$query.id]-pad, end(query)[ov$query.id]-pad),
                  strand = strand(query)[ov$query.id])

    
    #nov = query[ov$query.id] %-% (start(subject)[ov$subject.id] + round(width(subject)[ov$subject.id]/2) + round(width(query)[ov$query.id]/2))
    values(nov) = values(ov)

    if (verbose)
    {
        message('Flipping ranges relative to subject peak')
    }
    flip = ifelse(strand(subject)[ov$subject.id] == '+', 1, -1)

    
    tmp.dt = as.data.table(cbind(start(nov)*flip, end(nov)*flip))[, rowid := 1:.N]
    tmp.dt = melt.data.table(tmp.dt, id.vars = "rowid")
    setkeyv(tmp.dt, c("rowid", "value"))
    tmp.dt[, colid := rep(1:2, length(nov))]
    tmp = matrix(NA, ncol = 2, nrow = length(nov))
    tmp[cbind(tmp.dt$rowid, tmp.dt$colid)] = tmp.dt$value

##   tmp = t(apply(cbind(start(nov)*flip, end(nov)*flip), 1, sort))}) ## 150X slower old version for posterity

    if (verbose)
    {
        message('Creating final output')
    }
    out = GRanges(seqname,  IRanges(tmp[,1], tmp[,2]))
    values(out)$subject.id = ov$subject.id
    values(out)$query.id = ov$query.id
    if (include.values){

        if (ncol(values(query)))
            values(out) = cbind(as.data.table(values(out)), as.data.table(values(query))[ov$query.id,, drop = FALSE])

        if (ncol(values(subject))>0)
            values(out) = cbind(as.data.table(values(out)), as.data.table(values(subject))[ov$subject.id, ,drop = FALSE])
    }

    return(out)
}


#' @name gr.quantile
#' @title gr.quantile
#' @description
#'
#' Computes width weighted quantiles of some GRanges metadata field
#' "field" 
#'
#' @author Marcin Imielinski
#' @import GenomicRanges
#' @param gr GRanges
#' @param qs quantiles ie numeric between 0 and 1
#' larger than 1, both boundary will be considered individual breakpoints
#' @param field metadata field of gr to query (default is the first metadata field)
#' @return numeric vector of length qs corresponding to the width weighted quantiles specified by qs
#' @export
gr.quantile = function(gr, qs, field = values(gr)[1], na.rm = FALSE)
{
  if (any(qs>1 | qs<0))
    stop('qs must be a vector of values between 0 and 1')
  
  if (na.rm)    
    gr = gr[!is.na(values(gr)[[field]])]
  else if (any(is.na(values(gr)[[field]])))
    return(rep(NA, length(qs)))

  x = values(gr)[[field]]
  sw = sum(as.numeric(width(gr))) ## to prevent integer overflow
  width2 = ceiling(width(gr)/sw*1e8) ## round so that sum is <1e8
  ix = order(x)
  this_rle = Rle(x[ix], width2[ix]) ## base by base copy number to get quantiles
  len = length(this_rle)
  res = this_rle[ceiling(qs*len)] %>% as.numeric
  return(res)
}


#' @name gr.breaks
#' @title gr.breaks
#' @description
#'
#' Break GRanges at given breakpoints into disjoint gr
#'
#' @author Xiaotong Yao
#' @import GenomicRanges
#' @param bps \code{GRanges} of width 1, locations of the bp; if any element width
#' larger than 1, both boundary will be considered individual breakpoints
#' @param query a disjoint \code{GRanges} object to be broken
#' @return \code{GRanges} disjoint object at least the same length as query,
#' with a metadata column \code{qid} indicating input index where new segment is from
#' @export
gr.breaks = function(bps=NULL, query=NULL){
   ## ALERT: big change! input parameter shuffled!
   ## if bps not provided, return back-traced disjoin wrapper
   if (is.null(bps)) {
       message("Argument 'bps' not provided")
       return(query)
   } else {
       ## only when bps is given do we care about what query is

       if (is.null(query)){
           ## message("Trying chromosomes 1-22 and X, Y.")
           ## query = hg_seqlengths()
           ## if (is.null(query)){
           ##     message("Default BSgenome not found, let's hardcode it.")
           ##     cs = system.file("extdata",
           ##                      "hg19.regularChr.chrom.sizes", package = "gUtils")
           ##     query = read.delim(cs, header=FALSE, sep="\t")
           ##     query = setNames(sl$V2, sl$V1)
           ## }
           query = gr.stripstrand(si2gr(seqinfo(bps)))
       }

       ## in case query is not a GRanges
       if (!is(query, "GRanges")){
           stop("Error: 'query' must be a GRanges object.")
       }

       ## ## preprocess query
       ## if (!isDisjoint(query)){
       ##     warning("Warning: Query GRanges not disjoint.")
       ##     queryDj = disjoin(query)
       ##     queryDj$qid = queryDj %N% query ## only retain the first occurence
       ##     values(queryDj) = cbind(values(queryDj),
       ##                             as.data.table(values(query))[queryDj$qid])
       ##     query = queryDj
       ## } else {
       ##     if ("qid" %in% colnames(values(query))){
       ##         warning("Warning: 'qid' col in query overwritten.")
       ##     }
       query$qid = seq_along(query)
       ## }

       ## preprocess bps
       ## having meta fields? remove them!
       bps = bps[, c()]

       ## remove things outside of ref
       oo.seqlength = which(start(bps)<1 | end(bps)>seqlengths(bps)[as.character(seqnames(bps))])
       
       if (length(oo.seqlength)>0){
           warning("Warning: Some breakpoints out of chr lengths. Removing.")
           bps = bps[-oo.seqlength]
       }

       if (any(!is.null(names(bps)))){
           warning("Warning: Removing row names from bps.")
           names(bps) = NULL
       }

       ## having strand info? remove it!
       if (any(strand(bps)!="*")){
           warning("Warning: Some breakpoints have strand info. Force to '*'.")
           bps = gr.stripstrand(bps)
       }

       ## solve three edge cases
       if (any(w.0 <- (width(bps)<1))){
           warning("Warning: Some breakpoint width==0. Discard.")
           ## right bound smaller coor
           ## and there's no negative width GR allowed
           ## bps[which(w.0)] = gr.start(bps[which(w.0)]) %-% 1
           bps = bps[-which(w.0)]
       }

       if (any(w.2 <- (width(bps)==2))){
           warning("Warning: Some breakpoint width>2. Will tear them apart and treat as two breakpoints.")
           ## this is seen as breakpoint by spanning two bases
           bps = c(bps[-which(w.2)],
                   gr.start(bps[which(w.2)]),
                   gr.end(bps[which(w.2)]))
       }

       if (any(w.l <- (width(bps)>2))){
           ## some not a point? turn it into a point
           warning("Warning: Some breakpoint width>2. Treat them as segmentations.")
           rbps = gr.end(bps[which(w.l)])
           lbps = gr.start(bps[which(w.l)])
           start(lbps) = pmax(start(lbps)-1, 1)
           bps = c(bps[which(!w.l)], streduce(c(lbps, rbps)))
       }

       bps$inQuery = bps %^% query
       if (any(bps$inQuery==FALSE)){
           warning("Warning: Some breakpoint not within query ranges.")
       }

       ## label and only consider breakpoints not already at the boundary of query
       bps$inner = bps$inQuery ## out of query automatically FALSE
       bps$inner[which(bps %^% gr.start(query) | bps %^% gr.end(query))]=FALSE
       
       ## maybe no inner bp at all, then no need to proceed
       if (!any(bps$inner)){
           return(query)
       }
       bpsInner = bps %Q% (inner==T)

       ## map query and inner breakpoints
       qbMap = gr.findoverlaps(query, bpsInner)
       mappedQ = seq_along(query) %in% qbMap$query.id
       ## raw coors to construct ranges from
       tmpRange = data.table(qid2 = qbMap$query.id,
                             startFrom = start(query[qbMap$query.id]),
                             breakAt = start(bpsInner[qbMap$subject.id]),
                             upTo = end(query[qbMap$query.id]))
       tmpCoor = tmpRange[, .(pos=sort(unique(c(startFrom, breakAt, upTo)))), by=qid2]

       ## construct new ranges
       newRange = tmpCoor[, .(tmp.start=pos[-which.max(pos)],
                              end=pos[-which.min(pos)]),
                          by=qid2]
       newRange[, ":="(seqnames = as.vector(seqnames(query)[qid2]),
                       strand = as.vector(strand(query)[qid2]))]
       newRange[, start := ifelse(tmp.start==min(tmp.start), tmp.start, tmp.start+1), by=qid2]

       ## put together the mapped and broken
       newGr = GRanges(newRange, seqinfo = seqinfo(query))
       values(newGr) = values(query)[newGr$qid2, , drop=F] ## preserve the input metacol
       ## with the intact not mapped part of query
       output = sort(c(newGr, query[!mappedQ]))
       ## %Q% (order(strand, seqnames, start))
       ## browser()
       return(output)
   }
}



#' @name ra.duplicated
#' @title ra.duplicated
#' @description
#'
#' Show if junctions are Deduplicated
#'
#' Determines overlaps between two or more piles of rearrangement junctions (as named or numbered arguments) +/- padding
#' and will merge those that overlap into single junctions in the output, and then keep track for each output junction which
#' of the input junctions it was "seen in" using logical flag  meta data fields prefixed by "seen.by." and then the argument name
#' (or "seen.by.ra" and the argument number)
#'
#' @author Xiaotong Yao
#' @param grl GRangesList representing rearrangements to be merged
#' @param pad non-negative integer specifying padding
#' @param ignore.strand whether to ignore strand (implies all strand information will be ignored, use at your own risk)
#' @return \code{GRangesList} of merged junctions with meta data fields specifying which of the inputs each outputted junction was "seen.by"
#' @name ra.duplicated
#' @examples
#'
#' @export
ra.duplicated = function(grl, pad=500, ignore.strand=FALSE){

   if (!is(grl, "GRangesList")){
       stop("Error: Input must be GRangesList!")
   }

   ##if (any(elementNROWS(grl)!=2)){
   ##    stop("Error: Each element must be length 2!")
   ##}

   if (length(grl)==0){
       return(logical(0))
   }

   if (length(grl)==1){
       return(FALSE)
   }

   if (length(grl)>1){

       ix.pair = as.data.table(ra.overlaps(grl, grl, pad=pad, ignore.strand = ignore.strand))[ra1.ix!=ra2.ix]

       if (nrow(ix.pair)==0){
           return(rep(FALSE, length(grl)))
       }
       else {
         ##           dup.ix = unique(rowMax(as.matrix(ix.pair)))
         dup.ix = unique(apply(as.matrix(ix.pair), 1, max))
         return(seq_along(grl) %in% dup.ix)
       }
   }
}



#' @name ra.overlaps
#' @title ra.overlaps
#' @description
#'
#' Determines overlaps between two piles of rearrangement junctions ra1 and ra2 (each GRangesLists of signed locus pairs)
#' against each other, returning a sparseMatrix that is T at entry ij if junction i overlaps junction j.
#'
#' if argument pad = 0 (default) then only perfect overlap will validate, otherwise if pad>0 is given, then
#' padded overlap is allowed
#'
#' strand matters, though we test overlap of both ra1[i] vs ra2[j] and gr.flipstrand(ra2[j])
#'
#' @param ra1 \code{GRangesList} with rearrangement set 1
#' @param ra2 \code{GRangesList} with rearrangement set 2
#' @param pad Amount to pad the overlaps by. Larger is more permissive. Default is exact (0)
#' @param arr.ind Default TRUE
#' @param ignore.strand Ignore rearrangement orientation when doing overlaps. Default FALSE
#' @param ... params to be sent to \code{\link{gr.findoverlaps}}
#' @name ra.overlaps
#' @export
ra.overlaps = function(ra1, ra2, pad = 0, arr.ind = TRUE, ignore.strand=FALSE, ...)
{
    bp1 = grl.unlist(ra1) + pad
    bp2 = grl.unlist(ra2) + pad
    ix = gr.findoverlaps(bp1, bp2, ignore.strand = ignore.strand, ...)

    .make_matches = function(ix, bp1, bp2)
    {
        if (length(ix) == 0){
            return(NULL)
        }
        tmp.match = cbind(bp1$grl.ix[ix$query.id], bp1$grl.iix[ix$query.id], bp2$grl.ix[ix$subject.id], bp2$grl.iix[ix$subject.id])
        tmp.match.l = lapply(split(1:nrow(tmp.match), paste(tmp.match[,1], tmp.match[,3])), function(x) tmp.match[x, , drop = F])

        ## match only occurs if each range in a ra1 junction matches a different range in the ra2 junction
        matched.l = sapply(tmp.match.l, function(x) all(c('11','22') %in% paste(x[,2], x[,4], sep = '')) | all(c('12','21') %in% paste(x[,2], x[,4], sep = '')))

        return(do.call('rbind', lapply(tmp.match.l[matched.l], function(x) cbind(x[,1], x[,3])[!duplicated(paste(x[,1], x[,3])), , drop = F])))
    }

    tmp = .make_matches(ix, bp1, bp2)

    if (is.null(tmp)){
        if (arr.ind){

            return(as.matrix(data.table(ra1.ix = as.numeric(NA), ra2.ix = as.numeric(NA))))
        }
        else{
            return(Matrix::sparseMatrix(length(ra1), length(ra2), x = 0))
        }
    }

    rownames(tmp) = NULL

    colnames(tmp) = c('ra1.ix', 'ra2.ix')

    if (arr.ind) {
        ro = tmp[order(tmp[,1], tmp[,2]), , drop = FALSE]
        if (inherits(ro, "integer")) {
            ro <- matrix(ro, ncol=2, nrow=1, dimnames=list(c(), c('ra1.ix', 'ra2.ix')))
        }
        return(ro)
    } else {
        ro = Matrix::sparseMatrix(tmp[,1], tmp[,2], x = 1, dims = c(length(ra1), length(ra2)))
        return(ro)
    }
}



#' @name grl.expand
#' @title grl.expand
#' @description
#'
#' Function wrapping around the `+` operator
#' for GRanges objects to work on GRangesLists.
#' Expands window of element GRanges within GRangesList
#'
#' @param grl \code{GRangesList} 
#' @param expand_win \code{integer}
#' @return GRangesList with added window
#' @author Kevin Hadi
#' @export
grl.expand = function(grl, expand_win) {
    tmp_vals = mcols(grl)
    tmp_gr = unlist(grl, use.names = FALSE)
    tmp_gr = tmp_gr + expand_win
    new_grl = relist(tmp_gr, grl)
    mcols(new_grl) = tmp_vals
    return(new_grl)
}

#' @name grl.shrink
#' @title grl.shrink
#' @description
#'
#' Function wrapping around the `-` operator
#' for GRanges objects to work on GRangesLists.
#' Shrnks window of element GRanges within GRangesList
#'
#' @param grl \code{GRangesList} 
#' @param shrink_win \code{integer}
#' @return GRangesList with shrunken windows
#' @author Kevin Hadi
#' @export
grl.shrink = function(grl, shrink_win) {
    tmp_vals = mcols(grl)
    tmp_gr = unlist(grl, use.names = FALSE)
    tmp_gr = tmp_gr - shrink_win
    new_grl = relist(tmp_gr, grl)
    mcols(new_grl) = tmp_vals
    return(new_grl)
}



#' @name +
#' @title adding functionality to GRangesLists
#' @description
#'
#' Adding operator to GRangesList 
#'
#' @return extended grangeslist
#' @rdname grl.expand shortcut
#' @exportMethod +
#' @aliases +, GRangesList method 
#' @author Kevin Hadi
#' @export
setMethod(`+`, 'GRangesList', function(e1, e2) {
    return(grl.expand(e1, e2))
})


#' @name -
#' @title adding functionality to GRangesLists
#' @description
#'
#' Adding operator to GRangesList 
#'
#' @return extended grangeslist
#' @rdname grl.shrink shortcut
#' @exportMethod -
#' @aliases -, GRangesList method 
#' @author Kevin Hadi
#' @export
setMethod(`-`, 'GRangesList', function(e1, e2) {
    return(grl.shrink(e1, e2))
})

#' @name grl.start
#' @title Get the left ends of a \code{GRangesList}
#' @description
#'
#' Function wrapping around gr.start to work
#' on GRangesLists
#'
#' @param grl \code{GRangesList} object containing GRanges to operate on 
#' @param width integer Specify subranges of greater width including the start of the range. (default = 1)
#' @param force boolean Allows returned \code{GRangesList} to have ranges outside of its \code{Seqinfo} bounds. (default = FALSE)
#' @param ignore.strand boolean If set to \code{FALSE}, will extend '-' strands from the other direction (default = TRUE)
#' @param clip boolean Trims returned \code{GRangesList} so that it does not extend beyond bounds of the input \code{GRanges} (default = TRUE)
#' @return \code{GRangesList} object of width 1 ranges representing start of each genomic range in the input.
#' @importFrom GenomicRanges GRangesList
#' @return \code{GRangesList} object of width 1+ ranges representing start of each genomic range in the input.
#' @export
grl.start = function(grl, width = 1, force = FALSE, ignore.strand = TRUE, clip = TRUE) {
    tmp_vals = mcols(grl)
    tmp_gr = unlist(grl, use.names = FALSE)
    tmp_gr = gr.start(tmp_gr, width, force, ignore.strand, clip)
    new_grl = relist(tmp_gr, grl)
    mcols(new_grl) = tmp_vals
    return(new_grl)
}


#' @name grl.end
#' @title Get the right ends of a \code{GRangesList}
#' @description
#'
#' Function wrapping around gr.end to work on
#' GRangesLists
#' 
#' @param x \code{GRangesList} object containing GRanges to operate on
#' @param width integer Specify subranges of greater width including the start of the range. (default = 1)
#' @param force boolean Allows returned \code{GRangesList} to have ranges outside of its \code{Seqinfo} bounds. (default = FALSE)
#' @param ignore.strand boolean If set to \code{FALSE}, will extend '-' strands from the other direction. (default = TRUE)
#' @param clip boolean Trims returned \code{GRangesList} so that it does not extend beyond bounds of the input (default = TRUE)
#' @return GRangesList object of width = \code{width} ranges representing end of each genomic range in the input.
#' @importFrom GenomeInfoDb seqlengths
#' @importFrom GenomicRanges strand seqnames values<- values
#' @author Kevin Hadi
#' @export
grl.end = function(grl, width = 1, force = FALSE, ignore.strand = TRUE, clip = TRUE) {
    tmp_vals = mcols(grl)
    tmp_gr = unlist(grl, use.names = FALSE)
    tmp_gr = gr.end(tmp_gr, width, force, ignore.strand, clilp)
    new_grl = relist(tmp_gr, grl)
    mcols(new_grl) = tmp_vals
    return(new_grl)
}


#' @name rebin
#' @title reaggregate WGS bins around a new target value
#' @description 
#' 
#' Given GRanges of bins will aggregate around new bin width. 
#' 
#' @param cov GRanges of binned genome-wide coverage
#' @param binwidth new binwidth
#' @return GRanges of binned genome-wide coverage at new bin
#' @export
#' @author Marcin Imielinski
rebin = function(cov, binwidth, field = names(values(cov))[1], FUN = mean, na.rm = TRUE)
{
  tmp = as.data.table(cov[, c()])
  tmp$value =  values(cov)[[field]]
  outdt = tmp[
    , FUN(value, na.rm = na.rm),
      by = .(seqnames, start = floor(start/binwidth)*binwidth + 1)]
  ## xtYao update:
  ## by = .(seqnames, start = ceiling(start/binwidth)*binwidth)]
  outdt[, end := start + binwidth -1]
  out = dt2gr(outdt)
  names(values(out)) = field
  return(out)
}


