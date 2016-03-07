setwd(file.path(Sys.getenv("GIT_HOME"), "gUtils", "rtdocs"))
HI=150
WI=600

.plot_gr <- function(gro, grn) {

  gro_p <- gro[strand(gro) %in% c("+","*")]
  gro_m <- gro[strand(gro) %in% c("-")]
  grn_p <- grn[strand(grn) %in% c("+","*")]
  grn_m <- grn[strand(grn) %in% c("-")]

  alpha = 0.5
  g <- ggplot()
  if (length(gro_p))
    g <- g + geom_rect(data=as.data.frame(gro_p), aes(xmin=start, xmax=end+1, ymin=4, ymax=4.9), fill='blue', color='black', alpha=alpha)
  if (length(gro_m))
    g <- g + geom_rect(data=as.data.frame(gro_m), aes(xmin=start, xmax=end+1, ymin=3, ymax=3.9), fill='blue', color='black', alpha=alpha)

  if (length(grn_p))
    g <- g + geom_rect(data=as.data.frame(grn_p), aes(xmin=start, xmax=end+1, ymin=1, ymax=1.9), fill='red', color='black', alpha=alpha)
  if (length(grn_m))
    g <- g + geom_rectct(data=as.data.frame(grn_m), aes(xmin=start, xmax=end+1, ymin=0, ymax=0.9), fill='red', color='black', alpha=alpha)

  g <- g + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_line(colour="gray65", size=0.5, linetype = 'dashed'), axis.ticks.x=element_blank())
  g <- g + scale_x_continuous(breaks=seq(from=1.5, to=20.5, by=1), minor_breaks=seq(from=1,to=20, by=1), labels=1:20) + coord_fixed()
  g <- g + scale_y_continuous(breaks = NULL)
  g <- g + geom_text(data=data.frame(x=c(0,0,0,0), y=c(3.5, 4.5, 0.5, 1.5), text=c("gr(-)", "gr(+)", "(-)","(+)")), aes(x=x,y=y,label=text)) + xlab("") + ylab("")
  g <- g + geom_line(data=data.frame(x=c(1,1), y=c(0,5)), aes(x=x,y=y), color='black')

  .lab_start <- function(grr) { return( (end(grr) - start(grr))/2 + start(grr) + 0.5)}
  ## add labels to boxes
  xs <- c(.lab_start(gro_p), .lab_start(gro_m), .lab_start(grn_p), .lab_start(grn_m))
  ys <- c(rep(4.5, length(gro_p)), rep(3.5, length(gro_m)), rep(1.5, length(grn_p)), rep(0.5, length(grn_m)))
  tx <- c(gro$name, grn$name)

  g <- g + geom_text(data=data.frame(x=xs, y=ys, text=tx), aes(x=x, y=y, label=text))
  return(g)
}

#HI=2
#WI=9
## make the granges
gr <- GRanges(1, IRanges(c(3,7,13), c(5,9,16)), strand=c('+','-','-'), seqinfo=Seqinfo("1", 25), name=c("A","B","C"))
gr2 <- GRanges(1, IRanges(c(1,9), c(6,14)), seqinfo=Seqinfo("1", 20))
dt <- data.table(seqnames=1, start=c(2,5,10), end=c(3,8,15))

## gr.start
png("rtdocs/figures/gr.start.png", height=HI, width=WI);
print(.plot_gr(gr, gr.start(gr, width=2)) + ggtitle("gr.start(gr, width=2)"));
dev.off()

## gr.start with strand
png("rtdocs/figures/gr.start_wstrand.png", height=HI, width=WI);
print(.plot_gr(gr, gr.start(gr, width=2, ignore.strand=FALSE)) + ggtitle("gr.start(gr, width=2, ignore.strand=FALSE)"));
dev.off()

## gr.end
png("rtdocs/figures/gr.end.png", height=HI, width=WI);
print(.plot_gr(gr, gr.end(gr, width=2)) + ggtitle("gr.end(gr, width=2)"));
dev.off()

## gr.end
png("rtdocs/figures/gr.end_wstrand.png", height=HI, width=WI);
print(.plot_gr(gr, gr.end(gr, width=2, ignore.strand=FALSE)) + ggtitle("gr.end(gr, width=2, ignore.strand=FALSE)"));
dev.off()

## gr.mid
png("rtdocs/figures/gr.mid.png", height=HI, width=WI);
print(.plot_gr(gr, gr.mid(gr)) + ggtitle("gr.mid(gr)"));
dev.off()

## shift
png("rtdocs/figures/shift.png", height=HI, width=WI);
print(.plot_gr(gr, GenomicRanges::shift(gr,2)) + ggtitle("GenomicRanges::shift(gr,2)"));
dev.off()

## flank
png("rtdocs/figures/flank.png", height=HI, width=WI);
print(.plot_gr(gr, GenomicRanges::flank(gr,width=1)) + ggtitle("GenomicRanges::flank(gr,width=1)"));
dev.off()

## flank w/start
png("rtdocs/figures/flank_start.png", height=HI, width=WI);
print(.plot_gr(gr, GenomicRanges::flank(gr,start=FALSE, width=1)) + ggtitle("GenomicRanges::flank(gr,start=FALSE, width=1)"));
dev.off()

## streduce
png("rtdocs/figures/streduce.png", height=HI, width=WI)
grn = streduce(gr+2)
grn$name <- c("ABC")
print(.plot_gr(gr+2, grn) + ggtitle("streduce(gr+2)"));
dev.off()

## gr.sample
set.seed(137)
png("rtdocs/figures/gr.sample.png", height=HI, width=WI)
print(.plot_gr(gr, gr.sample(gr, 3, len=2, replace=TRUE)) + ggtitle("streduce(gr+2)"));
dev.off()

g <- ggplot() +
  geom_rect(data=as.data.frame(gr), aes(xmin=start, xmax=end, ymin=1, ymax=2), fill='blue') +
  geom_rect(data=as.data.frame(gr.sample(gr, 3, len=2, replace=TRUE)), aes(xmin=start, xmax=end, ymin=0, ymax=0.95), fill='red') +
  theme_bw() + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) + scale_x_continuous(breaks=0:20, labels=0:20) + coord_fixed()
pdf("figures/gr.sample.pdf", height=HI, width=WI); print(g); dev.off()
system.call("convert figures/gr.sample.pdf figures/gr.sample.png")

## gr.rand
set.seed(137)
g <- ggplot() +
  geom_rect(data=as.data.frame(gr), aes(xmin=start, xmax=end, ymin=1, ymax=2), fill='blue') +
  geom_rect(data=as.data.frame(gr.rand(c(2,5,3), seqinfo(gr))), aes(xmin=start, xmax=end, ymin=0, ymax=0.95), fill='red') +
  theme_bw() + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) + scale_x_continuous(breaks=0:20, labels=0:20) + coord_fixed()
pdf("figures/gr.rand.pdf", height=HI, width=WI); print(g); dev.off()
system.call("convert figures/gr.rand.pdf figures/gr.rand.png")

## gr.simplify
gr3 <- GRanges(1, IRanges(c(1,4,8), c(4,8,10)))
df <- data.frame(start=c(1,4,8), end=c(3.95, 7.96, 10))
g <- ggplot() +
  geom_rect(data=df, aes(xmin=start, xmax=end, ymin=1, ymax=2), fill='blue') +
  geom_rect(data=as.data.frame(gr.simplify(gr3)), aes(xmin=start, xmax=end, ymin=0, ymax=0.95), fill='red') +
  theme_bw() + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) + scale_x_continuous(breaks=0:20, labels=0:20) + coord_fixed()
pdf("figures/gr.simplify.pdf", height=HI, width=WI); print(g); dev.off()
system.call("convert figures/gr.simplify.pdf figures/gr.simplify.png")

## gr.tile
df1 <- data.frame(xmin=c(0,3,6), xmax=c(4,7,10), ymin=c(0,1,2), ymax=c(1,2,3))
df2 <- data.frame(xmin=1, xmax=9, ymin=3, ymax=4)
g <- ggplot() +
  geom_rect(data=df1, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill='blue') +
  geom_rect(data=df2, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill='red') +
  theme_bw() + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) + scale_x_continuous(breaks=0:20, labels=0:20) + coord_fixed()
pdf("figures/gr.tile.pdf", height=HI, width=WI); print(g); dev.off()
system.call("convert figures/gr.tile.pdf figures/gr.tile.png")

## gr.refactor
gr3 <- GRanges(1, IRanges(c(4,11), c(6, 20)))
g <- ggplot() + geom_rect(data=as.data.frame(gr3), aes(xmin=start, xmax=end, ymin=1, ymax=2), fill='blue') + geom_rect(data=as.data.frame(gr.refactor(gr3, "3")), aes(xmin=start, xmax=end, ymin=0, ymax=0.9), fill='red') + theme_bw() + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) + scale_x_continuous(breaks=0:20, labels=0:20) + coord_fixed()
pdf("figures/gr.refactor.pdf", height=HI, width=WI); print(g); dev.off()
system.call("convert figures/gr.refactor.pdf figures/gr.refactor.png")

## flank start=FALSE
g <- ggplot() + geom_rect(data=as.data.frame(gr), aes(xmin=start, xmax=end, ymin=1, ymax=2), fill='blue') + geom_rect(data=as.data.frame(flank(gr,width=2, start=FALSE)), aes(xmin=start, xmax=end, ymin=0, ymax=0.9), fill='red') + theme_bw() + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) + scale_x_continuous(breaks=0:20, labels=0:20) + coord_fixed()

## gr.round
g <- ggplot() + geom_rect(data=as.data.frame(gr), aes(xmin=start, xmax=end, ymin=1, ymax=2), fill='blue') + geom_rect(data=as.data.frame(gr.round(gr, gr)), aes(xmin=start, xmax=end, ymin=0, ymax=0.9), fill='red') + theme_bw() + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) + scale_x_continuous(breaks=0:20, labels=0:20) + coord_fixed()


