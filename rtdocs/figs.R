setwd(file.path(Sys.getenv("GIT_HOME"), "gUtils", "rtdocs"))
HI=100
WI=450

HI=2
WI=9
## make the granges
gr <- GRanges(1, IRanges(c(2,5,10), c(4,9,16)), seqinfo=Seqinfo("1", 20))
gr2 <- GRanges(1, IRanges(c(1,9), c(6,14)), seqinfo=Seqinfo("1", 20))
dt <- data.table(seqnames=1, start=c(2,5,10), end=c(3,8,15))

## gr.start
g <- ggplot() + geom_rect(data=as.data.frame(gr), aes(xmin=start, xmax=end, ymin=1, ymax=2), fill='blue') + geom_rect(data=as.data.frame(gr.start(gr, width=3)), aes(xmin=start, xmax=end, ymin=0, ymax=0.9), fill='red') + theme_bw() + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) + scale_x_continuous(breaks=1:20, labels=1:20) + coord_fixed()
pdf("figures/gr.start.pdf", height=HI, width=WI); print(g); dev.off()
system.call("convert figures/gr.start.pdf figures/gr.start.png")

## gr.end
g <- ggplot() + geom_rect(data=as.data.frame(gr), aes(xmin=start, xmax=end, ymin=1, ymax=2), fill='blue') + geom_rect(data=as.data.frame(gr.end(gr, width=3)), aes(xmin=start, xmax=end, ymin=0, ymax=0.9), fill='red') + theme_bw() + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) + scale_x_continuous(breaks=1:20, labels=1:20) + coord_fixed()
pdf("figures/gr.end.pdf", height=HI, width=WI); print(g); dev.off()
system.call("convert figures/gr.end.pdf figures/gr.end.png")

## gr.mid
g <- ggplot() + geom_rect(data=as.data.frame(gr), aes(xmin=start, xmax=end, ymin=1, ymax=2), fill='blue') + geom_rect(data=as.data.frame(gr.mid(gr)+1), aes(xmin=start, xmax=end, ymin=0, ymax=0.9), fill='red') + theme_bw() + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) + scale_x_continuous(breaks=0:20, labels=0:20) + coord_fixed()
pdf("figures/gr.mid.pdf", height=HI, width=WI); print(g); dev.off()
system.call("convert figures/gr.mid.pdf figures/gr.mid.png")

## shift
g <- ggplot() + geom_rect(data=as.data.frame(gr), aes(xmin=start, xmax=end, ymin=1, ymax=2), fill='blue') + geom_rect(data=as.data.frame(shift(gr, 2)), aes(xmin=start, xmax=end, ymin=0, ymax=0.9), fill='red') + theme_bw() + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) + scale_x_continuous(breaks=0:20, labels=0:20) + coord_fixed()
pdf("figures/shift.pdf", height=HI, width=WI); print(g); dev.off()
system.call("convert figures/shift.pdf figures/shift.png")

## flank
g <- ggplot() + geom_rect(data=as.data.frame(gr), aes(xmin=start, xmax=end, ymin=1, ymax=2), fill='blue') + geom_rect(data=as.data.frame(flank(gr,width=2)), aes(xmin=start, xmax=end, ymin=0, ymax=0.9), fill='red') + theme_bw() + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) + scale_x_continuous(breaks=0:20, labels=0:20) + coord_fixed()
pdf("figures/flank.pdf", height=HI, width=WI); print(g); dev.off()
system.call("convert figures/flank.pdf figures/flank.png")

## streduce
g <- ggplot() +
  geom_rect(data=as.data.frame(gr), aes(xmin=start, xmax=end, ymin=2.05, ymax=3), fill='blue') +
  geom_rect(data=as.data.frame(gr2),aes(xmin=start, xmax=end, ymin=1, ymax=1.95), fill='blue') +
  geom_rect(data=as.data.frame(streduce(c(gr, gr2))), aes(xmin=start, xmax=end, ymin=0, ymax=0.95), fill='red') +
  theme_bw() + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) + scale_x_continuous(breaks=0:20, labels=0:20) + coord_fixed()
pdf("figures/streduce.pdf", height=HI, width=WI); print(g); dev.off()
system.call("convert figures/streduce.pdf figures/streduce.png")

## gr.sample
set.seed(137)
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


