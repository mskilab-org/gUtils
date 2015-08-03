setwd(file.path(Sys.getenv("GIT_HOME"), "gUtils", "rtdocs"))
HI=100
WI=450

HI=2
WI=9
## make the granges
gr <- GRanges(1, IRanges(c(2,5,10), c(4,9,16)), seqinfo=Seqinfo("1", 20))
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

## flank start=FALSE
g <- ggplot() + geom_rect(data=as.data.frame(gr), aes(xmin=start, xmax=end, ymin=1, ymax=2), fill='blue') + geom_rect(data=as.data.frame(flank(gr,width=2, start=FALSE)), aes(xmin=start, xmax=end, ymin=0, ymax=0.9), fill='red') + theme_bw() + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) + scale_x_continuous(breaks=0:20, labels=0:20) + coord_fixed()

## gr.round
g <- ggplot() + geom_rect(data=as.data.frame(gr), aes(xmin=start, xmax=end, ymin=1, ymax=2), fill='blue') + geom_rect(data=as.data.frame(gr.round(gr, gr)), aes(xmin=start, xmax=end, ymin=0, ymax=0.9), fill='red') + theme_bw() + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) + scale_x_continuous(breaks=0:20, labels=0:20) + coord_fixed()


