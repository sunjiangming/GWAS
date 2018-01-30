#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
ppiFile=args[1]
#interP bait pi pb
genomewideline=args[2]
suggestiveline=args[3]

ppi=read.table(ppiFile,header=T,sep="\t")
interP=as.character(ppi$interP)
bait=as.character(ppi$bait)

uniqInterP <- ppi[!duplicated(ppi[, 1]),c(1,3)]
colnames(uniqInterP)=c("bait","pb")
uniqBait <- ppi[!duplicated(ppi[, 2] ),c(2,4)]
p=rbind(uniqInterP,uniqBait)[,2]
logP=-log10(p)
minX=floor(min(logP))
maxX=ceiling(max(logP))

### make stat on x span, e.g. number of proteins with -logP range at 5-6
br = seq(minX,maxX,1)
freq=hist(logP, breaks=br, include.lowest=TRUE, plot=FALSE)
spanLenBait = max(max(freq$counts)/(length(uniqBait[,1]) + 1), 1)
uniqBait$x=-log10(uniqBait$pb)
uniqBait$y=NULL
for( i in length(uniqBait[,1]):1){
  uniqBait$y[i]=2*(length(uniqBait[,1])-i+1)*spanLenBait
}

ppi$y2=NULL
for( i in 1:dim(ppi)[1]){
  ppi$y2[i] = uniqBait$y[ ppi$bait[i]==uniqBait[,1] ]
}
dymean = aggregate(ppi$y2, list(ppi$interP), mean)
dymax = aggregate(ppi$y2, list(ppi$interP), max)
###--------------------------------------------------
y1=dymean[,2]
t=1:dim(dymax)[1]
y1=sapply(t, function(i) { ifelse( any(dymean[i,2] == dymax[i,2]), runif(1,dymean[i,2]-1*spanLenBait, dymean[i,2]+3*spanLenBait),dymean[i,2]+runif(1,-1,4)) } )

dymax$x = -log10(aggregate(ppi$pi, list(ppi$interP), mean)[,2])
dymax$y=y1

ppi$x1 = -log10(ppi$pi)
ppi$x2 = -log10(ppi$pb)

ppi$y1=NULL
for( i in 1:dim(ppi)[1]){
  ppi$y1[i] = y1[ ppi$interP[i]==dymean[,1] ]
}

library(ggplot2)
library(ggrepel)

# Simple version.
p1 = ggplot(ppi,aes(x=x1, xend=x2, y=y1, yend=y2,colour=bait)) +
     geom_segment(size=0.2) +
     geom_point(aes(colour=ppi$bait)) + 
     geom_text(data=ppi,aes(label=interP, x=x1, y=y1-0.25),size=3) + 
     geom_text(data=ppi,aes(label=bait, x=x2, y=y2-0.25)) + 
     scale_x_continuous(breaks=c(5:25), minor_breaks=NULL) +
     labs(x="-log10P",y="") +
     theme(legend.position="none")


geom_text(data=uniqBait,aes(label=uniqBait$bait, x=uniqBait$x, y=uniqBait$y,colour=uniqBait$bait)) +

p1 + geom_text_repel(data=uniqBait,aes(label=bait, x=x, y=y)) +
p1 + geom_text(data=dymax,aes(label=dymax[,1], x=dymax$x, y=dymax$y))

ggsave(plot=p1, filename="plot_1.png", height=3.5, width=6)
ggsave(plot=p1, filename="plot_2.png")

