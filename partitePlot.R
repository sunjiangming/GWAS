#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
ppiFile=args[1]
#interP bait pi pb
genomewideline=args[2]
suggestiveline=args[3]

ppi=read.table(ppiFile,header=T,sep="\t")
interP=as.character(ppi$interP)
bait=as.character(ppi$bait)

uniqInterP <- ppi[!duplicated(ppi[ , 1]),c(1,3)]
colnames(uniqInterP)=c("bait","pb")
uniqBait <- ppi[!duplicated(ppi[, 2] ),c(2,4)]
p=rbind(uniqInterP,uniqBait)[,2]
logP=-log10(p)
minX=floor(min(logP))
maxX=ceiling(max(logP))

br = seq(minX,maxX,1)
freq=hist(logP, breaks=br, include.lowest=TRUE, plot=FALSE)
spanLenBait = max(max(freq$counts)/(length(uniqBait[,1]) + 1), 1)
uniqBait$y=NULL
for( i in length(uniqBait[,1]):1){
  uniqBait$y[i]=(length(uniqBait[,1])-i+1)*spanLenBait
}

ppi$yp=NULL
for( i in 1:dim(ppi)[1]){
  ppi$yp[i] = uniqBait$y[ ppi$bait[i]==uniqBait[,1] ]
}
dymean = aggregate(ppi$yp, list(ppi$interP), mean)
dymax = aggregate(ppi$yp, list(ppi$interP), max)
y1=dymean
xy1=sapply(1:dim(dymax[1]), function(i) { if ( dymean[i,2] == dymax[i,2]) y1[i,2]=runif(1,dymean[i,2]-spanLenBait/2, dymean[i,2]+spanLenBait/2) } )


#maxY=freq$counts
x1 = log10(uniqInterP[,2])
x2 = log10(uniqBait[,2])
y2 = uniqBait$y

library(ggplot2)
# Simple version.
p1 = ggplot(dat, aes(x=x1, xend=x2, y=y1, yend=y2)) +
     geom_segment(size=1.2)

ggsave(plot=p1, filename="plot_1.png", height=3.5, width=6)


