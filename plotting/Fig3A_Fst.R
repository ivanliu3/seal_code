
### 4_A, Fst matrix ###
t1 <- read.csv("~/Downloads/Fst/collect.world.fst.tsv",sep="\t",header=F)[-1]

m1<-as.matrix(t1[-1,])
pop <- sapply(t1[1,],as.character)
names(pop) <- NULL

m1[which(m1=="Not found")] <- "NA"
m1

m2<-as.matrix( apply(m1,2,as.numeric) )
colnames(m2) <-pop
rownames(m2) <- pop
m2


t2 <- read.table("~/Downloads/Fst/fst.info.rmPlateandRelate.csv",sep=",",header=F)
t2
pop.ord<-sapply(t2$V2, function(x) which(pop ==x))

m3<-m2[pop.ord,pop.ord]
colnames(m3) <- pop[pop.ord]
rownames(m3) <- pop[pop.ord]

m4 <- m3
m5 <- m3
m4[lower.tri(m4)] <- NA
m5[upper.tri(m5)] <- NA
library(scico)
library(gplots)
library("RColorBrewer")
pal <-scico(100, palette = 'hawaii')
colo <- as.integer(t2$V3)


par(mar=c(5.1,2.1,4.1,2.1))


heatmap.2(m5,dendrogram='none',Rowv=F,
          Colv = F,col=pal,tracecol=NA ,offsetRow =-50,colRow = 'black',colCol = 'black',key.title = "",lwid = c(1,2),
          keysize=1,cexRow = 1.5,cexCol = 1.5) 
## Fig size 900 X 787
#locator(1)
#text(0.02,0.85,"P_NE",col=color [5],cex=1.5)
#text(0.02,0.77,"P_NW",col=color[4],cex=1.5)
#text(0.02,0.627,"A_NW",col=color[3],cex=1.5)
#text(0.02,0.48,"A_NC",col=color[1],cex=1.5)
#text(0.02, 0.194,"A_NE",col=color[2],cex=1.5)

