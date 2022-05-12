library(RColorBrewer) 
palette(c(brewer.pal(12, "Paired"),"#1B9E77","#999999","#8DD3C7"))

colors.ori <- c(brewer.pal(12, "Paired"),"#1B9E77","#999999","#8DD3C7")
colors = colors.ori[c(6,4,1,2,3,5,7:14)]
palette(colors)
##### functions from admixFun.R

getRMSD<-function(Qest,Qtrue){
  require(gtools)
  
  npop<-nrow(Qtrue)
  perms<-permutations(npop,npop)
  Nperm<-dim(perms)[1]
  
  theMin<-c()
  for(g in 1:nrow(Qest)){
    cat("N perm",Nperm,"\n")
    RMSEs<-rep(0,Nperm)
    Qnew<-Qest[-g,]
    for(i in 1:Nperm){
      RMSEs[i]<-sqrt(sum((Qtrue-Qnew[perms[i,],])^2))
    }
    theMin[g]<-min(RMSEs)
  }
  
  w<-which.min(theMin)
  Qnew<-Qest[-w,]
  
  RMSEs<-rep(0,Nperm)
  for(i in 1:Nperm){
    RMSEs[i]<-sqrt(sum((Qtrue-Qnew[perms[i,],])^2))
  }
  #print(RMSEs)
  Qnew<-rbind(Qnew[perms[which.min(RMSEs),],],Qest[w,])
  Qnew
}

## returns the Q matrix with minimum difference to the old Q(that has 1 less population)
getFast<-function(Q,Qold){
  npop<-nrow(Qold)
  res<-c()
  for(g in 1:nrow(Q)){
    w<-rowSums((rep(Q[g,],each=npop)-Qold)^2)
    res<-rbind(res,c(which.min(w),min(w)))
  }
  dub <- duplicated(res[,1])
  dd<-res[dub,1]
  ww<-which.max(res[res[,1]==dd,2])
  res[which(res[,1]==dd)[ww],1]<-npop+1
  Q[order(res[,1]),]
}

colorFun<-function(x,Q){
  o<-Q[,x]
  if(length(x)==1)
    o<-matrix(o,ncol=1)
  
  K<-nrow(o)
  res<-matrix(0,ncol=ncol(o),nrow=K^2)
  most<-order(rowSums(o),decreasing=T)
  o<-o[,order(o[most[1],],decreasing=T)]
  if(length(x)==1)
    o<-matrix(o,ncol=1)
  
  for(k in 1:K)
    res[ (k-1)*K+most[k],]<-o[most[k],]
  res
}

mkOrd<-function(p,res,popOrd){
  x<-res[,popOrd==p]
  if(is.matrix(x)){
    w<-which.max(rowSums(x))
    m<-(1:length(popOrd))[popOrd==p][1]-1
    order(x[w,])+m
  }else{
    m<-(1:length(popOrd))[popOrd==p][1]-1
    c(1)+m
  }
}



#### global admixture ####
t<-read.table('~/Downloads/PCA/world.rmLow3andWSpecies.rmPlate1.v2.bamlist.info.csv',sep=",",header=T)
t2 <-t

# set the order
ord  <- order( t2$GlobalDataOrder)
region <- t2$RegionID[ord]
pop <- t2$PopID[ord]
pop2 = as.character(t2$RegionID[ord])
levels(pop2) = c(levels(pop2),"Atlantic Arctic")
pop2[which(pop2=="AtlanticNC")] = "Atlantic Arctic"
Kall <- 2:13

l<-list.files("~/Downloads/Admixture/World/",full=T,pattern="K[0-9].best.1.qopt")
files<-sort(l)

l2<-list.files("~/Downloads/Admixture/World/",full=T,pattern="K1[0-4].best.1.qopt")
files2<-sort(l2)

## revision ##
#l<-list.files("~/Documents/Project/Seal/popgen/revision/admixture/",full=T,pattern="K[0-9].best.1.qopt")
#files<-sort(l)

#l2<-list.files("~/Documents/Project/Seal/popgen/revision/admixture/",full=T,pattern="K1[0-4].best.1.qopt")
#files2<-sort(l2)
#files <- append(files,files2)


allQ <- list()
for (K in Kall) {
  allQ[[K]]<-t(read.table(files[K-min(Kall)+1]))[,ord]
}



### plot here
n<-length(pop)
reorder = 1
fast = T
lty =2 
lwd =2 

if(reorder==1){
  res<-allQ[[length(allQ)]]
  ord<-order(pop)
  res <-res[,ord]    
  popOrd<-pop[ord]
  u <-unlist(lapply(unique(pop), mkOrd,res=res,popOrd=popOrd))
  ordd <- ord[u]
}

pop<-pop[ordd]
Qold<-NA
par(mfrow=c(length(Kall)+2,1))

### information needed for modifying the plot
K3K4.ice = which(pop=="ICE")
K3K4.gre = which(pop=="GRE")

K9K10K11.bic = which(pop=="BIC")
K9K10K11.nfl = which(pop=="NFL")

for(K in Kall){
  par(mar=c(.1,5.1,.6,2.1))
  
  cat(K," ")
  Q<-allQ[[K]]
  Q<-Q[,ordd]
  if(K!=Kall[1]){
    if(fast)
      Q<-getFast(Q,Qold)
    else
      Q<-getRMSD(Q,Qold)
    
  }
  Qold<-Q
  
  Q<-do.call(cbind,tapply(1:length(pop),match(pop,unique(pop)),colorFun,Q=Q))
  
  #if (K==3 ) {
  #  Q[4,K3K4.ice] = Q[1,K3K4.ice]
  #  Q[3,K3K4.ice] = Q[6,K3K4.ice]
  #  Q[6,K3K4.ice] = 0
  #  Q[1,K3K4.ice] = 0
  #}
  
  if (K==3 || K==4) {
    row.ice.high = order(Q[,K3K4.ice[1]],decreasing = T)[1]
    row.ice.low = order(Q[,K3K4.ice[1]],decreasing = T)[2]
    
    row.gre.high = order(Q[,K3K4.gre[1]],decreasing = T)[1]
    row.gre.low = order(Q[,K3K4.gre[1]],decreasing = T)[2]
    
    Q[row.gre.high,K3K4.ice] = Q[row.ice.low,K3K4.ice]
    Q[row.gre.low,K3K4.ice] = Q[row.ice.high,K3K4.ice]
    Q[row.ice.low,K3K4.ice] = 0
    Q[row.ice.high ,K3K4.ice] = 0
  }
  
  if (K==9||K==10||K==11) {
    row.bic.high = order(Q[,K9K10K11.bic[3]],decreasing = T)[1]
    row.bic.low = order(Q[,K9K10K11.bic[3]],decreasing = T)[2]
    
    row.nfl.high = order(Q[,K9K10K11.nfl[3]],decreasing = T)[1]
    row.nfl.low = order(Q[,K9K10K11.nfl[3]],decreasing = T)[2]
    #col.index.new = which.max(Q[,K3K4.ref])
    #col.index.old = which.max(Q[,K3K4.ice[1]])
    
    Q[row.nfl.high,K9K10K11.bic] = Q[row.bic.low,K9K10K11.bic]
    Q[row.nfl.low,K9K10K11.bic] = Q[row.bic.high,K9K10K11.bic]
    Q[row.bic.low,K9K10K11.bic] = 0
    Q[row.bic.high ,K9K10K11.bic] = 0
  }
  
  ta<-tapply(pop,pop,length)
  small<- names(ta)[ta==1]
  #h<- barplot(Q,border=NA,col=1:K,space=0,ylab="Admixture proportion",main=paste("K = ",K,sep=""))
  h<- barplot(Q,border=NA,col=1:K,space=0,ylab="",yaxt='n')
  axis(side=2,las=1,at=c(0,1),c("0","1"),cex.axis=1.2)

  legend("right", inset = c(-0.2,0), legend= paste("K = ",K,sep="") )
  abline(v=tapply(h,pop,max)+0.5,col="black",lwd=lwd,lty=lty)
  med<-tapply(h,pop,median)
  
  med2 <- tapply(h,pop2,median)
  
}

par(mar=c(2.9,5.1,.1,2.1))
h<- barplot(Q,border=NA,col="transparent",space=0, axes=F)
text(med,rep(0.5,length(unique(pop))),names(med),xpd=TRUE,srt = 90,cex=1.5,adj=1)
par(mar=c(3,5.1,.1,2.1))
h<- barplot(Q,border=NA,col="transparent",space=0, axes=F)
text(med2,rep(0.5,length(unique(pop2))),names(med2),xpd=TRUE,srt = 30,cex=1.5,adj=1)





#### modify here, for certain K and pop, make the color order (vertical) consistent with the previous pop ####
allQ2 = allQ
K3K4.ice=which(pop=="ICE")
allQ2[[3]][2,K3K4.ice]  = allQ[[3]][3,K3K4.ice] 
allQ2[[3]][3,K3K4.ice] = allQ[[3]][2,K3K4.ice] 
allQ2[[4]][1,K3K4.ice]  = allQ[[4]][3,K3K4.ice] 
allQ2[[4]][3,K3K4.ice] = allQ[[4]][1,K3K4.ice] 

K9K10K11.bic = which(pop=="BIC")
allQ2[[9]][3,K9K10K11.bic] = allQ[[9]][4,K9K10K11.bic] 
allQ2[[9]][4,K9K10K11.bic] = allQ[[9]][3,K9K10K11.bic]
allQ2[[10]][1,K9K10K11.bic] = allQ[[10]][9,K9K10K11.bic] 
allQ2[[10]][9,K9K10K11.bic] = allQ[[10]][1,K9K10K11.bic] 
allQ2[[11]][4,K9K10K11.bic] = allQ[[11]][9,K9K10K11.bic] 
allQ2[[11]][9,K9K10K11.bic] = allQ[[11]][4,K9K10K11.bic]

plotMulti(allQ2,2:13,pop=pop,pop2=pop2,reorder=1,fast=T) 

dev.off()

#### Atlantic Structure ####
t<-read.table('~/Downloads/PCA/atlantic.rmLow3andWSpecies.rmPlate1.v2.bamlist.info.csv',sep=",",header=T)
t2 <-t

# set the order
ord  <- order( t2$GlobalDataOrder)
region <- t2$RegionID[ord]
pop <- t2$PopID[ord]


Kall <- 2:15

l<-list.files("~/Downloads/Admixture/Atlantic//",full=T,pattern="K[0-9].round.best.1.qopt")
files<-sort(l)

l2<-list.files("~/Downloads/Admixture/Atlantic/",full=T,pattern="K1[0-5].round.best.1.qopt")
files2<-sort(l2)

files <- append(files,files2)

allQ <- list()
for (K in Kall) {
  allQ[[K]]<-t(read.table(files[K-min(Kall)+1]))[,ord]
}

plotMulti(allQ,reorder=1,pop,fast=T) 


#### Pacific Strucutre ####
t<-read.table('~/Downloads/PCA/pacific.rmLow3andWSpecies.rmPlate1.v2.bamlist.info.csv',sep=",",header=T)
t2 <-t

# set the order
ord  <- order( t2$GlobalDataOrder)
region <- t2$RegionID[ord]
pop <- t2$PopID[ord]

Kall <- 2:4

l<-list.files("~/Downloads/Admixture/Pacific///",full=T,pattern="K[0-4].best.1.qopt")
files<-sort(l)
allQ <- list()
for (K in Kall) {
  allQ[[K]]<-t(read.table(files[K-min(Kall)+1]))[,ord]
}

plotMulti(allQ,reorder=1,pop,fast=T) 
