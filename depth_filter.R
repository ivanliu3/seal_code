
#### final
t <- read.table('/home/users/xiaodong/Documents/Project/Seal/depth/bed/FilterSample.depth.bg',header=F,sep="\t")
t2 <- t[,4:356]


means <- apply(t2,1,mean) 
stdv <- apply(t2,1,sd)
means2 <- apply(t2,1,function(x) sum(x)/length(which(x !=0)) )
stdv2 <- apply(t2,1,function(x) sd(x[which(x!=0)]))
no_zero<- apply( t2,1,function(x) length(which(x==0)) )
### missingness at locus
pdf('missingness_atlocus.pdf')
hist(no_zero,xlim=c(0,400),breaks=100,xlab='Number of indiviuals of no reads',ylab='Frequency of loci')
abline(v=300,col='red',lty=2)
abline(v=180,lty=2)
dev.off()

i_miss <- which(no_zero<300)
length(i_miss)
##[1] 4340984

#### the average depth per sample and standard variance at loci, using overall mean
par(mfrow=c(2,1))

hist(means,xlim=c(0,70),breaks=5000,ylim=c(0,800000),main='Average depth at locus without filter',ylab='frequency of loci')
abline(v=2.5,col='red',lty=2)
abline(v=25,col='red',lty=2)
abline(v=50,col='red',lty=2)

hist(means[i_miss],xlim=c(0,70),breaks=5000,ylim=c(0,800000),main='Average depth at locus after removing loci with high datamissing',ylab='frequency of loci')

lapply(list(means,means[i_miss]),function(x) quantile(x,c(0.05,0.95,0.975)))
## without removing loci of high data missing 
##      5%      95%    97.5% 
## 0.00000 38.96307 44.93182 

## with removing loci of high data missing
##        5%        95%      97.5% 
## 0.4829545 43.6051136 48.8863636 
imiss_mean25 <-intersect(i_miss,which(means > 2.5 & means <48.886))
imiss_mean2 <-intersect(i_miss,which(means > 2 & means <48.886))
imiss_mean15 <-intersect(i_miss,which(means > 1.5 & means <48.886))
imiss_mean1 <-intersect(i_miss,which(means > 1 & means <48.886))

### standard deviation of depth at loci

hist(stdv[imiss_mean1],xlim=c(0,40),breaks=80,main='sd of depths after filter missingness and average depth',xlab='sd of depth',ylab='frequency of loci')

imiss_mean1_sd <-intersect(imiss_mean1, which(stdv<=40))
imiss_mean15_sd <-intersect(imiss_mean15, which(stdv<=40))
imiss_mean2_sd <-intersect(imiss_mean2, which(stdv<=40))
imiss_mean25_sd <-intersect(imiss_mean25, which(stdv<=40))                                                      


#### what used below ####
#### using mean of sequenced ind, means2 ####
j_miss <- which(no_zero<180) # half of the individuals
length(j_miss)
#3868686

pdf('Average depth of sequenced individuals at locus.pdf') 
par(mfrow=c(2,1))
hist(means2,xlim=c(0,70),breaks=5000,ylim=c(0,800000),main='Average depth at locus without filter',ylab='frequency of loci')
abline(v=2.5,col='red',lty=2)
abline(v=25,col='red',lty=2)
abline(v=50,col='red',lty=2)

hist(means2[j_miss],xlim=c(0,70),breaks=5000,ylim=c(0,800000),main='Average depth at locus after removing loci with high datamissing',ylab='frequency of loci')
dev.off()

lapply(list(means2,means2[j_miss]),function(x) quantile(x,c(0.05,0.95,0.975),na.rm=T))
## without removing loci of high data missing 
##     5%     95%   97.5% 
## 1.0000 40.0085 45.7762 

## with removing loci of high data missing, go with this one
##       5%       95%     97.5% 
## 2.370911 44.504249 49.623229 

jmiss_meantwo25 <-intersect(j_miss,which(means2 > 2.5 & means <50))
jmiss_meantwo2 <-intersect(j_miss,which(means2 > 2 & means < 50)) # this what we used,3696232
jmiss_meantwo15 <-intersect(j_miss,which(means2 > 1.5 & means < 50))
jmiss_meantwo1 <-intersect(j_miss,which(means2 > 1 & means <50))

### standard deviation of depth at loci
pdf('SD depth of sequenced individuals at locus.pdf')
hist(stdv2[jmiss_meantwo2],xlim=c(0,40),breaks=80,main='sd of depths after filter missingness and average depth',xlab='sd of depth',ylab='frequency of loci')
dev.off()

jmiss_meantwo2_sd <-intersect(jmiss_meantwo2, which(stdv2<=30)) # 3691847




write.table(t[jmiss_meantwo2_sd,1:3],file="int_meantwo2_sd.bed",col.names=F,row.names=F,quote=FALSE,sep="\t")

####
bed.t<-read.table('bed.include.sample.list',header=F)
samples <-factor( gsub('.cov.bg','',gsub('/home/users/xiaodong/Documents/Project/Seal/depth/bed/','',bed.t$V1)))

info <- read.table('../../sample_list/World.new.rmLow3andWSpecies.info.csv',sep=",",header=T)

order<-match(samples,info$id)
info2<-info[order,]

index <-imiss_mean1_sd[which(! imiss_mean1_sd %in% jmiss_meantwo1_sd)]
pop.miss <- matrix(NA,ncol=25,nrow=length(index))
for (i in 1:25) {
    pop.i <- which(info2$GlobalDataOrder ==i )
    pop.miss[,i] <-  apply( t2[index,pop.i],1,function(x) length(which(x==0))/length(pop.i) )
}

par(mfrow=c(5,5))
for (i in 1:25) {
    hist(pop.miss[,i],ylim=c(0,2000000),main=i)
}



#### 20201226 ####
setwd("/home/users/xiaodong/Documents/Project/Seal/depth_v2/20201227_bed")
t <- read.table('./FilterSample.rmPlate1a.depth.bg',header=F,sep="\t")
t2<-t[,4:297]
means <- apply(t2,1,mean)
stdv <- apply(t2,1,sd)
means2 <- apply(t2,1,function(x) sum(x)/length(which(x !=0)) )
stdv2 <- apply(t2,1,function(x) sd(x[which(x!=0)]))
no_zero<- apply( t2,1,function(x) length(which(x==0)) )
length(no_zero)

### missingness at locus
pdf('missingness_atlocus.pdf')
hist(no_zero,xlim=c(0,300),breaks=100,xlab='Number of indiviuals of no reads',ylab='Frequency of loci')
abline(v=240,col='red',lty=2)
abline(v=150,lty=2)
dev.off()

j_miss <- which(no_zero<150) # half of the individuals
length(j_miss) #3298006


hist(means2,xlim=c(0,70),breaks=5000,ylim=c(0,800000),main='Average depth at locus without filter',ylab='frequency of loci')
abline(v=2.5,col='red',lty=2)
abline(v=25,col='red',lty=2)
abline(v=50,col='red',lty=2)


### mean of sequenced ind, means2 ###

pdf('Average depth of sequenced individuals at locus.pdf')
par(mfrow=c(2,1))
hist(means2,xlim=c(0,60),breaks=5000,ylim=c(0,800000),main='Average depth at locus without filter',ylab='frequency of loci')
abline(v=2.5,col='red',lty=2)
abline(v=25,col='red',lty=2)
abline(v=50,col='red',lty=2)

hist(means2[j_miss],xlim=c(0,60),breaks=5000,ylim=c(0,800000),main='Average depth at locus after removing loci with high datamissing',ylab='frequency of loci')
dev.off()


> lapply(list(means2,means2[j_miss]),function(x) quantile(x,c(0.05,0.95,0.975),na.rm=T))
#[[1]]
#      5%      95%    97.5% 
# 1.00000 40.48639 46.77211 

#[[2]]
#       5%       95%     97.5% 
# 2.269841 45.445578 51.129252 


jmiss_meantwo25 <-intersect(j_miss,which(means2 > 2.5 & means <50))
jmiss_meantwo2 <-intersect(j_miss,which(means2 > 2 & means < 50)) # this is what we used, #3112270
jmiss_meantwo15 <-intersect(j_miss,which(means2 > 1.5 & means < 50))
jmiss_meantwo1 <-intersect(j_miss,which(means2 > 1 & means <50))


### standard deviation of depth at loci 
pdf('SD depth of sequenced individuals at locus.pdf')
hist(stdv2[jmiss_meantwo2],xlim=c(0,40),breaks=80,main='sd of depths after filter missingness and average depth',xlab='sd of depth',ylab='frequency of loci')
abline(v=30,col="red")
dev.off()

jmiss_meantwo2_sd <-intersect(jmiss_meantwo2, which(stdv2<=30)) # 3107291

write.table(t[jmiss_meantwo2_sd,1:3],file="int_meantwo2_sd.new.bed",col.names=F,row.names=F,quote=FALSE,sep="\t")


#### 20201231 ####
setwd("/home/users/xiaodong/Documents/Project/Seal/depth_v2/20201227_bed_r2")
t <- read.table('./FilterSample.rmPlate1andRelate.depth.bg',header=F,sep="\t")
dim(t) # 6268433     289
t2<-t[,4:289]
means <- apply(t2,1,mean)
stdv <- apply(t2,1,sd)
means2 <- apply(t2,1,function(x) sum(x)/length(which(x !=0)) )
stdv2 <- apply(t2,1,function(x) sd(x[which(x!=0)]))
no_zero<- apply( t2,1,function(x) length(which(x==0)) )
length(no_zero)

### missingness at locus
pdf('missingness_atlocus.pdf')
hist(no_zero,xlim=c(0,300),breaks=100,xlab='Number of indiviuals of no reads',ylab='Frequency of loci')
abline(v=230,col='red',lty=2)
abline(v=145,lty=2) # half of the individuals
dev.off()

j_miss <- which(no_zero<145) # half of the individuals
length(j_miss) #3253859


hist(means2,xlim=c(0,70),breaks=5000,ylim=c(0,800000),main='Average depth at locus without filter',ylab='frequency of loci')
abline(v=2.5,col='red',lty=2)
abline(v=25,col='red',lty=2)
abline(v=50,col='red',lty=2)


### mean of sequenced ind, means2 ###

pdf('Average depth of sequenced individuals at locus.pdf')
par(mfrow=c(2,1))
hist(means2,xlim=c(0,60),breaks=5000,ylim=c(0,800000),main='Average depth at locus without filter',ylab='frequency of loci')
abline(v=2.5,col='red',lty=2)
abline(v=25,col='red',lty=2)
abline(v=50,col='red',lty=2)

hist(means2[j_miss],xlim=c(0,60),breaks=5000,ylim=c(0,800000),main='Average depth at locus after removing loci with high datamissing',ylab='frequency of loci')
dev.off()


lapply(list(means2,means2[j_miss]),function(x) quantile(x,c(0.05,0.95,0.975),na.rm=T))
#[[1]]
#      5%      95%    97.5% 
# 1.00000 40.60839 46.84266 

#[[2]]
#       5%       95%     97.5% 
# 2.291457 45.524476 51.174825 


jmiss_meantwo25 <-intersect(j_miss,which(means2 > 2.5 & means <50))
jmiss_meantwo2 <-intersect(j_miss,which(means2 > 2 & means < 50)) # this is what we used, #3075067
jmiss_meantwo15 <-intersect(j_miss,which(means2 > 1.5 & means < 50))
jmiss_meantwo1 <-intersect(j_miss,which(means2 > 1 & means <50))


### standard deviation of depth at loci 
pdf('SD depth of sequenced individuals at locus.pdf')
hist(stdv2[jmiss_meantwo2],xlim=c(0,40),breaks=80,main='sd of depths after filter missingness and average depth',xlab='sd of depth',ylab='frequency of loci')
abline(v=30,col="red")
dev.off()

jmiss_meantwo2_sd <-intersect(jmiss_meantwo2, which(stdv2<=30)) # 3069708

write.table(t[jmiss_meantwo2_sd,1:3],file="int_meantwo2_sd.new.r2.bed",col.names=F,row.names=F,quote=FALSE,sep="\t")



#### 20210120 r2 ####
setwd("/home/users/xiaodong/Documents/Project/Seal/depth_v2/202101_bed_r2" )
t <- read.table('FilterSample.depth.r2.bg',header=F,sep="\t")
dim(t) #  5980690     289
t2<-t[,4:289]
means <- apply(t2,1,mean)
stdv <- apply(t2,1,sd)
means2 <- apply(t2,1,function(x) sum(x)/length(which(x !=0)) )
stdv2 <- apply(t2,1,function(x) sd(x[which(x!=0)]))
no_zero<- apply( t2,1,function(x) length(which(x==0)) )
length(no_zero) #  5980690

### missingness at locus
pdf('missingness_atlocus.pdf')
hist(no_zero,xlim=c(0,300),breaks=100,xlab='Number of indiviuals of no reads',ylab='Frequency of loci')
abline(v=230,col='red',lty=2)
abline(v=145,lty=2) # half of the individuals
dev.off()

j_miss <- which(no_zero<145) # half of the individuals
length(j_miss) # 3211253



### mean of sequenced ind, means2 ###

pdf('Average depth of sequenced individuals at locus.pdf')
par(mfrow=c(2,1))
hist(means2,xlim=c(0,60),breaks=5000,ylim=c(0,800000),main='Average depth at locus without filter',ylab='frequency of loci')
abline(v=2.5,col='red',lty=2)
abline(v=25,col='red',lty=2)
abline(v=50,col='red',lty=2)

hist(means2[j_miss],xlim=c(0,60),breaks=5000,ylim=c(0,800000),main='Average depth at locus after removing loci with high datamissing',ylab='frequency of loci')
dev.off()

lapply(list(means2,means2[j_miss]),function(x) quantile(x,c(0.05,0.95,0.975),na.rm=T))
#[[1]]
#      5%      95%    97.5% 
# 1.00000 41.11228 47.26573 

#[[2]]
#       5%       95%     97.5% 
# 2.286432 45.617544 51.237762 


jmiss_meantwo25 <-intersect(j_miss,which(means2 > 2.5 & means <50))
jmiss_meantwo2 <-intersect(j_miss,which(means2 > 2 & means < 50)) # this is what we used, # 3033193
jmiss_meantwo15 <-intersect(j_miss,which(means2 > 1.5 & means < 50))
jmiss_meantwo1 <-intersect(j_miss,which(means2 > 1 & means <50))


### standard deviation of depth at loci 
pdf('SD depth of sequenced individuals at locus.pdf')
hist(stdv2[jmiss_meantwo2],xlim=c(0,40),breaks=80,main='sd of depths after filter missingness and average depth',xlab='sd of depth',ylab='frequency of loci')
abline(v=30,col="red")
dev.off()

jmiss_meantwo2_sd <-intersect(jmiss_meantwo2, which(stdv2 <30)) #3027883

write.table(t[jmiss_meantwo2_sd,1:3],file="int_meantwo2_sd.new.r2.bed",col.names=F,row.names=F,quote=FALSE,sep="\t")
