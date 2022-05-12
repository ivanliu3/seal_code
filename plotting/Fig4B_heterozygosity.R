library(RColorBrewer)
library(dplyr)
### 2021 01 25 ####
N<-11
mypalette <- brewer.pal(N,"RdBu")
mypalette <- append(mypalette,brewer.pal(5,"Spectral")[2] )
#world.palette <- mypalette[c(3,1,5,9,11)]
world.palette <- mypalette[c(11,9,12,3,1)]
color= world.palette
setwd('~/Downloads/heterozygosity/')
info <- t<-read.table('~/Downloads/PCA/world.rmLow3andWSpecies.rmPlate1.v2.bamlist.info.csv',sep=",",header=T)

files = list.files('./',pattern=".*fs$")
pops<- unique( info$PopID[order(info$GlobalDataOrder)] )

listA <- list()
for (pop in pops) {
  samples <-info$id[which(info$PopID  %in% pop)]
  #pop.files <- files[grep(pop, files)]
  pop.files <- paste(samples,'.fs',sep="")
  m <- matrix( rep(NA, 2*length(pop.files)),ncol = 2, byrow = T)
  
  
  for (i in 1:length(pop.files)) {
    t <- scan(pop.files[i])
    m[i,1] <- gsub('.fs','',pop.files[i])
    m[i,2] <- as.numeric( t[2]/sum(t) )
  }
  listA[[pop]] <-m
  print(pop)
}




library(reshape2)
library(ggplot2)
#library(ggpubr)
listB <- lapply(names(listA),function(x)cbind(name=x,listA[[x]]))

df2 <- data.frame(Reduce(rbind,listB))
names(df2) <- c('pop','ind','value')
df2$value <- as.numeric(as.character(df2$value))
df2$pop <- factor(df2$pop,levels=pops)
df2$color <- info$RegionID[order(info$GlobalDataOrder)]
df2$color <- factor(df2$color,levels= unique(df2$color))
df2$color <- color [as.numeric(df2$color)]

new_order <- c('HOK','KOD','END','CAL',
               'NFL','BIC', 'MET', 'NEW', 'GRE', 
               'ICE', 'SVA','ORK', 'LIS', 'ROL', 
               'SOT', 'TRO', 'WAS', 'WNL',
               'LIE', 'ANH', 'MAA', 'ROD')
               



df2$order <- apply(df2,1,function(x) which(new_order==x[1]))
df3 <-df2[order(df2$order),]
levels(df3$pop) <- as.factor(c('HOK','KOD','END','CAL',
                                'NFL','BIC', 'MET', 'NEW', 'GRE', 
                                'ICE', 'SVA','ORK', 'LIS', 'ROL', 
                                'SOT', 'TRO', 'WAS', 'WNL',
                                'LIE', 'ANH', 'MAA', 'ROD'))


p1 <- ggplot(df3, aes(x=pop, y=value)) + 
  geom_boxplot(alpha = 0.50) +
  geom_point( size = 2,col=as.factor(df3$color) )

mean2 <- df2 %>%
  group_by(color) %>%
  dplyr::summarise(Mean=mean(value,na.rm = T)) %>%
  ungroup %>%
  as.data.frame() 

mean3 <- df2 %>%
  group_by(pop) %>%
  dplyr::summarise(Mean=mean(value,na.rm = T)) %>%
  ungroup %>%
  as.data.frame() 

het_plot <-p1 + xlab('Populations') + ylab('Heterozygosity') + theme_bw()

het_plot+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


het_plot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                 axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
                 panel.background = element_blank(), axis.line = element_line(colour = "black"))


het_plot + geom_segment(data=df2, mapping=aes(x=0.5,xend=1.5,y=0.0006883950,yend=0.0006883950),col="#053061",size=2,linetype='twodash') + 
  geom_segment(data=df2, mapping=aes(x=1.5,xend=4.5,y=0.0012507392,yend=0.0012507392),col="#4393C3",size=2,linetype='twodash') +
  geom_segment(data=df2, mapping=aes(x=4.5,xend=8.5,y=0.0006342250,yend=0.0006342250),col="#FDDBC7",size=2,linetype='twodash') + 
  geom_segment(data=df2, mapping=aes(x=8.5,xend=11.5,y=0.0005020388,yend=0.0005020388),col="#D6604D",size=2,linetype='twodash') + 
  geom_segment(data=df2, mapping=aes(x=11.5,xend=26.5,y=0.0003642373,yend=0.0003642373),col="#67001F",size=2,linetype='twodash')
##########################

het_plot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                 axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
                 panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14,face="bold"))
