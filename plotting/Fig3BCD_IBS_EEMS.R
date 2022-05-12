## set the plot arragnment here ##
layout.matrix <- matrix (c(1,2,3,4),nrow=2,ncol=2,byrow = T)
layout.matrix
layout(mat = layout.matrix,
       heights = c(1, 1), # Heights of the two rows
       widths = c(1,1)) # Widths of the three columns

layout.show(4)

### 3A , Fst, empty
plot(0,type='n',axes=FALSE,ann=FALSE)
### end of 3A ###

### 3B, IBS water distance ###
load("~/Documents/Project/Seal/popgen/Figure_final/IBS_waterdistance.final.RData")
final3

atlantic.i <- which(final3$type =='atlantic')
pacific.i <- which(final3$type == 'pacific')
across.i <- which(final3$type == 'across')
par(mar=c(5,6,1,2))
plot(final3$value.y[across.i],final3$value.x[across.i],xlab="Geographical distance (km)",
     ylab=expression(paste( italic(F[ST])," value")), pch=1,ylim=c(0,1),
     xlim=c(0,13000000),col="#B385FF",
     yaxt = "n",xaxt="n",cex.lab=1.5,cex=1.5)

axis(side=2,las=2,mgp = c(3, 0.75, 0),cex.axis=1.5)
axis(side=1, at=seq(0,12000000,3000000), labels=c('0','3000','6000','9000','12000'),cex.axis=1.5)
points(final3$value.y[atlantic.i],final3$value.x[atlantic.i],col="#F8766D",pch=1,cex=1.5)
points(final3$value.y[pacific.i],final3$value.x[pacific.i],col="#00A6FF",pch=1,cex=1.5)
legend('topleft',c("Within Pacific","Within Atlantic","Across Oceans"),col=c("#00A6FF","#F8766D","#B385FF"),pch=1,bty="n",cex=1.5)
legend('topright',c("Pacific with HOK", "Across Ocean with HOK"),col=c("#00A6FF","#B385FF"),pch=15,bty="n",cex=1.5)


###
index1 <- which(grepl('HOK',final3$name) & final3$type=='pacific')
points(final3$value.y[index1],final3$value.x[index1],col="#00A6FF",pch=15,cex=1.5)

#index2 <- which(!grepl('HOK',final3$name) & final3$type=='pacific')
#points(final3$value.y[index2],final3$value.x[index2],col="pink",pch=16)


index3 <- which(grepl('HOK',final3$name) & final3$type=='across')
points(final3$value.y[index3],final3$value.x[index3],col="#B385FF",pch=15,cex=1.5)

## figure size 800 X 700
