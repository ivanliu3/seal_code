library(rworldmap)
library(dplyr)
library(ggplot2)
library(geosphere)
library(gpclib)
library(ggrepel)
###
# https://egallic.fr/en/maps-with-r/
###
# World map
worldMap <- getMap()
world.points <- fortify(worldMap)
world.points$region <- world.points$id

world.df <- world.points[,c("long","lat","group", "region")]
Nsphere.df <- world.df[world.df$lat>=40,]
worldmap <- ggplot() + 
  geom_polygon(data = world.df, aes(x = long, y = lat, group = group)) +
  scale_y_continuous(breaks = (-2:2) * 30) +
  scale_x_continuous(breaks = (-4:4) * 45)

worldmap

## change the projection
worldmap <- ggplot() + 
  geom_polygon(data = Nsphere.df, fill="grey",aes(x = long, y = lat, group = group)) +
  scale_y_continuous(breaks = (-2:2) * 30) +
  scale_x_continuous(breaks = (-4:4) * 45) +
  coord_map("ortho", orientation=c(75, 0, 0))+
      theme_light() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) 
  
worldmap 

# Given the long/lat coordinates of an origin (x) and a radius (radius) in km,
# returns the coordinates of 360 points on the circle of center x and radius radius km.
# --
# x (numeric vector) : coordinates of the origin of the circle
# radius (numeric) : radius of the circle
# --
distantCircle <- function(x, radius) {
  # Creation de 360 points distincts sur le cercle de centre
  # x et de rayon radius
  resul <- do.call("rbind", lapply(0:360, function(bearing) {
    res <- destPoint(p = x, b = bearing, d = radius)
    rownames(res) <- NULL
    return(data.frame(res))
  }))
  resul$dist <- radius / 1000
  return(resul)
}
circle.5600 <- distantCircle(x = c(0.0000001,89.9999999), radius = 5600*1000)
circles <- rbind(circle.5600)
circles$dist <- factor(circles$dist)



library('dplyr')

#t <- read.table('~/Downloads/world.rmLow3andWSpecies.bamlist.info.csv',sep=",",header=T)
#pops <- t %>% group_by(GlobalDataOrder) %>% summarise(first=first(PopID)) %>% as.data.frame()
#pops


#pops$logitude <- c(144.8635, -152.4072, -133.5895, -119.4179, -70.6620, -68.8608, -68.0154, -57.6605, -40.6043, -19.0208,
#                   20.9752, -5.5263, -2.9605,  6.0143,  10.3487, 18.9553, 0.2500, 6.9495,  8.5139,  8.9514,
#                   8.0182,  11.0674, 11.5547,  12.8134, 17.0275, 16.3307)

#pops$latitude <- c(43.2203,  57.7900,   57.8569,   36.7783,   42.6159,  48.3994,  48.6666, 53.1355,  65.7069,   64.9631,
#                   77.8750, 56.5006, 58.9809,  59.1490, 63.0137, 69.6492, 52.9167, 53.5787, 55.1389, 56.8937,
#                   58.1599, 58.5774, 56.7120,  55.3593, 69.1339, 56.5911)

t2 <- read.table("~/Documents/Project/Seal/popgen/LocalityCoordinates.seal.new.csv",sep = ",",header=T)
library(RColorBrewer)
N <- 11
mypalette <- brewer.pal(N,"RdBu")
world.palette <- mypalette[c(11,9,5,3,1)]

t2$col <- c(world.palette[1],rep(world.palette[2],3),rep(world.palette[3],4),rep(world.palette[4],3),rep(world.palette[5],11))
(sites <- data.frame(longitude = t2$E,
                     latitude =  t2$N,
                     name = t2$Locality
                     ) )
P <- worldmap + geom_path(data = circles, aes(x = lon, y = lat), linetype = 1,lwd=1,color="black") 
P

P +  
  geom_point(data = sites, aes(x = longitude, y = latitude), size = 4, 
             shape = 16, fill = "darkred",col=(t2$col)) +
  labs(x="",y="") +
  theme(
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank()) +
  annotate("text",x=-30,y=55,label="North Atlantic Ocean",fontface =2,cex=6) +
  annotate("text",x=0,y=130,label="North Pacific Ocean",fontface=2,cex=6) + 
  geom_text_repel(data=sites, label=sites$name,
                   mapping = aes(x = longitude, y = latitude),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50',size=6) 



### the version below not used, Morten doesn't want many lines ###
P +  
  geom_point(data = sites, aes(x = longitude, y = latitude), size = 4, 
             shape = 16, fill = "darkred",col=(t2$col)) +
  labs(x="",y="") +
  theme(
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  annotate("text",x=-30,y=55,label="North Atlantic Ocean",fontface =2,cex=6) +
  annotate("text",x=0,y=130,label="North Pacific Ocean",fontface=2,cex=6) + 
  geom_text_repel(data          = sites[sites$longitude<0,], label=sites[sites$longitude<0,]$name,
                  size =6,
                  mapping = aes(x = longitude, y = latitude),
                  direction     = "y",
                  force         = 50,
                  nudge_x =  0- subset(sites,  longitude<0)$longitude ) +
                
  geom_text_repel(data          = sites[sites$longitude>0,], label=sites[sites$longitude>0,]$name,
                  size =6,
                  mapping = aes(x = longitude, y = latitude),
                  direction     = "y",
                  force         = 50,
                  nudge_x =  0+ subset(sites,  longitude>0)$longitude ) 

