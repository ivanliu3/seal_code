
library(reemsplots2)
library(dplyr)
mcmcpath = paste0("/Users/pls394/Documents/Project/Seal/popgen/EEMS/result_world_202106/world_202106/barrier-nIndiv286-nSites13540-nDeme600-chain",
  2:10,"-restart")
plots <- make_eems_plots(mcmcpath, longlat = TRUE)

m_surface <-plots$mrates01

#### prepare map and project here ####

require(maps)
require(ggplot2)
world_map = data.frame(map(plot=FALSE)[c("x","y")])
names(world_map) = c("lon","lat")
world_map = within(world_map, {
  lon = ifelse(lon < 0, lon + 360, lon)
})
ggplot(aes(x = lon, y = lat), data = world_map) + geom_path()

world_map2 = within(world_map, {
  lon = ifelse(lon < 0, lon + 360, lon)
  lon = ifelse((lon < 1) | (lon > 359), NA, lon)
})
ggplot(aes(x = lon, y = lat), data = world_map2) + geom_path()

world_map3 = within(world_map2,{
  lon = ifelse(lon < 50, lon + 360, lon)
})
ggplot(aes(x = lon, y = lat), data = world_map3) + geom_path() + coord_map()


ggworld <- geom_polygon(data=world_map3, aes(x=lon,  y=lat),
                        fill=NA, colour="black")

m_surface2 <- m_surface + ggworld

### add points ###

### not used ###
sites <- read.table('~/Documents/Project/Seal/popgen/LocalityCoordinates.seal.new.csv',sep=",",header=T) 
head(sites)
names(sites) = c("site","lat","lon")


sites2 = within(sites, {
  lon = ifelse(lon < 0, lon + 360, lon)
})

sites3 = within(sites, {
  lon = ifelse(lon < 50, lon + 360, lon)
})

### end of not used ####

sites <- read.table("~/Documents/Project/Seal/popgen/EEMS/result_world_202106/world_202106/barrier-nIndiv286-nSites13540.coord",
                    sep=" ",header=F)

names(sites) <- c("lon","lat")

sites2 <- sites %>% 
   group_by(lon,lat) %>% summarize(n=n()) %>%
    as.data.frame()


m_surface2 +  geom_point(shape=20,color="darkviolet" ,aes(x=sites2$lon,y=sites2$lat,size=n/100),data = sites2, show.legend = F)
  
save.image("~/Documents/Project/Seal/popgen/EEMS/result_world_202106/eems_202106_ggPlot.RData")
load("~/Documents/Project/Seal/popgen/EEMS/result_world_202106/eems_202106_ggPlot.RData")
