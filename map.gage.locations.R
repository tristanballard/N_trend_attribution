## Creates map of the GAGESII basins, HUC8 watersheds, and HUC2 watersheds
library(ggplot2); library(rgdal); library(raster)

#-------------------------------------------------------------------
## Read in selected GAGESII station names
sta.id=read.table('./data/GAGESII/gagesII_station_ids.txt', colClasses="character", head=F)[,1]
id.remove=c('09498500', '09504000') #remove Arizona stations
sta.id=sta.id[-which(sta.id %in% id.remove)]; length(sta.id)

## Get lat and lon for each station
sites.lat.lon=read.csv('./data/GAGESII/gagesII_lat_lon.csv', head=T, colClasses = "character") #3 columns
sites.idx=which(sites.lat.lon[,1] %in% sta.id)
sites.lat.lon=sites.lat.lon[sites.idx,]
lat=as.numeric(sites.lat.lon[,2]); lon=as.numeric(sites.lat.lon[,3])
sta.id=cbind.data.frame(sta.id, lat, lon)

## Read in GAGESII basin shapefiles, n=123
bas.all=readOGR('./data/shapefiles/', layer='GAGESII_basins')

## Prepare shapefiles for ploting in ggplot2
bas.all@data$id=rownames(bas.all@data)
bas.all.points=fortify(bas.all, region="id")
bas.all.df=merge(bas.all.points, bas.all@data, by="id", all=TRUE)

## Read in locations of original 70 watersheds used in the paper
sites.original=read.table('./data/GAGESII/gagesII_station_ids_Sinha2016.txt', colClasses="character", head=F)[,1]
bas.original=bas.all.df$GAGE_ID %in% sites.original

## Read in HUC2 watershed boundaries
shape.huc2=readRDS('./data/shapefiles/shapefile.HUC2.rds')

## Read in HUC8 watershed boundaries
shape=shapefile("./data/shapefiles/WBDHU8_Reg_1_18.shp")
shape.huc8=fortify(shape) #transforms into format ggplot understands

## Read in shapefile for USA boundaries
usa=map_data('usa')

#-------------------------------------------------------------------
## Plot the basins, takes a few minutes 
## Watersheds with points and w/ HUC2 and HUC8 overlaid
p5=ggplot() +
  geom_polygon(data=usa, aes(x=long,y=lat,group=group),fill='white',color='black', size=.18) +
  geom_polygon(data=shape.huc8, aes(x=long, y=lat, group=group), fill=NA, color='grey62', size=.03) +
  geom_polygon(data=shape.huc2, aes(x=long, y=lat, group=group), fill=NA, color='black', size=.32) +
  geom_path(data=bas.all.df[bas.original,], aes(x=long,y=lat,group=group), color='steelblue2', size=.6) + #'navy'
  geom_path(data=bas.all.df[!bas.original,], aes(x=long,y=lat,group=group), color='olivedrab3', size=.65) + #'gold2' #'olivedrab3' #'steelblue1'
  geom_point(data=sta.id[,],aes(x=lon, y=lat),color='black',size=.7) +
    theme_classic() + theme(axis.line=element_blank(),axis.text.x=element_blank(), axis.text.y=element_blank(),axis.ticks=element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank())
ggsave(p5, file = "./figures/map.gage.locations.png", width=7.4, height=5, type = "cairo-png") #takes a few minutes

