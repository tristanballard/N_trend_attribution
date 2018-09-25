## This script matches up each basin watersheds with the fractional coverage of any HUC8's it overlaps
## It then takes the HUC8 flux estimates and aggregates it up for that basin

library(raster); library(ggplot2); library(sf); library(lwgeom)
## HUC8 shapefile
shape.huc8=shapefile("./data/shapefiles/WBDHU8_Reg_1_18.shp")
shape.huc8=spTransform(shape.huc8, CRS("+proj=longlat +datum=WGS84"))

## HUC2 shapefile used for plotting
shape.huc2=readRDS('./data/shapefiles/shapefile.HUC2.rds')

## HUC8 TN estimates
TN.huc8=readRDS('./data/HUC8/HUC8.TN.rds')
shape.idx=readRDS('./data/HUC8/HUC8.shapefile.idx.rds')
TN.huc8=TN.huc8[,shape.idx]
n.years=nrow(TN.huc8) #26 years

## Basin attributes
basin.meta=readRDS('./data/NASQAN/nasqan.metadata.rds')
n.basins=nrow(basin.meta); n.basins

extract.TN.estimates=function(shape.target, shape.target.area, shape.huc8=shape.huc8, TN.huc8=TN.huc8){
  a=st_as_sf(shape.huc8)
  b=st_as_sf(shape.target)
  
  overlap=st_intersection(a, b)
  fractional.coverage = as.numeric(st_area(overlap)) / 1000000 / overlap$AreaSqKm
  
  ## Minor quality control
  pred.area = as.numeric(sum(st_area(overlap))/1000000) #km2
  comp.area = pred.area / shape.target.area * 100
  print(paste0('Overlapped region is ',round(comp.area,2),'% of actual shapefile area')) # Won't be exact, but should be close
  
  ## Extract relevant HUC8 TN
  if(length(unique(overlap$HUC8))<length(overlap$HUC8)){print('Stop: HUC8 Error')}
  overlap.idx=shape.huc8$HUC8 %in% overlap$HUC8
  TN.overlap=TN.huc8[,overlap.idx] 
  TN.target = t(t(TN.overlap) * shape.huc8$AreaSqKm[overlap.idx] * (fractional.coverage)) #n.years x n.huc8 in target
  TN.target = apply(TN.target, 1, sum) / pred.area
  return(TN.target)
}

## Now apply the function to each basin and combine TN estimates into a matrix (basin.TN), takes awhile
## Also plot the shapefile of the basin on a map while the shapefile is loaded
basin.TN=matrix(rep(NA, n.years*n.basins), ncol=n.basins)
for (i in 1:n.basins){
  print(i)
  filename=paste0('./data/shapefiles/NASQAN/b',basin.meta$name[i],'.shp')
  shape.target=shapefile(filename)
  shape.target=spTransform(shape.target, CRS("+proj=longlat +datum=WGS84"))
  basin.TN[,i]=extract.TN.estimates(shape.target=shape.target, shape.target.area=basin.meta$Area_km2[i], shape.huc8=shape.huc8, TN.huc8=TN.huc8)

  ## Plot the basin for good measure
  s1=fortify(shape.target); rm(shape.target)
  map = ggplot() + geom_polygon(data=shape.huc2, aes(x=long, y=lat, group=group), fill=NA, color='black', size=.25) +
        geom_polygon(data=s1, aes(x=long, y=lat, group=group), fill='pink1', color='black', size=.23) +
        geom_point(data=basin.meta[i,], aes(x=lon, y=lat), color='black') +
        labs(title = expression(paste('', sep='')), fill = "") +
        theme_classic() + theme(axis.line=element_blank(),axis.text.x=element_blank(), axis.text.y=element_blank(),axis.ticks=element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank())
  ggsave(map, file = paste0("./data/NASQAN/b",basin.meta$name[i],".map.png"), width=7.5, height=5, type = "cairo-png")
  rm(s1)
}
rownames(basin.TN)=1987:2012
colnames(basin.TN)=basin.meta$name
#saveRDS(basin.TN, './data/nasqan.TN.rds')

