## Calculate HUC2 trends in each NANI component and compile results into barplot
library(raster); library(ggplot2); library(RColorBrewer)

#-------------------------------------------------------------------
## Read in HUC8 NANI components
dat.nani=read.table('./data/HUC8/HUC8.nani.components.txt', head=T, colClasses=c("HUC8"="character"))
sites.nani=unique(dat.nani$HUC8) #2107 sites
nani=matrix(dat.nani$NANI, ncol=length(sites.nani)) #26 x 2107; 1987-2012
fert=matrix(dat.nani$NANI_FERT, ncol=length(sites.nani))
nfix=matrix(dat.nani$NANI_NFIX, ncol=length(sites.nani))
ndep=matrix(dat.nani$NDEP, ncol=length(sites.nani))
nfcx=matrix(dat.nani$NANI_NFCX, ncol=length(sites.nani))
fdfd=matrix(dat.nani$NANI_FDFD, ncol=length(sites.nani))

#-------------------------------------------------------------------
## Read in shapefiles for the HUC8 watersheds
shape=shapefile("./data/shapefiles/WBDHU8_Reg_1_18.shp")
shape.df=fortify(shape) #transforms into format ggplot understands
huc2=shape$HUC2
huc2.area=aggregate(shape$AreaSqKm, by=list(huc2), FUN=sum)$x
huc2.names=as.character(aggregate(shape$Name_HUC2, by=list(huc2), FUN=unique)$x)
shape.idx=c() #figures out which idx in nani data corresponds to values in the shapefile
for (i in 1:2107){
  j=which(sites.nani==shape$HUC8[i])
  shape.idx=c(shape.idx, j)
}
#saveRDS(shape.idx, './Data/HUC8/HUC8.shapefile.idx.rds')

#-------------------------------------------------------------------
## Aggregate NANI fluxes to HUC2 scale
aggregate.to.huc2=function(dat, shape.idx, shape, huc2){
  dat.adj=dat[,shape.idx]
  dat.total=dat.adj*matrix(rep(shape$AreaSqKm, each=nrow(dat)), nrow=nrow(dat)) #26 x 2107
  dat.huc2=t(unlist(apply(dat.total, 1, function(x) aggregate(x, by=list(huc2), FUN=sum)[,2]))) #26 x 18 regions
  return(dat.huc2)
}
nani.huc2=aggregate.to.huc2(nani, shape.idx, shape, huc2)
fert.huc2=aggregate.to.huc2(fert, shape.idx, shape, huc2)
nfix.huc2=aggregate.to.huc2(nfix, shape.idx, shape, huc2)
nfcx.huc2=aggregate.to.huc2(nfcx, shape.idx, shape, huc2)
ndep.huc2=aggregate.to.huc2(ndep, shape.idx, shape, huc2)
fdfd.huc2=aggregate.to.huc2(fdfd, shape.idx, shape, huc2)

#-------------------------------------------------------------------
## Calculate trends in each HUC2 NANI component
trend.lm=function(data, return.pval=F){ #extracts the trend
  year=1:length(data)
  fit=coef(summary(lm(data~year)))
  if (return.pval==T){
    return(fit[2,4])
  } else {
    return(fit[2,1])
  }
}
trends.nani=apply(nani.huc2, 2, trend.lm) / huc2.area
trends.fert=apply(fert.huc2, 2, trend.lm) / huc2.area 
trends.nfix=apply(nfix.huc2, 2, trend.lm) / huc2.area  
trends.nfcx=apply(nfcx.huc2, 2, trend.lm) / huc2.area
trends.ndep=apply(ndep.huc2, 2, trend.lm) / huc2.area
trends.fdfd=apply(fdfd.huc2, 2, trend.lm) / huc2.area
pvals.nani=apply(nani.huc2, 2, trend.lm, return.pval=T)
nani.huc2.mean=apply(nani.huc2, 2, mean) / huc2.area


#-------------------------------------------------------------------
## Component barplots 
pdf(file='./figures/nani.component.trends.boxplot.pdf', width=7, height=6, useDingbats=F)
  mm = c(trends.fert, trends.nfix, trends.fdfd, trends.ndep, trends.nfcx) / 100 #divide by 100 to convert to hectares
  nn = c(rep('a', 18), rep('b', 18), rep('c', 18), rep('d', 18), rep('e', 18))
  oo = rep(huc2.names, 5)
  pp = rep(trends.nani, 5) / 100 #divide by 100 to convert to hectares
  qq=rep(nani.huc2.mean, 5) / 100 #divide by 100 to convert to hectares
  trends.nani.sig=trends.nani / 100; trends.nani.sig[pvals.nani>=.1]=NA  #divide by 100 to convert to hectares
  ggdat = data.frame(region=oo, Component=nn, Trend=mm, NANI_trend=pp, NANI_trend_sig=trends.nani.sig, NANI_mean=qq)
  order.idx=sort(trends.nani, index.return=T, decreasing=T)$ix
  
  p=ggplot() + geom_bar(data=ggdat, aes(x=region, y=Trend*10, fill=Component, color=Component), size=.6, stat='identity') +
    scale_x_discrete(limits=ggdat$region[order.idx]) +
    scale_color_manual(values=c(rev(brewer.pal(5, 'Accent'))), labels=c('Fertilizer','N Fixation','Food and Feed Import','Atmospheric Deposition','Non-Food Crop Export')) +
    scale_fill_manual(values=c(rev(brewer.pal(5, 'Accent'))), labels=c('Fertilizer','N Fixation','Food and Feed Import','Atmospheric Deposition','Non-Food Crop Export')) +
    geom_hline(yintercept=0, size=1) +
    geom_point(data=ggdat, aes(x=region, y=NANI_trend*10), color='orangered2', size=2.3, shape=22) +
    geom_point(data=ggdat, aes(x=region, y=NANI_trend_sig*10), color='orangered2', fill='orangered2', size=2, shape=22) +
    labs(fill='Component Trends', y='kg/ha/decade') +
    geom_point(data=ggdat, aes(x=(region), y=NANI_mean), color='black', size=1.3) +
    guides(color=F, size=F) +
    coord_flip() +theme_classic() + theme(legend.position='none',axis.title.y=element_blank(),axis.text.y=element_text(size=10), axis.title.x=element_text(size=11))
  p
dev.off()

