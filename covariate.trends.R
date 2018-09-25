## Plots 1987-2012 HUC8 decadal trends in annual precip, Mar-May extreme precip, Mar-May temp, and NANI.

library(raster); library(ggplot2); library(dplyr)
#-------------------------------------------------------------------
## Read in HUC8 NANI covariate
nani=readRDS('./data/HUC8/HUC8.nani.rds')
sites.nani=colnames(nani)
fnani=log(nani/2 + ((nani/2)^2+1)^.5) #asinh tranform
n.years=dim(nani)[1]
n.sites=ncol(fnani)
dim(fnani)

#-------------------------------------------------------------------
## Read in precipitation covariates (precip.mam and precip (annual))
precip.annual=readRDS('./data/HUC8/HUC8.precip.rds')
identical(colnames(precip.annual), sites.nani)

precip.mam=readRDS('./data/HUC8/HUC8.precip.mam.rds') #extreme precip
identical(colnames(precip.mam), sites.nani)

dim(precip.mam); dim(precip.annual)

#-------------------------------------------------------------------
## Read in temperature covariate
temp.mam=readRDS('./data/HUC8/HUC8.temp.mam.rds') 
identical(colnames(temp.mam), sites.nani)
dim(temp.mam)  

#-------------------------------------------------------------------
## Read in shapefiles
shape=shapefile('./data/shapefiles/WBDHU8_Reg_1_18.shp') #HUC8
shape.df=fortify(shape) #transforms into format ggplot understands
shape.huc2=readRDS('./data/shapefiles/shapefile.HUC2.rds') 
shape.idx=readRDS('./data/HUC8/HUC8.shapefile.idx.rds') #order of shapefile HUC8 watersheds corresponding to covariate HUC8 ordering

#-------------------------------------------------------------------
## Create function to estimate trend
trend.lm=function(data, return.pval=F){
  year=1:length(data)
  fit=coef(summary(lm(data~year)))
  if (return.pval==T){
    return(fit[2,4])
  } else {
    return(fit[2,1])
  }
}

#-------------------------------------------------------------------
## Estimate covariate trends
trend.precip.annual=apply(precip.annual, 2, trend.lm)
pvals.precip.annual=apply(precip.annual, 2, trend.lm, return.pval=T)
trend.precip.mam=apply(precip.mam, 2, trend.lm)
pvals.precip.mam=apply(precip.mam, 2, trend.lm, return.pval=T)
trend.temp.mam=apply(temp.mam, 2, trend.lm)
pvals.temp.mam=apply(temp.mam, 2, trend.lm, return.pval=T)
trend.nani=apply(nani, 2, trend.lm)
pvals.nani=apply(nani, 2, trend.lm, return.pval=T)

#-------------------------------------------------------------------
## Figure out lat/lon to put dots (center of HUC8 watersheds)
midpoint=function(x){ #calculates midpoint of a vector, for plotting a point over the geographic center of watersheds w/ significant trends
  return((max(x)-min(x))/2+min(x))
}
lat.huc8=aggregate(shape.df$lat, by=list(as.numeric(shape.df$id)), FUN=midpoint)[,2]
lon.huc8=aggregate(shape.df$lon, by=list(as.numeric(shape.df$id)), FUN=midpoint)[,2]
lat.lon.huc8=data.frame(lat=lat.huc8, lon=lon.huc8)

#-------------------------------------------------------------------
## plot plot plot plot plot plot plot plot plot plot plot plot plot
source('plot.huc8.fxn.R') ## Read in plotting fxn 'plot.huc8()'
colorbar=c('#357A7E','#74A4A7','#A8C7C8','#DEE9E8','#feedde','#F9C79E','#F5A25F','#F3801F') #teal->orange v5

## Precip trends (multiply by 10 to plot decadal trends)
cut.seq=seq(-80, 80, by=20)
plot.huc8(plotme=trend.precip.annual*10, filename='./figures/precip.annual.trend.huc8.png', shape.df=shape.df, shape.huc2=shape.huc2, 
          cut.seq=cut.seq, colorbar=colorbar, shape.idx=shape.idx, lat.lon.huc8=lat.lon.huc8, huc8.pvals=pvals.precip.annual)
plot.huc8(plotme=trend.precip.mam*10, filename='./figures/precip.mam.trend.huc8.png', shape.df=shape.df, shape.huc2=shape.huc2, 
          cut.seq=cut.seq, colorbar=colorbar, shape.idx=shape.idx, lat.lon.huc8=lat.lon.huc8, huc8.pvals=pvals.precip.mam)

## Sprimg Temperature trends
cut.seq=seq(-.6, .6, by=.15)
plot.huc8(plotme=trend.temp.mam*10, filename='./figures/temp.mam.trend.huc8.png', shape.df=shape.df, shape.huc2=shape.huc2, 
          cut.seq=cut.seq, colorbar=colorbar, shape.idx=shape.idx, lat.lon.huc8=lat.lon.huc8, huc8.pvals=pvals.temp.mam)

## NANI trends
cut.seq=seq(-400, 400, by=100)
plot.huc8(plotme=trend.nani*10, filename='./figures/nani.trend.huc8.v2.png', shape.df=shape.df, shape.huc2=shape.huc2, 
          cut.seq=cut.seq, colorbar=colorbar, shape.idx=shape.idx, lat.lon.huc8=lat.lon.huc8, huc8.pvals=pvals.nani)


