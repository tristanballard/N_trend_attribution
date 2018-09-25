## This reads in the estimated 1987-2012 HUC8 TN fluxes using the updated model 
## and using the original Sinha et al model and plots maps of the mean, standard 
## deviation, and coefficient of variation for each, as well as difference maps 
## in raw value and in percent.

library(ggplot2); library(raster); library(dplyr); library(RColorBrewer)

#-------------------------------------------------------------------
## Read in 1987-2012 TN fluxes estimated using updated and original (Sinha2016) model
TN=readRDS('./data/HUC8/HUC8.TN.rds')
TN.2016=readRDS('./data/HUC8/HUC8.TN.Sinha2016.rds')

## Read in shapefiles
shape=shapefile('./data/shapefiles/WBDHU8_Reg_1_18.shp') #HUC8
shape.df=fortify(shape) #transforms into format ggplot understands
shape.idx=readRDS('./data/HUC8/HUC8.shapefile.idx.rds') #order of shapefile HUC8 watersheds corresponding to covariate HUC8 ordering
shape.huc2=readRDS('./data/shapefiles/shapefile.HUC2.rds') 

#-------------------------------------------------------------------
## Calculate HUC8 mean, standard deviation, coefficient of variation
tn.mean=apply(TN, 2, mean)
tn.mean.2016=apply(TN.2016, 2, mean)
tn.sd=apply(TN, 2, sd)
tn.sd.2016=apply(TN.2016, 2, sd)
tn.cov=tn.sd/tn.mean
tn.cov.2016=tn.sd.2016/tn.mean.2016

#-------------------------------------------------------------------
## plot plot plot plot plot plot plot plot plot plot plot plot plot 
source('plot.huc8.fxn.R') ## Read in plotting fxn 'plot.huc8()'

## Mean and standard deviation
cut.seq=c(0,25,50,100,200,400,600,1000,1800,2600,3000) #break points for the colorbar
colorbar=c('#fdfce0','#fcf9ac','#fde078','#feb92a','#f88608','#e85c14','#c0391e','#852514','#521d08','#424242') # light yellow to darkred/grey, from Sinha2016
plot.huc8(plotme=tn.mean, filename='./figures/tn.huc8.mean.png', shape.df=shape.df, shape.huc2=shape.huc2, shape.idx=shape.idx, cut.seq=cut.seq, colorbar=colorbar)
plot.huc8(plotme=tn.mean.2016, filename='./figures/tn.huc8.mean.sinha2016.png', shape.df=shape.df, shape.huc2=shape.huc2, shape.idx=shape.idx, cut.seq=cut.seq, colorbar=colorbar)
plot.huc8(plotme=tn.sd, filename='./figures/tn.huc8.sd.png', shape.df=shape.df, shape.huc2=shape.huc2, shape.idx=shape.idx, cut.seq=cut.seq, colorbar=colorbar)
plot.huc8(plotme=tn.sd.2016, filename='./figures/tn.huc8.sd.sinha2016.png', shape.df=shape.df, shape.huc2=shape.huc2, shape.idx=shape.idx, cut.seq=cut.seq, colorbar=colorbar)

## Coefficient of Variation (sd / mean)
cut.seq=c(0,.25,.5,1,4)
colorbar=c('#C1C1C1','#878888', '#535454','#252525') #same greys from Sinha2016
plot.huc8(plotme=tn.cov, filename='./figures/tn.huc8.cov.png', shape.df=shape.df, shape.huc2=shape.huc2, shape.idx=shape.idx, cut.seq=cut.seq, colorbar=colorbar)
plot.huc8(plotme=tn.cov.2016, filename='./figures/tn.huc8.cov.sinha2016.png', shape.df=shape.df, shape.huc2=shape.huc2, shape.idx=shape.idx, cut.seq=cut.seq, colorbar=colorbar)

## Difference
cut.seq=c(-2600,-1800,-1000,-600,-400,-200,-100,-50,-25,0,25,50,100,200,400,600,1000,1800,2600,3000)
colorbar=c(rev(blues9),'#fdfce0','#fcf9ac','#fde078','#feb92a','#f88608','#e85c14','#c0391e','#852514','#521d08','#424242') #blues to yellow to darkred/grey
plot.huc8(plotme=tn.mean-tn.mean.2016, filename='./figures/tn.huc8.mean.difference.png', shape.df=shape.df, shape.huc2=shape.huc2, shape.idx=shape.idx, cut.seq=cut.seq, colorbar=colorbar)
plot.huc8(plotme=tn.sd-tn.sd.2016, filename='./figures/tn.huc8.sd.difference.png', shape.df=shape.df, shape.huc2=shape.huc2, shape.idx=shape.idx, cut.seq=cut.seq, colorbar=colorbar)
plot.huc8(plotme=tn.cov-tn.cov.2016, filename='./figures/tn.huc8.cov.difference.png', shape.df=shape.df, shape.huc2=shape.huc2, shape.idx=shape.idx, cut.seq=cut.seq, colorbar=colorbar)

## Percent difference
cut.seq=seq(-100,100, length=11)
colorbar=c(rev(brewer.pal(5, 'Blues')), brewer.pal(5, 'Reds'))
plot.huc8(plotme=(tn.mean-tn.mean.2016)/tn.mean.2016*100, filename='./figures/tn.huc8.mean.difference.percent.png', shape.df=shape.df, shape.huc2=shape.huc2, shape.idx=shape.idx, cut.seq=cut.seq, colorbar=colorbar)
plot.huc8(plotme=(tn.sd-tn.sd.2016)/tn.sd.2016*100, filename='./figures/tn.huc8.sd.difference.percent.png', shape.df=shape.df, shape.huc2=shape.huc2, shape.idx=shape.idx, cut.seq=cut.seq, colorbar=colorbar)
plot.huc8(plotme=(tn.cov-tn.cov.2016)/tn.cov.2016*100, filename='./figures/tn.huc8.cov.difference.percent.png', shape.df=shape.df, shape.huc2=shape.huc2, shape.idx=shape.idx, cut.seq=cut.seq, colorbar=colorbar)

