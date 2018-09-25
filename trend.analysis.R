## Estimates HUC8 TN, aggregates to HUC2 scale, and then calculates and plots the HUC2 trends. 
## Runs various trend attribution simulations with detrended covariates and plots the resulting HUC2 trends.
library(ggplot2); library(raster); library(dplyr)

#-------------------------------------------------------------------
## Read in the model obect
model=readRDS('TN_model.rds')
model.sd=sigma(model)

#-------------------------------------------------------------------
## Read in HUC8 NANI
nani=readRDS('./data/HUC8/HUC8.nani.rds')
sites.nani=colnames(nani)
fnani=log(nani/2 + ((nani/2)^2+1)^.5) #asinh tranform
n.years=dim(nani)[1]
n.sites=ncol(fnani)
dim(fnani)

#-------------------------------------------------------------------
## Read in precipitation predictors (precip.mam and precip (annual))
precip.annual=readRDS('./data/HUC8/HUC8.precip.rds')
identical(colnames(precip.annual), sites.nani)

precip.mam=readRDS('./data/HUC8/HUC8.precip.mam.rds') #extreme precip
identical(colnames(precip.mam), sites.nani)

dim(precip.mam); dim(precip.annual)

#-------------------------------------------------------------------
## Read in temperature predictors
temp.mam=readRDS('./data/HUC8/HUC8.temp.mam.rds') 
identical(colnames(temp.mam), sites.nani)
dim(temp.mam)  

#-------------------------------------------------------------------
## Read in land use predictors
dat.lu=read.csv('./data/HUC8/HUC8_landuse.csv', head=T, colClasses = c('HUC8'='character'))
sites.lu=dat.lu[,1]
luw=cbind(dat.lu$WOODYWETNLCD06 + dat.lu$EMERGWETNLCD06)
luf=cbind(dat.lu$FORESTNLCD06 + dat.lu$SHRUBNLCD06 + dat.lu$GRASSNLCD06)
idx.keep=which(sites.lu %in% sites.nani) #remove HUC8s not needed
sites.lu=sites.lu[idx.keep]; luw=luw[idx.keep]; luf=luf[idx.keep]
sort.idx=sort(sites.lu, index.return=T)$ix
luw=luw[sort.idx]; luf=luf[sort.idx] #sort land use into correct order

#-------------------------------------------------------------------
## Read in shapefiles
shape=shapefile('./data/shapefiles/WBDHU8_Reg_1_18.shp') #HUC8
shape.df=fortify(shape) #transforms into format ggplot understands

shape.huc2=readRDS('./data/shapefiles/shapefile.HUC2.rds') 
huc2=shape$HUC2
huc2.area=aggregate(shape$AreaSqKm, by=list(huc2), FUN=sum)$x
n.huc2=length(huc2.area)

shape.idx=readRDS('./data/HUC8/HUC8.shapefile.idx.rds') #order of shapefile HUC8 watersheds corresponding to covariate HUC8 ordering
huc2.names=as.character(aggregate(shape$Name_HUC2, by=list(huc2), FUN=unique)$x)

#-------------------------------------------------------------------
## Create function to detrend predictors
detrend.lm=function(data){ #linear regression detrend
      year=1:length(data)
      fit=lm(data~year)
      return(fit$residuals+mean(data))
}

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
## Detrend individual covariates
fnani.detrended=apply(fnani, 2, detrend.lm)
precip.mam.detrended=apply(precip.mam, 2, detrend.lm)
precip.annual.detrended=apply(precip.annual, 2, detrend.lm)
temp.mam.detrended=apply(temp.mam, 2, detrend.lm)

#-------------------------------------------------------------------
## Individual trend contributions

## Function to output the decadal trend (per km2) for each HUC2:
trend.huc2=function(TN=TN, shape.idx=shape.idx, shape=shape, huc2=huc2, return.pvals=F){
    ## Aggregate HUC8 estimates to HUC2 level
    TN.adj=TN[,shape.idx]
    TN.total=TN.adj*matrix(rep(shape$AreaSqKm, each=nrow(TN)), nrow=nrow(TN)) #26 x 2107
    huc2.area=aggregate(shape$AreaSqKm, by=list(huc2), FUN=sum)$x
    TN.huc2=t(unlist(apply(TN.total, 1, function(x) aggregate(x, by=list(huc2), FUN=sum)[,2]))) #26 x 18regions
    ## Calculate HUC2 trends
    trends.tn<-pvals.tn<-rep(NA, n.huc2)
    for (i in 1:n.huc2){
      trends.tn[i]=trend.lm(TN.huc2[,i])
      pvals.tn[i]=trend.lm(TN.huc2[,i], return.pval=T)
    }
    if (return.pvals==T){
      return(pvals.tn)
    }
    return(trends.tn=trends.tn*10/huc2.area) # multiply by 10 to get decadal trends
}


## Calculate predicted HUC2 TN trend
df.all=data.frame(fnani=as.vector(as.matrix(fnani)),
              precip=as.vector(as.matrix(precip.annual)),
              precip.mam=as.vector(as.matrix(precip.mam)),
              precipsq=as.vector(as.matrix(precip.annual^2)),
              temp.mam=as.vector(as.matrix(temp.mam)),
              luw=as.vector(matrix(rep(luw, each=n.years))),
              luf=as.vector(matrix(rep(luf, each=n.years))))

yhat=predict(model, df.all) #output is logTN 
yhat=matrix(yhat, nrow=n.years)
TN=exp(yhat)*exp(model.sd^2/2)
    TN.adj=TN[,shape.idx]
    TN.total=TN.adj*matrix(rep(shape$AreaSqKm, each=nrow(TN)), nrow=nrow(TN)) #26 x 2107
    huc2.area=aggregate(shape$AreaSqKm, by=list(huc2), FUN=sum)$x
    TN.huc2=t(unlist(apply(TN.total, 1, function(x) aggregate(x, by=list(huc2), FUN=sum)[,2]))) #26 x 18regions
    mean.tn=apply(TN.huc2, 2, mean) / huc2.area #HUC2 mean TN
    mean.tn.huc8=apply(TN, 2, mean) #HUC8 mean TN
trends.tn=trend.huc2(TN=TN, shape.idx=shape.idx, shape=shape, huc2=huc2)
pvals.tn=trend.huc2(TN=TN, shape.idx=shape.idx, shape=shape, huc2=huc2, return.pvals=T)
trends.tn.huc8=apply(TN, 2, trend.lm)
pvals.tn.huc8=apply(TN, 2, trend.lm, return.pval=T)

## Hold everything detrended except temp.mam  
df.temp.mam=data.frame(fnani=as.vector(as.matrix(fnani.detrended)),
              precip=as.vector(as.matrix(precip.annual.detrended)),
              precip.mam=as.vector(as.matrix(precip.mam.detrended)),
              precipsq=as.vector(as.matrix(precip.annual.detrended^2)),
              temp.mam=as.vector(as.matrix(temp.mam)),
              luw=as.vector(matrix(rep(luw, each=n.years))),
              luf=as.vector(matrix(rep(luf, each=n.years))))

yhat=predict(model, df.temp.mam) #output is logTN 
yhat=matrix(yhat, nrow=n.years)
TN=exp(yhat)*exp(model.sd^2/2)
trends.tn.t.mam=trend.huc2(TN=TN, shape.idx=shape.idx, shape=shape, huc2=huc2)
pvals.tn.t.mam=trend.huc2(TN=TN, shape.idx=shape.idx, shape=shape, huc2=huc2, return.pvals=T)
trends.tn.t.mam.huc8=apply(TN, 2, trend.lm)
pvals.tn.t.mam.huc8=apply(TN, 2, trend.lm, return.pval=T)

## Hold everything detrended except precip.annual  
df.precip.annual=data.frame(fnani=as.vector(as.matrix(fnani.detrended)),
              precip=as.vector(as.matrix(precip.annual)),
              precip.mam=as.vector(as.matrix(precip.mam.detrended)),
              precipsq=as.vector(as.matrix(precip.annual^2)),
              temp.mam=as.vector(as.matrix(temp.mam.detrended)),
              luw=as.vector(matrix(rep(luw, each=n.years))),
              luf=as.vector(matrix(rep(luf, each=n.years))))

yhat=predict(model, df.precip.annual) #output is logTN 
yhat=matrix(yhat, nrow=n.years)
TN=exp(yhat)*exp(model.sd^2/2)
trends.tn.p.annual=trend.huc2(TN=TN, shape.idx=shape.idx, shape=shape, huc2=huc2)
pvals.tn.p.annual=trend.huc2(TN=TN, shape.idx=shape.idx, shape=shape, huc2=huc2, return.pvals=T)
trends.tn.p.annual.huc8=apply(TN, 2, trend.lm)
pvals.tn.p.annual.huc8=apply(TN, 2, trend.lm, return.pval=T)

## Hold everything detrended except precip.mam  
df.precip.mam=data.frame(fnani=as.vector(as.matrix(fnani.detrended)),
              precip=as.vector(as.matrix(precip.annual.detrended)),
              precip.mam=as.vector(as.matrix(precip.mam)),
              precipsq=as.vector(as.matrix(precip.annual.detrended^2)),
              temp.mam=as.vector(as.matrix(temp.mam.detrended)),
              luw=as.vector(matrix(rep(luw, each=n.years))),
              luf=as.vector(matrix(rep(luf, each=n.years))))

yhat=predict(model, df.precip.mam) #output is logTN 
yhat=matrix(yhat, nrow=n.years)
TN=exp(yhat)*exp(model.sd^2/2)
trends.tn.p.mam=trend.huc2(TN=TN, shape.idx=shape.idx, shape=shape, huc2=huc2)
pvals.tn.p.mam=trend.huc2(TN=TN, shape.idx=shape.idx, shape=shape, huc2=huc2, return.pvals=T)
trends.tn.p.mam.huc8=apply(TN, 2, trend.lm)
pvals.tn.p.mam.huc8=apply(TN, 2, trend.lm, return.pval=T)

## Hold everything detrended except both precip variables
df.precip=data.frame(fnani=as.vector(as.matrix(fnani.detrended)),
             precip=as.vector(as.matrix(precip.annual)),
             precip.mam=as.vector(as.matrix(precip.mam)),
             precipsq=as.vector(as.matrix(precip.annual^2)),
             temp.mam=as.vector(as.matrix(temp.mam.detrended)),
             luw=as.vector(matrix(rep(luw, each=n.years))),
             luf=as.vector(matrix(rep(luf, each=n.years))))

yhat=predict(model, df.precip) #output is logTN 
yhat=matrix(yhat, nrow=n.years)
TN=exp(yhat)*exp(model.sd^2/2)
trends.tn.p=trend.huc2(TN=TN, shape.idx=shape.idx, shape=shape, huc2=huc2)
pvals.tn.p= trend.huc2(TN=TN, shape.idx=shape.idx, shape=shape, huc2=huc2, return.pvals=T)
trends.tn.p.huc8=apply(TN, 2, trend.lm)
pvals.tn.p.huc8=apply(TN, 2, trend.lm, return.pval=T)

## Hold everything detrended except fNANI 
df.nani=data.frame(fnani=as.vector(as.matrix(fnani)),
              precip=as.vector(as.matrix(precip.annual.detrended)),
              precip.mam=as.vector(as.matrix(precip.mam.detrended)),
              precipsq=as.vector(as.matrix(precip.annual.detrended^2)),
              temp.mam=as.vector(as.matrix(temp.mam.detrended)),
              luw=as.vector(matrix(rep(luw, each=n.years))),
              luf=as.vector(matrix(rep(luf, each=n.years))))

yhat=predict(model, df.nani) #output is logTN 
yhat=matrix(yhat, nrow=n.years)
TN=exp(yhat)*exp(model.sd^2/2)
trends.tn.fnani=trend.huc2(TN=TN, shape.idx=shape.idx, shape=shape, huc2=huc2)
pvals.tn.fnani=trend.huc2(TN=TN, shape.idx=shape.idx, shape=shape, huc2=huc2, return.pvals=T)
trends.tn.fnani.huc8=apply(TN, 2, trend.lm)
pvals.tn.fnani.huc8=apply(TN, 2, trend.lm, return.pval=T)

## Hold everything detrended except all climate variables (detrend only fnani)
df.temp.clim=data.frame(fnani=as.vector(as.matrix(fnani.detrended)),
              precip=as.vector(as.matrix(precip.annual)),
              precip.mam=as.vector(as.matrix(precip.mam)),
              precipsq=as.vector(as.matrix(precip.annual^2)),
              temp.mam=as.vector(as.matrix(temp.mam)),
              luw=as.vector(matrix(rep(luw, each=n.years))),
              luf=as.vector(matrix(rep(luf, each=n.years))))

yhat=predict(model, df.temp.clim) #output is logTN 
yhat=matrix(yhat, nrow=n.years)
TN=exp(yhat)*exp(model.sd^2/2)
trends.tn.clim=trend.huc2(TN=TN, shape.idx=shape.idx, shape=shape, huc2=huc2)
pvals.tn.clim=trend.huc2(TN=TN, shape.idx=shape.idx, shape=shape, huc2=huc2, return.pvals=T)
trends.tn.clim.huc8=apply(TN, 2, trend.lm)
pvals.tn.clim.huc8=apply(TN, 2, trend.lm, return.pval=T)

  
#-------------------------------------------------------------------
## Function to plot trends at HUC2 level
plot.huc2=function(trendz, file.out){
  TN.plot=rep(NA, n.sites)
  for (i in 1:n.huc2){
    TN.plot[as.numeric(huc2)==i]=trendz[i]
  }
  cut.seq=c(-100,-75,-50,-25,-5,0,5,25,50,75,100)
  colorbar=c('#357A7E','#74A4A7','#A8C7C8','#DEE9E8','white','white','#feedde','#F9D3B0','#F5A25F','#F3801F') #teal->orange 
  legend.min=min(cut.seq); legend.max=max(cut.seq);  n.col=length(cut.seq)-1; 
  legend.midpoints=diff(cut.seq)/2+cut.seq[-length(cut.seq)] 
  TN.plot[TN.plot<=legend.min]=legend.min+.01; TN.plot[TN.plot>=legend.max]=legend.max-.01

  dat=cbind((0:(n.sites-1)), TN.plot); colnames(dat)=c("id", "TN")
  dat=as.data.frame(dat); dat[,1]=as.character(dat[,1])
  ss=left_join(shape.df, dat) #merges datasets by the 'id' variable
  ss.append=as.data.frame(tail(ss,n.col)); ss.append$group=3000; ss.append$TN=legend.midpoints
  ss=rbind(ss, ss.append) #warning "In `[<-.factor`(`*tmp*`, ri, value = c(3000,..." is OK. Ensures that the plotter gets to use all the colors you specify
  map = ggplot() + geom_polygon(data=ss, aes(x=long, y=lat, group=group, fill=cut(TN, cut.seq))) +
      scale_fill_manual(values=colorbar, name=expression(paste('kg/km'^'2','/decade', sep='')), na.value='orange', guide=FALSE) +#guide = guide_legend(reverse = TRUE)) +
      labs(title = '', fill = "") +
      geom_polygon(data=shape.huc2, aes(x=long, y=lat, group=group), fill=NA, color='black', size=.25) +
      theme_classic() + theme(axis.line=element_blank(),axis.text.x=element_blank(), axis.text.y=element_blank(),axis.ticks=element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank())
  ggsave(map, file = file.out, width=7.4, height=5, type = "cairo-png")
}

#-------------------------------------------------------------------
#-------------------------------------------------------------------
#-------------------------------------------------------------------
## plot plot plot plot plot plot plot plot plot plot plot plot plot 
source('plot.huc8.fxn.R') ## Read in plotting fxn 'plot.huc8()'

## Figure out lat/lon to put dots (center of HUC8 watersheds)
midpoint=function(x){ #calculates midpoint of a vector, for plotting a point over the geographic center of watersheds w/ significant trends
  return((max(x)-min(x))/2+min(x))
}
lat.huc8=aggregate(shape.df$lat, by=list(as.numeric(shape.df$id)), FUN=midpoint)[,2]
lon.huc8=aggregate(shape.df$lon, by=list(as.numeric(shape.df$id)), FUN=midpoint)[,2]
lat.lon.huc8=data.frame(lat=lat.huc8, lon=lon.huc8)

## Plot huc8 tn trend
cut.seq=c(-100,-75,-50,-25,-5,0,5,25,50,75,100)
colorbar=c('#357A7E','#74A4A7','#A8C7C8','#DEE9E8','white','white','#feedde','#F9D3B0','#F5A25F','#F3801F') #teal->orange 
plot.huc8(plotme=trends.tn.huc8*10, filename='./figures/tn.trend.huc8.png', shape.df=shape.df, shape.huc2=shape.huc2, 
          shape.idx=shape.idx, cut.seq=cut.seq, colorbar=colorbar, lat.lon.huc8=lat.lon.huc8, huc8.pvals=pvals.tn.huc8)
# plot.huc8(plotme=trends.tn.fnani.huc8*10, filename='./figures/tn.trend.from.fnani.huc8.png', shape.df=shape.df, shape.huc2=shape.huc2, 
#           shape.idx=shape.idx, cut.seq=cut.seq, colorbar=colorbar, lat.lon.huc8=lat.lon.huc8, huc8.pvals=pvals.tn.fnani.huc8)
# plot.huc8(plotme=trends.tn.p.annual.huc8*10, filename='./figures/tn.trend.from.pannual.huc8.png', shape.df=shape.df, shape.huc2=shape.huc2, 
#           shape.idx=shape.idx, cut.seq=cut.seq, colorbar=colorbar, lat.lon.huc8=lat.lon.huc8, huc8.pvals=pvals.tn.p.annual.huc8)
# plot.huc8(plotme=trends.tn.p.mam.huc8*10, filename='./figures/tn.trend.from.pmam.huc8.png', shape.df=shape.df, shape.huc2=shape.huc2, 
#           shape.idx=shape.idx, cut.seq=cut.seq, colorbar=colorbar, lat.lon.huc8=lat.lon.huc8, huc8.pvals=pvals.tn.p.mam.huc8)
# plot.huc8(plotme=trends.tn.t.mam.huc8*10, filename='./figures/tn.trend.from.tmam.huc8.png', shape.df=shape.df, shape.huc2=shape.huc2, 
#           shape.idx=shape.idx, cut.seq=cut.seq, colorbar=colorbar, lat.lon.huc8=lat.lon.huc8, huc8.pvals=pvals.tn.t.mam.huc8)
# plot.huc8(plotme=trends.tn.p.huc8*10, filename='./figures/tn.trend.from.both.precip.huc8.png', shape.df=shape.df, shape.huc2=shape.huc2, 
#           shape.idx=shape.idx, cut.seq=cut.seq, colorbar=colorbar, lat.lon.huc8=lat.lon.huc8, huc8.pvals=pvals.tn.p.huc8)
# plot.huc8(plotme=trends.tn.clim.huc8*10, filename='./figures/tn.trend.from.climate.huc8.png', shape.df=shape.df, shape.huc2=shape.huc2, 
#           shape.idx=shape.idx, cut.seq=cut.seq, colorbar=colorbar, lat.lon.huc8=lat.lon.huc8, huc8.pvals=pvals.tn.clim.huc8)

#-------------------------------------------------------------------
## Plot various HUC2 trends
plot.huc2(trends.tn, file.out="./figures/tn.trend.huc2.png")
plot.huc2(trends.tn.fnani, file.out="./figures/tn.trend.from.fnani.huc2.png")
plot.huc2(trends.tn.p.annual, file.out="./figures/tn.trend.from.pannual.huc2.png")
plot.huc2(trends.tn.p.mam, file.out="./figures/tn.trend.from.pmam.huc2.png")
plot.huc2(trends.tn.t.mam, file.out="./figures/tn.trend.from.tmam.huc2.png")
plot.huc2(trends.tn.p, file.out="./figures/tn.trend.from.both.precip.huc2.png")
plot.huc2(trends.tn.clim, file.out="./figures/tn.trend.from.climate.huc2.png")

