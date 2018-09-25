## Fits the ultimate six-variable model previously identified through a BIC optimization to observed GAGESII TN fluxes. 
## Includes quality of fit checks (residuals plot, partial residuals plots)

library(car)
min.n=6 #minimum number of TN observations in a given year to qualify for annual flux estimates inclusion
id.remove=c('09498500', '09504000') #Disregard both GAGESII stations in Arizona, large negative NANI due to negative FDFD, not in the original 2016 paper
yr.start=1987; yr.end=2012; n.years=yr.end-yr.start+1

#-------------------------------------------------------------------
## Read in TN station attributes, subsetting by the 123 stations in sites.tn
sites.tn=read.table('./data/GAGESII/gagesII_station_ids.txt', colClasses="character", head=F)[,1]
sites.tn=sites.tn[-which(sites.tn %in% id.remove)] #comment this out of looking at 125 sites
n.sites=length(sites.tn); n.sites

gages.info=read.csv('./data/GAGESII/gagesII_station_info.csv', head=T, colClasses=c("STAID"="character"))
sort.idx=sort(gages.info$STAID, index.return=T)$ix #not actually necessary since basins are already in order, but good just in case
gages.info=gages.info[sort.idx,] #not actually necessary since basins are already in order, but good just in case
gages.info=gages.info[which(gages.info$STAID %in% sites.tn),] #subset just the 123 basins
head(gages.info)

## Read in WRTDS estiamted annual TN flux data, several have missing values for certain years, as well as data pre-1987 and post 2007
tn=matrix(rep(NA, n.years*n.sites), nrow=n.years)
for (i in 1:n.sites){
  site=sites.tn[i]
  dat.tn=read.table(paste0('./data/GAGESII/WRTDS_estimates/',site,'_Annual_TN_Estimated_flux.txt'), head=T)
  idx=dat.tn$Year>=yr.start & dat.tn$Year<=yr.end
  dat.tn=dat.tn[idx,] #remove any years outside of 1987:2012
  area=gages.info$DRAIN_SQKM[gages.info$STAID==site]
  low.n.idx=dat.tn$Count_obs<min.n #only use annual TN values if there are at least 6 obs that year
  dat.tn$Flux[low.n.idx]=NA #only use annual TN values if there are at least 6 obs that year
  tn.site=rep(NA, n.years) 
  tn.site[(yr.start:yr.end %in% dat.tn$Year)]=dat.tn$Flux/area
  tn[,i]=tn.site #store cleaned tn values of i'th station in the 'tn' dataframe
}
colnames(tn)=sites.tn; rownames(tn)=yr.start:yr.end
tn=as.data.frame(tn)
dim(tn) # 26 years x 123 stations

#-------------------------------------------------------------------
## Read in NANI data
nani=readRDS('./data/GAGESII/gagesII_NANI_1987_2012.rds')
nani=nani[,-which(names(nani) %in% id.remove)] #remove two AZ stations, comment out if including
sites.nani=names(nani)

fnani=log(nani/2 + ((nani/2)^2+1)^.5) #inverse hyperbolic sine (approx) of NANI
dim(fnani)

#-------------------------------------------------------------------
## Read in precipitation data, which contains monthly values from 1950-2015, and convert to annual 1987-2007/2012 as needed
precip.annual=matrix(rep(NA, n.years*n.sites), nrow=n.years); precip.mam=precip.annual
for (i in 1:n.sites){
  site=sites.tn[i]
  dat.precip=read.table(paste('./data/GAGESII/precip_temperature/',site,'.precip', sep=""), head=T)
  year=as.numeric(gsub('.{2}$', '', as.character(dat.precip$Date)))
  idx=year>=yr.start & year<=yr.end
  dat.precip=dat.precip[idx,] #remove years outside yr.start:yr.end
  dat.precip=matrix(dat.precip$Precip, nrow=12)
  precip=apply(dat.precip, 2, sum) 
  precip.annual[,i]=precip
  
  ## read in extreme MAM precip 1981-2015
  filename=paste('./data/GAGESII/precip_temperature/',site,'.precip.mam.p95.rds', sep="")
  if (file.exists(filename)){
    dat.pmam=readRDS(filename)
    precip.mam[,i]=dat.pmam$pmam95[dat.pmam$year %in% yr.start:yr.end] #extract only the values from yr.start:yr.end
  }
}  
colnames(precip.annual)=sites.tn; rownames(precip.annual)=yr.start:yr.end
precip=as.data.frame(precip.annual)
colnames(precip.mam)=sites.tn; rownames(precip.mam)=yr.start:yr.end
precip.mam=as.data.frame(precip.mam)

## Quick data check
id.precip.mam=apply(precip.mam, 2, function(x) sum(is.na(x))) #sum of NA values for each station
sum(id.precip.mam==0) #n.stations w/ 0 NA values (should be equal to n.sites)
dim(precip); dim(precip.mam)

#-------------------------------------------------------------------
## Read in temperature data
temp.mam=matrix(rep(NA, 26*n.sites), nrow=26)
colnames(temp.mam)=sites.tn
for (i in 1:n.sites){
  dat.temp=read.table(paste0('./data/GAGESII/precip_temperature/',sites.tn[i],'.tmean'), head=T, row.names=NULL)
  dat.temp=dat.temp[dat.temp$Date %in% c(1987:2012),]
  dat.temp=dat.temp[dat.temp$row.names %in% c('Mar','Apr','May'),]
  dat.mam=colMeans(matrix(dat.temp[,3], nrow=3))
  temp.mam[,i]=dat.mam
}  
dim(temp.mam)  


#-------------------------------------------------------------------
## Read in land use data (luw and luf)
## luw = % classified as wetland; luf = % classifed as forest, shrubland, or herbaceous
lu.dat=read.table('./data/GAGESII/gagesII_landcover.txt', head=T, sep=",", colClasses=c('STAID'="character"))
lu.dat=lu.dat[which(lu.dat$STAID %in% sites.tn),]
## I verified that stationID is in numerical order, as is sites.tn, so that there are not mismatches by mistake with above command
luw=cbind(lu.dat$WOODYWETNLCD06 + lu.dat$EMERGWETNLCD06)
luf=cbind(lu.dat$FORESTNLCD06 + lu.dat$SHRUBNLCD06 + lu.dat$GRASSNLCD06)
length(luw); length(luf)


#-------------------------------------------------------------------
## Now fit model to all 123 watersheds in years 1987,1992,1997,2002,2007,2012
## Subset by those years because NANI values are not technically available in the other years
yr.id=c(1,6,11,16,21,26) #subset for 1987,1992,1997,2002,2007,2012
tn.vector=as.vector(as.matrix(tn[yr.id,]))
fnani.vector=as.vector(as.matrix(fnani[yr.id,]))
precip.vector=as.vector(as.matrix(precip[yr.id,]))
precip.mam.vector=as.vector(as.matrix(precip.mam[yr.id,]))
temp.mam.vector=as.vector(as.matrix(temp.mam[yr.id,]))
luf.vector=as.vector(matrix(rep(luf[], each=length(yr.id)), nrow=length(yr.id)))
luw.vector=as.vector(matrix(rep(luw[], each=length(yr.id)), nrow=length(yr.id)))
df.5yr=data.frame(tn=tn.vector, fnani=fnani.vector, precip=precip.vector, precip.mam=precip.mam.vector, 
                  luf=luf.vector, luw=luw.vector, temp.mam=temp.mam.vector, precipsq=precip.vector^2)
sum(!is.na(tn.vector)) #sample size
fit.87.12=lm(log(tn)~fnani + precip + precip.mam + temp.mam + luw + luf, data=df.5yr, na.action=na.exclude)
summary(fit.87.12)
#saveRDS(fit.87.12, 'TN_model.rds')

#-------------------------------------------------------------------
## Now do some model fit checks
## Basic residuals plot
plot(fit.87.12$fitted.values, fit.87.12$residuals, cex=.7, las=1, ylab='Residuals',
     xlab='Fitted Values (log space)', pch=16, col='orange', ylim=c(-3.2,3.3)) 
abline(h=0, lty=2)

## Basic predicted vs. observed plots in native log space and non-log space
#pdf('./figures/pred.vs.observed.pdf', width=9, height=5, useDingbats=F)
  par(mfrow=c(1,2))
  plot(log(df.5yr$tn[!is.na(df.5yr$tn)]), fit.87.12$fitted.values, cex=.6, las=1, ylab='Predicted TN [log(kg/km2)]',
       xlab='Observed TN [log(kg/km2)]', pch=16, col='orange', xlim=c(0,10), ylim=c(0,10)) 
  abline(0,1, lty=2)
  cor(log(df.5yr$tn[!is.na(df.5yr$tn)]), fit.87.12$fitted.values)^2 #r2
  
  ## Now in non-log space:
  fit.pred=exp(predict.lm(fit.87.12, df.5yr))*exp(sigma(fit.87.12)^2/2) #backtransform
  plot(df.5yr$tn, fit.pred, cex=.6, las=1, ylab='Predicted TN [kg/km2]', 
       xlab='Observed TN [kg/km2]', pch=16, col='orange', xlim=c(0,8000), ylim=c(0,8000))
  abline(0,1, lty=2)
  cor(df.5yr$tn, fit.pred, use='na.or.complete')^2 #r2
#dev.off()

## Partial residual plots to check for potential nonlinear relationships between variables and response, among other things
## General interpretation of each is the effect of each variable on the response after accounting for the impact
## of the other variables. Requires 'car' package.
## Note to get the cex argument below to make a difference you neeed to change the source code
## for crPlots by doing trace(car:::crPlot.lm, edit=T) and adding ",..." into line 69. 
#pdf('./figures/partial.resid.plots.pdf', width=9, height=5, useDingbats=F)
  crPlots(fit.87.12, smooth=F, id=F, lwd=1, ylab='', main='', cex=.7, las=1, layout=c(2,3), pch=16, col='grey45', col.lines='black')
#dev.off()


#-------------------------------------------------------------------
## Predictive accuracy on the out of sample (non AG-census years) TN values
## Results: R2 of .74 on the test set (n=1493)
yr.id.test=c(1:26)[-yr.id] #1988,89,90,91,93,...
tn.vector=as.vector(as.matrix(tn[yr.id.test,]))
fnani.vector=as.vector(as.matrix(fnani[yr.id.test,]))
precip.vector=as.vector(as.matrix(precip[yr.id.test,]))
precip.mam.vector=as.vector(as.matrix(precip.mam[yr.id.test,]))
temp.mam.vector=as.vector(as.matrix(temp.mam[yr.id.test,]))
luf.vector=as.vector(matrix(rep(luf[], each=length(yr.id.test)), nrow=length(yr.id.test)))
luw.vector=as.vector(matrix(rep(luw[], each=length(yr.id.test)), nrow=length(yr.id.test)))
df.test=data.frame(tn=tn.vector, fnani=fnani.vector, precip=precip.vector, precip.mam=precip.mam.vector, 
                  luf=luf.vector, luw=luw.vector, temp.mam=temp.mam.vector, precipsq=precip.vector^2)
sum(!is.na(tn.vector)) #n=1493
#pdf('./figures/pred.vs.observed.testset.pdf', width=9, height=5, useDingbats=F)
  par(mfrow=c(1,2))
  test.pred=predict.lm(fit.87.12, df.test)
  plot(log(df.test$tn), test.pred, cex=.45, las=1, ylab='Predicted TN [log(kg/km2)]',
         xlab='Observed TN [log(kg/km2)]', pch=16, col='orange', xlim=c(0,10), ylim=c(0,10)) 
    abline(0,1, lty=2)
  test.r2=cor(log(df.test$tn), test.pred, use="na.or.complete")^2; test.r2 # R2=.73
  
  ## Now in non-log space:
  test.pred=exp(test.pred)*exp(sigma(fit.87.12)^2/2) #backtransform
  plot(df.test$tn, test.pred, cex=.45, las=1, ylab='Predicted TN [kg/km2]',
         xlab='Observed TN [log(kg/km2)]', pch=16, col='orange', xlim=c(0,8000), ylim=c(0,8000)) 
    abline(0,1, lty=2)
  test.r2=cor(df.test$tn, test.pred, use="na.or.complete")^2; test.r2 # R2=.59
#dev.off()
  
  