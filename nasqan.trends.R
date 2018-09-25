## This begins by reading in modeled and observed TN fluxes for the NASQAN stations and sees if
## there is a significant difference between modeled and observed annual trends at each basin.
## It then plots the modeled and observed trends as a time series for each basin.

library(png)
## Read in observed Basin fluxes
basin.TN.obs=readRDS('./data/NASQAN/nasqan.TN.observed.rds')
basin.meta=readRDS('./data/NASQAN/nasqan.metadata.rds')
basin.names=basin.meta$name
n.basins=nrow(basin.meta); n.basins

## Read in estimated Basin fluxes
basin.TN.model=readRDS('./data/NASQAN/nasqan.TN.rds')
dim(basin.TN.obs); dim(basin.TN.model)
identical(basin.names, colnames(basin.TN.model)) #TRUE = good

## Change estimated fluxes to have NA's corresponding to NA observations
basin.TN.model[is.na(basin.TN.obs)]=NA

## Test for similarity of trends between modeled and observed 
## Test for significnce of interaction term for group(observed v. modeled) against year
## Interaction term asks: Does the effect of year on TN (i.e. the TN trend) differ between groups 0 (observed) and 1 (modeled)
n.years=26
out.p=rep(NA, n.basins)
for (i in 1:n.basins){
  y=c(basin.TN.model[,i],basin.TN.obs[,i])
  year=rep(1:n.years, 2)
  group=c(rep(1,n.years), rep(0,n.years))
  fit=lm(y~year + group + year*group)
  out.p[i]=coef(summary(fit))[4,4] #store p-values of the interaction term
}
sum(out.p<.1) #number of interaction terms w/ p-values less than 0.1 (n=0 of 10)
basin.names[which(out.p<.1)] #not applicable since all have p-values above 0.1
## When running multiple tests you should do a multiple comparisons correction as below, but not in this case necessary
## because none are significant anyways, and the correction will only shift p-values higher.
#p.mcc=p.adjust(out.p, method='BH') #Benjamini-Hochberg correction for multiple comparisons
#sum(p.mcc<.1)

## Extract linear trend function 
trend.lm=function(data, return.pval=F){ #extracts the trend
  year=1:length(data)
  fit=coef(summary(lm(data~year)))
  if (return.pval==T){
    return(fit[2,4])
  } else {
    return(fit[2,1])
  }
}


## Plotting function
## A lot of extra plotting options have been commented out
plot.two.ts=function(obs.ts, model.ts, ylab=expression('kg/km'^'2'), legend.on=TRUE){
    yr.lab=c(1987, rep(NA, 4),1992, rep(NA, 4),1997, rep(NA, 4),2002, rep(NA, 4),2007, rep(NA, 4), 2012)
    plot.ts(model.ts, ylab=ylab, xlab='', xaxt='n', las=1, lwd=2, col='deepskyblue4', bty='n',
    ylim=c(0, max(c(model.ts, obs.ts), na.rm=T)*1.4))
    axis(1, at=1:26, labels=yr.lab)
    lines(obs.ts, col='deepskyblue', lwd=2)
    abline(lm(obs.ts~c(1:26)), lty=2, col='deepskyblue') #'cadetblue3'
    abline(lm(model.ts~c(1:26)), lty=2, col='deepskyblue4') #'cadetblue4'

    r2=round(cor(model.ts,obs.ts, use='na.or.complete')^2, 2)
    if(legend.on==TRUE){
    #  legend('topleft', c('','Modeled','Observed', 'Trend'),lty=c(0,1,1,2), lwd=c(1,2,2,1), seg.len=1.75, col=c('','deepskyblue4', 'deepskyblue', 'black'), bty='n', cex=.8)
      legend('bottomright', legend=c(as.expression(bquote(R^2 == .(r2)))), lty=0, lwd=1, bty='n', cex=.9)
    }
}

pdf(file='./figures/nasqan.obs.modeled.comparison.pdf', width=7.7, height=35)
  par(mfrow=c(n.basins,2))
  for (i in 1:n.basins){
    img=readPNG(paste0('./data/NASQAN/b',basin.names[i],'.map.png')) 
    plot(0, type='n', xlim=c(0,.8), ylim=c(0,.8), bty='n', axes=F, xlab='', ylab='', main=basin.names[i])
    rasterImage(img, 0, 0, .8, .8) # map of the basin
    plot.two.ts(basin.TN.obs[,i], basin.TN.model[,i]) # time series + trend lines of the basin
  }
dev.off()

