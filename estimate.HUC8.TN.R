## Reads in fitted model and applies it to the HUC8 watersheds to get 1987-2012 annual 
## TN flux estimates for all continental US HUC8 watersheds. It also similarly predicts 
## fluxes using the Sinha et al. 2016 model.

#-------------------------------------------------------------------
## Read in the model obect
model=readRDS('TN_model.rds')
model.sd=sigma(model)
model.2016=readRDS('TN_model_Sinha2016.rds') #slightly different names used for variables
model.2016.sd=sigma(model.2016)

#-------------------------------------------------------------------
## Read in HUC8 NANI, NDEP, FERT
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
idx.keep=which(sites.lu %in% sites.nani)
sites.lu=sites.lu[idx.keep]; luw=luw[idx.keep]; luf=luf[idx.keep]
sort.idx=sort(sites.lu, index.return=T)$ix
luw=luw[sort.idx]; luf=luf[sort.idx]

#-------------------------------------------------------------------
## Estimate HUC8 TN 
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
colnames(TN)=sites.nani
rownames(TN)=1987:2012
#saveRDS(TN, './data/HUC8/HUC8.TN.rds')

## Now make predictions on same years using the Sinha et al model, used for comparison in 'new.old.model.comparison.R'
yhat.2016=predict(model.2016, newdata=list(precip=df.all$precip, R95pTOT_MAM_1=df.all$precip.mam,
                                           wetland=df.all$luw, forest_shrubgrass=df.all$luf,
                                           asinh_NANI_0=df.all$fnani)) #output is logTN 
yhat.2016=matrix(yhat.2016, nrow=n.years)
TN.2016=exp(yhat.2016)*exp(model.2016.sd^2/2)
colnames(TN.2016)=sites.nani
rownames(TN.2016)=1987:2012
#saveRDS(TN.2016, './data/HUC8/HUC8.TN.Sinha2016.rds')
  
