SCRIPTS:
All scripts run assuming you have set the working directory (see setwd() command in R) as the one containing the scripts. Change as needed. Some packages, mentioned at the top of each script, may need to be installed. Code implemented with R version 3.5.1.

covariate.trends.R: Plots 1987-2012 HUC8 decadal trends in annual precip, Mar-May extreme precip, Mar-May temp, and NANI.

estimate.HUC8.TN.R: Reads in fitted model and applies it to the HUC8 watersheds to get 1987-2012 annual TN flux estimates for all continental US HUC8 watersheds. It also similarly predicts fluxes using the Sinha et al. 2016 model, with the TN output comparison conducted in new.old.model.comparison.R.

estimate.nasqan.TN.R: Reads in HUC8 TN estimates and aggregates them to the NASQAN basin level to get annual TN estimates for each NASQAN basin.

fit.model.R: Fits the six variable model (fit.87.12) identified through BIC and explores various aspects of its fit (e.g. residuals, predicted vs. observed plots). This model is used later on to estimate fluxes at HUC8 watersheds.

map.gage.locations.R: Creates a map of the 123 gage locations used in the study along with their basin boundaries. It also overlays the HUC2 and HUC8 watershed boundaries.

nani.component.trend.barplots.R: Estimates HUC2 trends in each of the 5 NANI components and summarizes the results in a barplot.

nasqan.trends.R: Tests for a significant difference between modeled and observed annual TN trends at each basin and plots the associated time series.

new.old.model.comparison.R: Reads in the estimated 1987-2012 HUC8 TN fluxes using the updated model and using the original Sinha et al model and plots maps of the mean, standard deviation, and coefficient of variation for each, as well as difference maps in raw value and in percent.

trend.analysis.R: Estimates HUC8 TN, aggregates to HUC2 scale, and then calculates and plots the HUC2 trends. Runs various trend attribution simulations with detrended covariates and plots the resulting HUC2 and HUC8 trends.


FUNCTIONS/OBJECTS:
plot.huc8.rxn.R: Function for plotting HUC8 values across the US, with HUC2 boundaries overlaid. 

TN_model.rds: Fitted model R object created in fit.model.R.

TN_model_Sinha2016.rds: Fitted model R object corresponding to Sinha et al. 2016 final model.



DATA: 
Note that many of the station IDs in the files begin with '0', e.g. '01111500', which can accidentally be converted to '1111500' if converting formats/copy and pasting. Blame USGS.


DATA/GAGESII:
gagesII_landcover.txt: land use classifications [%] for each GAGESII watershed from the NCLD 2006 dataset. 

gagesII_lat_lon.csv: lat and lon coordinates of all GAGES-II outlets/monitoring locations. Same information is found in gagesII_station_info.csv.

gagesII_NANI_1987_2012.rds: 1987-2012 NANI values for the 125 GAGES-II stations.

gagesII_station_info.csv: various attributes of the GAGES-II stations such as station name, drainage basin area, lat, and lon. 

gagesII_station_ids.txt: GAGES-II station ids for the 125 sites used in this analysis. Two of the sites in Arizona were ultimately removed due to anomalous NANI values.

gagesII_station_ids_Sinha2016.txt: GAGES-II station id's for the 70 sites used in Sinha et al. 2016.

DATA/GAGESII/WRTDS_ESTIMATES:
This folder contains the WRTDS TN flux estimates for each GAGES-II station in kg (not kg/km2) along with a column denoting the actual number of TN measurements taken in each year.

DATA/GAGESII/PRECIP_TEMPERATURE:
This folder contains the GAGESII PRISM monthly mean temperature [degC], monthly mean precip [mm], and annual amount of extreme precip [mm], defined as the sum of precip on days where the precip exceeded the  95th percentile of historical (1981-2010) daily precipitation with the percentile threshold calculated only with days having at least 1mm of precip to remove 'dry' days. All of the values have been aggregated over each GAGESII basin area. 



DATA/HUC8:
HUC8.landuse.csv: land use classifications for each HUC8 watershed from the 2006 NCLD database.

HUC8.nani.rds: 1987-2012 annual NANI [kg/km2] for each HUC8 watershed

HUC8.precip.rds: 1987-2012 annual total precipitation [mm] for each HUC8 watershed

HUC8.precip.mam.rds: 1987-2012 annual extreme Mar-May precipitation [mm] for each HUC8 watershed.

HUC8.NANI.components.txt: 1987-2012 annual values for NANI and its five components [kg/km2] for each HUC8 watershed in the continental U.S. See manuscript for a description of each of component.

HUC8.shapefile.idx.rds: vector indicating the order of HUC8 watersheds in the shapefile corresponding to the ascending order of the HUC8 sites used elsewhere in the analysis. Important for plotting correctly. Created in nani.component.trend.barplots.R

HUC8.TN.rds: 1987-2012 annual TN fluxes [kg/km2] estimated by applying the TN_model.rds model to annual HUC8 covariate inputs. 

HUC8.TN.Sinha2016.rds: 1987-2012 annual TN fluxes [kg/km2] estimated by applying the model from Sinha et al. 2016 to annual HUC8 covariate inputs. 



DATA/NASQAN:
b____.png (and others): Map showing the location of this NASQAN station. Created in estimate.nasqan.TN.R script.

nasqan.metadata.rds: Attributes of the 10 NASQAN watersheds (e.g. area, long name)

nasqan.TN.rds: 1987-2012 annual TN fluxes [kg/km2] estimated by aggregating HUC8 TN estimates from HUC8.TN.rds for those watersheds making up each basin.

nasqan.TN.observed.rds: 1987-2012 observed (WRTDS-estimated) annual TN fluxes [kg/km2].


DATA/SHAPEFILES:
GAGESII_basins.shp: shapefiles for the 123 GAGES sites.

shapefile.HUC2.rds: 18 continental US HUC2 shapefiles. 

WBDHU8_Reg_1_18.shp: 2,107 HUC8 watershed boundary shapefiles and various info about each watershed (e.g. basin area, full name).

DATA/SHAPEFILES/NASQAN: shapefiles for various NASQAN watersheds.




