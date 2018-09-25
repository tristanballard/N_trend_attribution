## Function for plotting HUC8 values across the US, with HUC2 boundaries overlaid. 
## Saves as a .png to working directory
## Takes 2-3 minutes on a standard macbook due to the large number of points in each shapefile.

## 'plotme' is a vector of values to plat
## 'filename' is the desired output filename (e.g. plot.png)
## 'shape.df' is the dataframe version of HUC8 watersheds
## 'shape.huc2' is the dataframe version of HUC2 watersheds
## 'shape.idx' is a numeric vector indicating the order of 'plotme' values associated with watershed order in 'shape.df' (shape.df is not in ascending HUC8 order)
## 'cut.seq' is the numeric vector of cut points for the colorbar
## 'colorbar' is the character vector of colors
## 'main' is the plot title
## 'legend.on' is TRUE/FALSE whether to include the legend in the plot
## 'legend.lab' is the legend label
## 'lat.lon.huc8' is a 2 column dataframe of lat/lon for the centers of each HUC8 watershed, in same order as in shape.df
## 'huc8.pvals' is a vector of pvalues. Only those less than 0.1 are plotted as points
## You will receive, and can ignore, the following warning: In `[<-.factor`(`*tmp*`, ri, value = c(3000, 3000, 3000, 3000, 3000,:invalid factor level, NA generated

plot.huc8=function(plotme, filename='', shape.df=shape.df, shape.huc2='', shape.idx='', cut.seq='', 
                   colorbar='', main='', legend.on=F, legend.lab='', lat.lon.huc8=NA, huc8.pvals=NA){
    n.sites=length(plotme)
    ## Set some plotting parameters
    legend.min=min(cut.seq); legend.max=max(cut.seq);  n.col=length(cut.seq)-1; 
    legend.midpoints=diff(cut.seq)/2+cut.seq[-length(cut.seq)]
    TN.plot=plotme[shape.idx]
    TN.plot[TN.plot<=legend.min]=legend.min+.01; TN.plot[TN.plot>=legend.max]=legend.max-.01
    
    ## Clean up dataframes for plotting w/ ggplot
    dat=cbind((0:(n.sites-1)), TN.plot); colnames(dat)=c("id", "TN")
    dat=as.data.frame(dat); dat[,1]=as.character(dat[,1]); ss=left_join(shape.df, dat) #merges datasets by the 'id' variable
    ss.append=as.data.frame(tail(ss,n.col)); ss.append$group=3000; ss.append$TN=legend.midpoints
    ss=rbind(ss, ss.append) #warnings OK. ensures that the plotter gets to use all the colors you specify
    if(legend.on==F){
      map = ggplot() + geom_polygon(data=ss, aes(x=long, y=lat, group=group, fill=cut(TN, cut.seq)), color='grey', size=.05) +
          scale_fill_manual(values=colorbar, name=legend.lab, na.value='purple', guide = FALSE) +
          geom_polygon(data=shape.huc2, aes(x=long, y=lat, group=group), fill=NA, color='black', size=.25) +
          theme_classic() + theme(axis.line=element_blank(),axis.text.x=element_blank(), axis.text.y=element_blank(),axis.ticks=element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank())
      } else {
      map = ggplot() + geom_polygon(data=ss, aes(x=long, y=lat, group=group, fill=cut(TN, cut.seq)), color='grey', size=.05) +
          scale_fill_manual(values=colorbar, name=legend.lab, na.value='purple', guide = guide_legend(reverse = TRUE)) +
          labs(title = main, fill = "") +
          geom_polygon(data=shape.huc2, aes(x=long, y=lat, group=group), fill=NA, color='black', size=.25) +
          theme_classic() + theme(axis.line=element_blank(),axis.text.x=element_blank(), axis.text.y=element_blank(),axis.ticks=element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank())
      }
     if(sum(!is.na(lat.lon.huc8))>0){
      huc8.pvals=p.adjust(huc8.pvals, method='BH') # Multiple comparisons correction (Benjamini-Hochberg)
      lat.lon.sig=lat.lon.huc8[huc8.pvals[shape.idx]<0.1,] # Identify which have p-values less than 0.1
      map = map + geom_point(data=lat.lon.sig, aes(x=lon, y=lat), size=.03) #adds dots over HUC8 watersheds w/ significant trends
     }
    
    ggsave(map, file=filename, width=7.4, height=5, type = "cairo-png")
}
