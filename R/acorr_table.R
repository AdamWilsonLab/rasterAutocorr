acor_table=function(x,maxdist=1500){
  ## run the autocorrelation function and write out the output raster
  wrapglobe=!(extent(x)@xmin==-180&extent(x)@xmax==180)  #does this region wrap the globe?
  ac=acorr(x,padlongitude=wrapglobe)
  ## build the table of values to construct the correlograms
  ftd=rbind.data.frame(
    data.frame(values=values(ac[["acor"]])/10,dist=values(ac[["dist"]]),n=values(ac[["nobs"]]))
  )
  ## filter to a reasonable distance, signal gets noisy above a few thousand km due to sparse measurements
  ftd <- filter(ftd, dist <= maxdist)
  ## normalize the covariogram to a correlogram by dividing value at lag=0
  ftd$values=ftd$values/max(ftd$values[which.min(ftd$dist)])
  ## round to approximate resolution of raster (in km)
  rasres=round(rasterRes(x))
  rasbins=c(-1,seq(0,max(ftd$dist)+rasres,by=rasres))
  ftd$dist2=as.numeric(as.character((cut(ftd$dist,rasbins,labels=rasres/2+rasbins[-length(rasbins)]))))
  ## take mean by distance bin
  ftd2 <- group_by(ftd, dist2,type)
  ftd2 <- summarise(ftd2,
                    min = min(values, na.rm = TRUE),
                    max = max(values, na.rm = TRUE),
                    sd = sd(values, na.rm = TRUE),
                    mean = mean(values, na.rm = TRUE)
  )
  return(ftd2)
}