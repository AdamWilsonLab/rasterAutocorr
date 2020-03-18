#' Calculates correlogram values on a raster using fft.
#'
#' @description Applies the Wiener-Khinchin theorem to extract spatial autocorrelation using Fast Fourier Transform techniques.  This results in an extremely fast way to calculate a complete correlogram (correlation as a function of distance) for a raster image.  
#' @param x A raster* object. Missing values are indicated by NA.
#' @param maxdist Maximum distance (in km) to include in correlogram table.  All possible distances are calculated using FFT, then trimmed to this value.
#' @return The table of distances
#' @references \url{en.wikipedia.org/wiki/WienerKhinchin_theorem}
#' @references Xianlin Ma, Tingting Yao, A program for 2D modeling (cross) correlogram tables using fast Fourier transform, Computers & Geosciences, Volume 27, Issue 7, August 2001, Pages 763-774, ISSN 0098-3004, \url{http://dx.doi.org/10.1016/S0098-3004(01)00007-3}.
#' @references \url{http://www.johnloomis.org/ece563/notes/freq/autoself/autoself.htm}
#' @references \url{http://www.seas.upenn.edu/~ese502/NOTEBOOK/Part_II/3_Spatially_Dependent_Random_Effects.pdf}
#' @export
#' @import raster dplyr



acor_table=function(x,maxdist=1500000,verbose=F){
  ## run the autocorrelation function and write out the output raster
  wrapglobe=!(extent(x)@xmin==-180&extent(x)@xmax==180)  #does this region wrap the globe?
  ac=acorr(x,padlongitude=wrapglobe,verbose=verbose)
  ## build the table of values to construct the correlograms
  if(verbose) print("Extracting autocorrelation values for table and filtering to maxdist")
  ftd=rbind.data.frame(
    data.frame(values=values(ac[["acor"]])/10,
               dist=values(ac[["dist"]]),
               n=values(ac[["nobs"]]))
  )
  ## filter to a reasonable distance, signal gets noisy above a few thousand km due to sparse measurements
  ftd <- filter(ftd, dist <= maxdist)
  ## normalize the covariogram to a correlogram by dividing value at lag=0
  ftd$values=ftd$values/max(ftd$values[which.min(ftd$dist)])
  ## round to approximate resolution of raster (in km)
  rasres=round(rasterRes(x))
  rasbins=c(-1,seq(0,max(ftd$dist)+rasres,by=rasres))
  ftd$dist=as.numeric(as.character((cut(ftd$dist,rasbins,labels=rasres/2+rasbins[-length(rasbins)]))))
  ## take mean by distance bin
  if(verbose) print("Summarizing by distance bin")
  ftd2 <- group_by(ftd, dist)
  ftd2 <- summarise(ftd2,
                    min = min(values, na.rm = TRUE),
                    max = max(values, na.rm = TRUE),
                    sd = sd(values, na.rm = TRUE),
                    mean = mean(values, na.rm = TRUE)
  )
  return(ftd2)
}