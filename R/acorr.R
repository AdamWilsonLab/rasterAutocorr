
#' Calculates the full spatial autocorrelation on a raster using fft.
#'
#' @description Applies the Wiener-Khinchin theorem to extract spatial autocorrelation using Fast Fourier Transform techniques.  This results in an extremely fast way to calculate a complete correlogram (correlation as a function of distance) for a raster image.  
#' @param x A raster* object
##' @param fmean function to use to calculate the overall mean of the raster.  For small objects, "mean" is fine, but for larger rasters it is recommended to use "CellStats"
#' @param file File to write results to as in writeRaster.  If NULL a temporary file is written as in the raster package.
#' @param gain number to multiply the values by prior to writing to disk.  If NULL, the original values in [0,1] will be written.  Note that the output file's metadata does not store this gain value.
#' @return The spatial autocorrelation matrix
#' @example examples/examples.R
#' @references http://en.wikipedia.org/wiki/Wiener%E2%80%93Khinchin_theorem
#' @references Xianlin Ma, Tingting Yao, A program for 2D modeling (cross) correlogram tables using fast Fourier transform, Computers & Geosciences, Volume 27, Issue 7, August 2001, Pages 763-774, ISSN 0098-3004, \url{http://dx.doi.org/10.1016/S0098-3004(01)00007-3}.
#' @references http://www.johnloomis.org/ece563/notes/freq/autoself/autoself.htm
#' @references http://en.wikipedia.org/wiki/Wiener%E2%80%93Khinchin_theorem



acorr=function(x,gain=NULL,...){
  ## convert to matrix
  xm=as(x,"matrix")
  ## subtract the mean of the object
  tmean=mean(xm,na.rm=T)
  xm=xm-tmean
  ## fill missing values with mean
  xm[is.na(xm)]=tmean
  ## should we pad the image here?
  ## create padded version
  #xpad=as.matrix(extend(tcld, extent(c(-39,90,-20,30)), value=0))
  ## take the fft of the matrix
  fftx=fft(xm)
  fftx2=Re(fft(fftx* Conj(fftx), inverse=TRUE))
  ## shift the matrix to bring the low frequency autocorrelation to the center of the image
  acor1=fftshift2(fftx2)
  # Scale to get values in [0,1]
  tmax=max(acor1,na.rm=T)
  acor2=acor1/tmax
  ## create a raster object to fill with the new values
  dims=dim(acor2)
  ## multiply by gain to enable storing as integer
  if(!is.null(gain)) acor2=acor2*gain
  acor2=raster(acor2,xmn=-dims[2]/2,xmx=dims[2]/2,ymn=-dims[1]/2,ymx=dims[1]/2)
  if(exists("filename",inherits=F)) acor2=writeRaster(acor2,...)
  rm(xm,fftx,fftx2,acor1);gc()
  return(acor2)
}
