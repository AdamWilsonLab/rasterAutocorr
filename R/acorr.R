
#' Calculates the full spatial autocorrelation on a raster using fft.
#'
#' @param x A raster* object
#' @param fmean function to use to calculate the overall mean of the raster.  For small objects, "mean" is fine, but for larger rasters it is recommended to use "CellStats"
#' @return The spatial autocorrelation matrix
#' @example examples/examples.R


acorr=function(x,fmean=c("mean","CellStats")){
  ## convert to matrix
  xm=as(x,"matrix")
  ## subtract the mean of the object
  if(fmean=="mean") tmean=mean(xm)
  if(fmean=="CellStats") tmean=cellStats(x,mean)
  xm=xm-tmean
  ## should we pad the image here?

  ## take the fft of the matrix
  fftx=fft(xm)
  fftx2=Re(fft(fftx* Conj(fftx), inverse=TRUE))
  ## shift the matrix to bring the low frequency autocorrelation to the center of the image
  acor1=fftshift2(fftx2)
  # Scale to get values in [0,1]
  acor2=acor1/max(acor1)
  ## create a raster object to fill with the new values
  dims=dim(acor2)
  acor2=raster(acor2,xmn=-dims[2]/2,xmx=dims[2]/2,ymn=-dims[1]/2,ymx=dims[1]/2)
  return(acor2)
}
