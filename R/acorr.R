
#' Calculates the full spatial autocorrelation on a raster using fft.
#'
#' @param x A raster* object
#' @return The spatial autocorrelation matrix
#' @example examples/examples.R


acorr=function(x){
  ## convert to matrix
  xm=as.matrix(x)
  ## subtract the mean of the object
  xm=xm-mean(xm)
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
  acor2=raster(acor2,xmn=-dims[1]/2,xmx=dims[1]/2,ymn=-dims[2]/2,ymx=dims[2]/2)
  return(acor2)
}
