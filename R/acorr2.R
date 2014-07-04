
#' Calculates the full spatial autocorrelation on a raster using fft.
#'
#' @description Applies the Wiener-Khinchin theorem to extract spatial autocorrelation using Fast Fourier Transform techniques.  This results in an extremely fast way to calculate a complete correlogram (correlation as a function of distance) for a raster image.  
#' @param x A raster* object. Missing values are indicated by NA.
#' @param file File to write results to as in writeRaster.  If NULL a temporary file is written as in the raster package.
#' @return The spatial autocorrelation matrix
#' @references \url{en.wikipedia.org/wiki/WienerKhinchin_theorem}
#' @references Xianlin Ma, Tingting Yao, A program for 2D modeling (cross) correlogram tables using fast Fourier transform, Computers & Geosciences, Volume 27, Issue 7, August 2001, Pages 763-774, ISSN 0098-3004, \url{http://dx.doi.org/10.1016/S0098-3004(01)00007-3}.
#' @references \url{http://www.johnloomis.org/ece563/notes/freq/autoself/autoself.htm}
#' @references \url{http://www.seas.upenn.edu/~ese502/NOTEBOOK/Part_II/4_Variograms.pdf}
#' @example examples/examples.R

acorr2=function(x,padlongitude=T,...){
  # dimensions of raster
  nr <- nrow(x)
  nc <- ncol(x)
  ## Images must be padded to size 2N-1 by 2M-1
  # find the closest multiple of 8 to obtain a good compromise between
  # speed (a power of 2) and memory required
  nr2=ifelse(nr<5,5,ceiling((2*nr-1)/8)*8)
  nc2=ifelse(!padlongitude,nc,ifelse(nc<5,5,ceiling((2*nc-1)/8)*8))  
  ## create a new extent by padding to the right and below with 0s
  resx=res(x)
  extx=extent(x)
  extx2=extent(c(xmin=extx@xmin,xmax=extx@xmax+(resx[1]*(nc2-nc)),ymin=extx@ymin-(resx[2]*(nr2-nr)),ymax=extx@ymax))
  rx=extend(x,extx2,val=0)
  ## convert to matrix
  x1=as(rx,"matrix")
  # make an indicator matrix with 1's for all data values & O's for missing values
  xnull=matrix(0,nrow=nr2,ncol=nc2)
  xnull[1:nr,1:nc][as.matrix(!is.na(x))]=1
  # in data matrix, replace missing values by 0;
  x1[xnull==0]=0

  fx1=fft(x1)  # fourier transform of xl
  fx1_x1=fft(x1*x1)    # fourier transform of x1*x1
  
  fxnull=fft(xnull)  # fourier transform of the indicator matriJ
  # compute number of pairs at all lags
  nobs=round(Re(ifft(Conj(fxnull)*fxnull)))
  mnobs=nobs
  mnobs[mnobs<1]=1
  ## compute the correlogram
  m1=Re(ifft(Conj(fx1)*fxnull))/mnobs
  m2=Re(ifft(Conj(fxnull)*fx1))/mnobs
  g=Re(ifft(Conj(fx1)*fx1)/mnobs-m1*m2)
  
nobs2=fftshift2(nobs)      
g2=fftshift2(g)*10      

## get distances in km
d1=acorr_dist(rx)

# convert back to raster
  g3=d1;values(g3)=g2
  nobs3=d1;values(nobs3)=log10(nobs2)
  acor=stack(g3,nobs3,d1)
  names(acor)=c("acor","nobs","dist")
#  if(exists("filename",inherits=F)) acor2=writeRaster(acor,...)
  rm(x,rx,x1,fx1,fx1_x1,fxnull,m1,m2,g,nobs,g3,d1,nobs3);gc()
  return(acor)
}