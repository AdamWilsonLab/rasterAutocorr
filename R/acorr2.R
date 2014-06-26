
#' Calculates the full spatial autocorrelation on a raster using fft.
#'
#' @description Applies the Wiener-Khinchin theorem to extract spatial autocorrelation using Fast Fourier Transform techniques.  This results in an extremely fast way to calculate a complete correlogram (correlation as a function of distance) for a raster image.  
#' @param x A raster* object. Missing values are indicated by NA.
##' @param fmean function to use to calculate the overall mean of the raster.  For small objects, "mean" is fine, but for larger rasters it is recommended to use "CellStats"
#' @param file File to write results to as in writeRaster.  If NULL a temporary file is written as in the raster package.
#' @return The spatial autocorrelation matrix
#' @example examples/examples.R
#' @references \url{http://en.wikipedia.org/wiki/WienerKhinchin_theorem}
#' @references Xianlin Ma, Tingting Yao, A program for 2D modeling (cross) correlogram tables using fast Fourier transform, Computers & Geosciences, Volume 27, Issue 7, August 2001, Pages 763-774, ISSN 0098-3004, \url{http://dx.doi.org/10.1016/S0098-3004(01)00007-3}.
#' @references \url{http://www.johnloomis.org/ece563/notes/freq/autoself/autoself.htm}

#r2=r
#values(r2)[sample(1:ncell(r2),10000)]=NA
#x=r2
#x=m1

acorr2=function(x,gain=NULL,...){
  ## convert to matrix
  x1=as(x,"matrix")
  ## pad the data
  # dimensions of data matrix
  nr <- nrow(x)
  nc <- ncol(x)
  ## Images must be padded to size 2N-1 by 2M-1
  # find the closest multiple of 8 to obtain a good compromise between
  # speed (a power of 2) and memory required
  nr2=ceiling((2*nr-1)/8)*8
  nc2=ceiling((2*nc-1)/8)*8  
  x2=matrix(0,nrow=nr2,ncol=nc2)
  x2[1:nr,1:nc]=x1

  # form an indicator matrix:
  # 1's for all data values
  # O's for missing values
  # in data matrix, replace missing values by 0;
  xnull=matrix(0,nrow=nr2,ncol=nc2)
  xnull[1:nr,1:nc][!is.na(x1)]=1
  x2[xnull==0]=0

  fx1=fft(x2)  # fourier transform of xl
  fx1_x1=fft(x2*x2)    # fourier transform of x1*x1
  
  fxnull=fft(xnull)  # fourier transform of the indicator matriJ
  # compute number of pairs at all lags
  nobs=round(Re(ifft(Conj(fxnull)*fxnull)))

  # compute the autocorrelation
  g=Re(ifft(Conj(fxnull)*fx1_x1+Conj(fx1_x1)*fxnull-2*Conj(fx1)*fx1))
  g=g/max(nobs,1)/2

  #rm(fx1_x1,fxnull,fx1)     ## clean up

  ## shift the matrix to bring the low frequency autocorrelation to the center of the image
  g2=fftshift2(g)
  nobs2=fftshift2(nobs)
#  rm(g,nobs);gc()
  ## crop to original dimensions
  nre=ifelse(nrow(g2)/2==round(nrow(g2)/2),1,2)
  nce=ifelse(ncol(g2)/2==round(ncol(g2)/2),1,2)

nobs2=rbind(
  cbind(nobs[1:nr,1:nc],          nobs[1:nr,(nc2-nc+2):nc2]),
  cbind(nobs[(nr2-nr+2):nr2,1:nc],nobs[(nr2-nr+2):nr2,(nc2-nc+2):nc2])
)
nobs3=fftshift2(nobs2)      

g2=rbind(
  cbind(g[1:nr,1:nc],          g[1:nr,(nc2-nc+2):nc2]),
  cbind(g[(nr2-nr+2):nr2,1:nc],g[(nr2-nr+2):nr2,(nc2-nc+2):nc2])
)
g3=fftshift2(g2)      

#g3=g2[(nre+(nr2-nr)/2):(nr+((nr2-nr)/2)),(nce+(nc2-nc)/2):(nc+(nc2-nc)/2)]  
#nobs3=nobs2[(nre+(nr2-nr)/2):(nr+((nr2-nr)/2)),(nce+(nc2-nc)/2):(nc+(nc2-nc)/2)]

## get distances in km
d1=acorr_dist(nobs3,nc=nc,nr=nr)

# convert to raster
  g4=raster(g3,xmn=-nc/2,xmx=nc/2,ymn=-nr/2,ymx=nr/2)
  nobs4=raster(nobs3,xmn=-nc/2,xmx=nc/2,ymn=-nr/2,ymx=nr/2)
  acor=stack(g4,nobs4,d1)
  names(acor)=c("acor","nobs","dist")
  if(exists("filename",inherits=F)) acor2=writeRaster(,...)
#  rm(x,xm);gc()
  return(acor)
}