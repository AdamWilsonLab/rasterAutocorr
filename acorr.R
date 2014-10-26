
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

acorr=function(x,padlongitude=T,verbose=T,...){
  # dimensions of raster
  nr <- nrow(x)
  nc <- ncol(x)
  ## FFT needs mulitple of 8
  ## Images must be padded to size 2N-1 by 2M-1
  # find the closest multiple of 8 to obtain a good compromise between
  # speed (a power of 2) and memory required
  nr2=ifelse(nr<5,5,ceiling((2*nr-1)/8)*8) #for elev 1072 rows instead of 536 #why doubling?..why not 1024 (2^10)
  nc2=ifelse(!padlongitude,nc,ifelse(nc<5,5,ceiling((2*nc-1)/8)*8))  
  ## create a new extent by padding to the right and below with 0s
  resx=res(x)
  extx=extent(x)
  #create a new extent object with new ajdusted rows and columns
  extx2=extent(c(xmin=extx@xmin,xmax=extx@xmax+(resx[1]*(nc2-nc)),ymin=extx@ymin-(resx[2]*(nr2-nr)),ymax=extx@ymax))
  if(verbose) print("Padding the array")
  rx=extend(x,extx2,val=0) #this replaces the expand function to  avoid conflict with matrix packgae
  ## convert to matrix
  x1=as(rx,"matrix") #make it a matrix,  this may become quite big...
  # make an indicator matrix with 1's for all data values & O's for missing values
  if(verbose) print("Identifying missing data")
  xnull=as(extend(!is.na(x),extx2,val=0),"matrix")
  # in data matrix, replace missing values by 0;
  x1[xnull==0]=0

  if(verbose) print("Running the initial FFTs")  
  fx1=fft(x1)  # fourier transform of x1, calling C code implementation of Fortran  Singleton   1979
  #fft is a 
  fx1_x1=fft(x1*x1)    # fourier transform of x1*x1, this is very fast! why is this here?
  
  fxnull=fft(xnull)  # fourier transform of the indicator matrix? ok 1,0
  # compute number of pairs at all lags
  if(verbose) print("Computing the number of observations at each lag")
  nobs=round(Re(ifft(Conj(fxnull)*fxnull))) 
  mnobs=nobs
  mnobs[mnobs<1]=1
  ## compute the correlogram
  m1=Re(ifft(Conj(fx1)*fxnull))/mnobs #is this the normalization?
  m2=Re(ifft(Conj(fxnull)*fx1))/mnobs # I'm confused...
  PSD <- Conj(fx1)*fx1 #compute hte square of amplitude/norm of FFT complex number which correspond o the amplitude of frequencies
  PSD_r <-raster(Re(PSD)) 
  plot(PSD_r) #plotting the image power spectrum....mm not working
  #sqrt(PSD_r)
  plot(sqrt((PSD_r)) #ok taking the amplitude
  test<-ifft(PSD) 
  plot(raster(Re(test))) #plotting the autocorrelation.. range is huge!!! 250,000 ah ok need to normalize?
  g=Re(ifft(Conj(fx1)*fx1)/mnobs-m1*m2) #Ok is this the calucation autocorrelation surface??
  g_r <-raster(g)
  plot(g_r)
  #http://faculty.olin.edu/bstorey/Notes/Fourier.pdf
  
  if(verbose) print("Shifting the FFT array") # why?
  
  nobs2=fftshift2(nobs)    #also in this package...  
  g2=fftshift2(g)*10  #why 10?    

  ## get distances in km
  if(verbose) print("Calculating distances")
  d1=acorr_dist(rx) #also created for this pack

  # convert back to raster
  if(verbose) print("Convert back to raster* format")
  g3=d1;values(g3)=g2
  nobs3=d1;values(nobs3)=log10(nobs2)
  acor=stack(g3,nobs3,d1)
  names(acor)=c("acor","nobs","dist")
#  if(exists("filename",inherits=F)) acor2=writeRaster(acor,...)
  if(verbose) print("Cleaning up")
  rm(x,rx,x1,fx1,fx1_x1,fxnull,m1,m2,g,nobs,g3,d1,nobs3);gc()
  return(acor)
}


## Is this from John Loomis

#where auto is the following Matlab file

#function f = auto(a)
#%
#% AUTO
#%   f = auto(a)
#%	 returns as f the auto-correlation of image a.

#f = fft2(a);
#f = f.*conj(f);
#f = ifft2(f);
#f = fftshift(f);
#f = real(f);
#f = f/max(max(f));

#http://www.johnloomis.org/ece563/notes/freq/autoself/autoself.htm