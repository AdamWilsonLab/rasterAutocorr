#' Calculates the full spatial autocorrelation on a raster using fft.
#'
#' @description Pads an image with zeros to facilitate accurate FFT.  
#' @param x A matrix object
#' @param padX logical, should the matrix be padded in the x (longitude) direction? Set to FALSE if matrix includes full -180:180 latitudes to incorporate wrapping.
#' @param padY logical, should the matrix be padded in the y (latitude) direction?
#' @return The padded matrix
#' @references Marcotte, Denis. 1996. “Fast Variogram Computation with FFT.” Computers & Geosciences 22 (10): 1175–86.  \url{http://dx.doi.org://10.1016/S0098-3004(96)00026-X}


pad=function(x,padX=T,padY=T){
  # dimensions of data matrix
  nr <- nrow(x)
  nc <- ncol(x)
  ## Images must be padded to size 2N-1 by 2M-1
  # find the closest multiple of 8 to obtain a good compromise between
  # speed (a power of 2) and memory required
  nr2=ceiling((2*nr-1)/8)*8
  nc2=ceiling((2*nc-1)/8)*8  
  if(padX&padY){
    ## make the new matrix
    nx=matrix(0,nrow=nr2,ncol=nc2)
    nx[c((ceiling((nr2-nr)/2)):(ceiling((nr2+nr)/2)-1)),c(ceiling((nc2-nc)/2)):(ceiling((nc2+nc)/2)-1)]=x
  }
  if(padX&!padY){
    ## make the new matrix
    nx=matrix(0,nrow=nr,ncol=nc2)
    nx[1:nr,c((nc/2):(nc+nc/2-1))]=x
  }
  if(!padX&padY){
    ## make the new matrix
    nx=matrix(0,nrow=nr2,ncol=nc)
    nx[c((nr/2):(nr+nr/2-1)),1:nc]=x
  }
  if(!padX&!padY){
  warning("No padding requested, returning original matrix")
  nx=x
  }  
    return(nx)
}