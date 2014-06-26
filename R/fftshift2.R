#' Rearranges FFT output to put zero-distance in center of image
#' @description Rearranges outputs from fft to move the zero-frequency component to the center of the matrix.  This is useful to visualize a Fourier transform with the zero-frequency component in the center of the image.
#' @param x A matrix returned from \code{fft()}
#' @return The transformed matrix
#' @references \url{http://www.mathworks.com/help/matlab/ref/fftshift.html}
#' @references \url{http://stackoverflow.com/questions/5735720/effcient-way-to-do-fft-shift-in-matlab-without-using-fftshift-function}
#' @references The waved package (\url{http://cran.r-project.org/web/packages/waved/index.html}) has a 1-dimensional fftshift function.
#' 

fftshift2=function(x){
  nr=nrow(x)
  nc=ncol(x)
  ## get values to add (accounts for odd and even rows)
  nre=ifelse(nr/2==round(nr/2),1,1.5)
  nce=ifelse(nc/2==round(nc/2),1,1.5)
  x2=x[
    c(((nr/2)+nre):(nr), 1:(1+nr/2)),
    c(((nc/2)+nce):(nc), 1:(1+nc/2))]
  return(x2)
}


