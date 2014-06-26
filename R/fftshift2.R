#' Rearranges FFT output to put zero-distance in center of image
#' @description Rearranges outputs from fft to move the zero-frequency component to the center of the matrix.  This is useful to visualize a Fourier transform with the zero-frequency component in the center of the image.
#' @param x A matrix returned from \code{fft()}
#' @return The transformed matrix
#' @references \url{http://www.mathworks.com/help/matlab/ref/fftshift.html}
#' @references \url{http://stackoverflow.com/questions/5735720/effcient-way-to-do-fft-shift-in-matlab-without-using-fftshift-function}
#' @references The waved package (\url{http://cran.r-project.org/web/packages/waved/index.html}) has a 1-dimensional fftshift function.
#' 

fftshift2=function(x){
  sz = floor(dim(x)/2)
       x2=x[
      c((sz[1]+1):(sz[1]*2), 1:sz[1]),
      c((sz[2]+1):(sz[2]*2), 1:sz[2])]
       return(x2)
}


function (y) 
{
  n = length(y)
  n2 = floor(n/2)
  ind1 = (1:n2)
  ind2 = ((n2 + 1):n)
  aux1 = (y[ind1])
  aux2 = y[ind2]
  y = c(aux2, aux1)
}