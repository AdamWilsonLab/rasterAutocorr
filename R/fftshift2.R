#' Rearranges FFT output to put zero-distance in center of image
#' @description Rearranges outputs from fft to move the zero-frequency component to the center of the matrix.  This is useful to visualize a Fourier transform with the zero-frequency component in the center of the image.
#' @param x A matrix returned from \code{fft()}
#' @return The transformed matrix
#' @references Adapted from GNU Octave's fftshift function (\url{http://octave.sourceforge.net/octave/function/fftshift.html} \url{http://hg.savannah.gnu.org/hgweb/octave/file/914c0b103a3d/scripts/signal/fftshift.m#l76})
#' @references \url{http://www.mathworks.com/help/matlab/ref/fftshift.html}
#' @references \url{http://stackoverflow.com/questions/5735720/effcient-way-to-do-fft-shift-in-matlab-without-using-fftshift-function}
#' @references The waved package (\url{http://cran.r-project.org/web/packages/waved/index.html}) has a 1-dimensional fftshift function.

fftshift2=function(x){
  nd = length(dim(x))
  sz = dim(x)
  sz2 = ceiling(sz/2);
  idx = list()
  for (i in 1:nd)  idx[[i]] = c((sz2[i]+1):sz[i], 1:sz2[i])
  retval = x[idx[[1]],idx[[2]]];
  return(retval)
}
