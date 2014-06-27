#' Use Octave's \code{xcorr2} to get spatial autocorrelation of complete image
#'
#' @param x A raster* object
#' @return A raster object of the autocorrelation of the image
#' @references \url{http://octave.sourceforge.net/signal/function/xcorr2.html}
#' 

acorr_xcorr2=function(x){
#  require(RcppOctave)
  o_eval("pkg load signal")
  t(.O$xcorr2(x,"coeff"))
}

