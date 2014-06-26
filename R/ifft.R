#' Calculates the normalized inverse fft of an array
#'
#' @description Convenience function to calculates the normalized inverse fft of an array (like in MatLab)
#' @param x An array as used by fft()
#' @return The normalized inverse fft of x
#' @references \url{http://stackoverflow.com/questions/8162562/is-there-a-package-in-r-that-gives-normalized-inverse-fft}


ifft <- function(x) { fft(x, inverse=TRUE ) / length(x) }
