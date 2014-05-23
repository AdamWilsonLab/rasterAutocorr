
#' Calculates the direction of each cell in a raster to the center pixel to faciltate comparison of results from \code{acorr} as a function of distance.
#'
#' @param x A raster* object
#' @return A raster object showing the direction (in degrees) from the center pixel
#' 

acorr_dir=function(x){
  x2=acorr_center(x)
  dir=direction(x2,degrees=T)
  return(dir)
}

