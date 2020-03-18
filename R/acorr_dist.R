
#' Calculates the distance of each cell in a raster to the center pixel to faciltate comparison of results from \code{acorr} as a function of distance.
#'
#' @param x A raster* object
#' @return A raster object showing the distances from the center pixel in the original units
#' 

acorr_dist=function(x){
  x2=acorr_center(x)
  dist=distance(x2)
  return(dist)
}

