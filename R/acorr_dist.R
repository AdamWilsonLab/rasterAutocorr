
#' Calculates the distance of each cell in a raster to the center pixel to faciltate comparison of results from \code{acorr} as a function of distance.
#'
#' @param x A raster* object
#' @return A raster object showing the distances from the center pixel
#' 

acorr_dist=function(x){
  x2=acorr_center(x)
  ## distance units need to be updated from raster somehow...
  dist=distance(x2)*res(x)[1]
  return(dist)
}

