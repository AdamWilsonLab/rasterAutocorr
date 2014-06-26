
#' Calculates the distance of each cell in a raster to the center pixel to faciltate comparison of results from \code{acorr} as a function of distance.
#'
#' @param x A raster* object
#' @return A raster object showing the distances from the center pixel
#' 

acorr_dist=function(x,...){
  x2=acorr_center(x,...)
  ## distance units need to be updated from raster somehow...
  ## currently it assumes x is in lat-lon and so distance returns meters
  ## divide by 1000 to km
  dist=distance(x2)/1000
  return(dist)
}

