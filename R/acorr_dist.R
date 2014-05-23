
#' Calculates the distance of each cell in a raster to the center pixel to faciltate comparison of results from \code{acorr} as a function of distance.
#'
#' @param x A raster* object
#' @return A raster object showing the distances from the center pixel
#' 

acorr_dist=function(x){
  center=c(round(nrow(x)/2),round(ncol(x)/2))
  dist=x
  dist[,]=NA
  dist[center[1],center[2]]=1
  dist=distance(dist)*res(r)[1]
  return(dist)
}

