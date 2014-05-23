
#' Finds the center pixel to faciltate comparison of results from \code{acorr}.
#'
#' @param x A raster* object
#' @return A raster object with all cells NA except the center pixel.
#' 

acorr_center=function(x){
  center = ceiling(dim(x)/2)
  center2=x
  center2[,]=NA
  center2[center[1]+1,center[2]+1]=1
  return(center2)
}

