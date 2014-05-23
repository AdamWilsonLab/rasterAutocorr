
#' Finds the center pixel to faciltate comparison of results from \code{acorr}.
#'
#'
#' @description Performs the following functions. 1) create an empty raster to hold the output, 2) fill the raster with NAs, 3) Set the center pixel equal to 1.
#' This faciliates using other functions in the \pkg{raster} package that calculate distance and direction to the closest non-NA pixel.  This is useful because the output of \code{\link{autocor}} returns the spatial autocorrelation for all directions simutaneously and thus needs to be linked to distance and direction to generate a correlogram.
#' @param x A raster* object
#' @return A raster object with all cells NA except the center pixel.
#' 


acorr_center=function(x){
  center = ceiling(dim(x)/2)
  # Create an empty raster to hold the output
  center2=raster(stack(x),layer=0)  
  # Fill the raster with NAs
  center2[,]=NA
  # Replace the center pixel with a one  
  center2[center[1]+1,center[2]+1]=1
  return(center2)
}

