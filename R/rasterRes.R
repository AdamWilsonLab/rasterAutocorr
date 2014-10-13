rasterRes=function(x){
  nr=nrow(x)
  nc=ncol(x)
##  extract middle portion of raster
  x2=x[round(nr/2-1):round(nr/2+1),round(nc/2-1):round(nc/2+1),drop=F]
## get distance matrix in km
  d1=as(acorr_dist(x2),"matrix")
## get mean in x & y directions
  return(mean(c(d1[2,1],d1[1,2],d1[3,2],d1[2,3])))
}