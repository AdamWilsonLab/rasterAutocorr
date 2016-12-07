# packages used for the data generation
require(raster)
require(fields)
library(dplyr)

## first show the examples in Marcotte, Denis. 1996. "Fast Variogram Computation with FFT." Computers & Geosciences 22 (10): 1175–86. doi:10.1016/S0098-3004(96)00026-X.
m1=raster(matrix(c(3,6,5,7,2,2,4,NA,0),ncol=3,byrow=T))
m2=raster(matrix(c(10,NA,5,NA,8,7,5,9,11),ncol=3,byrow=T))

ac=acorr(m1,padlongitude=T,verbose=T)

## confirm nobs == nh11 on top of page 1179
10^as.matrix(ac[["nobs"]])

## comfirm acor= gh11 on top of page 1179
round(as.matrix(ac[["acor"]])/10,3)

