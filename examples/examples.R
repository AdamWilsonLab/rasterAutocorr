\dontrun{
# packages used for the data generation
require(raster)
require(gstat)

## create a random raster with some spatial structure
nx=100  # size of raster
ny=nx
r=raster(nrows=ny, ncols=nx,vals=1,xmn=-nx/2, xmx=nx/2, ymn=-ny/2, ymx=ny/2)
names(r)="z"

# define a spatial model using gstat
range=nx*10
sill=10
smod = gstat(formula=~1, data=r, dummy=TRUE, 
             beta=1, model=vgm(psill=sill,model="Exp",range=range),nmax=100)

# predict a new surface using that model
s=interpolate(r,smod,nsim=5)

## check out the random image
## for some reason there are bands in the image, which I'm sure can be fixed
## but they don't matter for this purpose...
plot(s)


################################
### now let's make the correlogram using fft()

## fit the complete spatial autocorrelation using fft()
## this function is the thing I need to confirm is working correctly
a1=acorr(s)

## get distances for each shift to facilitate plotting of the correlogram
## I think this may be off by a pixel or two, need to check...
d1=acorr_dist(s)


## plot the autocorrlation and distance
  plot(a1,ylab="Shift in Y",xlab="Shift in X",main="Autocorrelation")
  plot(d1,ylab="Shift in Y",xlab="Shift in X",main="Distance in units of original raster")

## now we can combine the acorr values with distance
ftd=data.frame(cor=values(a1),dist=values(d1))

## this table has the correlation values and distance for each pixel
head(ftd)

## plot the correlogram
plot(cor~dist,data=ftd,ylab="Correlation",xlab="Distance",main="Correlogram",sub="Different pixels within a distance class correspond to shifts of different directions (north, south, etc.) from the origin")

xyplot(cor~dist,type=c("p","spline"),data=ftd,ylab="Correlation",xlab="Distance",main="Correlogram",sub="Different pixels within a distance class correspond to shifts of different directions (north, south, etc.) from the origin",col=c("black"),col.line="red",lwd=2)




## fit the correlogram with ncf using points
#library(ncf)
#pdata=cbind.data.frame(coordinates(s),z=values(s))
## sample cause it's slow
#samp=sample(1:nrow(pdata),1000)
#ncf.cor <- correlog(pdata$x[samp],pdata$y[samp],pdata$z[samp],increment=2)
#plot(ncf.cor)
}
