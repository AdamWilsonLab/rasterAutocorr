\dontrun{
# packages used for the data generation
require(raster)
require(fields)
require(latticeExtra)

## create a random raster with some spatial structure
nx=1000  # size of raster
ny=round(nx*1.2)
r=raster(nrows=ny, ncols=nx,vals=1,xmn=-nx/2, xmx=nx/2, ymn=-ny/2, ymx=ny/2)
names(r)="z"

#Simulate a Gaussian random field with an exponential covariance function
## Theta is the scale of an exponential decay function.  This controls degree of autocorrelation, values close to 1 are close to random while values near nx/4 have high autocorrelation
theta=50

grid=list(x=seq(xmin(r),xmax(r)-1,by=res(r)[1]),y=seq(ymin(r),ymax(r)-1,res(r)[2]))

obj<-Exp.image.cov(grid=grid, theta=theta, setup=TRUE)
look<- sim.rf( obj)

values(r)=t(look)

image(r,asp=1)

################################
### now let's make the correlogram using fft()

## fit the complete spatial autocorrelation using fft()
## this function is the thing I need to confirm is working correctly
a1=acorr(r,fmean="CellStats")

## get distances for each shift to facilitate plotting of the correlogram
d1=acorr_dist(r)
## get directions for each shift to facilitate plotting of the correlogram
d2=acorr_dir(r)


## plot the autocorrlation and distance
  plot(r,ylab="Y",xlab="X",main="Original (Simulated) Raster")
  plot(a1,ylab="Shift in Y",xlab="Shift in X",main="Autocorrelation")
  plot(d1,ylab="Shift in Y",xlab="Shift in X",main="Distance from center in units of original raster")
  plot(d2,ylab="Shift in Y",xlab="Shift in X",main="Direction from center (degrees)")

## now we can combine the acorr values with distance
ftd=data.frame(cor=values(a1),dist=values(d1)/1000,dir=values(d2))

## this table has the correlation values and distance for each pixel
head(ftd)

## plot the correlogram
xyplot(cor~dist,groups=cut(dir,5),type=c("p","smooth"),span=.3,data=ftd,ylab="Correlation",xlab="Distance",main="Correlogram",sub="Different pixels within a distance class correspond to shifts of different directions (north, south, etc.) from the origin.  There is one point on the plot for each pixel.",auto.key=T,lwd=2,pch=16,cex=.25)+
  layer(panel.abline(h=0))

bwplot(cor~cut(dist,pretty(ftd$dist,20)),type=c("p","smooth"),data=ftd,ylab="Correlation",xlab="Distance",main="Correlogram",sub="Different pixels within a distance class correspond to shifts of different directions (north, south, etc.) from the origin",col="black",fill="grey",scales=list(x=list(rot=45)))+
  layer(panel.abline(h=0))


## fit the correlogram with ncf using points
#library(ncf)
#pdata=cbind.data.frame(coordinates(s),z=values(s))
## sample cause it's slow
#samp=sample(1:nrow(pdata),1000)
#ncf.cor <- correlog(pdata$x[samp],pdata$y[samp],pdata$z[samp],increment=2)
#plot(ncf.cor)
}
