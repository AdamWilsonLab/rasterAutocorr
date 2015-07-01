\dontrun{
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


#######################################################
## create a random raster with some spatial structure
nx=500  # size of raster
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

################################
### now let's make the correlogram using fft()

## fit the complete spatial autocorrelation using fft()
## this function is the thing I need to confirm is working correctly
a1=acorr(r)

## get directions for each shift to facilitate plotting of the correlogram
d2=acorr_dir(r)


## plot the autocorrlation and distance
  plot(r,ylab="Y",xlab="X",main="Original (Simulated) Raster")
  plot(a1,ylab="Shift in Y",xlab="Shift in X",main="Autocorrelation")
  plot(d2,ylab="Shift in Y",xlab="Shift in X",main="Direction from center (degrees)")

## now we can combine the acorr values with distance
ftd=data.frame(cor=values(a1))#,dir=values(d2))

## round to  faciliate binning
ftd$dist=round(ftd$cor.dist,3)#,c(0:50,seq(51,1000,by=10)))
## take mean by km
ftd2 <- group_by(ftd, dist)
ftd2 <- summarise(ftd2,
                  min = min(cor.acor, na.rm = TRUE),
                  max = max(cor.acor, na.rm = TRUE),
                  sd = sd(cor.acor, na.rm = TRUE),
                  mean = mean(cor.acor, na.rm = TRUE)
)

## this table has the correlation values and distance for each pixel
head(ftd2)

## plot the correlogram
xyplot(mean/10~dist,type=c("point","l"),span=.3,data=ftd2,ylab="Correlation",xlab="Distance",main="Correlogram",sub="Different pixels within a distance class correspond to shifts of different directions (north, south, etc.) from the origin.\n  There is one point on the plot for each pixel.",auto.key=T,lwd=2,pch=16,cex=.25)+
  layer(panel.abline(h=0))

bwplot(cor.acor/10~cut(cor.dist,pretty(ftd$cor.dist,20)),data=ftd,ylab="Correlation",xlab="Distance",main="Correlogram",sub="Different pixels within a distance class correspond\nto shifts of different directions (north, south, etc.) from the origin",col="black",fill="grey",scales=list(x=list(rot=45)))+
  layer(panel.abline(h=0))


###
## Matlab/Octave script
#x1=[3 6 5 ; 7 2 2 ; 4 NaN 0]
#[n,p]=size(x1); 
#ncols=2*p-1;
#nrows=2*n-1;
#nr2=5
#nc2=5
#x1id=~isnan(x1);
#x1(~x1id)=zeros(sum(sum(~x1id),1);
#fx1=fft2(x1,nr2,nc2);
#fx1id=fft2(x1id,nr2,nc2);
#nh11=round(real(ifft2(conj(fx1id) .*fx1id)));
#m1=real(ifft2(conj(fx1) .*fx1id)) ./max(nh11,1);
#m2=real(ifft2(conj(fx1id) .*fx1)) ./max(nh11,1)
#gh11=real(ifft2(conj(fx1) .*fx1)); gh11=gh11./max(nh11,1)-m1.*m2;

}
