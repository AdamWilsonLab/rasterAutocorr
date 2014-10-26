####################################  FFT to Spatial correlogram  #######################################
############################  Analyses to assess the FFT algorithm  ##############################
#Testing rasterAutocor packages and generating 
#Analyses, figures, tables and data for the  paper are also produced in the script.
#AUTHOR: Benoit Parmentier 
#CREATED ON: 10/16/2014  
#MODIFIED ON: 10/26/2014            
#Version: 1
#PROJECT: Environmental Layers project                                     
#################################################################################################

#\dontrun{
# packages used for the data generation
require(raster)
require(fields)
require(latticeExtra)
library(dplyr)
library(devtools) 
#install_github("adammwilson/rasterAutocorr")
library(rasterAutocorr) #the new library

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

#r=raster(nrows=ny, ncols=nx,vals=1,xmn=-nx/2, xmx=nx/2, ymn=-ny/2, ymx=ny/2) #doe not run
r=raster(nrows=ny, ncols=nx,xmn=-nx/2, xmx=nx/2, ymn=-ny/2, ymx=ny/2)
r<- setValues(r,1)
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
debug(acorr)
## get directions for each shift to facilitate plotting of the correlogram
d2=acorr_dir(r)

## plot the autocorrlation and distance
plot(r,ylab="Y",xlab="X",main="Original (Simulated) Raster")
plot(a1,ylab="Shift in Y",xlab="Shift in X",main="Autocorrelation")
plot(d2,ylab="Shift in Y",xlab="Shift in X",main="Direction from center (degrees)")

## now we can combine the acorr values with distance
ftd=data.frame(cor=values(a1))#,dir=values(d2))

## round to  faciliate binning
ftd$dist2=round(ftd$cor.dist,3)#,c(0:50,seq(51,1000,by=10)))
## take mean by km
ftd2 <- group_by(ftd, dist2)
ftd2 <- summarise(ftd2,
                  min = min(cor.acor, na.rm = TRUE),
                  max = max(cor.acor, na.rm = TRUE),
                  sd = sd(cor.acor, na.rm = TRUE),
                  mean = mean(cor.acor, na.rm = TRUE)
)

## this table has the correlation values and distance for each pixel
head(ftd2)

## plot the correlogram
xyplot(mean/10~dist2,type=c("point","l"),span=.3,data=ftd2,ylab="Correlation",xlab="Distance",main="Correlogram",sub="Different pixels within a distance class correspond to shifts of different directions (north, south, etc.) from the origin.\n  There is one point on the plot for each pixel.",auto.key=T,lwd=2,pch=16,cex=.25)+
  layer(panel.abline(h=0))

bwplot(cor.acor/10~cut(cor.dist,pretty(ftd$cor.dist,20)),data=ftd,ylab="Correlation",xlab="Distance",main="Correlogram",sub="Different pixels within a distance class correspond\nto shifts of different directions (north, south, etc.) from the origin",col="black",fill="grey",scales=list(x=list(rot=45)))+
  layer(panel.abline(h=0))

######## PART II use new data ###############
#### get data from  interpolation project

load_obj <- function(f){
env <- new.env()
nm <- load(f, env)[1]
env[[nm]]
}

#Use elevation image

#Read in elevation map and LST map
#compute the spatial correlogram
#Make a small subset and show how it relates to the predifined functions
#Go through  every line of the code...

s_raster <- brick("/home/parmentier/Data/Spatial_correlogram_Autocorr/data/covariates__oregon_region_TMAX__OR_11032013.tif")
dim(s_raster)
#load_obj

#covar_obj_file_1<- list.files(path=in_dir1,pattern="covar_obj.*.RData")
# "covar_obj__365d_gam_daily_lst_comb5_11012013.RData"
covar_obj_file_1 <- "covar_obj__365d_gam_cai_lst_comb5_11032013.RData"
met_obj_file_1 <- "met_stations_outfiles_obj_gam_daily__365d_gam_daily_lst_comb5_11012013.RData"

#Load objects containing training, testing, models objects 
#met_stations_obj <- load_obj(file.path(in_dir1,met_obj_file_1))
covar_obj <-load_obj(file.path("/home/parmentier/Data/Spatial_correlogram_Autocorr/data",covar_obj_file_1)) #Reading covariates object for GAM daily method
infile_covariates <- covar_obj$infile_covariates
infile_reg_outline <- covar_obj$infile_reg_outline
covar_names<- covar_obj$covar_names
#####
s_raster <- brick(infile_covariates)
names(s_raster)<-covar_names
#write things out later so you don't have to load the RData?
elev <- subset(s_raster,"elev_s")


### now let's make the correlogram using fft()

## fit the complete spatial autocorrelation using fft()
## this function is the thing I need to confirm is working correctly
debug(acorr)

a1=acorr(elev)
## get directions for each shift to facilitate plotting of the correlogram
d2=acorr_dir(elev)


## plot the autocorrlation and distance
plot(elev,ylab="Y",xlab="X",main="Elevation in Oregon")
#check the scaling here!!
plot(a1,ylab="Shift in Y",xlab="Shift in X",main="Autocorrelation")
plot(d2,ylab="Shift in Y",xlab="Shift in X",main="Direction from center (degrees)")

## now we can combine the acorr values with distance
ftd=data.frame(cor=values(a1))#,dir=values(d2))

## round to  faciliate binning
ftd$dist2=round(ftd$cor.dist,3)#,c(0:50,seq(51,1000,by=10)))
## take mean by km
ftd2 <- group_by(ftd, dist2)
ftd2 <- summarise(ftd2,
                  min = min(cor.acor, na.rm = TRUE),
                  max = max(cor.acor, na.rm = TRUE),
                  sd = sd(cor.acor, na.rm = TRUE),
                  mean = mean(cor.acor, na.rm = TRUE)
)

## this table has the correlation values and distance for each pixel
head(ftd2)

## plot the correlogram
#mean by 10??
xyplot(mean/10~dist2,type=c("point","l"),span=.3,data=ftd2,ylab="Correlation",xlab="Distance",main="Correlogram",sub="Different pixels within a distance class correspond to shifts of different directions (north, south, etc.) from the origin.\n  There is one point on the plot for each pixel.",auto.key=T,lwd=2,pch=16,cex=.25)+
  layer(panel.abline(h=0))

bwplot(cor.acor/10~cut(cor.dist,pretty(ftd$cor.dist,20)),data=ftd,ylab="Correlation",xlab="Distance",main="Correlogram",sub="Different pixels within a distance class correspond\nto shifts of different directions (north, south, etc.) from the origin",col="black",fill="grey",scales=list(x=list(rot=45)))+
  layer(panel.abline(h=0))


#### Use interpolated data from Multi-time scale paper

### generate filter for Moran's I function in raster package
autocor_filter_fun <-function(no_lag=1,f_type="queen"){
  if(f_type=="queen"){
    no_rows <- 2*no_lag +1
    border_row <-rep(1,no_rows)
    other_row <- c(1,rep(0,no_rows-2),1)
    other_rows <- rep(other_row,no_rows-2)
    mat_data<- c(border_row,other_rows,border_row)
    autocor_filter<-matrix(mat_data,nrow=no_rows)
  }
  #if(f_type=="rook){} #add later
  return(autocor_filter)
}
#MODIFY: calculate for multiple dates and create averages...
#Now run Moran's I for raster image given a list of  filters for different lags and raster stack
moran_multiple_fun<-function(i,list_param){
  #Parameters:
  #list_filters: list of filters with different lags in the image
  #r_stack: stack of raster image, only the selected layer is used...
  list_filters <-list_param$list_filters
  r <- subset(list_param$r_stack,i)
  moran_list <- lapply(list_filters,FUN=Moran,x=r)
  moran_v <-as.data.frame(unlist(moran_list))
  names(moran_v)<-names(r)
  return(moran_v)
}

#Extract moran's I profile from list of images...the list may contain sublist!!! e.g. for diffeferent
#methods in interpolation
calculate_moranI_profile <- function(lf,nb_lag){
  list_filters<-lapply(1:nb_lag,FUN=autocor_filter_fun,f_type="queen") #generate lag 10 filters
  #moran_list <- lapply(list_filters,FUN=Moran,x=r)
  list_moran_df <- vector("list",length=length(lf))
  for (j in 1:length(lf)){
    r_stack <- stack(lf[[j]])
    list_param_moran <- list(list_filters=list_filters,r_stack=r_stack) #prepare parameters list for function
    #moran_r <-moran_multiple_fun(1,list_param=list_param_moran)
    nlayers(r_stack) 
    moran_I_df <-mclapply(1:nlayers(r_stack), list_param=list_param_moran, FUN=moran_multiple_fun,mc.preschedule=FALSE,mc.cores = 10) #This is the end bracket from mclapply(...) statement

    moran_df <- do.call(cbind,moran_I_df) #bind Moran's I value 10*nlayers data.frame
    moran_df$lag <-1:nrow(moran_df)
  
    list_moran_df[[j]] <- moran_df
  }
  names(list_moran_df) <- names(lf)
  return(list_moran_df)
}

plot_moranI_profile_fun <- function(i,list_param){
  #extract relevant parameters
  list_moran_df <- list_param$list_moran_df[[i]]
  title_plot <- list_param$list_title_plot[[i]]
  layout_m <- list_param$layout_m
  names_panel_plot <- list_param$names_panel_plot
  
  #date_selected <- list_param$list_date_selected[i]
  list_dd <- vector("list",length=length(list_moran_df))

  ## Reorganize data to use xyplot
  for(j in 1:length(list_moran_df)){
    method_name <- names(list_moran_df)[j]
    mydata <- list_moran_df[[j]]
    dd <- do.call(make.groups, mydata[,-ncol(mydata)]) 
    dd$lag <- mydata$lag
    dd$method_v <- method_name
    list_dd[[j]] <- dd
  }
  dd_combined<- do.call(rbind,list_dd)
  
  ## Set up plot
  #layout_m<-c(2,4) #one row two columns
  
  png(paste("Spatial_correlogram_prediction_models_levelplot_",i,"_",out_prefix,".png", sep=""),
      height=480*layout_m[1],width=480*layout_m[2])
  par(mfrow=layout_m)
  p<-xyplot(data ~ lag | which , data=dd_combined,group=method_v,type="b", as.table=TRUE,
            pch=1:3,,pch.cex=3,auto.key=list(columns=1,space="right",title="Method",cex=2.5,font=2), #Legend information
            main=title_plot,
            par.settings = list(
              superpose.symbol = list(lty=1,pch=1:3,col=1:3,cex=2.3),
              axis.text = list(font = 2, cex = 2),layout=layout_m,
              par.main.text=list(font=2,cex=2.5),strip.background=list(col="white")),
            par.strip.text=list(font=2,cex=2),
            strip=strip.custom(factor.levels=names_panel_plot),
            xlab=list(label="Spatial lag neighbor", cex=2.5,font=2),
            ylab=list(label="Moran's I", cex=2.5, font=2))
  
  print(p)
  
  dev.off()
  
  return(p)
}

#### Figure 5: Spatial lag profiles  
#This figure is generated to show the spatial Moran'I for 10 spatial 
#for Jan 1 and Sept 1 in 2010 for all models (1 to 7) and methods

index <- 1 #index corresponding to Jan 1 #For now create Moran's I for only one date...
lf_moran_list_date1 <-lapply(list_raster_obj_files[c("gam_daily","gam_CAI","gam_fss")],
                               FUN=function(x){x<-load_obj(x);x$method_mod_obj[[index]][[y_var_name]]})                           
index <- 244 #index corresponding to Sept 1 #For now create Moran's I for only one date...
lf_moran_list_date2 <-lapply(list_raster_obj_files[c("gam_daily","gam_CAI","gam_fss")],
                               FUN=function(x){x<-load_obj(x);x$method_mod_obj[[index]][[y_var_name]]})                           

#date_selected <- "20100101"
#methods_names <-c("gam","kriging","gwr")
methods_names <-c("gam_daily","gam_CAI","gam_FSS") #drop gam FSS later

names_layers<-methods_names
#Subset images to eliminate mod_kr
lf1 <- list(lf_moran_list_date1[[1]],lf_moran_list_date1[[2]][1:7]) #,lf_moran_list_date1[[3]][1:7])
lf2 <- list(lf_moran_list_date2[[1]],lf_moran_list_date2[[2]][1:7]) #,lf_moran_list_date2[[3]][1:7])
names(lf1)<-c("gam_daily","gam_CAI") #,"gam_FSS")
names(lf2)<-c("gam_daily","gam_CAI") #,"gam_FSS")

### Now extract Moran's I for a range of lag using a list of images

#set maximum lag range
nb_lag <-10
#Provide list of raster images:
list_lf <- list(lf1,lf2)

list_moran_df1 <- calculate_moranI_profile(list_lf[[1]],nb_lag) #for January 1
list_moran_df2 <- calculate_moranI_profile(list_lf[[2]],nb_lag) #for September 1
names(list_moran_df1)<-c("gam_daily","gam_CAI")#,"gam_FSS")
names(list_moran_df2)<-c("gam_daily","gam_CAI")#,"gam_FSS")

#Run accross two dates...
#list_moran lapply(list_lf,FUN=calculate_moranI_profile,nb_lag=nb_lag)

### Prepare to plot lag Moran's I profiles

#generate automatic title for exploration if necessary!!
list_title_plot<- list(c("Spatial lag profile on January 1, 2010"),
                  c("Spatial lag profile on September 1, 2010"))
#name used in the panel!!!
names_panel_plot <-c("mod1 = lat*long","mod2 = lat*long + LST","mod3 = lat*long + elev","mod4 = lat*long + elev + N_w*E_w",
                 "mod5 = lat*long + elev + DISTOC","mod6 = lat*long + elev + LST","mod7 = lat*long + elev + LST*FOR")
layout_m<-c(2,4) # works if set to this?? ok set the resolution...
list_moran_df <- list(list_moran_df1,list_moran_df2)
list_param_plot_moranI_profile_fun <- list(list_moran_df,list_title_plot,names_panel_plot,layout_m)
names(list_param_plot_moranI_profile_fun) <- c("list_moran_df","list_title_plot","names_panel_plot","layout_m")

#debug(plot_moranI_profile_fun)
#p<- plot_moranI_profile_fun(1,list_param=list_param_plot_moranI_profile_fun)
  
list_moranI_plots <- lapply(1:2,FUN=plot_moranI_profile_fun,list_param=list_param_plot_moranI_profile_fun)


#This function uses list moran_df object from calculate_moranI_profile function!!

#layout_m<-c(2,4) # works if set to this?? ok set the resolution...
#layout_m<-c(2*3,4) # works if set to this?? ok set the resolution...

png(paste("Figure5_paper_spatial_correlogram_prediction_models_levelplot_",out_prefix,".png", sep=""),
    height=480*layout_m[1],width=480*layout_m[2])
    #height=3*480*layout_m[1],width=2*480*layout_m[2])
    #height=480*6,width=480*4)
#png(paste("Figure11_paper_spatial_correlogram_prediction_models_levelplot_",out_prefix,".png", sep=""),
#    height=480,width=480)
    #height=480*6,width=480*4)

#p1 <- list_moranI_plots[[1]]
p2 <- list_moranI_plots[[2]]
#grid.arrange(p2,ncol=1)
print(p2)
#grid.arrange(p1,p2,ncol=1)
dev.off()

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
