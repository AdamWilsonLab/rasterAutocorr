#' Really Remove Raster Temp files.
#'
#' @description Really Remove Raster Temp files. This function checks the path of a raster object and removes it from disk as well as from the R workspace.
#' @param x A raster* object. 

rmr=function(x){                                                                                                                                                                                   
  ## function to truly delete raster and temporary files associated with them                                                                                                                    


  if(class(x)=="RasterLayer"&grepl("^/tmp",x@file@name)&fromDisk(x)==T){                                                                                                                     
    file.remove(x@file@name,sub("grd","gri",x@file@name))                                                                                                                                  
    rm(x)                                                                                                                                                                                  
  }                                                                                                                                                                                              
}     