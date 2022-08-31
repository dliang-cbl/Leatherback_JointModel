## intensity_2.R
## ver 2) allow more resolution to be shown
require(raster)
require(MBA)
intensity.rescale <- function(object,eps=NULL)
{
  if(!is.null(eps)){
    kfac <- max(ceiling(eps/res(object)))
    eraster1 <- object
    res(eraster1) <- eps
    if(kfac>1){
      object <- aggregate(object,fact=kfac)
    }
    output <- resample(object,eraster1,method="ngb")
  }
  else{
    output <- object
  }
  output
}
intensity <- function(object,xyz,eps=NULL,n=1,m=1,h=8,extend=FALSE,sp=FALSE,...)
{
  ## helper script to show irregular data on a raster by MBA package.
  ## object: a raster object
  ## xyz   : xyz matrix to visualize on raster
  ## eps   : resolution of the raster object
  ## Value: a raster object for plotting

  object <- intensity.rescale(object,eps=eps)  ## rescale the raster* object
  
  l.xyz <- rasterToPoints(object)      ## extract all non-missing cells
  l.id <- cellFromXY(object,l.xyz[,1:2]) ## the cell number in the raster 
  
  l.xyz.est <- mba.points(as.matrix(xyz),l.xyz[,1:2])[[1]] ## multilevel B-splines interpolation on raster
  
  r.object <- object
  r.object[l.id] <- l.xyz.est[,3]
  
  r.object
}