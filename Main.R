## R Vignette 
## Integrated Specie Distribution Modeling 
## D. Liang
## University of Maryland Center for Environmental Science

## Clean workspace and load R packages
rm(list=ls())
options(warn = 1)
library(raster); library(ctmcmove); library(mgcv);
options("rgdal_show_exportToProj4_warnings"="none")

## load simulated data
cat("loading simulated data \n")
## Fishery observations
loc <- read.csv("loc.csv")
str(loc)
## Long/Lat: observation location
## amon    : observation month (used to link to environmental variable)

## Tracking observations
track <- read.csv("track.csv")
str(track)
## x/y: tracking locations (same projection as loc)
## t  : time-stamp (e.g. days, used to calculate residence time)
## p  : observation (e.g. month, used to link environmental data)
## id : individual tracked

## Load seascape raster stacks
load(file="seascape.RData")
names(seascape)
## bathy: static bathymetry
## sst  : sea surface temp.
## sst2 : square of sea surface temp.
## fpi  : frontal probability index
## ssh  : sea surface height
dim(seascape$bathy)
## 60x60 horizontal raster layers
## 10   : monthly layers (used to link to tracks and fishery observations)
sapply(seascape,res)
## the same resolution 0.5 decimal degrees


## load tracking data codes
source("nmiss_2.R")
source("s2inla_6.R")
source("ctmc2glm_study_3.R")

# Process tracks into a data structure
# that allows generalized linear modeling of telemetry data
tmp_ <- split(track,track$id) ## process each individual
tmp2_ <- lapply(tmp_,function(elmt){
  line <-  with(elmt,ctmc2glm.study(
    xy=cbind(x,y),t=t,tmap=p,seascape,NULL))
  if(!is.null(line)){
    line$id <- elmt$id[1]  
  }
  line
})
glmobj <- do.call(rbind,tmp2_)
glmobj$id <- factor(glmobj$id)

## rename the environmental variables 
## consistent with the names of stacks in the seascape 
names(glmobj)[c(2:6)] <- c("bathy","sst","sst2","fpi","ssh")


## load point process codes
source("extract_1.R")
source("cellFromXY_nbr.R")
source("ppm_st_glm_6f.R")
source("ppm_st_utils.R")

## study area raster layer (NA=outside area, non-NA=in the area)
sp <- raster(seascape[["bathy"]],layer=1)
## plot(sp)

## build pseudo-absence data for dynamic point process modeling
quad <- ppm.st.glm(loc, window = sp, eps = 0.5,
                  covariates = seascape,
                  trend=~bathy+sst+sst2+fpi+ssh, 
                  coordinate=~long+lat+amon)

## Load ISDM modeling codes
library(abind)
source("term_utils.R")
source("sim2jam_study_1.R")
source("all_group_1.R")
source("predict_gam_term_2.R")
source("sdm_bc_5d.R")

## telemetry model: linking motility to seascape variables
term0 <- c("bathy","s(sst,bs='cr',k=3)","fpi","ssh")
for0 <- reformulate(term0,response="z")
## bias correction model
for1 <- ~s(long,lat,bs="ds",m=c(1,0.5),k=10)

## build scaling model list all subset for sdm terms
scale_terms <- all_group(
  reformulate(term0,response = "z"),verbose=F)
head(scale_terms)
## each row denotes a possible combination of terms
## value in the columns denote whether the scaling terms are equal
## for example 1,2,2,2 indicates that the last three terms
##  share the same scaling constant 

term1 <- scale_terms[15,] ## pick the most general scaling model
                          ## each environmental data has distinct scaling

iter <- 9 ## Number of Monte Carlo iterations: 
          ## kept small for demo only, at least 199 recommended.
## Integrated modeling
ts0 <- proc.time()
bsdm.obj <- sdm.bc(
  ## see sdm_bc_5d.R for argument definitions
  formula = for0,
  term = term1,
  bias = for1,
  quad = quad,
  glmobj = glmobj,
  wcoord = "wt", ## weight for pseudo-absence, 
  gcoord="id",tcoord = "tau",cterm="crw", ## terms defining individual, residence time and last movement angle
  d=-1,S=3/8, ## lognormal prior mean and standard deviation for scaler
  timefact = 30,  ## 30 days in a month assumed
  iter =iter,trace=T
)
ts1 <- proc.time()
cat("Joint analyses took",ts1-ts0,"\n")


## load prediction codes
source("predict_sdm_bc_5d.R")
source("intensity_2.R")

## prepare parallel nodes
library(doParallel)
if(FALSE){
  cl <- makeCluster(nodes)
  dump <- clusterEvalQ(cl,{
    source("ppm_st_glm_6f.R")
    source("extract_1.R")
    source("intensity_2.R")
    source("predict_gam_term_2.R")
    source("term_utils.R")
    source("predict_sdm_bc_5d.R")})
  registerDoParallel(cl)
}

newd <- ppm.st.glm.new(
  window=sp,newdata=seascape,
  coord=~long+lat+amon)


cat("predicting monthly intensity...\n")
(ts1 <- proc.time())
raw.pred <- predict_sdm_bc(bsdm.obj,newd,all=F,
                           verbose=F,batch=10)
(ts2 <- proc.time())


predMonth <- ppm.st.glm.fill(
  window = sp,
  newdata = raw.pred,
  formula = fit~long+lat+amon,
  verbose = T,
)


dim(predMonth)
## 10 layers of monthly predictions
output <- raster(predMonth,layer=1)

plot(output)
writeRaster(output,file=paste0("intensity.tif"),overwrite=T)
