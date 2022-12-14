## s2inla: helper function to write GAM within Bayesian framework
## D. Liang
## ver. 6: minor fix to allow interface to stan as well as INLA
## bug: no smooth interface.
## helper function to define smooth term fitting in INLA
## supporting offset through formula instead
require(mgcv)
## helper function to fix smooth terms in INLA
## formula: the formula containing only smooth terms
## data: the data to define smooth terms
## U  : standard deviations of T Normal priors on standard deviation of random effects
##    : a vector of two, the first entry is for the smooth with rank > 1
##    :                  the second entry is for the smooth with rank 1
##    : NOTE: these numbers appear work with single smooth, but the linear
##    combination part is tricky with z models
## rankone: logical, whether include rankone smooth or not.
## value: X: the Z matrices for INLA definition
##        formula: the INLA formula for the smooth terms
##        idx: the indices for Z models
s2inla <- function(formula,data,U=c(0.00001,10),
                   rankone=FALSE,family=gaussian)
{
  catch <- try(model.matrix(formula,data),silent = TRUE)
  if(class(catch)!="try-error"){
    ## the formula contain no s term
    ## this may require a more permanant fix in case model.matrix
    ## starts to support s term..
    frame <- model.frame(formula,data)
    offset <- model.offset(frame)
    if(is.null(offset)){
      offset <- rep(0,length(model.response(frame)))
    }
    return(list(y=model.response(frame),offset=offset,
                Design=catch,X=NULL,formula='',idx=NULL))
  }
  
  ## remove missing values
  data <- na.omit(data[,all.vars(formula)])
  
  jf <- tempfile()
  jd <- jagam(formula,data=data,file=jf,diagonalize = TRUE,
              na.action=na.pass,family=family)
  
  ## creat fixed design matrix
  nms <- names(jd$pregam$cmX) ## if this name is not null, it is a design matrix
  Design <- jd$pregam$X[,nchar(nms)>0]
  #browser()
  offset <- model.offset(model.frame(jd$pregam$terms,data))
  if(is.null(offset)){
    offset <- rep(0,nrow(data))
  }
  
  ## create list for random design matrix, formula and simply IDs
  X <- vector("list",length(jd$pregam$smooth))
  fout <- vector("list",length(X))
  varnames <- letters[1:length(X)]
  ones <- vector("list",length(X))
  
  for(i in 1:length(jd$pregam$smooth)){ ## loop thru smooth term
    ## create id matrices
    tmp <- matrix(1:nrow(data),nrow(data),length(jd$pregam$smooth[[i]]$rank))
    colnames(tmp) <- paste(varnames[i],1:length(jd$pregam$smooth[[i]]$rank),sep="")
    ones[[i]] <- as.data.frame(tmp)
    
    ## create design matrices
    X[[i]] <- vector("list",length(jd$pregam$smooth[[i]]$rank))
    start <- jd$pregam$smooth[[i]]$first.para
    ff <- rep("",length(jd$pregam$smooth[[i]]$rank))
    for(j in 1:length(jd$pregam$smooth[[i]]$rank)){
      idx <- start+seq(0,jd$pregam$smooth[[i]]$rank[j]-1)
      cat("idx=",idx,"\n")
      if(length(idx)>1){ ## smooth rank > 1
        l.U <- U[1]
      }                  ## end
      else{              ## smooth rank == 1
        l.U <- U[2]
      }                  ## end
      X[[i]][[j]] <- jd$jags.data$X[,idx,drop=FALSE]
      ## INLA formula
      tmpfor <- paste("f(",varnames[i],j,',model="z",',"Z=X[[",i,"]][[",j,"]],",
                      'hyper=list(prec=list(prior="logtnormal",param=c(0,',
                      l.U,'))))',
                      sep="")
      if(length(idx)>1){
        ff[j] <- tmpfor
      }
      else{
        if(rankone){
          ff[j] <- tmpfor
        }
      }
      start <- max(idx)+1
    }
    fout[[i]] <- paste(ff[nchar(ff)>0],collapse = "+")
  }
  list(y=jd$pregam$y,Design=Design,X=X,formula=do.call(c,fout),idx=do.call(cbind,ones),
       pregam=jd$pregam,offset=offset)
}
