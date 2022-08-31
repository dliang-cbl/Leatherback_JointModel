## predict from fitted object from sdm_bc
## ver 4): use cut workflow to reflect uncertainty in prediction
## sub ver h); use seasonal and longterm median instead of mean
##             when aggregating.
## sub ver g): simple parallel computing


## object: fitted object from sdm.bc (ver 3.d onward)
## formula : a formula, LHS a temporal variable to aggregate by
##           RHS spatial variables to aggregate by
##         if formula is NULL then no aggregation is done
## newd  : data frame containing new predicted values
## all   : T/F to include all terms in prediction
##         when False, only derived terms from tracking data used
## type  : link or response, whether report prediction on response 
##         scale default=response
## terms : similar definition to predict.gam include only these terms
##         in prediction 
##        NOTE: terms with all=T does not make sense and ignored
## batch: default all in one go, a positive number of batches to process
## value : a list returning fitted value and standard errors
##       : after aggregation using by variable in newd
predict_sdm_bc <- function(object,newd,formula=NULL,
                           all=F,type="response",terms=NULL,
                           verbose=FALSE,batch=NULL){
  
  stopifnot(type %in% c("response","link"))
  
  if(is.null(formula)){
    if(verbose) cat("generating predictive samples..\n")
    if(is.null(batch)){
      return(predict_sdm_bc_work(
        object=object,newd=newd,all = all,
        type=type,terms=terms,verbose=verbose
      ))
    }else{
      uid <- seq(1,nrow(newd))
      batch.ii <- rep(1:batch,length=nrow(newd))
      #sam.list <- vector("list",batch)
      sam.list <- foreach(i=1:batch) %dopar% {
        if(verbose) cat("processing batch ",i,"\n")
        bsel_ <- batch.ii==i 
        line_ <- predict_sdm_bc_work(
          object=object,newd=newd[bsel_,],
          all = all,type=type,terms=terms,verbose=verbose
        )
        cbind(uid=uid[bsel_],fit=line_$fit,sd=line_$sd)
      }
      sam_line <- as.data.frame(do.call(rbind,sam.list))
      return(cbind(newd,sam_line[order(sam_line[,"uid"]),-1]))
    }
  }
  else{
    if(is.null(batch)){
      return(predict_sdm_bc_work_by(
        object = object,newd=newd,formula=formula,
        all=all,type=type,terms=terms,verbose=verbose))
    }else{
      ## split data according to formula
      if(verbose) cat("preparing aggregation.\n")
      keys <- unique(newd[,all.vars(formula)])
      batch.ii <- rep(1:batch,length=nrow(keys))
      #sam.list <- vector("list",batch)
      sam.list <- foreach(i=1:batch) %dopar% {
        if(verbose) cat("processing batch ",i,"\n")
        bsel_ <- batch.ii==i
        newd_batch <- merge(newd,keys[bsel_,])
        predict_sdm_bc_work_by(
          object=object,newd=newd_batch,formula=formula,
          all = all,type=type,terms=terms,verbose=verbose
        )
      }
      return(do.call(rbind,sam.list))
    }
  }
}
predict_sdm_bc_work <- function(object,newd,all=T,type="response",terms=NULL,verbose=F){
  
  ## prepare design matrix of unbiased prediction
  newd$crw <- 0
  newd$tau <- 1
  newd$id <- object$ctds$model$id[1]
  Xp <- predict(object$ctds,newd,type="lpmatrix")
  
  if(!all){
    sel <- c("(Intercept)",paste0("E.",1:length(object$.args$term2)))
    sel <- which(names(object$value$theta[[1]]$coef) %in% sel)
  }else{
    if(!is.null(terms)){
      ## terms with all=TRUE does not make sense
      if(verbose) cat("ignoring terms =",terms,"\n")
      terms <- NULL 
    }
    sel <- seq(1,length(object$value$theta[[1]]$coef))
  }
  
  ## check input terms
  if(!is.null(terms)){
    stopifnot(all(terms %in% all.vars(object$ctds$formula)))
  }
  
  iter <- nrow(object$value$beta) ## number of iterations
  storage <- matrix(NA,nrow(Xp),iter) ## predictive samples
  for(i in 1:iter){
    
    if(verbose) cat("predictive distribution ",i,"\n")
    ## reuse the fitted ppm object
    ##1. unbiased and un-scaled prediction of log-intensity
    
    curr_beta <- object$value$beta[i,]
    # if(!is.null(terms)){
    #   ## zero out certain terms
    #   term_lst <- lapply(terms,function(x) grep(x,names(curr_beta)))
    #   term_ii <- do.call(c,term_lst)
    #   curr_beta[-term_ii] <- 0
    # }
    
    if(object$.args$grouped){
      E.lst.pred <- predict.gam.term(
        formula=object$.args$formula,
        x = object$.args$term2,
        beta = curr_beta[-object$.args$group.col],
        X = Xp[,-object$.args$group.col,drop=FALSE]
      )
    }else{
      E.lst.pred <- predict.gam.term(
        formula=object$ctds$formula,
        x=object$.arg$term2,
        beta=curr_beta,
        X=Xp
      )
    }


    ##2. prepare prediction matrix
    if(i==1){
      Zp <- predict(object$ppm,newdata=cbind(newd,E.lst.pred),
                    type="lpmatrix")
    }else{
      Zp[,names(E.lst.pred)] <- as.matrix(E.lst.pred)
    }
    
    ## 3.prediction
    curr_theta <- object$value$theta[[i]]$b
    if(!is.null(terms)){
      if(colnames(Zp)[1]=="(Intercept)") curr_theta[1] <- 0
    }
    eta <- as.numeric(Zp[,sel,drop=FALSE]%*%curr_theta[sel])
    if(type=="response"){
      storage[,i] <- exp(eta)
    }else{
      storage[,i] <- eta
    }
    
  }
  
  if(verbose) cat("calculating posterior summaries.\n")

  fit <- apply(storage,1,median)
  fit.var <- apply(storage,1,var)
  cbind(newd,fit=fit,sd=sqrt(fit.var))
}


predict_sdm_bc_work_by <- function(object,newd,formula,all=T,
                                   type="response",terms=NULL,verbose=F){
  ## prepare design matrix of unbiased prediction
  newd$crw <- 0
  newd$tau <- 1
  newd$id <- object$ctds$model$id[1]
  Xp <- predict(object$ctds,newd,type="lpmatrix")
  
  ## select terms to include.
  if(!all){
    sel <- c("(Intercept)",paste0("E.",1:length(object$.args$term2)))
    sel <- which(names(object$value$theta[[1]]$coef) %in% sel)
  }else{
    if(!is.null(terms)){
      if(verbose) cat("ignoring terms =",terms,"\n")
      terms <- NULL 
    }
    sel <- seq(1,length(object$value$theta[[1]]$coef))
  }

  ## check input terms
  if(!is.null(terms)){
    stopifnot(all(terms %in% all.vars(object$ctds$formula)))
  }
  
  iter <- nrow(object$value$beta) ## number of iterations
  storage <- vector("list",iter) ## predictive samples
  for(i in 1:iter){
    
    if(verbose) cat("predictive distribution ",i,"\n")
    ## reuse fitted ppm object
    ##3. unbiased and un-scaled prediction of log-intensity
    curr_beta <- object$value$beta[i,]
    # if(!is.null(terms)){
    #   ## zero out certain terms
    #   term_lst <- lapply(terms,function(x) grep(x,names(curr_beta)))
    #   term_ii <- do.call(c,term_lst)
    #   curr_beta[-term_ii] <- 0
    # }
    
    if(object$.args$grouped){
      E.lst.pred <- predict.gam.term(
        formula=object$.args$formula,
        x = object$.args$term2,
        beta = curr_beta[-object$.args$group.col],
        X = Xp[,-object$.args$group.col,drop=FALSE]
      )
    }else{
      E.lst.pred <- predict.gam.term(
        formula=object$ctds$formula,
        x=object$.arg$term2,
        beta=curr_beta,
        X=Xp
      )
    }
    
    ##4. prepare prediction matrix
    if(i==1){
      Zp <- predict(object$ppm,newdata=cbind(newd,E.lst.pred),
                    type="lpmatrix")
    }else{
      Zp[,names(E.lst.pred)] <- as.matrix(E.lst.pred)
    }
    

    curr_beta <- object$value$theta[[i]]$b
    if(!is.null(terms)){
      if(colnames(Zp)[1]=="(Intercept)") curr_beta[1] <- 0
    }
    
    ##5. generate prediction
    eta <- as.numeric(Zp[,sel,drop=FALSE]%*%curr_beta[sel])
    if(type=="response"){
      lambda_ <- exp(eta) 
    }else{
      lambda_ <- eta
    }
    
    ##6. aggregation
    storage[[i]] <- aggregate(
      reformulate(all.vars(formula),response="lambda_"),
      data=newd,FUN=median)
  }
  
  if(verbose) cat("calculating posterior summaries.\n")
  pred_ <- do.call(rbind,storage)
  ret_ <- aggregate(reformulate(all.vars(formula),response="lambda_"),
            data=pred_,FUN=function(x){
              c(fit=median(x),sd=sd(x))
            })
  cbind(ret_[,1:length(all.vars(formula))],ret_[,ncol(ret_)])
}  
