## loglik function
## x: the quasi-likelihood fit
## data: the observed and quadrature points
## value: log-likelihood of the fitted poisson process
## adjust the total area so it match
ppm.st.lik <- function(x,data){
  lambda <- fitted(x)
  l.lik <- with(subset(data,!is.na(z)),{
    wt * (z*log(lambda)-lambda)
  })
  sum(l.lik)
}

## value: residual measure of each quadrature point
ppm.st.residuals <- function(x,data){
  l.data <- subset(data,!is.na(z))
  lambda <- predict(x,newdata=l.data,type="response")
  l.res <- with(l.data,{
    wt*(z-lambda)*1/sqrt(lambda)
  })
  l.res
}


## x: predicted data with residuals
## formula: define residual on left and area on right and time period on right
## value:
## a residual plot
## plus the aggregated time
plot.ppm.st.residuals <- function(formula,x,cex=3,pch=21,bg="steelblue1",...)
{
  tmp1 <- aggregate(formula,data=x,FUN=sum)
  vars <- all.vars(formula)
  wt <- tmp1[,vars[1]]
  upr <- 2*sqrt(wt)
  lwr <- -2*sqrt(wt)
  e <- tmp1[,vars[2]]
  t_ <- tmp1[,vars[3]]
  plot(t_,e,ylim=range(e,upr,lwr),
       xlab="Time period",ylab="Person Residual",...)
  polygon(x=c(t_,rev(t_)),y=c(lwr,rev(upr)),col=grey(0.9))
  points(t_,e,cex=cex,pch=pch,bg=bg)
  lines(lowess(t_,e),col="red")
  tmp1
}

