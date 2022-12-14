# Created by Fabio M. Bayer (bayer@ufsm.br)
# Department of Statistics, UFSM, Brazil
# December 2022
#
# Some info:
#
# diag = 0 : without graphs (useful for simulations)
# diag = 1 : plot diagnostic graphs
# diga = 2 : make pdf with diagnostic graphs
#
# h : number of forecast steps
# 
# resid = 1 : Pearson residuals
# resid = 2 : deviance residuals
# resid = 3 : quantile residuals (the best one)
#
# Reference:
# BAYER, F.M.; PUMI, G.; PEREIRA, T.L.; SOUZA, T.C.
# Inflated Beta Autoregressive Moving Average Models. 
# Forthcoming
# 2023

ibarma<- function (y, ar=NA, ma=NA, link = "logit",diag=1,h=6,X=NA,X_hat=NA,resid=3)
{  
  source("ibarma.fit.r")
  
  if (min(y) < 0 || max(y) > 1)
    stop("OUT OF RANGE [0,1]!")
  
  if(is.ts(y)==T)
  {
    freq<-frequency(y)
  }else stop("data can be a time-series object")
  
  
  if(any(is.na(ar))==F) names_phi<-c(paste("phi",ar,sep=""))
  
  if(any(is.na(ma))==F) names_theta<-c(paste("theta",ma,sep=""))
  
  if(any(is.na(X))==F)
  {
    names_beta<-c(paste("beta",1:ncol( as.matrix(X) ),sep=""))
  }
  
  p <- max(ar)
  q <- max(ma)
  n <- length(y)
  
  m <- max(p,q,na.rm=T)
  
  p1 <- length(ar)
  q1 <- length(ma)
  
  linktemp <- substitute(link)
  if (!is.character(linktemp))
  {
    linktemp <- deparse(linktemp)
    if (linktemp == "link")
      linktemp <- eval(link)
  }
  if (any(linktemp == c("logit", "probit", "cloglog")))
    stats <- make.link(linktemp)
  else stop(paste(linktemp, "link not available, available links are \"logit\", ",
                  "\"probit\" and \"cloglog\""))
  
  link1 <- structure(list(link = linktemp, 
                          linkfun = stats$linkfun,
                          linkinv = stats$linkinv, 
                          mu.eta = stats$mu.eta, 
                          diflink = function(t) 1/(stats$mu.eta(stats$linkfun(t)))
  )
  )
  
  
  fit1 <- ibarma.fit(y, ar, ma, link1, names_phi, names_theta, names_beta, diag, h, X, X_hat,resid=resid) # model estimation
  
  return(fit1)
}

dibeta<-function(x,alpha0,alpha1,gama,phi)
{
  c<- 1- alpha0*(1-gama)-alpha1*gama
  mu<- gama*(1-alpha1)/c
  
  f0 = ifelse( (x == 0), alpha0*(1-gama), 0)
  f1 = ifelse( (x == 1), (alpha1*gama), 0)
  f2 =  ifelse(((x == 0) | (x == 1)), 0, c*suppressWarnings(dbeta(x, mu * phi, (1 - mu) * phi)))
  
  (f0 + f1 + f2)
}

