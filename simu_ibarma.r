# Created by Fabio M. Bayer (bayer@ufsm.br)
# Department of Statistics, UFSM, Brazil
# December 2022
#
# Reference:
# BAYER, F.M.; PUMI, G.; PEREIRA, T.L.; SOUZA, T.C.
# Inflated Beta Autoregressive Moving Average Models. 
# Forthcoming
# 2023

library(gamlss)

ribeta<- function(n,alpha0,alpha1,gama,prec)
{
  c<- 1- alpha0*(1-gama)-alpha1*gama
  
  mu<- gama*(1-alpha1)/c
  delta0<- alpha0 * (1-gama)
  delta1<- alpha1 * (gama)
  
  p2<- 1- delta0-delta1
  
  if( (alpha0 != 0) && (alpha1 != 0) ) 
    x<- rBEINF(n, mu = mu, sigma = sqrt(1/(prec+1)), nu = delta0/p2, tau = delta1/p2)
  
  if( (alpha0 == 0) && (alpha1 !=0))
    x<- rBEINF1(n, mu = mu, sigma = sqrt(1/(prec+1)), nu = delta1/p2)
  
  if( (alpha0 != 0) && (alpha1 ==0))
    x<- rBEINF0(n, mu = mu, sigma = sqrt(1/(prec+1)), nu = delta0/p2)
  
  if( (alpha0 == 0) && (alpha1 ==0))
    x<- rbeta(n, mu*prec, (1-mu)*prec)
    
  return(x)

}


simu.ibarma <- function(n,phi=NA,theta=NA,alpha0=0.05,alpha1=0.05,
                       alpha=0.0,prec=100,freq=12,link="logit",beta=NA,X=NA)
{
  ar<-NA
  ma<-NA
  
  if(any(is.na(phi)==F))
  {
    ar <- 1:length(phi)
  }
  
  if(any(is.na(theta)==F))
  {
    ma <- 1:length(theta)
  }
  
  linktemp <- substitute(link)
  if (!is.character(linktemp))
  {
    linktemp <- deparse(linktemp)
    if (linktemp == "link")
      linktemp <- eval(link)
  }
  if (any(linktemp == c("logit", "probit", "cloglog")))
  {  
    stats <- make.link(linktemp)
  }else{
    stop(paste(linktemp, "link not available, available links are \"logit\", ",
               "\"probit\" and \"cloglog\""))
  } 
  
  link <- structure(list(link = linktemp, 
                         linkfun = stats$linkfun,
                         linkinv = stats$linkinv
  )
  )
  
  linkfun <- link$linkfun
  linkinv <- link$linkinv
  
  ###### IBARMA model
  if(any(is.na(phi)==F) && any(is.na(theta)==F) && any(is.na(X)==T))
  {
    p <- max(ar)
    q <- max(ma)
    m <- 2*max(p,q)
    
    ynew <-rep(alpha,(n+m))
    gama <- linkinv(ynew)
    
    error<-rep(0,n+m) 
    eta<- y <- rep(0,n+m)
    
    for(i in (m+1):(n+m))
    {
      eta[i]  <- alpha + as.numeric(phi%*%y[i-ar]) + as.numeric(theta%*%error[i-ma])
      gama[i]   <- linkinv(eta[i])
      y[i]    <- ribeta(1,alpha0,alpha1,gama[i],prec)
      error[i] <- y[i]-gama[i]
    }
    
    return( ts(y[(m+1):(n+m)],frequency=freq) )
  } 
  
  ###### IBAR model
  if(any(is.na(phi)==F) && any(is.na(theta)==T) && any(is.na(X)==T))
  {
    p <- max(ar)
    m <- 2*p
    
    ynew <-rep(alpha,(n+m))
    gama <- linkinv(ynew)
    
    error<-rep(0,n+m) 
    eta<- y <- rep(0,n+m)
    
    for(i in (m+1):(n+m))
    {
      eta[i]  <- alpha + as.numeric(phi%*%y[i-ar]) 
      gama[i]   <- linkinv(eta[i])
      y[i]    <- ribeta(1,alpha0,alpha1,gama[i],prec)
      error[i] <- y[i]-gama[i]
    }
    
    return( ts(y[(m+1):(n+m)],frequency=freq) )
  } 
  
  ###### IBMA model
  if(any(is.na(phi)==T) && any(is.na(theta)==F) && any(is.na(X)==T))
  {
    q <- max(ma)
    m <- 2*q
    
    ynew <-rep(alpha,(n+m))
    gama <- linkinv(ynew)
    
    error<-rep(0,n+m) 
    eta<- y <- rep(0,n+m)
    
    for(i in (m+1):(n+m))
    {
      eta[i]  <- alpha + as.numeric(theta%*%error[i-ma])
      gama[i]   <- linkinv(eta[i])
      y[i]    <- ribeta(1,alpha0,alpha1,gama[i],prec)
      error[i] <- y[i]-gama[i]
    }
    
    return( ts(y[(m+1):(n+m)],frequency=freq) )
  } 
  
  ###### IBARMAX model
  if(any(is.na(phi)==F) && any(is.na(theta)==F) && any(is.na(X)==F))
  {
    p <- max(ar)
    q <- max(ma)
    m <- 2*max(p,q)
    
    ynew <-rep(alpha,(n+m))
    gama <- linkinv(ynew)
    
    error<-rep(0,n+m) 
    eta<- y <- rep(0,n+m)
    
    X<-as.matrix(X)

    for(i in (m+1):(n+m))
    {
      eta[i]  <- alpha + as.numeric(X[(i-m),]%*%as.matrix(beta)) + as.numeric(phi%*%y[i-ar]) + as.numeric(theta%*%error[i-ma])
      gama[i]   <- linkinv(eta[i])
      y[i]    <- ribeta(1,alpha0,alpha1,gama[i],prec)
      error[i] <- y[i]-gama[i]
    }
    
    return( ts(y[(m+1):(n+m)],frequency=freq) )
  } 
  
  ###### IBARX model
  if(any(is.na(phi)==F) && any(is.na(theta)==F) && any(is.na(X)==F))
  {
    p <- max(ar)
    m <- 2*p
    
    ynew <-rep(alpha,(n+m))
    gama <- linkinv(ynew)
    
    error<-rep(0,n+m) 
    eta<- y <- rep(0,n+m)
    
    X<-as.matrix(X)
    
    for(i in (m+1):(n+m))
    {
      eta[i]  <- alpha + as.numeric(X[(i-m),]%*%as.matrix(beta)) + as.numeric(phi%*%y[i-ar])
      gama[i]   <- linkinv(eta[i])
      y[i]    <- ribeta(1,alpha0,alpha1,gama[i],prec)
      error[i] <- y[i]-gama[i]
    }
    
    return( ts(y[(m+1):(n+m)],frequency=freq) )
  } 
  
  ###### IBMAX model
  if(any(is.na(phi)==F) && any(is.na(theta)==F) && any(is.na(X)==F))
  {
    q <- max(ma)
    m <- 2*q
    
    ynew <-rep(alpha,(n+m))
    gama <- linkinv(ynew)
    
    error<-rep(0,n+m) 
    eta<- y <- rep(0,n+m)
    
    X<-as.matrix(X)
    
    for(i in (m+1):(n+m))
    {
      eta[i]  <- alpha + as.numeric(X[(i-m),]%*%as.matrix(beta)) +  as.numeric(theta%*%error[i-ma])
      gama[i]   <- linkinv(eta[i])
      y[i]    <- ribeta(1,alpha0,alpha1,gama[i],prec)
      error[i] <- y[i]-gama[i]
    }
    
    return( ts(y[(m+1):(n+m)],frequency=freq) )
  } 
}
