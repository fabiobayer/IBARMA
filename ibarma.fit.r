# Created by Fabio M. Bayer (bayer@ufsm.br)
# Department of Statistics, UFSM, Brazil
# December 2022
#
# Reference:
# BAYER, F.M.; PUMI, G.; PEREIRA, T.L.; SOUZA, T.C.
# Inflated Beta Autoregressive Moving Average Models. 
# Forthcoming
# 2023

ibarma.fit<- function (y, ar, ma, link, names_phi,names_theta,names_beta,diag,h1,X,X_hat,resid)
{
  maxit1<-1000
  
  ###########################
  
  n<- length(y) # sample size
  n0<- sum(y==0) # number of zeros
  n1<- sum(y==1) # number of ones
  
  y_bar<- mean(y)
  
  delta0<- n0/n
  delta1<- n1/n
  
  a0<- delta0/(1-y_bar)
  a1<- delta1/(y_bar)
  
  ys01<- y
  if(any(y==0)) ys01<- ys01[-which(y==0)]
  if(any(y==1)) ys01<- ys01[-which(ys01==1)]
  
  y01<-y
  y01[which(y01==0)]<- min(ys01)
  y01[which(y01==1)]<- max(ys01)
  
  linkfun <- link$linkfun
  linkinv <- link$linkinv
  mu.eta <-  link$mu.eta
  diflink <- link$diflink
  
  ynew01 = linkfun(y01)
  ystar01 = log(y01/(1-y01))
  
  p <- max(ar)
  q <- max(ma)
  m <- max(p,q,na.rm=T)
  p1 <- length(ar)
  q1 <- length(ma)
  
  y1<-y[(m+1):n]
  y_prev01 <- c(rep(NA,(n+h1)))
  
  if(any(is.na(ar)==F)) 
  {
    P <- matrix(rep(NA,(n-m)*p1),ncol=p1)
    
    for(i in 1:(n-m))
    {
      P[i,] <- y[i+m-ar]
    }
    
    Z <- cbind(rep(1,(n-m)),P)
  }else{
    Z <- as.matrix(rep(1,(n-m)))
  }
  
  if(any(is.na(X)==T)) 
  {
    x <- as.matrix(Z)
    Y01 <- y01[(m+1):n]
    Ynew01 = linkfun(Y01)
    Ystar01 = log(Y01/(1-Y01))
    ajuste = lm.fit(x, Ynew01)
    mqo = c(ajuste$coef)
    k = length(mqo)
    ntemp = length(Y01)
    mean = fitted(ajuste)
    mean = linkinv(mean)
    dlink = diflink(mean)
    er = residuals(ajuste)
    sigma2 = sum(er^2)/((ntemp - k) * (dlink)^2)
    prec = 1/ntemp * sum(mean * (1 - mean)/sigma2 - 1)
  }else{ # com regressores
    X_hat<-as.matrix(X_hat)
    X<-as.matrix(X)
    x <- cbind(as.matrix(Z),X[(m+1):n,])
    Y01 <- y01[(m+1):n]
    Ynew01 = linkfun(Y01)
    Ystar01 = log(Y01/(1-Y01))
    ajuste = lm.fit(x, Ynew01)
    mqo = c(ajuste$coef)
    k = length(mqo)
    ntemp = length(Y01)
    mean = fitted(ajuste)
    mean = linkinv(mean)
    dlink = diflink(mean)
    er = residuals(ajuste)
    sigma2 = sum(er^2)/((ntemp - k) * (dlink)^2)
    prec = 1/ntemp * sum(mean * (1 - mean)/sigma2 - 1)
    
  }
  
  ######### ARMA model
  if(any(is.na(ar)==F) && any(is.na(ma)==F) && any(is.na(X)==T))
  { 
    reg <- c(rep(0,(p1+q1+1)), 1.2*prec, a0, a1) # initializing the parameter values
    
    loglik <- function(z) 
    {
      alpha <- z[1]
      phi = z[2:(p1+1)] 
      theta = z[(p1+2):(p1+q1+1)]
      prec <- z[p1+q1+2]
      alpha0 <- z[p1+q1+3]
      alpha1 <- z[p1+q1+4]
      
      error<-rep(0,n) 
      eta<-rep(NA,n)
      
      for(i in (m+1):n)
      {
        eta[i] <- alpha + (phi%*%y[i-ar]) + (theta%*%error[i-ma])
        error[i] <- y[i]-linkinv(eta[i])
      }
      gama <- linkinv(eta[(m+1):n])
      
      alpha0 = ifelse(n0==0,0,alpha0)
      alpha1 = ifelse(n1==0,0,alpha1)
      
      c<- 1- alpha0*(1-gama)-alpha1*gama
      mu<- gama*(1-alpha1)/c
      
      l0 = ifelse( (y1 == 0), log(alpha0*(1-gama)), 0)
      l1 = ifelse( (y1 == 1), log(alpha1*gama), 0)
      l2 =  ifelse(((y1 == 0) | (y1 == 1)), 0, log(c))
      l3 =  ifelse(((y1 == 0) | (y1 == 1)), 0, suppressWarnings(dbeta(y1, mu * prec, (1 - mu) * prec, log = TRUE)))
      
      sum(l0 + l1 + l2+ l3)
    } 
    
    escore <- function(z)
    {
      alpha <- z[1]
      phi = z[2:(p1+1)]
      theta = z[(p1+2):(p1+q1+1)]
      prec <- z[p1+q1+2] 
      alpha0 <- z[p1+q1+3]
      alpha1 <- z[p1+q1+4]
      
      error<-rep(0,n) 
      eta<-rep(0,n)
      
      for(i in (m+1):n)
      {
        eta[i] <- alpha + (phi%*%y[i-ar]) + (theta%*%error[i-ma])
        error[i] <- y[i]-linkinv(eta[i])
      }
      mu <- linkinv(eta[(m+1):n])
      y1<-y[(m+1):n]
      
      alpha0 = ifelse(n0==0,0,alpha0)
      alpha1 = ifelse(n1==0,0,alpha1)

      c <- 1- alpha0*(1-mu)-alpha1*mu
      nu <- mu*(1-alpha1)/c
      
      ystar  <- ifelse(((y1==0)|(y1==1)), 0, log(y1/(1-y1)) )
      mustar <- ifelse(((y1==0)|(y1==1)), 0, digamma(nu * prec) - digamma((1 - nu) * prec) )
      
      v1<- ifelse(((y1==0)|(y1==1)), 0, ((alpha0-alpha1)/c + prec*(ystar-mustar)*((alpha0-1)*(alpha1-1)/(c^2))))
      v2<- ifelse((y1==1), 1/mu, 0) 
      v3<- ifelse((y1==0), 1/(1-mu), 0) 
      v <- v1 + v2 - v3
      
      Q <- matrix(0,nrow=(n-m),ncol=q1)
      
      for(i in 1:(n-m))
      {
        Q[i,]<-error[i+m-ma]
      }

      deta.dalpha <- rep(0,n)
      deta.dphi <- matrix(0, ncol=p1,nrow=n)
      deta.dtheta<- matrix(0, ncol=q1,nrow=n)
      
      for(i in (m+1):n)
      {
        deta.dalpha[i]<- 1 - theta%*% (mu.eta(eta[i-ma]) * deta.dalpha[i-ma])
        deta.dphi[i,]<- P[(i-m),] - theta%*% (mu.eta(eta[i-ma]) *deta.dphi[i-ma,])
        deta.dtheta[i,]<- Q[(i-m),] - theta%*% (mu.eta(eta[i-ma]) *deta.dtheta[i-ma,])
      }
      
      a <- deta.dalpha[(m+1):n]
      rP <- deta.dphi[(m+1):n,]
      rQ <- deta.dtheta[(m+1):n,]
      
      mT <- diag(mu.eta(eta[(m+1):n]))
      
      Ualpha <- a %*% mT %*% v
      Uphi <-   t(rP) %*% mT %*% v
      Utheta <- t(rQ) %*% mT %*% v
      
      Ualpha01 = ifelse((y1 == 0), (1/alpha0), 0)
      Ualpha02 = ifelse(((y1==0)|(y1==1)),0,((mu-1)/c) + ((prec*(ystar-mustar)* mu*(mu-1)*(alpha1-1))/(c^2)))
      Ualpha0 = sum(Ualpha01 + Ualpha02)
      Ualpha0 = ifelse(n0==0,0,Ualpha0)
      
      Ualpha11 = ifelse((y1 == 1), (1/alpha1), 0)
      Ualpha12 = ifelse(((y1==0)|(y1==1)),0,(-mu/c)+(prec*(ystar-mustar)* mu*(1-mu)*(alpha0-1))/(c^2))
      Ualpha1 = sum(Ualpha11 + Ualpha12)
      Ualpha1 = ifelse(n1==0,0,Ualpha1)
      
      Uprec <-  sum(ifelse(((y1 == 0)|(y1 == 1)), 0, nu*(ystar-mustar)+log(1-y1)-digamma((1-nu)*prec)+digamma(prec) ) )

      rval <- c(Ualpha,Uphi,Utheta,Uprec,Ualpha0,Ualpha1)
    }
    
    names_par <- c("alpha",names_phi,names_theta,"precision","alpha0","alpha1")

    opt <- optim(reg, loglik, escore, 
                 method = "BFGS", 
                 control = list(fnscale = -1))
    
    z <- c()
 
    if (opt$conv != 0) 
      warning("FUNCTION DID NOT CONVERGE!")
    
    z$conv <- opt$conv
    coef <- (opt$par)
    names(coef)<-names_par
    z$coeff <- coef
    
    alpha <-coef[1]
    phi <- coef[2:(p1+1)]
    theta <- coef[(p1+2):(p1+q1+1)]
    prec <- coef[p1+q1+2] 
    alpha0 <- coef[p1+q1+3]
    alpha1 <- coef[p1+q1+4]
    
    z$alpha<- alpha
    z$phi <- phi
    z$theta <- theta
    z$prec <- prec
    z$alpha0 <- alpha0
    z$alpha1 <- alpha1
    
    errorhat<-rep(0,n)
    etahat<-rep(0,n)
    
    for(i in (m+1):n)
    {
      etahat[i]<-alpha + (phi%*%y[i-ar]) + (theta%*%errorhat[i-ma])
      errorhat[i] <- y[i]-linkinv(etahat[i]) 
    }
    muhat <- linkinv(etahat[(m+1):n])
    y1<-y[(m+1):n]
    
    z$fitted <- ts(c(rep(NA,m),muhat),start=start(y),frequency=frequency(y))
    z$etahat <- etahat
    z$errorhat <- errorhat
    
    yp <- c(y,rep(NA,h1))
    y_prev <- c(z$fitted,rep(NA,h1))
    eta_prev <- rep(NA,h1)
    
    for(i in 1:h1)
    {
      eta_prev[i] <- alpha + (phi%*%yp[n+i-ar]) + (theta%*%errorhat[n+i-ma])
      y_prev[n+i] <- linkinv(eta_prev[i])
      yp[n+i] <- y_prev[n+i]
      errorhat[n+i] <- 0 
    }
    y_prev <- ts(y_prev, start=start(y), frequency=frequency(y))
    
    c <- 1- alpha0*(1-muhat)-alpha1*muhat
    nu <- muhat*(1-alpha1)/c
    
    ystar  <- ifelse(((y1==0)|(y1==1)), 0, log(y1/(1-y1)) )
    mustar <- ifelse(((y1==0)|(y1==1)), 0, digamma(nu * prec) - digamma((1 - nu) * prec) )
    
    v1<- ifelse(((y1==0)|(y1==1)), 0, ((alpha0-alpha1)/c + prec*(ystar-mustar)*((alpha0-1)*(alpha1-1)/(c^2))))
    v2<- ifelse((y1==1), 1/muhat, 0)
    v3<- ifelse((y1==0), 1/(1-muhat), 0)
    v <- v1 + v2 - v3
    
    Q <- matrix(rep(0,(n-m)*q1),ncol=q1)
    for(i in 1:(n-m))
    {
      Q[i,] <- errorhat[i+m-ma]
    }

    deta.dalpha <- rep(0,n)
    deta.dphi <- matrix(0, ncol=p1,nrow=n)
    deta.dtheta<- matrix(0, ncol=q1,nrow=n)
    
    for(i in (m+1):n)
    {
      deta.dalpha[i]<- 1 - theta%*% (mu.eta(etahat[i-ma]) * deta.dalpha[i-ma])
      deta.dphi[i,]<- P[(i-m),] - theta%*% (mu.eta(etahat[i-ma]) * deta.dphi[i-ma,])
      deta.dtheta[i,]<- Q[(i-m),] - theta%*% (mu.eta(etahat[i-ma]) * deta.dtheta[i-ma,])
    }
    
    a <- deta.dalpha[(m+1):n]
    rP <- deta.dphi[(m+1):n,]
    rQ <- deta.dtheta[(m+1):n,]
    
    s0 <- (((prec^2) *(alpha0-1)*(alpha1-1)^2 * (1-muhat)*muhat)/(c^3))* (trigamma(nu*prec)+trigamma((1-nu)*prec) )+1+((muhat-1)*(alpha0-alpha1)/c)
    s1 <- -(((prec^2) *(alpha0-1)^2 *(alpha1-1) * (1-muhat)*muhat)/(c^3)) * (trigamma(nu*prec)+trigamma((1-nu)*prec) )-1+((-muhat)*(alpha0-alpha1)/c)
    
    mT <- diag(mu.eta(etahat[(m+1):n]))
    
    d <- (((1-alpha0)*(1-alpha1)*prec)/(c^2))*( (1-nu)*trigamma((1-nu)*prec) - nu*trigamma(nu*prec))
      
    dnuda0 <- (alpha1-1)*(muhat-1)*muhat/(c^2)
    dnuda1 <- -(alpha0-1)*(muhat-1)*muhat/(c^2)
    dcda0 <- muhat-1
    dcda1 <- -muhat
    
    dlda0a0 <- (muhat-1)/alpha0 -c*(prec^2)*(trigamma(nu*prec) + trigamma((1-nu)*prec) )*dnuda0*dnuda0-(1/c)*dcda0*dcda0
    dlda1a1 <- (-muhat)/alpha1 -c*(prec^2)*(trigamma(nu*prec) + trigamma((1-nu)*prec) )*dnuda1*dnuda1-(1/c)*dcda1*dcda1
    dlda0a1 <- -c*(prec^2)*(trigamma(nu*prec) + trigamma((1-nu)*prec) )*dnuda0*dnuda1-(1/c)*dcda0*dcda1
    
    dlda0prec <- ((prec*(alpha1-1)*(muhat-1)*muhat)/(c^2))*(-nu*trigamma(nu*prec)+(1-nu)*trigamma((1-nu)*prec))
    dlda1prec <- ((-prec*(alpha0-1)*(muhat-1)*muhat)/(c^2))*(-nu*trigamma(nu*prec)+(1-nu)*trigamma((1-nu)*prec))
    
    dldprec2 <- c*(trigamma(prec)-(nu^2)*trigamma(nu*prec)-((1-nu)^2) * trigamma((1-nu)*prec))
    
    V1 <- (alpha0*muhat + alpha1*(1-muhat))/(muhat*(muhat-1))
    V2 <- ((alpha0-alpha1)^2)/c
    V3 <- ((prec*(alpha0-1)*(alpha1-1))^2)/(c^3)
    V4 <- trigamma(nu*prec)+trigamma((1-nu)*prec)
    V <- diag(V1-V2-V3*V4)
    
    Ka0a0 <- sum(dlda0a0)
    Ka1a1 <- sum(dlda1a1)
    Ka0a1 <- sum(dlda0a1)
    Ka1a0 <- Ka0a1
    
    Ka0prec <- sum(dlda0prec)
    Kpreca0 <- Ka0prec 
    Ka1prec <- sum(dlda1prec)
    Kpreca1 <- Ka1prec 
    
    Ka0a <- t(a) %*% mT %*% s0
    Kaa0 <- Ka0a
    Ka1a <- t(a) %*% mT %*% s1
    Kaa1 <- Ka1a 
    
    Kpa0 <- t(rP) %*% mT %*% s0
    Ka0p <- t(Kpa0)
    
    Kpa1 <- t(rP) %*% mT %*% s1
    Ka1p <- t(Kpa1)
    
    Kta0 <- t(rQ) %*% mT %*% s0
    Ka0t <- t(Kta0)
    
    Kta1 <- t(rQ) %*% mT %*% s1
    Ka1t <- t(Kta1)
    
    Kprecprec <- sum(dldprec2)
    
    Kaprec <- t(a) %*% mT %*% d
    Kpreca <- Kaprec
    
    Kpprec <- t(rP) %*% mT %*% d
    Kprecp <- t(Kpprec)
    
    Ktprec <- t(rQ) %*% mT %*% d
    Kprect <- t(Ktprec)
    
    Kpp <- t(rP) %*% (mT^2) %*% V %*% rP
    
    Ktt <- t(rQ) %*% (mT^2) %*% V %*% rQ
    
    Kpt <- t(rP) %*% (mT^2) %*% V %*% rQ
    Ktp <- t(Kpt)
    
    Kaa <- t(a) %*% (mT^2) %*% V %*% a
    
    Kap <- t(a) %*% (mT^2) %*% V %*% rP
    Kpa <- t(Kap)
    
    Kat <- t(a) %*% (mT^2) %*% V %*% rQ
    Kta <- t(Kat)

    
    K <- -rbind(
      cbind(Kaa,Kap,Kat,Kaprec,Kaa0,Kaa1),
      cbind(Kpa,Kpp,Kpt,Kpprec,Kpa0,Kpa1),
      cbind(Kta,Ktp,Ktt,Ktprec,Kta0,Kta1),
      cbind(Kpreca,Kprecp,Kprect,Kprecprec,Kpreca0,Kpreca1),
      cbind(Ka0a,Ka0p,Ka0t,Ka0prec,Ka0a0,Ka0a1),
      cbind(Ka1a,Ka1p,Ka1t,Ka1prec,Ka1a0,Ka1a1)
    )
    z$K <- K 
  }  
  
  
  ##### only AR model
  if(any(is.na(ar)==F) && any(is.na(ma)==T) && any(is.na(X)==T))
  {
    print("IBAR model is not yet implemented",quote=F)
  }
  
  ##### only MA model
  if(any(is.na(ar)==T) && any(is.na(ma)==F) && any(is.na(X)==T))
  { 
    print("IBMA model is not yet implemented",quote=F)
  }
  
########### BARMAX model #######################
################################################
  
  if(any(is.na(ar)==F) && any(is.na(ma)==F) && any(is.na(X)==F))
  { 
    beta1<- mqo[(p1+2):length(mqo)]
    reg <- c(rep(0,(p1+q1+1)), 1.2*prec, a0, a1, beta1) 
    
    loglik <- function(z) 
    {
      alpha <- z[1]
      phi = z[2:(p1+1)] 
      theta = z[(p1+2):(p1+q1+1)]
      prec <- z[p1+q1+2] 
      alpha0 <- z[p1+q1+3]
      alpha1 <- z[p1+q1+4]
      beta <- z[(p1+q1+5):length(z)]
      
      error<-rep(0,n)
      eta<-rep(NA,n)
      
      for(i in (m+1):n)
      {
        eta[i] <- alpha + X[i,]%*%as.matrix(beta) + (phi%*%y[i-ar]) + (theta%*%error[i-ma])
        error[i] <- y[i]-linkinv(eta[i])
      }
      gama <- linkinv(eta[(m+1):n])
      
      alpha0 = ifelse(n0==0,0,alpha0)
      alpha1 = ifelse(n1==0,0,alpha1)
    
      c<- 1- alpha0*(1-gama)-alpha1*gama
      mu<- gama*(1-alpha1)/c
      
      l0 = ifelse( (y1 == 0), log(alpha0*(1-gama)), 0)
      l1 = ifelse( (y1 == 1), log(alpha1*gama), 0)
      l2 =  ifelse(((y1 == 0) | (y1 == 1)), 0, log(c))
      l3 =  ifelse(((y1 == 0) | (y1 == 1)), 0, suppressWarnings(dbeta(y1, mu * prec, (1 - mu) * prec, log = TRUE)))
      
      sum(l0 + l1 + l2+ l3)
      
    } 
    
    escore <- function(z)
    {
      alpha <- z[1]
      phi = z[2:(p1+1)]
      theta = z[(p1+2):(p1+q1+1)]
      prec <- z[p1+q1+2] 
      alpha0 <- z[p1+q1+3]
      alpha1 <- z[p1+q1+4]
      beta <- z[(p1+q1+5):length(z)]
      
      error<-rep(0,n) 
      eta<-rep(0,n)
      
      for(i in (m+1):n)
      {
        eta[i] <- alpha + X[i,]%*%as.matrix(beta) + (phi%*%y[i-ar]) + (theta%*%error[i-ma])
        error[i] <- y[i]-linkinv(eta[i])
      }
      mu <- linkinv(eta[(m+1):n])
      y1<-y[(m+1):n]
      
      alpha0 = ifelse(n0==0,0,alpha0)
      alpha1 = ifelse(n1==0,0,alpha1)
      
      c <- 1- alpha0*(1-mu)-alpha1*mu
      nu <- mu*(1-alpha1)/c
      
      ystar  <- ifelse(((y1==0)|(y1==1)), 0, log(y1/(1-y1)) )
      mustar <- ifelse(((y1==0)|(y1==1)), 0, digamma(nu * prec) - digamma((1 - nu) * prec) )
      
      v1<- ifelse(((y1==0)|(y1==1)), 0, ((alpha0-alpha1)/c + prec*(ystar-mustar)*((alpha0-1)*(alpha1-1)/(c^2))))
      v2<- ifelse((y1==1), 1/mu, 0) 
      v3<- ifelse((y1==0), 1/(1-mu), 0) 
      v <- v1 + v2 - v3
      
      Q <- matrix(0,nrow=(n-m),ncol=q1)
      
      for(i in 1:(n-m))
      {
        Q[i,]<-error[i+m-ma]
      }
      
      k1<- length(beta)
      R <- matrix(rep(NA,(n-m)*k1),ncol=k1)
      for(i in 1:(n-m))
      {
        for(j in 1:k1)
          R[i,j] <- X[i+m,j]
      }
      
      deta.dalpha <- rep(0,n)
      deta.dphi <- matrix(0, ncol=p1,nrow=n)
      deta.dtheta<- matrix(0, ncol=q1,nrow=n)
      deta.dbeta<- matrix(0, ncol=k1,nrow=n)
      
      for(i in (m+1):n)
      {
        deta.dalpha[i]<- 1 - theta%*% (mu.eta(eta[i-ma]) * deta.dalpha[i-ma])
        deta.dphi[i,]<- P[(i-m),] - theta%*% (mu.eta(eta[i-ma]) *deta.dphi[i-ma,])
        deta.dtheta[i,]<- Q[(i-m),] - theta%*% (mu.eta(eta[i-ma]) *deta.dtheta[i-ma,])
        deta.dbeta[i,]<- R[(i-m),] - theta%*% (mu.eta(eta[i-ma]) *deta.dbeta[i-ma,])
      }
      
      a <- deta.dalpha[(m+1):n]
      rP <- deta.dphi[(m+1):n,]
      rQ <- deta.dtheta[(m+1):n,]
      rR <- deta.dbeta[(m+1):n,]
      
      mT <- diag(mu.eta(eta[(m+1):n]))
      
      Ualpha <- a %*% mT %*% v
      Uphi <-   t(rP) %*% mT %*% v
      Utheta <- t(rQ) %*% mT %*% v
      
      Ualpha01 = ifelse((y1 == 0), (1/alpha0), 0)
      Ualpha02 = ifelse(((y1==0)|(y1==1)),0,((mu-1)/c) + ((prec*(ystar-mustar)* mu*(mu-1)*(alpha1-1))/(c^2)))
      Ualpha0 = sum(Ualpha01 + Ualpha02)
      Ualpha0 = ifelse(n0==0,0,Ualpha0)
      
      Ualpha11 = ifelse((y1 == 1), (1/alpha1), 0)
      Ualpha12 = ifelse(((y1==0)|(y1==1)),0,(-mu/c)+(prec*(ystar-mustar)* mu*(1-mu)*(alpha0-1))/(c^2))
      Ualpha1 = sum(Ualpha11 + Ualpha12)
      Ualpha1 = ifelse(n1==0,0,Ualpha1)
      
      Uprec <-  sum(ifelse(((y1 == 0)|(y1 == 1)), 0, nu*(ystar-mustar)+log(1-y1)-digamma((1-nu)*prec)+digamma(prec) ) )
      
      Ubeta <- t(rR) %*% mT %*% v
      
      rval <- c(Ualpha,Uphi,Utheta,Uprec,Ualpha0,Ualpha1,Ubeta)
    }
    
    names_par <- c("alpha",names_phi,names_theta,"precision","alpha0","alpha1",names_beta)
 
    opt <- optim(reg, loglik, escore, 
                 method = "BFGS", 
                 control = list(fnscale = -1))
    
    z <- c()
    
    if (opt$conv != 0) 
      warning("FUNCTION DID NOT CONVERGE!")
    
    
    z$conv <- opt$conv
    coef <- (opt$par)#[1:(p1+q1+2)]
    names(coef)<-names_par
    z$coeff <- coef
    
    alpha <-coef[1]
    phi <- coef[2:(p1+1)]
    theta <- coef[(p1+2):(p1+q1+1)]
    prec <- coef[p1+q1+2] 
    alpha0 <- coef[p1+q1+3]
    alpha1 <- coef[p1+q1+4]
    beta <- coef[(p1+q1+5):length(coef)]
    
    z$alpha<- alpha
    z$phi <- phi
    z$theta <- theta
    z$prec <- prec
    z$alpha0 <- alpha0
    z$alpha1 <- alpha1
    z$beta <- beta
    
    errorhat<-rep(0,n) 
    etahat<-rep(0,n)
    
    for(i in (m+1):n)
    {
      etahat[i]<-alpha + X[i,]%*%as.matrix(beta) + (phi%*%y[i-ar]) + (theta%*%errorhat[i-ma])
      errorhat[i] <- y[i]-linkinv(etahat[i])
    }
    muhat <- linkinv(etahat[(m+1):n])
    y1<-y[(m+1):n]

    z$fitted <- ts(c(rep(NA,m),muhat),start=start(y),frequency=frequency(y))
    z$etahat <- etahat
    z$errorhat <- errorhat

    yp <- c(y,rep(NA,h1))
    y_prev <- c(z$fitted,rep(NA,h1))
    eta_prev <- rep(NA,h1)

    for(i in 1:h1)
    {
      eta_prev[i] <- alpha + X_hat[i,]%*%as.matrix(beta) + (phi%*%yp[n+i-ar]) + (theta%*%errorhat[n+i-ma])
      y_prev[n+i] <- linkinv(eta_prev[i])
      yp[n+i] <- y_prev[n+i]
      errorhat[n+i] <- 0
    }
    y_prev <- ts(y_prev, start=start(y), frequency=frequency(y))

    c <- 1- alpha0*(1-muhat)-alpha1*muhat
    nu <- muhat*(1-alpha1)/c

    ystar  <- ifelse(((y1==0)|(y1==1)), 0, log(y1/(1-y1)) )
    mustar <- ifelse(((y1==0)|(y1==1)), 0, digamma(nu * prec) - digamma((1 - nu) * prec) )

    v1<- ifelse(((y1==0)|(y1==1)), 0, ((alpha0-alpha1)/c + prec*(ystar-mustar)*((alpha0-1)*(alpha1-1)/(c^2))))
    v2<- ifelse((y1==1), 1/muhat, 0)
    v3<- ifelse((y1==0), 1/(1-muhat), 0)
    v <- v1 + v2 - v3

    Q <- matrix(rep(0,(n-m)*q1),ncol=q1)
    for(i in 1:(n-m))
    {
      Q[i,] <- errorhat[i+m-ma]
    }
    
    k1<- length(beta)
    R <- matrix(rep(NA,(n-m)*k1),ncol=k1)
    for(i in 1:(n-m))
    {
      for(j in 1:k1)
        R[i,j] <- X[i+m,j]
    }

    deta.dalpha <- rep(0,n)
    deta.dphi <- matrix(0, ncol=p1,nrow=n)
    deta.dtheta<- matrix(0, ncol=q1,nrow=n)
    deta.dbeta<- matrix(0, ncol=k1,nrow=n)
 
    for(i in (m+1):n)
    {
      deta.dalpha[i]<- 1 - theta%*% (mu.eta(etahat[i-ma]) * deta.dalpha[i-ma])
      deta.dphi[i,]<- P[(i-m),] - theta%*% (mu.eta(etahat[i-ma]) * deta.dphi[i-ma,])
      deta.dtheta[i,]<- Q[(i-m),] - theta%*% (mu.eta(etahat[i-ma]) * deta.dtheta[i-ma,])
      deta.dbeta[i,]<- R[(i-m),] - theta%*% (mu.eta(etahat[i-ma]) *deta.dbeta[i-ma,])
    }

    a <- deta.dalpha[(m+1):n]
    rP <- deta.dphi[(m+1):n,]
    rQ <- deta.dtheta[(m+1):n,]
    rR <- deta.dbeta[(m+1):n,]

    s0 <- (((prec^2) *(alpha0-1)*(alpha1-1)^2 * (1-muhat)*muhat)/(c^3))* (trigamma(nu*prec)+trigamma((1-nu)*prec) )+1+((muhat-1)*(alpha0-alpha1)/c)
    s1 <- -(((prec^2) *(alpha0-1)^2 *(alpha1-1) * (1-muhat)*muhat)/(c^3)) * (trigamma(nu*prec)+trigamma((1-nu)*prec) )-1+((-muhat)*(alpha0-alpha1)/c)

    mT <- diag(mu.eta(etahat[(m+1):n]))

    d <- (((1-alpha0)*(1-alpha1)*prec)/(c^2))*( (1-nu)*trigamma((1-nu)*prec) - nu*trigamma(nu*prec))

    dnuda0 <- (alpha1-1)*(muhat-1)*muhat/(c^2)
    dnuda1 <- -(alpha0-1)*(muhat-1)*muhat/(c^2)
    dcda0 <- muhat-1
    dcda1 <- -muhat

    dlda0a0 <- (muhat-1)/alpha0 -c*(prec^2)*(trigamma(nu*prec) + trigamma((1-nu)*prec) )*dnuda0*dnuda0-(1/c)*dcda0*dcda0
    dlda1a1 <- (-muhat)/alpha1 -c*(prec^2)*(trigamma(nu*prec) + trigamma((1-nu)*prec) )*dnuda1*dnuda1-(1/c)*dcda1*dcda1
    dlda0a1 <- -c*(prec^2)*(trigamma(nu*prec) + trigamma((1-nu)*prec) )*dnuda0*dnuda1-(1/c)*dcda0*dcda1

    dlda0prec <- ((prec*(alpha1-1)*(muhat-1)*muhat)/(c^2))*(-nu*trigamma(nu*prec)+(1-nu)*trigamma((1-nu)*prec))
    dlda1prec <- ((-prec*(alpha0-1)*(muhat-1)*muhat)/(c^2))*(-nu*trigamma(nu*prec)+(1-nu)*trigamma((1-nu)*prec))

    dldprec2 <- c*(trigamma(prec)-(nu^2)*trigamma(nu*prec)-((1-nu)^2) * trigamma((1-nu)*prec))

    V1 <- (alpha0*muhat + alpha1*(1-muhat))/(muhat*(muhat-1))
    V2 <- ((alpha0-alpha1)^2)/c
    V3 <- ((prec*(alpha0-1)*(alpha1-1))^2)/(c^3)
    V4 <- trigamma(nu*prec)+trigamma((1-nu)*prec)
    V <- diag(V1-V2-V3*V4)

    Ka0a0 <- sum(dlda0a0)
    Ka1a1 <- sum(dlda1a1)
    Ka0a1 <- sum(dlda0a1)
    Ka1a0 <- Ka0a1

    Ka0prec <- sum(dlda0prec)
    Kpreca0 <- Ka0prec
    Ka1prec <- sum(dlda1prec)
    Kpreca1 <- Ka1prec

    Ka0a <- t(a) %*% mT %*% s0
    Kaa0 <- Ka0a
    Ka1a <- t(a) %*% mT %*% s1
    Kaa1 <- Ka1a

    Kpa0 <- t(rP) %*% mT %*% s0
    Ka0p <- t(Kpa0)

    Kpa1 <- t(rP) %*% mT %*% s1
    Ka1p <- t(Kpa1)

    Kta0 <- t(rQ) %*% mT %*% s0
    Ka0t <- t(Kta0)

    Kta1 <- t(rQ) %*% mT %*% s1
    Ka1t <- t(Kta1)

    Kprecprec <- sum(dldprec2)

    Kaprec <- t(a) %*% mT %*% d
    Kpreca <- Kaprec

    Kpprec <- t(rP) %*% mT %*% d
    Kprecp <- t(Kpprec)

    Ktprec <- t(rQ) %*% mT %*% d
    Kprect <- t(Ktprec)

    Kpp <- t(rP) %*% (mT^2) %*% V %*% rP

    Ktt <- t(rQ) %*% (mT^2) %*% V %*% rQ

    Kpt <- t(rP) %*% (mT^2) %*% V %*% rQ
    Ktp <- t(Kpt)

    Kaa <- t(a) %*% (mT^2) %*% V %*% a

    Kap <- t(a) %*% (mT^2) %*% V %*% rP
    Kpa <- t(Kap)

    Kat <- t(a) %*% (mT^2) %*% V %*% rQ
    Kta <- t(Kat)

    Kba0 <- t(rR) %*% mT %*% s0
    Ka0b <- t(Kba0)
    
    Kba1 <- t(rR) %*% mT %*% s1
    Ka1b <- t(Kba1)
    
    Kbprec <- t(rR) %*% mT %*% d
    Kprecb <- t(Kbprec)
    
    Kab <- t(a) %*% (mT^2) %*% V %*% rR
    Kba <- t(Kab)
    
    Kbb <- t(rR) %*% (mT^2) %*% V %*% rR
    
    Kpb <- t(rP) %*% (mT^2) %*% V %*% rR
    Kbp <- t(Kpb)
    
    Kbt <- t(rR) %*% (mT^2) %*% V %*% rQ
    Ktb <- t(Kbt)
    
    K <- -rbind(
      cbind(Kaa,Kap,Kat,Kaprec,Kaa0,Kaa1,Kab),
      cbind(Kpa,Kpp,Kpt,Kpprec,Kpa0,Kpa1,Kpb),
      cbind(Kta,Ktp,Ktt,Ktprec,Kta0,Kta1,Ktb),
      cbind(Kpreca,Kprecp,Kprect,Kprecprec,Kpreca0,Kpreca1,Kprecb),
      cbind(Ka0a,Ka0p,Ka0t,Ka0prec,Ka0a0,Ka0a1,Ka0b),
      cbind(Ka1a,Ka1p,Ka1t,Ka1prec,Ka1a0,Ka1a1,Ka1b),
      cbind(Kba,Kbp,Kbt,Kbprec,Kba0,Kba1,Kbb)
    )
    
    z$K <- K
    
  }  
  
  ######### BARX model
  if(any(is.na(ar)==F) && any(is.na(ma)==T) && any(is.na(X)==F))
  { 
    print("BARX model is not yet implemented",quote=F)
  }  
  
  ######### BMAX model
  if(any(is.na(ar)==T) && any(is.na(ma)==F) && any(is.na(X)==F))
  { 
    print("BMAX model is not yet implemented",quote=F)
  }  
  
  z$serie <- y
  z$barma <- names_par
  z$forecast <- y_prev[(n+1):(n+h1)]
  
  c<- 1-z$alpha0*(1-z$fitted[(m+1):n])-z$alpha1*z$fitted[(m+1):n]
  nu <- z$fitted[(m+1):n]*(1-z$alpha1)/c
  
  vary <- ((1+z$alpha1*z$prec)/(1+z$prec))*z$fitted[(m+1):n]+((((1-z$alpha1)^2)*z$prec)/(c*(1+z$prec))-1)*(z$fitted[(m+1):n]^2)
  
  resid1 = (y1 - z$fitted[(m+1):n])/sqrt(vary)
  
  resid2 <- sign(y1 - z$fitted[(m+1):n])*sqrt(2*abs(log(dibeta(y1,z$alpha0,z$alpha1,y1,z$prec))-log(dibeta(y1,z$alpha0,z$alpha1,z$fitted[(m+1):n],z$prec)) ) )

  FA = z$alpha0*(1-z$fitted[(m+1):n])
  
  FB = ifelse( (y1 == 1), z$alpha1*z$fitted[(m+1):n], 0)
  
  FBeta = c*pbeta(y1,nu*z$prec,(1-nu)*z$prec)
  
  FBinf = FA+FB+FBeta
  
  at0=0
  bt0=z$alpha0*(1-z$fitted[(m+1):n])
  ut0=ifelse((y1 == 0), runif((n-m),at0,bt0),0)
  
  at1=(1-z$alpha1)*z$fitted[(m+1):n]
  bt1=1
  ut1=ifelse((y1 == 1), runif((n-m),at1,bt1),0)
  
  ut01=ifelse( ((y1==0)|(y1==1)),0,FBinf)
  
  ut<-ut0+ut1+ut01
  
  resid3 = qnorm(ut)
  
  z$resid1 <- resid1
  z$resid2 <- resid2
  z$resid3 <- resid3
  
  if(resid==1) residc <- z$resid1
  if(resid==2) residc <- z$resid2
  if(resid==3) residc <- z$resid3

  z$y <- y
  z$m <- m

  Kchol<- tryCatch(chol(K), error = function(e) return("error"), warning = function(o) return("error"))

  if(Kchol[1]== "error")
  {
    z$KOK <- FALSE
    z$vcov <- K
    warning("We have problems with information matrix inversion!")
    
    stderror <- rep(1,length(z$coef))
    z$stderror <- stderror
    
  }else{
    z$KOK <- TRUE
    vcov <- try(chol2inv(Kchol))
    z$vcov <- vcov
    stderror <- sqrt(diag(vcov))
    z$stderror <- stderror
  }

  z$zstat <- z$coef/stderror
  z$pvalues <- 2*(1 - pnorm( abs(z$zstat) ) )

  z$loglik <- opt$value*(n/(n-m))
  z$counts <- as.numeric(opt$counts[1])
  z$aic <- -2*z$loglik+2*(length(z$coef))
  z$bic <- -2*z$loglik+log(n)*(length(z$coef))


  model_presentation <- cbind(round(z$coef,4),round(z$stderror,4),round(z$zstat,4),round(z$pvalues,4))
  colnames(model_presentation)<-c("Estimate","Std. Error","z value","Pr(>|z|)")

  z$model <- model_presentation
  
  z$link <- link
  
  z$rank <- length(z$coef)
  
  if(diag>0)
  {
    print(model_presentation)
    print(" ",quote=F)
    print(c("Log-likelihood:",round(z$loglik,4)),quote=F)
    print(c("Number of iterations in BFGS optim:",z$counts),quote=F)
    print(c("AIC:",round(z$aic,4)," BIC:",round(z$bic,4)),quote=F)
    
    print("Residuals:",quote=F)
    print(summary(residc))
    
    t<-seq(-5,n+6,by=1)
    
    par(mfrow=c(1,1))
    par(mar=c(2.8, 2.7, 1.2, 1)) 
    par(mgp=c(1.7, 0.45, 0))
    plot(residc,main=" ",xlab="index",ylab="residuals", pch = "+",ylim=c(-4,4))
    lines(t,rep(-3,n+12),lty=2,col=1)
    lines(t,rep(3,n+12),lty=2,col=1)
    lines(t,rep(-2,n+12),lty=3,col=1)
    lines(t,rep(2,n+12),lty=3,col=1)
    
    
    max_y<- max(c(z$fitted,y),na.rm=T)
    min_y<- min(c(z$fitted,y),na.rm=T)
    plot(as.vector(z$fitted), as.vector(y), main=" ", pch = "+",
         xlab="fitted values",ylab="actual values",
         xlim=c(0.95*min_y,max_y*1.05),
         ylim=c(0.95*min_y,max_y*1.05))
    lines(c(-0.2,1.2),c(-0.2,1.2),lty=2)
    
    densidade<-density(residc)
    plot(densidade,ylab="density",main=" ")
    lines(densidade$x,dnorm(densidade$x),lty=2)
    legend("topleft",c("Estimated density","Normal density"),#pch=vpch,
           pt.bg="white", lty=c(1,2), bty="n")
    
    acf(residc,ylab="ACF",xlab="lag") 
    
    pacf(residc,ylab="PACF",xlab="lag") 
    
    max_r<- max(residc,na.rm=T)
    min_r<- min(residc,na.rm=T)
    qqnorm(residc, pch = "+",
           xlim=c(0.95*min_r,max_r*1.05),
           ylim=c(0.95*min_r,max_r*1.05),
           main="",xlab="normal quantiles",ylab="empirical quantiles")
    lines(c(-10,10),c(-10,10),lty=2)
    
    par(mfrow=c(1,1))
    plot(y,type="l",ylab="data",xlab="time")
    lines(z$fitted,col="red")
    
    fim<-end(y)[1]+end(y)[2]/12
    
    par(mfrow=c(1,1))
    par(mar=c(2.8, 2.7, 1, 1)) 
    par(mgp=c(1.7, 0.45, 0))
    plot(y,type="l",ylab="data",xlab="time")
    lines(z$fitted,col="red")

    par(mfrow=c(1,1))
    par(mar=c(2.8, 2.7, 1, 1)) 
    par(mgp=c(1.7, 0.45, 0))
    plot(y_prev,type="l",col="red", ylim=c(min(y),max(y)),ylab="data",xlab="time")
    abline(v=fim,lty=2)
    lines(y)
    
    w1<-3
    h1<-3
    
    if(diag>1)
    {
      pdf(file = "resid_v_ind.pdf",width = w1, height = h1,family = "Times")
      {
        par(mfrow=c(1,1))
        par(mar=c(2.8, 2.7, 1, 1)) 
        par(mgp=c(1.7, 0.45, 0))
        plot(residc,main=" ",xlab="index",ylab="residuals", pch = "+",ylim=c(-4,4))
        lines(t,rep(-3,n+12),lty=2,col=1)
        lines(t,rep(3,n+12),lty=2,col=1)
        lines(t,rep(-2,n+12),lty=3,col=1)
        lines(t,rep(2,n+12),lty=3,col=1)
      }
      dev.off()
      
      pdf(file = "obs_v_fit.pdf",width = w1, height = h1,family = "Times")
      {
        par(mfrow=c(1,1))
        par(mar=c(2.8, 2.7, 1, 1))
        par(mgp=c(1.7, 0.45, 0))
        plot(as.vector(z$fitted), as.vector(y), main=" ", pch = "+",
             xlab="fitted values",ylab="actual values",
             xlim=c(0.95*min_y,max_y*1.05),
             ylim=c(0.95*min_y,max_y*1.05))
        lines(c(-0.2,1.2),c(-0.2,1.2),lty=2)
      }
      dev.off()
      
      pdf(file = "resid_density.pdf",width = w1, height = h1,family = "Times")
      {
        par(mfrow=c(1,1))
        par(mar=c(1.5, 2.7, 1, 1)) 
        par(mgp=c(1.7, 0.45, 0))
        
        plot(densidade,ylab="density",main=" ",xlab=" ",ylim=c(0,1.15*max(densidade$y)))
        lines(densidade$x,dnorm(densidade$x),lty=2)
        legend("topleft",c("estimated density","normal density"),#pch=vpch,
               pt.bg="white", lty=c(1,2), bty="n")
      }
      dev.off()
      
      pdf(file = "resid_FAC.pdf",width = w1, height = h1,family = "Times")
      {
        par(mfrow=c(1,1))
        par(mar=c(2.8, 2.7, 1, 1)) 
        par(mgp=c(1.7, 0.45, 0))
        acf(residc,ylab="ACF",xlab="lag") 
      }
      dev.off()
      
      pdf(file = "resid_FACP.pdf",width = w1, height = h1,family = "Times")
      {
        par(mfrow=c(1,1))
        par(mar=c(2.8, 2.7, 1, 1)) # 
        par(mgp=c(1.7, 0.45, 0))
        pacf(residc,ylab="PACF",xlab="lag") 
      }
      dev.off()
      
      pdf(file = "qq_plot.pdf",width = w1, height = h1,family = "Times")
      {  
        par(mfrow=c(1,1))
        par(mar=c(2.8, 2.7, 1, 1)) 
        par(mgp=c(1.7, 0.45, 0))
        qqnorm(residc, pch = "+",
               xlim=c(0.95*min_r,max_r*1.05),
               ylim=c(0.95*min_r,max_r*1.05),
               main="",xlab="normal quantiles",ylab="empirical quantiles")
        lines(c(-10,10),c(-10,10),lty=2)
      }
      dev.off()
      
      pdf(file = "adjusted.pdf",width = 6, height = 4,family = "Times")
      {
        par(mfrow=c(1,1))
        par(mar=c(2.8, 2.7, 1, 1)) 
        par(mgp=c(1.7, 0.45, 0))
        plot(y,type="l",ylab="data",xlab="time")
        lines(z$fitted,col="red")
      }
      dev.off()
      
      pdf(file = "forecast.pdf",width = 6, height = 4,family = "Times")
      {
        par(mfrow=c(1,1))
        par(mar=c(2.8, 2.7, 1, 1)) 
        par(mgp=c(1.7, 0.45, 0))
        plot(y_prev,type="l",col="red",lty=2, ylim=c(min(y),max(y)),ylab="data",xlab="time")
        abline(v=fim,lty=2)
        lines(y)
      }
      dev.off()
    }    
  }  
  return(z)
}