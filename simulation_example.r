# Created by Fabio M. Bayer (bayer@ufsm.br)
# Department of Statistics, UFSM, Brazil
# December 2022
#
# Description: Simulation study of IBARMA model.
#
# Reference:
# BAYER, F.M.; PUMI, G.; PEREIRA, T.L.; SOUZA, T.C.
# Inflated Beta Autoregressive Moving Average Models. 
# Forthcoming
# 2023

source('ibarma.r') # read the ibarma function
source('simu_ibarma.r') # read the simu.ibarma function

#set.seed(10)
R<-100 # number of replications
pa<-c(-0.5,1.5,-1,20,0.07,0.08,1) # parameter values
vn<-c(100,300,500) # sample sizes
h1<-12 # number of forecast steps

for(n in vn) # loop for sample sizes
{
  x_hat1<-runif(h1)
  X<-runif(n)
  TC<-rep(0,length(pa))
  coef_results<-c()
  results<-c()
  i<-1
  nc<-0
  while (i <= R) # Monte Carlo loop
  {
    y<-simu.ibarma(n,alpha=pa[1],phi=pa[2],theta=pa[3],alpha0=pa[5],alpha1=pa[6],prec=pa[4],beta=pa[7],X=X)
    fit1<-ibarma(y,ar=c(1),ma=c(1),diag=0,X=X,X_hat=x_hat1,h=h1)
    
    if(fit1$conv==0)
    {
      results<-rbind(results,fit1$pvalues)
      coef_results<-rbind(coef_results,fit1$coef)

      LI<- fit1$coef - qnorm(0.975)* fit1$stderror
      LS<- fit1$coef + qnorm(0.975)* fit1$stderror
      
      TC <- TC + ( (pa<LI) + (pa>LS))
      
      i<-i+1
    }else{
      nc<-nc+1
      print(c("non convergence",i,nc))
    }
  }

  m<-colMeans((coef_results),na.rm=T)
  sd<-apply(coef_results,2,sd)
  bias<-m-pa
  rb<- 100*bias/pa
  tc <- 1-TC/R
  
  M<-rbind(pa,m,sd,bias,rb,tc)
  row.names(M)<-c("Parameters","Mean","SE","Bias","RB","CR")
  
  print(c("n",n),quote=F)
  print(round(M,3))
}



