# Created by Fabio M. Bayer (bayer@ufsm.br)
# Department of Statistics, UFSM, Brazil
# December 2022
#
# Description: This is an IBARMA model application example. 
# In this application we consider a time series with 143 observations 
# of percentage of useful volume of the Samuel water reservoir in the 
# Rondonia State, Brazil.
#
# Reference:
# BAYER, F.M.; PUMI, G.; PEREIRA, T.L.; SOUZA, T.C.
# Inflated Beta Autoregressive Moving Average Models. 
# Forthcoming
# 2023

# read the ibarma function
source("ibarma.r")

# read the data
data<-read.table("samuel.txt") 
attach(data)
y1<-V1/100 

# number of forecast steps
h1<-12

N<-length(y1)
y2<-y1[1:(N-h1)]
yout<-y1[(N-(h1-1)):N]
y<-ts(y2,start=c(2011,1),frequency=12) # time series

n<-length(y) # sample size

sum(y==0) # total number of zeros
sum(y==0)/n # percentage of zeros in the sample

plot(y)
acf(y)
pacf(y)

## Example 1: fit the IBARMA(1,1) without covariates
fit1<-ibarma(y,ar=c(1),ma=c(1),h=h1,diag=1,resid=3) # fit the model
names(fit1) # you can see all the returns
fit1$model # it show the fitted model

# Ljung-Box test based on the quantile residuals
Box.test(fit1$resid3, lag = 20, type = c("Ljung-Box"), fitdf = 2)


## Example 2: fit the IBARMA(1,1) with covariates
# seasonal covariates (deterministic seasonality)
m<-5 # for phase adjusting
t <- (1+m):(n+m) # in-sample
t_hat <- (n+1+m):(n+h1+m) # out-of-sample

C<-cos(2*pi*t/12) # in-sample 
C_hat<-cos(2*pi*t_hat/12) # out-of-sample

S<-sin(2*pi*t/12) # in-sample
S_hat<-sin(2*pi*t_hat/12) # out-of-sample

mX<-cbind(S,C) # in-sample
mX_hat<-cbind(S_hat,C_hat) # out-of-sample

fit2<-ibarma(y,ar=c(1),ma=c(2),h=h1,diag=2,X=mX,X_hat=mX_hat,resid=3) # fit the model
fit2$model

# Ljung-Box test based on the quantile residuals
Box.test(fit2$resid3, lag = 20, type = c("Ljung-Box"), fitdf = 2)

