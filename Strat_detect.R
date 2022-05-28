

###Version History########

## 2022-MM-DD	Ver 1.0 has released (creator: Shunichiro Orihara, Yokohama City Univ.)

##########################


###READ ME################

## This is a function for the propensity score subclassification estimator (see Orihara and Hamada, 2021).
## Please prepare a propensity score estimator.

## TT: treatment variable, PS: propensity score, YY: outcome variable

## <Limitations>
## Treatment variables are only dichotomous.
## The maximum number of strata is necessary to use the method.

## <Caution>
## The programs has not validated correctly such as double programming (there is only self check).
## The programs are freely available; however, there are no responsibility.

##########################


###PROGRAMS###############

#Functions of selection methods proposed by Orihara and Hamada, 2021.
Strat_func <- function(nn,strat,max_strat){	
#Store of datasets
Data_mat <- cbind(TT,ps,YY)

#Sort by propensity score
Data_mat <- Data_mat[order(Data_mat[,2]),]

kai <- rep(0,strat)
out1 <- out0 <- out1_m <- out0_m <- 0
pw11 <- pw10 <- pw21 <- pw20 <- yy <- yy_m <- c()
ZZ <- ZZ_m <- matrix(0,nrow=nn,ncol=1)
			
for(jj in 1:strat){
if(jj<strat){
#Datasets of jth stratum
  D_group <- Data_mat[(ceiling(nn/strat)*jj-ceiling(nn/strat)+1):(ceiling(nn/strat)*jj),]
  ZZ[(ceiling(nn/strat)*jj-ceiling(nn/strat)+1):(ceiling(nn/strat)*jj),1] <- D_group[,1]
}else{
  D_group <- Data_mat[(ceiling(nn/strat)*strat-ceiling(nn/strat)+1):nn,]
  ZZ[(ceiling(nn/strat)*jj-ceiling(nn/strat)+1):nn,1] <- D_group[,1]
}
#Restore datasets
yy <- c(yy,D_group[,3])

kai[jj] <- mean(D_group[,1])
out1 <- out1+(ceiling(nn/strat)/nn)*(mean(D_group[,3][D_group[,1]==1]))
out0 <- out0+(ceiling(nn/strat)/nn)*(mean(D_group[,3][D_group[,1]==0]))
  		
#Weights for proposed criterion
pw11 <- c(pw11,D_group[,1]/D_group[,2]); pw10 <- c(pw10,(1-D_group[,1])/(1-D_group[,2]))
pw21 <- c(pw21,D_group[,1]/kai[jj]); pw20 <- c(pw20,(1-D_group[,1])/(1-kai[jj]))
}

#Estimate error terms
for(jj in 1:max_strat){
if(jj<max_strat){
  D_group <- Data_mat[(ceiling(nn/max_strat)*jj-ceiling(nn/max_strat)+1):(ceiling(nn/max_strat)*jj),]
  ZZ_m[(ceiling(nn/max_strat)*jj-ceiling(nn/max_strat)+1):(ceiling(nn/max_strat)*jj),1] <- D_group[,1]
}else{
  D_group <- Data_mat[(ceiling(nn/max_strat)*max_strat-ceiling(nn/max_strat)+1):nn,]
  ZZ_m[(ceiling(nn/max_strat)*jj-ceiling(nn/max_strat)+1):nn,1] <- D_group[,1]
}
yy_m <- c(yy_m,D_group[,3])

out1_m <- out1_m+(ceiling(nn/max_strat)/nn)*(mean(D_group[,3][D_group[,1]==1]))
out0_m <- out0_m+(ceiling(nn/max_strat)/nn)*(mean(D_group[,3][D_group[,1]==0]))
}

#Store necessary estimators
beta <- c(c(out0),c(out1-out0))
beta_m <- c(c(out0_m),c(out1_m-out0_m))
ZZ1 <- cbind(1,ZZ); ZZ1_m <- cbind(1,ZZ_m)
PW11 <- diag(pw11); PW10 <- diag(pw10)

#Estimation of error terms	
err1 <- pw11%*%(yy-ZZ1_m%*%beta_m)*pw21%*%(yy-ZZ1_m%*%beta_m)/nn^2
err0 <- pw10%*%(yy-ZZ1_m%*%beta_m)*pw20%*%(yy-ZZ1_m%*%beta_m)/nn^2

#Propose criterion and estimator
crt <- sqrt(t(yy-ZZ1%*%beta)%*%(PW11+PW10)%*%(yy-ZZ1%*%beta)/nn+2*(err1+err0))
out <- out1-out0

kekka <- list(out,crt)
names(kekka) <- c("Causal effects","Criterion")
return(kekka)
}

##########################


###EXAMPLES###############

#rm(list=ls(all.names=T))
library(MASS)

set.seed(1234)
nn <- 1000

#Covariates
mu <- c(0,0,0,0); Sigma <- 1*diag(4)
XX <- mvrnorm(nn,mu,Sigma)

#Treatment
beta_tx <- 1*c(1,-0.5,0.25,0.1,0)
ps <- (1+exp(-(beta_tx[1]*XX[,1]+beta_tx[2]*XX[,2]+beta_tx[3]*XX[,3]+beta_tx[4]*XX[,4]+beta_tx[5])))^(-1)
TT <- rbinom(nn,1,ps)

#Outcomes
beta <- 13.7*1

yy1 <- 210+2*beta*XX[,1]+beta*(XX[,2]+XX[,3]+XX[,4])+rnorm(nn)
yy0 <- 100+2*beta*XX[,1]+beta*(XX[,2]+XX[,3]+XX[,4])+rnorm(nn)

YY=TT*yy1+(1-TT)*yy0

for(kk in 5:11){
  kekka <- Strat_func(nn,kk,11)
  print(c(kk,kekka[[2]],kekka[[1]]))
}

##########################



