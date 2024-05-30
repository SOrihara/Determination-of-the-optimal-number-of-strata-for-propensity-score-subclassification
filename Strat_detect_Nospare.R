


#Functions of selection methods proposed by Orihara and Hamada, 2021.
Strat_func <- function(AA, ps, YY, strat, max_strat){
  nn <- length(AA)
  
  ## Storing datasets
  Data_mat <- cbind(AA,ps,YY)
    
  #Sort by propensity score
  Data_mat <- Data_mat[order(Data_mat[,2]),]
  
  
  out <- out_m <- 0
  var <- var_m <- 0
  pw11 <- pw10 <- pw21 <- pw20 <- yy <- yy_m <- c()

  for(jj in 1:strat){
    if(jj<strat){
      #Datasets of jth stratum
      D_group <- Data_mat[(ceiling(nn/strat)*jj-ceiling(nn/strat)+1):(ceiling(nn/strat)*jj),]
    }else{
      D_group <- Data_mat[(ceiling(nn/strat)*strat-ceiling(nn/strat)+1):nn,]
    }
    #Restore datasets
    kk <- summary(lm(D_group[,3] ~ D_group[,1]))
    kk <- kk$coef[2,1:2]
    
    out <- out + kk[1]
    var <- var + kk[2]
  }
  
  #Estimate error terms
  for(jj in 1:max_strat){
    if(jj<max_strat){
      #Datasets of jth stratum
      D_group <- Data_mat[(ceiling(nn/max_strat)*jj-ceiling(nn/max_strat)+1):(ceiling(nn/max_strat)*jj),]
    }else{
      D_group <- Data_mat[(ceiling(nn/max_strat)*strat-ceiling(nn/max_strat)+1):nn,]
    }
    #Restore datasets
    kk <- summary(lm(D_group[,3] ~ D_group[,1]))
    kk <- kk$coef[2,1:2]
    
    out_m <- out_m + kk[1]
    var_m <- var_m + kk[2]
  }
  
  
  #Propose criterion and estimator
  crt <- var/strat^2 + (out/strat-out_m/max_strat)^2
  
  kekka <- list(as.numeric(out/strat), as.numeric(var/strat^2), as.numeric(crt))
  names(kekka) <- c("Causal effects","Standard error","Criterion")
  return(kekka)
}







