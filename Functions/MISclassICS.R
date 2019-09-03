library(rootSolve)



### solve naive expected estimating equations
### Xki should contain intercept, true exposure Xk (and NAs for main sample), and correctly measured covariates
### Xk is the correctly observed exposure (and NAs for main sample)
fit_naive2stage_R <- function(Xki,Zxki,Nk_actual,Yki,Xk,IDk,WEIGHTk=NA,VALIDk,start_param=NA,nboot=100){

  ## make unit level
  VALIDki <- rep(VALIDk,Nk_actual)
  X1ki <- rep(Xk,Nk_actual)
  if(length(WEIGHTk)!=length(VALIDk)){
    WEIGHTki <- rep(1,length(VALIDki))
  }else{
    WEIGHTki <- rep(WEIGHTk,Nk_actual)
  }

  ## NEEE estimates
  est <- est_naive2stage_R(Xki=Xki,Zxki=Zxki,Yki=Yki,X1ki=X1ki,WEIGHTki=WEIGHTki,VALIDki=VALIDki)



  ## bootstrap
  if(!is.na(nboot) & nboot>0){
    IDki <- 1:nrow(Xki) ## identify rows
    IDS <- data.frame(IDki=IDki,IDk=rep(IDk,Nk_actual)) ## matrix to identify resampled rows

    boot_ests <- c()
    for(bb in 1:nboot){
      valid_bootIDk <- sample(unique(IDk[VALIDk==1]),replace=TRUE) ## resample validation clusters
      main_bootIDk <- sample(unique(IDk[VALIDk==0]),replace=TRUE)  ## resample main sample clusters
      valid_bootID <- data.frame(IDk=valid_bootIDk) ## make dataframe for merge
      main_bootID <- data.frame(IDk=main_bootIDk) ## make dataframe for merge

      # which rows in resample
      valid_bootIDki <- merge(valid_bootID, IDS, by="IDk")$IDki
      main_bootIDki <- merge(main_bootID, IDS, by="IDk")$IDki
      ## row IDs in bootstrap resample
      bootIDki <- c(valid_bootIDki,main_bootIDki)

      ## fit to sample
      b_est <- est_naive2stage_R(Xki=Xki[bootIDki,],Zxki=Zxki[bootIDki,],Yki=Yki[bootIDki],X1ki=X1ki[bootIDki],WEIGHTki=WEIGHTki[bootIDki],VALIDki=VALIDki[bootIDki])

      ## collect results
      boot_ests <- rbind(boot_ests,b_est)
    }
    SE <- apply(boot_ests,2,sd,na.rm=TRUE)
  }else{
    SE <- rep(NA,length(est))
  }

  return(list(ests=est,SEs=SE))

}



est_naive2stage_R <- function(Xki,Zxki,Yki,X1ki,WEIGHTki,VALIDki,start_param=NA){

  ## stage 1
  epsX <- glm(X1ki[VALIDki==1]~Zxki[VALIDki==1,]-1,family="binomial")$coef

  ## stage 2
  if(length(start_param)!=(ncol(Xki))){start_param <- rep(0,ncol(Xki))}
  return(multiroot(naive2stage_R,start=start_param,epsX=epsX,Xki=Xki,Zxki=Zxki,Yki=Yki,WEIGHTki=WEIGHTki,VALIDki=VALIDki)$root)

}



### function to compute naive expected estimating equations
naive2stage_R <- function(params,epsX,Xki,Zxki,Yki,WEIGHTki,VALIDki){


  ## covariate matrices
  Xki0 <- Xki1 <- Xki
  Xki0[,2] <- 0
  Xki1[,2] <- 1

  ### collect parameters
  beta <- params

  ### stage 1
  Px1ki <- expit(Zxki%*%epsX)

  ### stage 2

  ## complete data estimating functions for outcome model
  Ski_valid <- WEIGHTki*Xki*c(Yki-expit(Xki%*%beta))
  ## expected
  S0ki <- WEIGHTki*Xki0*c(Yki-expit(Xki0%*%beta))
  S1ki <- WEIGHTki*Xki1*c(Yki-expit(Xki1%*%beta))
  Ski_main <- S0ki*c(1-Px1ki)+S1ki*c(Px1ki)

  ## combine validation or main components
  Ski <- (VALIDki)*Ski_valid+(1-VALIDki)*Ski_main
  EEE <- c(apply(Ski,2,sum,na.rm=TRUE))

  ## return EEE
  return(EEE)


}

















### solve naive expected estimating equations
fit_naiveEEE_R <- function(Zki,Zxki,Nk_actual,Yki,Xk,IDk,WEIGHTk=NA,VALIDk,start_param=NA){

  if(length(start_param)!=(ncol(Zki)+1+ncol(Zxki))){start_param <- rep(0,ncol(Zki)+1+ncol(Zxki))}
  ## solve for estimates
  est <- multiroot(naiveEEE_R,start=start_param,Zki=Zki,Zxki=Zxki,Nk_actual=Nk_actual,Yki=Yki,Xk=Xk,IDk=IDk,WEIGHTk=WEIGHTk,VALIDk=VALIDk)$root
  ## compute inverse of hessian
  bread <- tryCatch(solve(jacobian.full(y = est, func = naiveEEE_hessR,Zki=Zki,Zxki=Zxki,Nk_actual=Nk_actual,Yki=Yki,Xk=Xk,IDk=IDk,WEIGHTk=WEIGHTk,VALIDk=VALIDk)),error=function(cond) NA)
  if(is.na(bread)){
    SE <- rep(NA,length(est))
  }else{
    ## compute middle of sandwich
    cheese <-naiveEEE_cheeseR(est,Zki=Zki,Zxki=Zxki,Nk_actual=Nk_actual,Yki=Yki,Xk=Xk,IDk=IDk,WEIGHTk=WEIGHTk,VALIDk=VALIDk)
    ## sandwich form variance-covariance matrix
    VAR <- bread%*%cheese%*%bread
    SE <- sqrt(diag(VAR))
  }

  ## return estimates and SEs
  return(list(est=est,SE=SE))


}

### function to compute naive expected estimating equations
naiveEEE_R <- function(params,Zki,Zxki,Nk_actual,Yki,Xk,IDk,WEIGHTk=NA,VALIDk){

  ## make unit level
  VALIDki <- rep(VALIDk,Nk_actual)
  if(length(WEIGHTk)!=length(VALIDk)){
    WEIGHTki <- rep(1,length(VALIDki))
  }else{
    WEIGHTki <- rep(WEIGHTk,Nk_actual)
  }

  ## covariate matrices
  Xki <- cbind(Zki[,1],rep(Xk,Nk_actual),Zki[,2:ncol(Zki)])
  X0ki <- cbind(Zki[,1],0,Zki[,2:ncol(Zki)])
  X1ki <- cbind(Zki[,1],1,Zki[,2:ncol(Zki)])

  ### collect parameters
  beta <- params[1:(ncol(Zki)+1)] ## outcome model parameters
  alphaX <- params[(ncol(Zki)+2):length(params)] ## exposure model parameters (second column of Zki must be W)

  ## complete data estimating functions for outcome model
  Ski_valid <- WEIGHTki*Xki*c(Yki-expit(Xki%*%beta))
  ## expected
  S0ki <- WEIGHTki*X0ki*c(Yki-expit(X0ki%*%beta))
  S1ki <- WEIGHTki*X1ki*c(Yki-expit(X1ki%*%beta))
  Px1ki <- expit(Zxki%*%alphaX)
  Ski_main <- S0ki*c(1-Px1ki)+S1ki*c(Px1ki)

  ## compete data estimating functions for exposure model (reduces to 0 in main sample)
  Sxki_valid <- Zxki*c(Xki[,2]-expit(Zxki%*%alphaX))

  ## combine validation or main components
  Ski <- (VALIDki)*Ski_valid+(1-VALIDki)*Ski_main
  Sxki <- (VALIDki)*Sxki_valid
  EEE <- c(apply(Ski,2,sum,na.rm=TRUE),
           apply(Sxki,2,sum,na.rm=TRUE))

  ## return EEE
  return(EEE)


}


### function to compute naive expected estimating equations (but reparametrized for hessian)
naiveEEE_hessR <- function(time=0,params,Zki,Zxki,Nk_actual,Yki,Xk,IDk,WEIGHTk=NA,VALIDk,parms=NULL){

  return(naiveEEE_R(params,Zki,Zxki,Nk_actual,Yki,Xk,IDk,WEIGHTk=NA,VALIDk))

}

### function to compute CHEESE for naive EEE
naiveEEE_cheeseR <- function(params,Zki,Zxki,Nk_actual,Yki,Xk,IDk,WEIGHTk=NA,VALIDk){

  ## make unit level
  VALIDki <- rep(VALIDk,Nk_actual)
  if(length(WEIGHTk)!=length(VALIDk)){
    WEIGHTki <- rep(1,length(VALIDki))
  }else{
    WEIGHTki <- rep(WEIGHTk,Nk_actual)
  }

  ## covariate matrices
  Xki <- cbind(Zki[,1],rep(Xk,Nk_actual),Zki[,2:ncol(Zki)])
  X0ki <- cbind(Zki[,1],0,Zki[,2:ncol(Zki)])
  X1ki <- cbind(Zki[,1],1,Zki[,2:ncol(Zki)])

  ### collect parameters
  beta <- params[1:(ncol(Zki)+1)] ## outcome model parameters
  alphaX <- params[(ncol(Zki)+2):length(params)] ## exposure model parameters (second column of Zki must be W)

  ## complete data estimating functions for outcome model
  Ski_valid <- WEIGHTki*Xki*c(Yki-expit(Xki%*%beta)); Ski_valid[is.na(Ski_valid)] <- 0
  ## expected
  S0ki <- WEIGHTki*X0ki*c(Yki-expit(X0ki%*%beta))
  S1ki <- WEIGHTki*X1ki*c(Yki-expit(X1ki%*%beta))
  Px1ki <- expit(Zxki%*%alphaX)
  Ski_main <- S0ki*c(1-Px1ki)+S1ki*c(Px1ki)

  ## compete data estimating functions for exposure model (reduces to 0 in main sample)
  Sxki_valid <- Zxki*c(Xki[,2]-expit(Zxki%*%alphaX)); Sxki_valid[is.na(Sxki_valid)] <- 0

  ## combine validation or main components
  Ski <- (VALIDki)*Ski_valid+(1-VALIDki)*Ski_main
  Sxki <- (VALIDki)*Sxki_valid

  EEEk <- cbind(Ski,Sxki)
  return(t(EEEk)%*%EEEk)

}

# ### solve naive expected estimating equations
# fit_naive2stage_R <- function(Zki,Zxki,Nk_actual,Yki,Xk,IDk,WEIGHTk=NA,VALIDk,start_param=NA){
#
#   ## make unit level
#   VALIDki <- rep(VALIDk,Nk_actual)
#   if(length(WEIGHTk)!=length(VALIDk)){
#     WEIGHTki <- rep(1,length(VALIDki))
#   }else{
#     WEIGHTki <- rep(WEIGHTk,Nk_actual)
#   }
#
#   ## covariate matrices
#   Xki <- cbind(Zki[,1],rep(Xk,Nk_actual),Zki[,2:ncol(Zki)])
#   X0ki <- cbind(Zki[,1],0,Zki[,2:ncol(Zki)])
#   X1ki <- cbind(Zki[,1],1,Zki[,2:ncol(Zki)])
#
#   ### collect parameters
#   beta <- params[1:(ncol(Zki)+1)] ## outcome model parameters
#   alphaX <- params[(ncol(Zki)+2):length(params)] ## exposure model parameters (second column of Zki must be W)
#
#   ## stage 1
#   alphaX <- glm(rep(Xk,Nk_actual)[VALIDki==1]~Zxki[VALIDki==1,]-1,family="binomial")$coef
#
#
#   if(length(start_param)!=(ncol(Zki)+1)){start_param <- rep(0,ncol(Zki)+1)}
#
#   return(multiroot(naive2stage_R,start=start_param,alphaX=alphaX,Zki=Zki,Zxki=Zxki,Nk_actual=Nk_actual,Yki=Yki,Xk=Xk,IDk=IDk,WEIGHTk=WEIGHTk,VALIDk=VALIDk))
#
#
# }
#
# ### function to compute naive expected estimating equations
# naive2stage_R <- function(params,alphaX,Zki,Zxki,Nk_actual,Yki,Xk,IDk,WEIGHTk=NA,VALIDk){
#
#   ## make unit level
#   VALIDki <- rep(VALIDk,Nk_actual)
#   if(length(WEIGHTk)!=length(VALIDk)){
#     WEIGHTki <- rep(1,length(VALIDki))
#   }else{
#     WEIGHTki <- rep(WEIGHTk,Nk_actual)
#   }
#
#   ## covariate matrices
#   Xki <- cbind(Zki[,1],rep(Xk,Nk_actual),Zki[,2:ncol(Zki)])
#   X0ki <- cbind(Zki[,1],0,Zki[,2:ncol(Zki)])
#   X1ki <- cbind(Zki[,1],1,Zki[,2:ncol(Zki)])
#
#   ### collect parameters
#   # beta <- params[1:(ncol(Zki)+1)] ## outcome model parameters
#   # alphaX <- params[(ncol(Zki)+2):length(params)] ## exposure model parameters (second column of Zki must be W)
#   beta <- params
#
#   ### stage 1
#   Px1ki <- expit(Zxki%*%alphaX)
#
#   ### stage 2
#
#   ## complete data estimating functions for outcome model
#   Ski_valid <- WEIGHTki*Xki*c(Yki-expit(Xki%*%beta))
#   ## expected
#   S0ki <- WEIGHTki*X0ki*c(Yki-expit(X0ki%*%beta))
#   S1ki <- WEIGHTki*X1ki*c(Yki-expit(X1ki%*%beta))
#   Ski_main <- S0ki*c(1-Px1ki)+S1ki*c(Px1ki)
#
#   ## combine validation or main components
#   Ski <- (VALIDki)*Ski_valid+(1-VALIDki)*Ski_main
#   EEE <- c(apply(Ski,2,sum,na.rm=TRUE))
#
#   ## return EEE
#   return(EEE)
#
#
# }

### solve naive expected estimating equations
fit_partialEEE_R <- function(Zki,Zwki,Znki,Zxki,Nk_actual,Yki,Xk,Wk,IDk,minNk=0,WEIGHTk=NA,VALIDk,start_param=NA){

  if(length(start_param)!=((ncol(Zki)+1)+(ncol(Zwki)+1)+(ncol(Znki)+1)+ncol(Zxki))){start_param <- rep(0,((ncol(Zki)+1)+(ncol(Zwki)+1)+(ncol(Znki)+1)+ncol(Zxki)))}

  return(multiroot(partialEEE_R,start=start_param,Zki=Zki,Zwki=Zwki,Znki=Znki,Zxki=Zxki,Nk_actual=Nk_actual,Yki=Yki,Xk=Xk,Wk=Wk,IDk=IDk,minNk=minNk,WEIGHTk=WEIGHTk,VALIDk=VALIDk))


}


### function to compute partial expected estimating equations
partialEEE_R <- function(params,Zki,Zwki,Znki,Zxki,Nk_actual,Yki,Xk,Wk,IDk,minNk=0,WEIGHTk=NA,VALIDk){

  ## make unit level
  VALIDki <- rep(VALIDk,Nk_actual)
  Nki <- rep(Nk_actual,Nk_actual)-minNk
  Xki_true <- rep(Xk,Nk_actual)
  Wki_mis <- rep(Wk,Nk_actual)
  if(length(WEIGHTk)!=length(VALIDk)){
    WEIGHTki <- rep(1,length(VALIDki))
  }else{
    WEIGHTki <- rep(WEIGHTk,Nk_actual)
  }

  ## covariate matrices -- outcome model
  Xki <- cbind(Zki[,1],rep(Xk,Nk_actual));if(ncol(Zki)>1) {Xki <- cbind(Xki,Zki[,2:ncol(Zki)])}
  X0ki <- cbind(Zki[,1],0);               if(ncol(Zki)>1) {X0ki <- cbind(X0ki,Zki[,2:ncol(Zki)])}
  X1ki <- cbind(Zki[,1],1);               if(ncol(Zki)>1) {X1ki <- cbind(X1ki,Zki[,2:ncol(Zki)])}
  ## covariate matrices -- misclassification model
  Xwki <- cbind(Zwki[,1],rep(Xk,Nk_actual));if(ncol(Zwki)>1) {Xwki <- cbind(Xwki,Zwki[,2:ncol(Zwki)])}
  Xw0ki <- cbind(Zwki[,1],0);               if(ncol(Zwki)>1) {Xw0ki <- cbind(Xw0ki,Zwki[,2:ncol(Zwki)])}
  Xw1ki <- cbind(Zwki[,1],1);               if(ncol(Zwki)>1) {Xw1ki <- cbind(Xw1ki,Zwki[,2:ncol(Zwki)])}
  ## covariate matrices -- size model
  Xnki <- cbind(Znki[,1],rep(Xk,Nk_actual));if(ncol(Znki)>1) {Xnki <- cbind(Xnki,Znki[,2:ncol(Znki)])}
  Xn0ki <- cbind(Znki[,1],0);               if(ncol(Znki)>1) {Xn0ki <- cbind(Xn0ki,Znki[,2:ncol(Znki)])}
  Xn1ki <- cbind(Znki[,1],1);               if(ncol(Znki)>1) {Xn1ki <- cbind(Xn1ki,Znki[,2:ncol(Znki)])}
  ## covariate matrices -- exposure model
  Xxki <- Zxki

  ### collect parameters
  beta <- params[1:(ncol(Zki)+1)] ## outcome model parameters
  alphaW <- params[(ncol(Zki)+1+1):((ncol(Zki)+1)+(ncol(Zwki)+1))] ## misclassification parameters
  alphaN <- params[((ncol(Zki)+1)+(ncol(Zwki)+1)+1):((ncol(Zki)+1)+(ncol(Zwki)+1)+(ncol(Znki)+1))] ## size model parameters
  alphaX <- params[((ncol(Zki)+1)+(ncol(Zwki)+1)+(ncol(Znki)+1)+1):length(params)] ## exposure model parameters (second column of Zki must be W)

  ## complete data estimating functions
  Ski_valid <- WEIGHTki*Xki*c(Yki-expit(Xki%*%beta)) ## outcome model
  Swki_valid <- Xwki*c(Wki_mis-expit(Xwki%*%alphaW)) ## outcome model
  Snki_valid <- Xnki*c(Nki-exp(Xnki%*%alphaN)) ## outcome model
  Sxki_valid <- Xxki*c(Xki_true-expit(Xxki%*%alphaX)) ## outcome model

  ## components for expected estimating functions
  S0ki <- WEIGHTki*X0ki*c(Yki-expit(X0ki%*%beta));  S1ki <- WEIGHTki*X1ki*c(Yki-expit(X1ki%*%beta))
  Sw0ki <- Xw0ki*c(Wki_mis-expit(Xw0ki%*%alphaW));  Sw1ki <- Xw1ki*c(Wki_mis-expit(Xw1ki%*%alphaW));
  Sn0ki <- Xn0ki*c(Nki-exp(Xn0ki%*%alphaN));        Sn1ki <- Xn1ki*c(Nki-exp(Xn1ki%*%alphaN));
  Sx0ki <- Xxki*c(0-expit(Xxki%*%alphaX));          Sx1ki <- Xxki*c(1-expit(Xxki%*%alphaX));
  ## Distribution P(X|Y,W,N,X)
  P1ki <- (Yki*expit(X1ki%*%beta)+(1-Yki)*(1-expit(X1ki%*%beta)))*
          (Wki_mis*expit(Xw1ki%*%alphaW)+(1-Wki_mis)*(1-expit(Xw1ki%*%alphaW)))*
          (exp(-exp(Xn1ki%*%alphaN))*exp(Nki*Xn1ki%*%alphaN)/factorial(Nki))*
          (expit(Xxki%*%alphaX))
  P0ki <- (Yki*expit(X0ki%*%beta)+(1-Yki)*(1-expit(X0ki%*%beta)))*
          (Wki_mis*expit(Xw0ki%*%alphaW)+(1-Wki_mis)*(1-expit(Xw0ki%*%alphaW)))*
          (exp(-exp(Xn0ki%*%alphaN))*exp(Nki*Xn0ki%*%alphaN)/factorial(Nki))*
          (1-expit(Xxki%*%alphaX))
  P1ki <- P1ki/(P1ki+P0ki) ## normalizing (to get conditional)
  ## marginalizing
  Ski_main <- S0ki*c(1-P1ki)+S1ki*c(P1ki)
  Swki_main <- Sw0ki*c(1-P1ki)+Sw1ki*c(P1ki)
  Snki_main <- Sn0ki*c(1-P1ki)+Sn1ki*c(P1ki)
  Sxki_main <- Sx0ki*c(1-P1ki)+Sx1ki*c(P1ki)

  ## combine validation or main components
  ## outcome model
  Ski <- Ski_main
  Ski[VALIDki==1,] <- Ski_valid[VALIDki==1,]
  Ski <- apply(Ski,2,sum)
  ## misclassification model
  Swki <- Swki_main
  Swki[VALIDki==1,] <- Swki_valid[VALIDki==1,]
  Swki <- apply(Swki,2,sum)
  ## size model
  Snki <- Snki_main
  Snki[VALIDki==1,] <- Snki_valid[VALIDki==1,]
  Snki <- apply(Snki,2,sum)
  ## exposure model
  Sxki <- Sxki_main
  Sxki[VALIDki==1,] <- Sxki_valid[VALIDki==1,]
  Sxki <- apply(Sxki,2,sum)

  ## combine and return
  EEE <- c(Ski,Swki,Snki,Sxki)
  return(EEE)


}
