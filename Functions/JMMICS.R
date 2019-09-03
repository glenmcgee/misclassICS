############################################################
## Wrapper function for joint marginalized models for ICS ##
############################################################
## Permits Poisson, Negative-Binomial, or ZIP size models (does not combine zero inflation with negative binomial)
## notation differs slightly from paper: Nk0 is 1(Nk==0), Z0k is Wk (ie covariates for zero model), and theta is Wk epsilon

## APRIL, 2019 EDIT: ## NEWER VERSION THAN IN JMMICSpack (paper2)
                     ## includes startpos=FALSE; TRUE sets starting gamma value to 1 (rather than -1)

## random slopes parametrization: sigma0 is the sd of intercepts among unexposed, sigma1 among exposed

require(gaussquad)
require(Rcpp)
if(runLOCAL==TRUE){
  install.packages("../Functions/JMMICSpack_0.1.0.tar.gz",repos=NULL)
}
library(JMMICSpack)


## get conditionally specified parameters (size model)
get_Omega <- function(gamma,sigma,zeta,condSize=FALSE,w,z){
  ## get Omega (implicitly defined conditional linear predictor for size model)
  if (condSize == FALSE){
    return(compute_Omega(gamma,sigma,zeta,w,z))
  } else{
    return(zeta)
  }
}

## get conditionally specified parameters (outcome model)
get_Delta <- function(sigma,eta,condOut,w,z){
  ## get Delta (implicitly defined conditional linear predictor for outcome model)
  if (condOut == FALSE){
    return(compute_Delta(sigma,eta,w,z))
  } else{
    return(eta)
  }
}

extract_params <- function(params,Nk,Nk_actual,Nk0,Zk,Yki,Xki,IDk,Z0k,Wki,minNk=0,NegBin=FALSE,ZIP=FALSE,slopes=FALSE,slope_col=2,condSize=FALSE,condOut=FALSE,joint=TRUE,w,z){


  if(joint==TRUE){
    if(slopes==FALSE){ ## random intercepts

      alpha <- params[1:ncol(Zk)]
      gamma <- params[ncol(Zk)+1]
      beta <- params[(ncol(Zk)+2):(ncol(Zk)+ncol(Xki)+1)]
      sigma <- exp(params[(ncol(Zk)+ncol(Xki)+2)]) # maximize wrt logsigma
      if(NegBin==TRUE){
        tau <- exp(params[length(params)]) # maximize wrt logtau
      } else{
        tau <- NA ## not used
      }
      if(ZIP==TRUE){
        epsilon <- params[(ncol(Zk)+ncol(Xki)+3):(length(params))]
      }else{
        epsilon <- NA ## not used
      }
      ## convert gamma and sigma to vectors
      gamma_vec <- rep(gamma,length(Nk))
      sigma_vec <- rep(sigma,length(Nk))

      rand_indic <- NA
    }else{ ## random slopes

      alpha <- params[1:ncol(Zk)]
      gamma0 <- params[ncol(Zk)+1]
      gamma1 <- params[ncol(Zk)+2]
      beta <- params[(ncol(Zk)+3):(ncol(Zk)+ncol(Xki)+2)]
      sigma0 <- exp(params[(ncol(Zk)+ncol(Xki)+3)]) # maximize wrt logsigma
      sigma1 <- exp(params[(ncol(Zk)+ncol(Xki)+4)]) # maximize wrt logsigma
      if(NegBin==TRUE){
        tau <- exp(params[length(params)]) # maximize wrt logtau
      } else{
        tau <- NA ## not used
      }
      if(ZIP==TRUE){
        epsilon <- params[(ncol(Zk)+ncol(Xki)+5):(length(params))]
      }else{
        epsilon <- NA ## not used
      }
      ## random slopes indicator (ie binary covariate for random slope)
      rand_indic <- Zk[,slope_col]

      ## convert gamma and sigma to vectors (for random slopes)
      gamma_vec <- rep(gamma0,length(Nk))
      gamma_vec[rand_indic==1] <- gamma1
      sigma_vec <- rep(sigma0,length(Nk))
      sigma_vec[rand_indic==1] <- sigma1
    }



    ## compute eta, zeta, possibly theta
    zeta <- Zk%*%alpha
    eta <- Xki%*%beta
    if(ZIP==TRUE){
      theta <- Z0k%*%epsilon
    }else{theta=NA}

    ## get implicitly-defined conditional params
    Omega <- get_Omega(gamma_vec,sigma_vec,zeta,condSize,w,z)
    #Delta <- get_Delta(sigma_vec,eta,condOut,w,z)
    Delta <- get_Delta(rep(sigma_vec,Nk_actual),eta,condOut,w,z) ## new version needs separate sigma for each cluster,
                                                               ## but Delta is observation-wise, so we need a separate sigma for each obs
    if(ZIP==TRUE){
      Xi <- theta #get_Xi(nu,sigma_vec,theta,condZero,w,z)
    }else{Xi=NA}

    return(list(alpha=alpha,gamma=gamma_vec,beta=beta,sigma=sigma_vec,tau=tau,epsilon=epsilon,
                zeta=zeta,eta=eta,theta=theta,Omega=Omega,Delta=Delta,Xi=Xi,
                rand_indic=rand_indic))

  }else{ ## outcome-only model
    if(slopes==FALSE){ ## random intercepts

      beta <- params[1:ncol(Xki)]
      sigma <- exp(params[(ncol(Xki)+1)]) ## maximize wrt logsigma
      sigma_vec <- rep(sigma,length(Nk)) ## convert sigma to vector
      rand_indic <- NA ## no random slope

    }else{ ## random slopes

      beta <- params[1:ncol(Xki)]
      sigma0 <- exp(params[(ncol(Xki)+1)]) ## maximize wrt logsigma
      sigma1 <- exp(params[(ncol(Xki)+2)]) ## maximize wrt logsigma
      rand_indic <- Zk[,slope_col]       ## random slopes indicator (ie binary covariate for random slope)
      sigma_vec <- rep(sigma0,length(Nk))
      sigma_vec[rand_indic==1] <- sigma1 ## convert sigma to vector (for random slopes)

    }

    ## compute eta
    eta <- Xki%*%beta
    ## get implicitly-defined conditional params
    Delta <- get_Delta(rep(sigma_vec,Nk_actual),eta,condOut,w,z) ## new version needs separate sigma for each cluster,
                                                                 ## but Delta is observation-wise, so we need a separate sigma for each obs

    return(list(beta=beta,sigma=sigma_vec,
                eta=eta,Delta=Delta,
                rand_indic=rand_indic))

  }

}


### function to compute log-likelihood
logLik_JMMICS <- function(params,Nk,Nk_actual,Nk0,Zk,Yki,Xki,IDk,Z0k,Wki,minNk=0,NegBin=FALSE,ZIP=FALSE,slopes=FALSE,slope_col=2,condSize=FALSE,condOut=FALSE,joint=TRUE,w,z){

  comp <- extract_params(params=params,Nk=Nk,Nk_actual=Nk_actual,Nk0=Nk0,Zk=Zk,Yki=Yki,Xki=Xki,IDk=IDk,Z0k=Z0k,Wki=Wki,minNk=minNk,NegBin=NegBin,ZIP=ZIP,slopes=slopes,slope_col=slope_col,condSize=condSize,condOut=condOut,joint=joint,w=w,z=z)

  if(joint==TRUE){ ## joint model

    if(NegBin==TRUE){
      return(-compute_logLik_JMMICS_NegBin(comp$gamma,comp$sigma,comp$tau,comp$Omega,comp$Delta,Nk,Nk_actual,Yki,IDk,Wki,minNk,w,z))
    }else if(ZIP==TRUE){
      return(-compute_logLik_JMMICS_ZIP(comp$gamma,comp$sigma,comp$Xi,comp$Omega,comp$Delta,Nk,Nk_actual,Nk0,Yki,IDk,Wki,minNk,w,z))
    }else{
      return(-compute_logLik_JMMICS(comp$gamma,comp$sigma,comp$Omega,comp$Delta,Nk,Nk_actual,Yki,IDk,Wki,minNk,w,z))
    }
  }else{ ## outcome-only

    return(-compute_logLik_GLMM(comp$sigma,comp$Delta,Nk_actual,Yki,IDk,Wki,w,z))
  }
}


### function to compute score
score_JMMICS <- function(params,Nk,Nk_actual,Nk0,Zk,Yki,Xki,IDk,Z0k,Wki,minNk=0,NegBin=FALSE,ZIP=FALSE,slopes=FALSE,slope_col=2,condSize=FALSE,condOut=FALSE,joint=TRUE,w,z){

  comp <- extract_params(params=params,Nk=Nk,Nk_actual=Nk_actual,Nk0=Nk0,Zk=Zk,Yki=Yki,Xki=Xki,IDk=IDk,Z0k=Z0k,Wki=Wki,minNk=minNk,NegBin=NegBin,ZIP=ZIP,slopes=slopes,slope_col=slope_col,condSize=condSize,condOut=condOut,joint=joint,w=w,z=z)

  if(joint==TRUE){ ## joint model

    if(slopes==FALSE){ ## random intercepts
      if(NegBin==TRUE){
        return(-compute_score_JMMICS_NegBin(comp$gamma,comp$sigma,tau,DG_tau=digamma(comp$tau),comp$Omega,comp$Delta,comp$zeta,comp$eta,Nk,Nk_actual,Zk,Yki,Xki,IDk,Wki,minNk,condSize,condOut,w,z,length(comp$alpha),length(comp$beta)))
      }else if(ZIP==TRUE){
        return(-compute_score_JMMICS_ZIP(comp$gamma,comp$sigma,comp$Xi,comp$Omega,comp$Delta,comp$theta,comp$zeta,comp$eta,Nk,Nk_actual,Nk0,Z0k,Zk,Yki,Xki,IDk,Wki,minNk,condSize,condOut,w,z,length(comp$epsilon),length(comp$alpha),length(comp$beta)))
      }else{
        return(-compute_score_JMMICS(comp$gamma,comp$sigma,comp$Omega,comp$Delta,comp$zeta,comp$eta,Nk,Nk_actual,Zk,Yki,Xki,IDk,Wki,minNk,condSize,condOut,w,z,length(comp$alpha),length(comp$beta)))
      }
    }else{ ## random slopes
      if(NegBin==TRUE){
        return(-compute_score_JMMICS_NegBin_slopes(comp$gamma,comp$sigma,tau,DG_tau=digamma(comp$tau),comp$Omega,comp$Delta,comp$zeta,comp$eta,comp$rand_indic,Nk,Nk_actual,Zk,Yki,Xki,IDk,Wki,minNk,condSize,condOut,w,z,length(comp$alpha),length(comp$beta)))
      }else if(ZIP==TRUE){
        return(-compute_score_JMMICS_ZIP_slopes(comp$gamma,comp$sigma,comp$Xi,comp$Omega,comp$Delta,comp$theta,comp$zeta,comp$eta,comp$rand_indic,Nk,Nk_actual,Nk0,Z0k,Zk,Yki,Xki,IDk,Wki,minNk,condSize,condOut,w,z,length(comp$epsilon),length(comp$alpha),length(comp$beta)))
      }else{
        return(-compute_score_JMMICS_slopes(comp$gamma,comp$sigma,comp$Omega,comp$Delta,comp$zeta,comp$eta,comp$rand_indic,Nk,Nk_actual,Zk,Yki,Xki,IDk,Wki,minNk,condSize,condOut,w,z,length(comp$alpha),length(comp$beta)))
      }
    }
  }else{ ## outcome only model

    if(slopes==FALSE){ ## random intercepts
      return(-compute_score_GLMM(comp$sigma,comp$Delta,comp$eta,Nk_actual,Yki,Xki,IDk,Wki,condOut,w,z,length(comp$beta)))
    }else{ ## random slopes
      return(-compute_score_GLMM_slopes(comp$sigma,comp$Delta,comp$eta,comp$rand_indic,Nk_actual,Yki,Xki,IDk,Wki,condOut,w,z,length(comp$beta)))
    }
  }

}

### function to compute cheese
cheese_JMMICS <- function(params,Nk,Nk_actual,Nk0,Zk,Yki,Xki,IDk,Z0k,Wki,minNk=0,NegBin=FALSE,ZIP=FALSE,slopes=FALSE,slope_col=2,condSize=FALSE,condOut=FALSE,joint=TRUE,w,z){

  comp <- extract_params(params=params,Nk=Nk,Nk_actual=Nk_actual,Nk0=Nk0,Zk=Zk,Yki=Yki,Xki=Xki,IDk=IDk,Z0k=Z0k,Wki=Wki,minNk=minNk,NegBin=NegBin,ZIP=ZIP,slopes=slopes,slope_col=slope_col,condSize=condSize,condOut=condOut,joint=joint,w=w,z=z)

  if(joint==TRUE){ ## joint model
    if(slopes==FALSE){ ## random intecepts
      if(NegBin==TRUE){
        return(compute_cheese_JMMICS_NegBin(comp$gamma,comp$sigma,comp$tau,DG_tau=digamma(comp$tau),comp$Omega,comp$Delta,comp$zeta,comp$eta,Nk,Nk_actual,Zk,Yki,Xki,IDk,Wki,minNk,condSize,condOut,w,z,length(comp$alpha),length(comp$beta)))
      }else if(ZIP==TRUE){
        return(compute_cheese_JMMICS_ZIP(comp$gamma,comp$sigma,comp$Xi,comp$Omega,comp$Delta,comp$theta,comp$zeta,comp$eta,Nk,Nk_actual,Nk0,Z0k,Zk,Yki,Xki,IDk,Wki,minNk,condSize,condOut,w,z,length(comp$epsilon),length(comp$alpha),length(comp$beta)))
      }else{
        return(compute_cheese_JMMICS(comp$gamma,comp$sigma,comp$Omega,comp$Delta,comp$zeta,comp$eta,Nk,Nk_actual,Zk,Yki,Xki,IDk,Wki,minNk,condSize,condOut,w,z,length(comp$alpha),length(comp$beta)))
      }
    }else{ ## random slopes
      if(NegBin==TRUE){
        return(compute_cheese_JMMICS_NegBin_slopes(comp$gamma,comp$sigma,comp$tau,DG_tau=digamma(comp$tau),comp$Omega,comp$Delta,comp$zeta,comp$eta,comp$rand_indic,Nk,Nk_actual,Zk,Yki,Xki,IDk,Wki,minNk,condSize,condOut,w,z,length(comp$alpha),length(comp$beta)))
      }else if(ZIP==TRUE){
        return(compute_cheese_JMMICS_ZIP_slopes(comp$gamma,comp$sigma,comp$Xi,comp$Omega,comp$Delta,comp$theta,comp$zeta,comp$eta,comp$rand_indic,Nk,Nk_actual,Nk0,Z0k,Zk,Yki,Xki,IDk,Wki,minNk,condSize,condOut,w,z,length(comp$epsilon),length(comp$alpha),length(comp$beta)))
      }else{
        return(compute_cheese_JMMICS_slopes(comp$gamma,comp$sigma,comp$Omega,comp$Delta,comp$zeta,comp$eta,comp$rand_indic,Nk,Nk_actual,Zk,Yki,Xki,IDk,Wki,minNk,condSize,condOut,w,z,length(comp$alpha),length(comp$beta)))
      }
    }
  }else{ ## outcome-only
    if(slopes==FALSE){ ## random intercepts
      return(compute_cheese_GLMM(comp$sigma,comp$Delta,comp$eta,Nk_actual,Yki,Xki,IDk,Wki,condOut,w,z,length(comp$beta)))
    }else{ ## random slopes
      return(compute_cheese_GLMM_slopes(comp$sigma,comp$Delta,comp$eta,comp$rand_indic,Nk_actual,Yki,Xki,IDk,Wki,condOut,w,z,length(comp$beta)))
    }

  }
}


### function to compute bks
EB_JMMICS <- function(params,Nk,Nk_actual,Nk0,Zk,Yki,Xki,IDk,Z0k,Wki,minNk=0,NegBin=FALSE,ZIP=FALSE,slopes=FALSE,slope_col=2,condSize=FALSE,condOut=FALSE,joint=TRUE,w,z){

  comp <- extract_params(params=params,Nk=Nk,Nk_actual=Nk_actual,Nk0=Nk0,Zk=Zk,Yki=Yki,Xki=Xki,IDk=IDk,Z0k=Z0k,Wki=Wki,minNk=minNk,NegBin=NegBin,ZIP=ZIP,slopes=slopes,slope_col=slope_col,condSize=condSize,condOut=condOut,joint=joint,w=w,z=z)

  if(joint==TRUE){ ## joint model

    if(NegBin==TRUE){
      return(compute_EB_JMMICS_NegBin(comp$gamma,comp$sigma,comp$tau,comp$Omega,comp$Delta,Nk,Nk_actual,Yki,IDk,Wki,minNk,w,z))
    }else if(ZIP==TRUE){
      return(compute_EB_JMMICS_ZIP(comp$gamma,comp$sigma,comp$Xi,comp$Omega,comp$Delta,Nk,Nk_actual,Nk0,Yki,IDk,Wki,minNk,w,z))
    }else{
      return(compute_EB_JMMICS(comp$gamma,comp$sigma,comp$Omega,comp$Delta,Nk,Nk_actual,Yki,IDk,Wki,minNk,w,z))
    }
  }else{ ## outcome-only

    return(compute_EB_GLMM(comp$sigma,comp$Delta,Nk_actual,Yki,IDk,Wki,w,z))
  }
}

## JMMICS wrapper function
## slopes==TRUE is for random slopes
  ## slope_col indicates the column of Zk that refers to the binary covariate for random slopes
JMMICS_fit <- function(Nk,Zk,Yki,Xki,IDk,Z0k=NULL,weights=NA,minNk=0,NegBin=FALSE,ZIP=FALSE,slopes=FALSE,slope_col=2,condSize=FALSE,condOut=FALSE,joint=TRUE,startpos=FALSE,nquad=25){
  ## turn GH weights/zeros into appropriately scaled terms
  w <- as.matrix(hermite.h.quadrature.rules(nquad)[[nquad]])[,2]; w <- w/sqrt(pi)
  z <- as.matrix(hermite.h.quadrature.rules(nquad)[[nquad]])[,1]; z <- sqrt(2)*z

  ## handle empty clusters
  Nk0 <- as.numeric(Nk==0) ## indicator of emptiness
  if(length(weights)<length(Yki)){
    Wki <- rep(1,length(Yki))
    Nk_actual <- Nk
  }else{
    Wki <- weights
    Nk_obs <- c(table(IDk))
    ## below added in case not all clusters are sampled (like in CS or RS)
    Nk_actual <- rep(0,length(Nk))
    Nk_actual[as.numeric(names(Nk_obs))] <- Nk_obs
  }

  ## starting value for gamma
  if(startpos==TRUE){
    startgam <- 1
  }else{startgam <- -1}

  if(joint==TRUE){ ## usual joint model
    ## alpha +gamma +beta +sigma
    ##  if negbin --> +possible tau
    ##  if ZIP -->  +possible epsilon +possible nu
    if(slopes==FALSE){ ## random intercepts
      startval <- c(rep(1,ncol(Zk)),startgam,rep(0,ncol(Xki)),log(0.5))
      resnames <- c(paste0(rep("alpha",ncol(Zk)),seq(0,ncol(Zk)-1)),"gamma",paste0(rep("beta",ncol(Xki)),seq(0,ncol(Xki)-1)),"sigma") ## names of parameters to print
      sig_id <- ncol(Zk)+1+ncol(Xki)+1 ## index for logsigma
    }else{ ## random slopes (for binary covariate)
      startval <- c(rep(1,ncol(Zk)),startgam,startgam,rep(0,ncol(Xki)),log(0.5),log(0.5))
      resnames <- c(paste0(rep("alpha",ncol(Zk)),seq(0,ncol(Zk)-1)),"gamma0","gamma1",paste0(rep("beta",ncol(Xki)),seq(0,ncol(Xki)-1)),"sigma0","sigma1") ## names of parameters to print
      sig_id <- ncol(Zk)+2+ncol(Xki)+1:2  ## index for logsigma0 and logsigma1
    }
    if(NegBin==TRUE){
      startval <- c(startval,0)
      resnames <- c(resnames,"tau")
    }
    if(ZIP==TRUE){
      startval <- c(startval,rep(0,ncol(Z0k)))
      resnames <- c(resnames,paste0(rep("epsilon",ncol(Z0k)),seq(0,ncol(Z0k)-1)))
    }
  }else{ ## outcome-only model
    NegBin <- FALSE
    ZIP <- FALSE
    if(slopes==FALSE){ ## random intercepts
      startval <- c(rep(0,ncol(Xki)),log(0.5))
      resnames <- c(paste0(rep("beta",ncol(Xki)),seq(0,ncol(Xki)-1)),"sigma") ## names of parameters to print
      sig_id <- ncol(Xki)+1 ## index for logsigma
    }else{ ## random slopes (for binary covariate)
      startval <- c(rep(0,ncol(Xki)),log(0.5),log(0.5))
      resnames <- c(paste0(rep("beta",ncol(Xki)),seq(0,ncol(Xki)-1)),"sigma0","sigma1") ## names of parameters to print
      sig_id <- ncol(Xki)+1:2  ## index for logsigma0 and logsigma1
    }
  }

    ## optimization
    opt <- nlminb(start=startval,logLik_JMMICS,gradient=score_JMMICS,Nk=Nk,Nk_actual=Nk_actual,Nk0=Nk0,Zk=Zk,Yki=Yki,Xki=Xki,IDk=IDk,Z0k=Z0k,Wki=Wki,minNk=minNk,NegBin=NegBin,ZIP=ZIP,slopes=slopes,slope_col=slope_col,condSize=condSize,condOut=condOut,joint=joint,w=w,z=z,control=list(eval.max=2000,iter.max=1500))$par
    est <- opt; ## ML estimates
    est[sig_id] <- exp(est[sig_id]) ## transform logsigma to sigma


    ## Standard Errors (with Delta Method)
    invI <- try(solve(optimHess(par=opt,logLik_JMMICS,gr=score_JMMICS,Nk=Nk,Nk_actual=Nk_actual,Nk0=Nk0,Zk=Zk,Yki=Yki,Xki=Xki,IDk=IDk,Z0k=Z0k,Wki=Wki,minNk=minNk,NegBin=NegBin,ZIP=ZIP,slopes=slopes,slope_col=slope_col,condSize=condSize,condOut=condOut,joint=joint,w=w,z=z)),silent=T)
    if(NegBin==TRUE){
      est[length(est)] <- exp(est[length(est)])
      if (length(invI)>1 & length(weights)<length(Yki)){ ## ML version
        SE <- c(sqrt(diag(invI)[1:(min(sig_id)-1)]),sqrt(diag(invI)[sig_id]*((exp(opt[sig_id]))^2) ),sqrt(diag(invI)[length(opt)]*((exp(opt[length(opt)]))^2) ))
      } else if(length(invI)>1){ ## weighted method requires sandwich form
        cheese <- cheese_JMMICS(params=opt,Nk=Nk,Nk_actual=Nk_actual,Nk0=Nk0,Zk=Zk,Yki=Yki,Xki=Xki,IDk=IDk,Z0k=Z0k,Wki=Wki,minNk=minNk,NegBin=NegBin,ZIP=ZIP,slopes=slopes,slope_col=slope_col,condSize=condSize,condOut=condOut,joint=joint,w=w,z=z)
        preDelta_vcov <- diag(invI%*%cheese%*%invI) # just to save us from matrix multiplying twice
        SE <- c(sqrt(preDelta_vcov[1:(min(sig_id)-1)]),sqrt(preDelta_vcov[sig_id]*((exp(opt[sig_id]))^2) ),sqrt(preDelta_vcov[length(opt)]*((exp(opt[length(opt)]))^2) ))
      } else { ## in case the matrix is numerically singular
        SE <- c(rep(NA,length(est)))
      }
    }else if(ZIP==TRUE){
      if (length(invI)>1 & length(weights)<length(Yki)){ ## ML version
        SE <- c(sqrt(diag(invI)[1:(min(sig_id)-1)]),sqrt(diag(invI)[sig_id]*((exp(opt[sig_id]))^2) ),sqrt(diag(invI)[(max(sig_id)+1):(length(opt))]))
      } else if(length(invI)>1){ ## weighted method requires sandwich form
        cheese <- cheese_JMMICS(params=opt,Nk=Nk,Nk_actual=Nk_actual,Nk0=Nk0,Zk=Zk,Yki=Yki,Xki=Xki,IDk=IDk,Z0k=Z0k,Wki=Wki,minNk=minNk,NegBin=NegBin,ZIP=ZIP,slopes=slopes,slope_col=slope_col,condSize=condSize,condOut=condOut,joint=joint,w=w,z=z)
        preDelta_vcov <- diag(invI%*%cheese%*%invI) # just to save us from matrix multiplying twice
        SE <- c(sqrt(preDelta_vcov[1:(min(sig_id)-1)]),sqrt(preDelta_vcov[sig_id]*((exp(opt[sig_id]))^2) ),sqrt(preDelta_vcov[(max(sig_id)+1):(length(opt))]))
              #c(sqrt(preDelta_vcov[1:(length(opt)-2)]),sqrt(preDelta_vcov[length(opt)-1]*((exp(opt[length(opt)-1]))^2) ),sqrt(preDelta_vcov[length(opt)]*((exp(opt[length(opt)]))^2) ))
      } else { ## in case the matrix is numerically singular
        SE <- c(rep(NA,length(est)))
      }
    }else{
      if (length(invI)>1 & length(weights)<length(Yki)){ ## ML version (for joint model and for outcome-only)
        SE <- c(sqrt(diag(invI)[1:(min(sig_id)-1)]),sqrt(diag(invI)[sig_id]*((exp(opt[sig_id]))^2) ))
      } else if(length(invI)>1){ ## weighted method requires sandwich form
        cheese <- cheese_JMMICS(params=opt,Nk=Nk,Nk_actual=Nk_actual,Nk0=Nk0,Zk=Zk,Yki=Yki,Xki=Xki,IDk=IDk,Z0k=Z0k,Wki=Wki,minNk=minNk,NegBin=NegBin,ZIP=ZIP,slopes=slopes,slope_col=slope_col,condSize=condSize,condOut=condOut,joint=joint,w=w,z=z)
        preDelta_vcov <- diag(invI%*%cheese%*%invI) # just to save us from matrix multiplying twice
        SE <- c(sqrt(preDelta_vcov[1:(min(sig_id)-1)]),sqrt(preDelta_vcov[sig_id]*((exp(opt[sig_id]))^2) ))
      } else { ## in case the matrix is numerically singular
        SE <- c(rep(NA,length(est)))
      }
    }

    ranef <- EB_JMMICS(opt,Nk=Nk,Nk_actual=Nk_actual,Nk0=Nk0,Zk=Zk,Yki=Yki,Xki=Xki,IDk=IDk,Z0k=Z0k,Wki=Wki,minNk=minNk,NegBin=NegBin,ZIP=ZIP,slopes=slopes,slope_col=slope_col,condSize=condSize,condOut=condOut,joint=joint,w=w,z=z)



  ## collect results
  results <- list(ests=est,SEs=SE,bks=ranef)
  ## print names of parameters
  names(results$ests) <- c(resnames)
  names(results$SEs) <- c(resnames)
  return(results)
}



