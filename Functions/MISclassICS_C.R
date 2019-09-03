############################################################
## Wrapper function for joint marginalized models for ICS ##
############################################################
## Permits Poisson, Negative-Binomial, or ZIP size models (does not combine zero inflation with negative binomial)
## notation differs slightly from paper: Nk0 is 1(Nk==0), Z0k is Wk (ie covariates for zero model), and theta is Wk epsilon

## May 27, 2018 EDIT: ## random slopes (for binary, cluster-level covariate) now permitted
## sigma and gamma are translated into vectors of length K (to acommodate random slopes)
runLOCAL=TRUE
## random slopes parametrization: sigma0 is the sd of intercepts among unexposed, sigma1 among exposed

require(gaussquad)
require(Rcpp)
if(runLOCAL==TRUE){
  install.packages("../Functions/MISclassICSpack_0.1.0.tar.gz",repos=NULL)
}
library(MISclassICSpack)


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

extract_params_MIS <- function(params,Nk,Nk_actual,Nk0,Zk,Yki,Xki,IDk,Z0k,Wk,Zwk,Xk,Zxk,minNk=0,NegBin=FALSE,ZIP=FALSE,slopes=FALSE,slope_col=2,condSize=FALSE,condOut=FALSE,joint=TRUE,interac=FALSE,inter_id=4,main_id=3,w,z){

  ## NegBin/ZIP not implemented yet
  theta=NA
  tau=NA
  epsilon=NA
  Xi=NA

  ## plug in 0 and 1
  Xki0 <- Xki1 <- Xki ## outcome model
  Xki0[,2] <- 0; Xki1[,2] <- 1
  Zk0 <- Zk1 <- Zk ## size model
  Zk0[,2] <- 0; Zk1[,2] <- 1
  Zwk0 <- Zwk1 <- Zwk ## misclassification model
  Zwk0[,2] <- 0; Zwk1[,2] <- 1
  if(interac==TRUE){ ## if interaction
    Zwk0[,inter_id] <- 0 ## set interaction to 0*maineffect
    Zwk1[,inter_id] <- Zwk1[,main_id] ## set interaction to 1*maineffect
  }

  if(joint==TRUE){
    if(slopes==FALSE){ ## random intercepts

      alpha <- params[1:ncol(Zk)]
      gamma <- params[ncol(Zk)+1]
      gamma0 <- gamma1 <- gamma ## for later use
      beta <- params[(ncol(Zk)+2):(ncol(Zk)+ncol(Xki)+1)]
      sigma <- exp(params[(ncol(Zk)+ncol(Xki)+2)]) # maximize wrt logsigma
      sigma0 <- sigma1 <- sigma ## for later use
      # if(NegBin==TRUE){
      #   tau <- exp(params[length(params)]) # maximize wrt logtau
      # } else{
      #   tau <- NA ## not used
      # }
      # if(ZIP==TRUE){
      #   epsilon <- params[(ncol(Zk)+ncol(Xki)+3):(length(params))]
      # }else{
      #   epsilon <- NA ## not used
      # }
      ## else
      epsW <- params[(ncol(Zk)+ncol(Xki)+3):(ncol(Zk)+ncol(Xki)+2+ncol(Zwk))]
      epsX <- params[(ncol(Zk)+ncol(Xki)+2+ncol(Zwk)+1):(ncol(Zk)+ncol(Xki)+2+ncol(Zwk)+ncol(Zxk))]

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
      # if(NegBin==TRUE){
      #   tau <- exp(params[length(params)]) # maximize wrt logtau
      # } else{
      #   tau <- NA ## not used
      # }
      # if(ZIP==TRUE){
      #   epsilon <- params[(ncol(Zk)+ncol(Xki)+5):(length(params))]
      # }else{
      #   epsilon <- NA ## not used
      # }
      ## else
      epsW <- params[(ncol(Zk)+ncol(Xki)+4+1):(ncol(Zk)+ncol(Xki)+4+ncol(Zwk))]
      epsX <- params[(ncol(Zk)+ncol(Xki)+4+ncol(Zwk)+1):(ncol(Zk)+ncol(Xki)+4+ncol(Zwk)+ncol(Zxk))]

      ## random slopes indicator (ie binary covariate for random slope)
      rand_indic <- Zk[,slope_col]

      ## convert gamma and sigma to vectors (for random slopes)
      gamma_vec <- rep(gamma0,length(Nk))
      gamma_vec[rand_indic==1] <- gamma1
      sigma_vec <- rep(sigma0,length(Nk))
      sigma_vec[rand_indic==1] <- sigma1
    }



    ## compute eta, zeta, possibly theta
    zeta <- Zk%*%alpha; zeta0 <- Zk0%*%alpha; zeta1 <- Zk1%*%alpha
    eta <- Xki%*%beta;  eta0 <- Xki0%*%beta;  eta1 <- Xki1%*%beta
    etaW <- Zwk%*%epsW; etaW0 <- Zwk0%*%epsW; etaW1 <- Zwk1%*%epsW
    etaX <- Zxk%*%epsX ## only one version
    # if(ZIP==TRUE){
    #   theta <- Z0k%*%epsilon
    # }else{theta=NA}

    ## get implicitly-defined conditional params
    Omega <- get_Omega(gamma_vec,sigma_vec,zeta,condSize,w,z)
    Omega0 <- get_Omega(rep(gamma0,length(Nk)),rep(sigma0,length(Nk)),zeta0,condSize,w,z)
    Omega1 <- get_Omega(rep(gamma1,length(Nk)),rep(sigma1,length(Nk)),zeta1,condSize,w,z)
    #Delta <- get_Delta(sigma_vec,eta,condOut,w,z)
    Delta <- get_Delta(rep(sigma_vec,Nk_actual),eta,condOut,w,z) ## new version needs separate sigma for each cluster, but Delta is observation-wise, so we need a separate sigma for each obs
    Delta0 <- get_Delta(rep(rep(sigma0,length(Nk)),Nk_actual),eta0,condOut,w,z)
    Delta1 <- get_Delta(rep(rep(sigma1,length(Nk)),Nk_actual),eta1,condOut,w,z)
    # if(ZIP==TRUE){
    #   Xi <- theta #get_Xi(nu,sigma_vec,theta,condZero,w,z)
    # }else{
    #   Xi=NA
    # }
    sigma_vec0 <- rep(sigma0,length(sigma_vec)); sigma_vec1 <- rep(sigma1,length(sigma_vec))
    gamma_vec0 <- rep(gamma0,length(gamma_vec)); gamma_vec1 <- rep(gamma1,length(gamma_vec))
    return(list(alpha=alpha,beta=beta,epsW=epsW,epsX=epsX,
                gamma=gamma_vec,gamma0=gamma_vec0,gamma1=gamma_vec1,
                sigma=sigma_vec,sigma0=sigma_vec0,sigma1=sigma_vec1,
                tau=tau,epsilon=epsilon,theta=theta,Xi=Xi,## not implemented
                zeta=zeta,eta=eta,Omega=Omega,Delta=Delta,etaW=etaW,etaX=etaX,
                zeta0=zeta0,eta0=eta0,Omega0=Omega0,Delta0=Delta0,etaW0=etaW0,
                zeta1=zeta1,eta1=eta1,Omega1=Omega1,Delta1=Delta1,etaW1=etaW1,
                Xki0=Xki0,Xki1=Xki1,Zk0=Zk0,Zk1=Zk1,Zwk0=Zwk0,Zwk1=Zwk1,
                rand_indic=rand_indic))

  }else{ ## outcome-only model
    if(slopes==FALSE){ ## random intercepts

      beta <- params[1:ncol(Xki)]
      sigma <- exp(params[(ncol(Xki)+1)]) ## maximize wrt logsigma
      sigma0 <- sigma1 <- sigma ## for later use
      epsW <- params[(ncol(Xki)+1+1):(ncol(Xki)+1+ncol(Zwk))]
      epsX <- params[(ncol(Xki)+1+ncol(Zwk)+1):(ncol(Xki)+1+ncol(Zwk)+ncol(Zxk))]

      sigma_vec <- rep(sigma,length(Nk)) ## convert sigma to vector
      rand_indic <- NA ## no random slope

    }else{ ## random slopes

      beta <- params[1:ncol(Xki)]
      sigma0 <- exp(params[(ncol(Xki)+1)]) ## maximize wrt logsigma
      sigma1 <- exp(params[(ncol(Xki)+2)]) ## maximize wrt logsigma
      epsW <- params[(ncol(Xki)+2+1):(ncol(Xki)+2+ncol(Zwk))]
      epsX <- params[(ncol(Xki)+2+ncol(Zwk)+1):(ncol(Xki)+2+ncol(Zwk)+ncol(Zxk))]
      rand_indic <- Zk[,slope_col]       ## random slopes indicator (ie binary covariate for random slope)
      sigma_vec <- rep(sigma0,length(Nk))
      sigma_vec[rand_indic==1] <- sigma1 ## convert sigma to vector (for random slopes)

    }

    ## compute eta
    eta <- Xki%*%beta; eta0 <- Xki0%*%beta; eta1 <- Xki1%*%beta
    etaW <- Zwk%*%epsW; etaW0 <- Zwk0%*%epsW; etaW1 <- Zwk1%*%epsW
    etaX <- Zxk%*%epsX ## only one version
    ## get implicitly-defined conditional params
    Delta <- get_Delta(rep(sigma_vec,Nk_actual),eta,condOut,w,z) ## new version needs separate sigma for each cluster, but Delta is observation-wise, so we need a separate sigma for each obs
    Delta0 <- get_Delta(rep(rep(sigma0,length(Nk)),Nk_actual),eta0,condOut,w,z)
    Delta1 <- get_Delta(rep(rep(sigma1,length(Nk)),Nk_actual),eta1,condOut,w,z)

    sigma_vec0 <- rep(sigma0,length(sigma_vec)); sigma_vec1 <- rep(sigma1,length(sigma_vec))
    return(list(beta=beta,epsW=epsW,epsX=epsX,
                sigma=sigma_vec,sigma0=sigma_vec0,sigma1=sigma_vec1,
                eta=eta,Delta=Delta,etaW=etaW,etaX=etaX,
                eta0=eta0,Delta0=Delta0,etaW0=etaW0,
                eta1=eta1,Delta1=Delta1,etaW1=etaW1,
                Xki0=Xki0,Xki1=Xki1,Zk0=Zk0,Zk1=Zk1,Zwk0=Zwk0,Zwk1=Zwk1,
                rand_indic=rand_indic))

  }

}








### function to compute log-likelihood
logLik_MISclassICS <- function(params,Nk,Nk_actual,Nk0,Zk,Yki,Xki,IDk,Z0k,Wk,Zwk,Xk,Zxk,VALIDk,minNk=0,NegBin=FALSE,ZIP=FALSE,slopes=FALSE,slope_col=2,condSize=FALSE,condOut=FALSE,joint=TRUE,interac=FALSE,inter_id=4,main_id=3,w,z){

  comp <- extract_params_MIS(params=params,Nk=Nk,Nk_actual=Nk_actual,Nk0=Nk0,Zk=Zk,Yki=Yki,Xki=Xki,IDk=IDk,Z0k=Z0k,Wk=Wk,Zwk=Zwk,Xk=Xk,Zxk=Zxk,minNk=minNk,NegBin=NegBin,ZIP=ZIP,slopes=slopes,slope_col=slope_col,condSize=condSize,condOut=condOut,joint=joint,interac=interac,inter_id=inter_id,main_id=main_id,w=w,z=z)

  if(joint==TRUE){ ## joint model

    if(NegBin==TRUE){
      # return(-compute_logLik_JMMICS_NegBin(comp$gamma,comp$sigma,comp$tau,comp$Omega,comp$Delta,Nk,Nk_actual,Yki,IDk,minNk,w,z))
    }else if(ZIP==TRUE){
      # return(-compute_logLik_JMMICS_ZIP(comp$gamma,comp$sigma,comp$Xi,comp$Omega,comp$Delta,Nk,Nk_actual,Nk0,Yki,IDk,minNk,w,z))
    }else{
      return(-compute_logLik_MISclassICS(comp$gamma,comp$gamma0,comp$gamma1,comp$sigma,comp$sigma0,comp$sigma1,comp$Omega,comp$Omega0,comp$Omega1,comp$Delta,comp$Delta0,comp$Delta1,comp$etaW,comp$etaW0,comp$etaW1,comp$etaX,Nk,Nk_actual,Yki,Wk,Xk,IDk,c(VALIDk==1),minNk,w,z))
    }
  }else{ ## outcome-only

    return(-compute_logLik_MISclassGLMM(comp$sigma,comp$sigma0,comp$sigma1,comp$Delta,comp$Delta0,comp$Delta1,comp$etaW,comp$etaW0,comp$etaW1,comp$etaX,Nk_actual,Yki,Wk,Xk,IDk,c(VALIDk==1),w,z))
  }
}


### function to compute score
score_MISclassICS <- function(params,Nk,Nk_actual,Nk0,Zk,Yki,Xki,IDk,Z0k,Wk,Zwk,Xk,Zxk,VALIDk,minNk=0,NegBin=FALSE,ZIP=FALSE,slopes=FALSE,slope_col=2,condSize=FALSE,condOut=FALSE,joint=TRUE,interac=FALSE,inter_id=4,main_id=3,w,z){

  comp <- extract_params_MIS(params=params,Nk=Nk,Nk_actual=Nk_actual,Nk0=Nk0,Zk=Zk,Yki=Yki,Xki=Xki,IDk=IDk,Z0k=Z0k,Wk=Wk,Zwk=Zwk,Xk=Xk,Zxk=Zxk,minNk=minNk,NegBin=NegBin,ZIP=ZIP,slopes=slopes,slope_col=slope_col,condSize=condSize,condOut=condOut,joint=joint,interac=interac,inter_id=inter_id,main_id=main_id,w=w,z=z)

  if(joint==TRUE){ ## joint model

    if(slopes==FALSE){ ## random intercepts
      if(NegBin==TRUE){
        # return(-compute_score_JMMICS_NegBin(comp$gamma,comp$sigma,tau,DG_tau=digamma(comp$tau),comp$Omega,comp$Delta,comp$zeta,comp$eta,Nk,Nk_actual,Zk,Yki,Xki,IDk,minNk,condSize,condOut,w,z,length(comp$alpha),length(comp$beta)))
      }else if(ZIP==TRUE){
        # return(-compute_score_JMMICS_ZIP(comp$gamma,comp$sigma,comp$Xi,comp$Omega,comp$Delta,comp$theta,comp$zeta,comp$eta,Nk,Nk_actual,Nk0,Z0k,Zk,Yki,Xki,IDk,minNk,condSize,condOut,w,z,length(comp$epsilon),length(comp$alpha),length(comp$beta)))
      }else{
        return(-compute_score_MISclassICS(comp$gamma,comp$gamma0,comp$gamma1,comp$sigma,comp$sigma0,comp$sigma1,comp$Omega,comp$Omega0,comp$Omega1,comp$Delta,comp$Delta0,comp$Delta1,comp$zeta,comp$zeta0,comp$zeta1,comp$eta,comp$eta0,comp$eta1,comp$etaW,comp$etaW0,comp$etaW1,comp$etaX,Nk,Nk_actual,Zk,comp$Zk0,comp$Zk1,Yki,Xki,comp$Xki0,comp$Xki1,Wk,Zwk,comp$Zwk0,comp$Zwk1,Xk,Zxk,IDk,c(VALIDk==1),minNk,condSize,condOut,w,z,length(comp$alpha),length(comp$beta),length(comp$epsW),length(comp$epsX)))
      }
    }else{ ## random slopes
      if(NegBin==TRUE){
        # return(-compute_score_JMMICS_NegBin_slopes(comp$gamma,comp$sigma,tau,DG_tau=digamma(comp$tau),comp$Omega,comp$Delta,comp$zeta,comp$eta,comp$rand_indic,Nk,Nk_actual,Zk,Yki,Xki,IDk,minNk,condSize,condOut,w,z,length(comp$alpha),length(comp$beta)))
      }else if(ZIP==TRUE){
        # return(-compute_score_JMMICS_ZIP_slopes(comp$gamma,comp$sigma,comp$Xi,comp$Omega,comp$Delta,comp$theta,comp$zeta,comp$eta,comp$rand_indic,Nk,Nk_actual,Nk0,Z0k,Zk,Yki,Xki,IDk,minNk,condSize,condOut,w,z,length(comp$epsilon),length(comp$alpha),length(comp$beta)))
      }else{
        return(-compute_score_MISclassICS_slopes(comp$gamma,comp$gamma0,comp$gamma1,comp$sigma,comp$sigma0,comp$sigma1,comp$Omega,comp$Omega0,comp$Omega1,comp$Delta,comp$Delta0,comp$Delta1,comp$zeta,comp$zeta0,comp$zeta1,comp$eta,comp$eta0,comp$eta1,comp$etaW,comp$etaW0,comp$etaW1,comp$etaX,comp$rand_indic,Nk,Nk_actual,Zk,comp$Zk0,comp$Zk1,Yki,Xki,comp$Xki0,comp$Xki1,Wk,Zwk,comp$Zwk0,comp$Zwk1,Xk,Zxk,IDk,c(VALIDk==1),minNk,condSize,condOut,w,z,length(comp$alpha),length(comp$beta),length(comp$epsW),length(comp$epsX)))
      }
    }
  }else{ ## outcome only model

    if(slopes==FALSE){ ## random intercepts
      return(-compute_score_MISclassGLMM(comp$sigma,comp$sigma0,comp$sigma1,comp$Delta,comp$Delta0,comp$Delta1,comp$eta,comp$eta0,comp$eta1,comp$etaW,comp$etaW0,comp$etaW1,comp$etaX,Nk_actual,Yki,Xki,comp$Xki0,comp$Xki1,Wk,Zwk,comp$Zwk0,comp$Zwk1,Xk,Zxk,IDk,c(VALIDk==1),condOut,w,z,length(comp$beta),length(comp$epsW),length(comp$epsX)))
    }else{ ## random slopes
      return(-compute_score_MISclassGLMM_slopes(comp$sigma,comp$sigma0,comp$sigma1,comp$Delta,comp$Delta0,comp$Delta1,comp$eta,comp$eta0,comp$eta1,comp$etaW,comp$etaW0,comp$etaW1,comp$etaX,comp$rand_indic,Nk_actual,Yki,Xki,comp$Xki0,comp$Xki1,Wk,Zwk,comp$Zwk0,comp$Zwk1,Xk,Zxk,IDk,c(VALIDk==1),condOut,w,z,length(comp$beta),length(comp$epsW),length(comp$epsX)))
    }
  }

}




## MISclassICS wrapper function
## slopes==TRUE is for random slopes
## slope_col indicates the column of Zk that refers to the binary covariate for random slopes
## include true exposure where available, and NA when not available
MISclassICS_fit <- function(Nk,Zk,Yki,Xki,Wk,Zwk,Xk,Zxk,VALIDk,IDk,Z0k=NULL,minNk=0,NegBin=FALSE,ZIP=FALSE,slopes=FALSE,slope_col=2,condSize=FALSE,condOut=FALSE,joint=TRUE,interac=FALSE,inter_id=4,main_id=3,startpos=FALSE,nquad=25){
  ## turn GH weights/zeros into appropriately scaled terms
  w <- as.matrix(hermite.h.quadrature.rules(nquad)[[nquad]])[,2]; w <- w/sqrt(pi)
  z <- as.matrix(hermite.h.quadrature.rules(nquad)[[nquad]])[,1]; z <- sqrt(2)*z

  ## handle NAs
  Zk[VALIDk==0,2] <- 0
  Xki[rep(VALIDk,Nk)==0,2] <- 0
  Zwk[VALIDk==0,2] <- 0
  if(interac==TRUE){Zwk[VALIDk==0,inter_id] <- 0}


  ## handle empty clusters
  Nk0 <- as.numeric(Nk==0) ## indicator of emptiness
  Nk_actual <- Nk

  ## starting value for gamma
  if(startpos==TRUE){
    startgam <- 1
  }else{startgam <- -1 }

  if(joint==TRUE){ ## usual joint model
    ## alpha +gamma +beta +sigma
    ##  if negbin --> +possible tau
    ##  if ZIP -->  +possible epsilon +possible nu
    if(slopes==FALSE){ ## random intercepts
      startval <- c(rep(1,ncol(Zk)),startgam,rep(0,ncol(Xki)),log(0.5),rep(1,ncol(Zwk)),rep(1,ncol(Zxk)))
      resnames <- c(paste0(rep("alpha",ncol(Zk)),seq(0,ncol(Zk)-1)),"gamma",paste0(rep("beta",ncol(Xki)),seq(0,ncol(Xki)-1)),"sigma",paste0(rep("epsW",ncol(Zwk)),seq(0,ncol(Zwk)-1)),paste0(rep("epsX",ncol(Zxk)),seq(0,ncol(Zxk)-1))) ## names of parameters to print
      sig_id <- ncol(Zk)+1+ncol(Xki)+1 ## index for logsigma
    }else{ ## random slopes (for binary covariate)
      startval <- c(rep(1,ncol(Zk)),startgam,startgam,rep(0,ncol(Xki)),log(0.5),log(0.5),rep(1,ncol(Zwk)),rep(1,ncol(Zxk)))
      resnames <- c(paste0(rep("alpha",ncol(Zk)),seq(0,ncol(Zk)-1)),"gamma0","gamma1",paste0(rep("beta",ncol(Xki)),seq(0,ncol(Xki)-1)),"sigma0","sigma1",paste0(rep("epsW",ncol(Zwk)),seq(0,ncol(Zwk)-1)),paste0(rep("epsX",ncol(Zxk)),seq(0,ncol(Zxk)-1))) ## names of parameters to print
      sig_id <- ncol(Zk)+2+ncol(Xki)+1:2  ## index for logsigma0 and logsigma1
    }
    if(NegBin==TRUE){
      # startval <- c(startval,0)
      # resnames <- c(resnames,"tau")
    }
    if(ZIP==TRUE){
      # startval <- c(startval,rep(0,ncol(Z0k)))
      # resnames <- c(resnames,paste0(rep("epsilon",ncol(Z0k)),seq(0,ncol(Z0k)-1)))
    }
  }else{ ## outcome-only model
    NegBin <- FALSE
    ZIP <- FALSE
    if(slopes==FALSE){ ## random intercepts
      startval <- c(rep(0,ncol(Xki)),log(0.5),rep(1,ncol(Zwk)),rep(1,ncol(Zxk)))
      resnames <- c(paste0(rep("beta",ncol(Xki)),seq(0,ncol(Xki)-1)),"sigma",paste0(rep("epsW",ncol(Zwk)),seq(0,ncol(Zwk)-1)),paste0(rep("epsX",ncol(Zxk)),seq(0,ncol(Zxk)-1))) ## names of parameters to print
      sig_id <- ncol(Xki)+1 ## index for logsigma
    }else{ ## random slopes (for binary covariate)
      startval <- c(rep(0,ncol(Xki)),log(0.5),log(0.5),rep(1,ncol(Zwk)),rep(1,ncol(Zxk)))
      resnames <- c(paste0(rep("beta",ncol(Xki)),seq(0,ncol(Xki)-1)),"sigma0","sigma1",paste0(rep("epsW",ncol(Zwk)),seq(0,ncol(Zwk)-1)),paste0(rep("epsX",ncol(Zxk)),seq(0,ncol(Zxk)-1))) ## names of parameters to print
      sig_id <- ncol(Xki)+1:2  ## index for logsigma0 and logsigma1
    }
  }

  ## optimization
  opt <- nlminb(start=startval,logLik_MISclassICS,gradient=score_MISclassICS,Nk=Nk,Nk_actual=Nk_actual,Nk0=Nk0,Zk=Zk,Yki=Yki,Xki=Xki,IDk=IDk,Z0k=Z0k,Wk=Wk,Zwk=Zwk,Xk=Xk,Zxk=Zxk,VALIDk=VALIDk,minNk=minNk,NegBin=NegBin,ZIP=ZIP,slopes=slopes,slope_col=slope_col,condSize=condSize,condOut=condOut,joint=joint,interac=interac,inter_id=inter_id,main_id=main_id,w=w,z=z,control=list(eval.max=2000,iter.max=1500))$par
  est <- opt; ## ML estimates
  est[sig_id] <- exp(est[sig_id]) ## transform logsigma to sigma


  ## Standard Errors (with Delta Method)
  invI <- try(solve(optimHess(par=opt,logLik_MISclassICS,gr=score_MISclassICS,Nk=Nk,Nk_actual=Nk_actual,Nk0=Nk0,Zk=Zk,Yki=Yki,Xki=Xki,IDk=IDk,Z0k=Z0k,Wk=Wk,Zwk=Zwk,Xk=Xk,Zxk=Zxk,VALIDk=VALIDk,minNk=minNk,NegBin=NegBin,ZIP=ZIP,slopes=slopes,slope_col=slope_col,condSize=condSize,condOut=condOut,joint=joint,interac=interac,inter_id=inter_id,main_id=main_id,w=w,z=z)),silent=T)
  if(NegBin==TRUE){
    # est[length(est)] <- exp(est[length(est)])
    # if (length(invI)>1 & length(weights)<length(Yki)){ ## ML version
    #   SE <- c(sqrt(diag(invI)[1:(min(sig_id)-1)]),sqrt(diag(invI)[sig_id]*((exp(opt[sig_id]))^2) ),sqrt(diag(invI)[length(opt)]*((exp(opt[length(opt)]))^2) ))
    # } else if(length(invI)>1){ ## weighted method requires sandwich form
    #   cheese <- cheese_JMMICS(params=opt,Nk=Nk,Nk_actual=Nk_actual,Nk0=Nk0,Zk=Zk,Yki=Yki,Xki=Xki,IDk=IDk,Z0k=Z0k,minNk=minNk,NegBin=NegBin,ZIP=ZIP,slopes=slopes,slope_col=slope_col,condSize=condSize,condOut=condOut,joint=joint,w=w,z=z)
    #   preDelta_vcov <- diag(invI%*%cheese%*%invI) # just to save us from matrix multiplying twice
    #   SE <- c(sqrt(preDelta_vcov[1:(min(sig_id)-1)]),sqrt(preDelta_vcov[sig_id]*((exp(opt[sig_id]))^2) ),sqrt(preDelta_vcov[length(opt)]*((exp(opt[length(opt)]))^2) ))
    # } else { ## in case the matrix is numerically singular
    #   SE <- c(rep(NA,length(est)))
    # }
  }else if(ZIP==TRUE){
    # if (length(invI)>1 & length(weights)<length(Yki)){ ## ML version
    #   SE <- c(sqrt(diag(invI)[1:(min(sig_id)-1)]),sqrt(diag(invI)[sig_id]*((exp(opt[sig_id]))^2) ),sqrt(diag(invI)[(max(sig_id)+1):(length(opt))]))
    # } else if(length(invI)>1){ ## weighted method requires sandwich form
    #   cheese <- cheese_JMMICS(params=opt,Nk=Nk,Nk_actual=Nk_actual,Nk0=Nk0,Zk=Zk,Yki=Yki,Xki=Xki,IDk=IDk,Z0k=Z0k,minNk=minNk,NegBin=NegBin,ZIP=ZIP,slopes=slopes,slope_col=slope_col,condSize=condSize,condOut=condOut,joint=joint,w=w,z=z)
    #   preDelta_vcov <- diag(invI%*%cheese%*%invI) # just to save us from matrix multiplying twice
    #   SE <- c(sqrt(preDelta_vcov[1:(min(sig_id)-1)]),sqrt(preDelta_vcov[sig_id]*((exp(opt[sig_id]))^2) ),sqrt(preDelta_vcov[(max(sig_id)+1):(length(opt))]))
    #   #c(sqrt(preDelta_vcov[1:(length(opt)-2)]),sqrt(preDelta_vcov[length(opt)-1]*((exp(opt[length(opt)-1]))^2) ),sqrt(preDelta_vcov[length(opt)]*((exp(opt[length(opt)]))^2) ))
    # } else { ## in case the matrix is numerically singular
    #   SE <- c(rep(NA,length(est)))
    # }
  }else{
    if (length(invI)>1 & sum(diag(invI)<0)==0 ){ ## ML version (for joint model and for outcome-only)
      SE <- c(sqrt(diag(invI)[1:(min(sig_id)-1)]),sqrt(diag(invI)[sig_id]*((exp(opt[sig_id]))^2)),sqrt(diag(invI)[(max(sig_id)+1):length(opt)]))
    } else { ## in case the matrix is numerically singular
      SE <- c(rep(NA,length(est)))
    }
  }



  ## collect results
  results <- list(ests=est,SEs=SE)
  ## print names of parameters
  names(results$ests) <- c(resnames)
  names(results$SEs) <- c(resnames)
  return(results)
}






#
# ### solve naive expected estimating equations
# fit_naiveEEE <- function(start_param,Zki,Zxki,Nk_actual,Yki,Xk,IDk,WEIGHTk=NA,VALIDk){
#
#
#   return(multiroot(naiveEEE,start=start_param,Zki=Zki,Zxki=Zxki,Nk_actual=Nk_actual,Yki=Yki,Xk=Xk,IDk=IDk,WEIGHTk=WEIGHTk,VALIDk=VALIDk))
#
#
# }
#
# ### function to compute naive expected estimating equations
# naiveEEE <- function(params,Zki,Zxki,Nk_actual,Yki,Xk,IDk,WEIGHTk=NA,VALIDk){
#
#   ## make unit level
#   Xki <- rep(Xk,Nk_actual)
#   VALIDki <- rep(VALIDk,Nk_actual)
#   if(length(WEIGHTk)!=length(VALIDk)){
#     WEIGHTki <- rep(1,length(VALIDki))
#   }else{
#     WEIGHTki <- rep(WEIGHTk,Nk_actual)
#   }
#
#
#   ### collect parameters
#   beta <- params[1:(ncol(Zki)+1)] ## outcome model parameters
#   betaX <- beta[2] ## X parameter
#   betaZ <- beta[-2] ## Z parameters
#   alphaX <- params[(ncol(Zki)+2):length(params)] ## exposure model parameters (second column of Zki must be W)
#
#   ## compute linear predictors for outcome model
#   muki <- Xki*betaX+Zki%*%betaZ
#   mu0ki <- 0*betaX+Zki%*%betaZ
#   mu1ki <- 1*betaX+Zki%*%betaZ
#   ## linear predictors for exposure model
#   muXki <- Zxki%*%alphaX
#
#   ## compute EEE
#   # return(compute_naiveEEE(muki,mu0ki,mu1ki,muXki,Xki,Nk_actual,Yki,IDk,WEIGHTki,VALIDki))
#   return(compute_naiveEEE(muki,mu0ki,mu1ki,muXki,Zki,Zxki,length(beta),length(alphaX),Xki,Nk_actual,Yki,IDk,WEIGHTki,VALIDki))
# }

