##############################
## Simulate Data with ICS   ##
##############################
## Updated: 11/21/2018
## Code generates dataset with ICS
## Reports 3 objects:
  ## data --> obs-level dataframe
  ## data_clstlvl --> cluster-level dataframe
  ## true_randeff --> true random intercepts used to generate data
    ## Currently reports K true random effects (not excluding the empty clusters)
    ## be wary if you test the MSEP, need to exclude the ones with Nk=0 (but we also report Nk=0 so should be easy)


## load packages
require(glmmML)
require(geepack)
require(MMLB)
require(gaussquad)
source("../Functions/JMMICS.R")

## set functions
expit <- function(y){  return(exp(y)/(1+exp(y)))}
logit <- function(y){  return(log(y/(1-y)))}


genICS <- function(K=1000, ## No. of clusters
                   Xrate=0.25, ## cluster-level exposure rate
                   alpha0=0.6,alpha1=-0.2, ## volume model coefficients
                   beta0=-4,beta1=0.5,## outcome model coefficients
                   beta2=0.2, ## obs-level coefficient--> set to 0 for NO observation-level covariate
                   gamma_bk=-0.1, ## scaling factor for random intercept in volume model
                   sig.bk=2, ## standard deviation of random intercepts
                   min_Nk=0, ## minimum cluster size added to every (potentially 0) cluster size
                   cond=TRUE, ## TRUE--> specify parameters conditionally
                   ZIP=FALSE, ## TRUE--> ZIP model
                     eps0=-2,eps1=0.2, ## ZIP model coefficients
                     nu=0, ## DEPRECATED/REMOVE: for informative zeroes in ZIP (overparametrized model)
                   slopes=FALSE, ## TRUE--> random slope for binary covariate
                     sig1.bk=2, ## sd for rand eff when X1k=1
                     gamma1_bk=-0.1, ## scaling factor when X1k=1
                   nquad=50){

  ## turn GH weights/zeros into appropriately scaled terms
  w <- as.matrix(hermite.h.quadrature.rules(nquad)[[nquad]])[,2]; w <- w/sqrt(pi)
  z <- as.matrix(hermite.h.quadrature.rules(nquad)[[nquad]])[,1]; z <- sqrt(2)*z


  

  ###########################
  ## Simulate Generation 1 ##
  ###########################

  ## Generate Cluster id (Gen1)
  IDk <- seq(1,K)

  ## Generate True Exposures (Gen1)
  X1k <- rbinom(K,1,Xrate)
  
  ###  Random effects specification  ###
  if(slopes==FALSE){ ## random intercepts
    sig1.bk <- sig.bk ## set both SDs to be equal
    gamma1_bk <- gamma_bk ## set both gammas to be equal
  }
  sigma_vec <- sig.bk*(1-X1k)+sig1.bk*X1k ## defined at the cluster level
  gamma_vec <- gamma_bk*(1-X1k)+gamma1_bk*X1k ## defined at the cluster level
  
  
  ## SHARED RANDOM EFFECT
  bk <- rnorm(K,0,sigma_vec)   ## simplifies to bk <- rnorm(K,0,sig.bk) for rand int


  ## Generate Cluster sizes (Gen2)
  zeta <- alpha0+alpha1*X1k
  omega <- get_Omega(gamma_vec,sigma_vec,zeta,condSize=cond,w,z)
  Nk <- min_Nk+rpois(K,exp(omega+gamma_vec*bk))
  Nk[Nk>50] <- 50 ## failsafe

  ## zero-inflation:
  if(ZIP==TRUE){
    zips <- rbinom(K,1,expit(eps0+eps1*X1k+nu*bk)) ## nu=0
    Nk[zips==1] <- 0
  }

  ## inverse cluster size weights
  #  icw <- 1/Nk

  ###########################
  ## Simulate Generation 2 ##
  ###########################

  ## cluster id (Gen2)
  IDki <- rep(IDk,times=Nk)

  ## random intercepts
  bki <- rep(bk,times=Nk)

  ## true exposure (Gen2)
  X1ki <- rep(X1k,times=Nk)
  
  ## obs-level covariate
  X2ki <- rnorm(length(IDki),0,1)

  ## inverse cluster size weights
  ICWki <- 1/rep(Nk,times=Nk) ## obviously excludes 0s
  Nki <- rep(Nk,times=Nk)

  ## generate outcome
  eta <- beta0+beta1*X1ki+beta2*X2ki ## Conditional linear predictor
  delta <- get_Delta(rep(sigma_vec,times=Nk),eta,condOut=cond,w,z)
  Yki <- rbinom(sum(Nk),1,expit(delta+bki))



  ## outcome-level datagrame
  df <- data.frame(cbind(Yki,X1ki,X2ki,IDki,ICWki,Nki))

  ## cluster-level dataframe
  df_cluster <- data.frame(cbind(Nk,X1k,IDk))






  return(list(data=df,data_clstlvl=df_cluster,true_randeff=bk))
}



# ## Testing the distribution of cluster sizes
# ## median is 2, 25%ile is 1, 75th is 3, max is between 7-12
# test <- c()
# for(ii in 1:100){
#   dd <- genICS()$data_clstlvl
#   test <- rbind(test,quantile(dd$Nk,c(0,0.25,0.5,0.75,1)))
# }
# test
