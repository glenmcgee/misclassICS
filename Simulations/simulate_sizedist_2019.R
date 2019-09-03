######################################################
##      Simulate Distribution of Cluster Sizes      ##
######################################################
## Updated: 04/07/2019
## generate dataset with k=1million for each setting and report distribution of cluster sizes


runLOCAL=TRUE ## set to TRUE to run locally ## FALSE is on cluster
set.seed(1000)

#####################
## setting parameters
#####################
RR <- 2000 ## No. of iterations
KK <- 1000000 ## No. of clusters
Vpercent <- 0.20 ## Validation sampling fraction
# VV <- 1000  ## No. of clusters in (internal) validation sample
JJ <- 500 ## Number of jobs in array (only matters on cluster)
Xprev <- 0.25 ## Exposure prevalence
a0 <- 0.6; a1 <- -0.2 ## volume model coefficients
b0 <- -4; b1 <- 0.5; b2 <- 0.2 ## outcome model coefficients
#gamma <- -0.1 ## scaling factor for random intercept in volume model
sigma0 <- 2 ## SD of rand eff for X=0
sigma1 <- 1.5 ## SD of rand eff for X=1

###############################################
## set folder to save data and load functions
if(runLOCAL==TRUE){
  path <- "../Simulations/Results/" ## path for results
  funpath <- "../Functions/" ## path for functions
  outpath <- "../Simulations/Plots/"
} else{
  # path <- "../Simulations/Results/" ## path for results
  # funpath <- "../Functions/" ## path for functions
  path <- "~/Paper3/Code/Simulations/Results/" ## path for results
  funpath <- "~/Paper3/Code/Functions/"  ## path for functions
}


####################
## load packages
source(paste0(funpath,"JMMICS.R"))
source(paste0(funpath,"MISclassICS_C.R"))
source(paste0(funpath,"MISclassICS.R"))
source(paste0(funpath,"genICS.R"))
source(paste0(funpath,"genMIS.R"))
require(lme4)
require(geepack)
require(JMMICSpack)
require(xtable)
require(ggplot2)


##############################################################
## function to generate conditional and marginal distributions
get_sizeDist <- function(gam=0,ZIP=FALSE,slopes=TRUE){
  # cdat <- genICS(gamma_bk=gam,gamma1_bk=gam,min_Nk=1,slopes=slopes,cond=TRUE,ZIP=ZIP,   K=KK,Xrate=Xprev,alpha0=a0,alpha1=a1,beta0=b0,beta1=b1,beta2=b2,sig.bk=sigma0,sig1.bk=sigma1) ## conditional
  mdat <- genICS(gamma_bk=gam,gamma1_bk=gam,min_Nk=1,slopes=slopes,cond=FALSE,ZIP=ZIP,  K=KK,Xrate=Xprev,alpha0=a0,alpha1=a1,beta0=b0,beta1=b1,beta2=b2,sig.bk=sigma0,sig1.bk=sigma1) ## marginal

  # cdist <- table(cdat$data_clstlvl$Nk) ## get no. of clusters
  mdist <- table(mdat$data_clstlvl$Nk);mdist0 <- table(mdat$data_clstlvl$Nk[mdat$data_clstlvl$X1k==0]);mdist1 <- table(mdat$data_clstlvl$Nk[mdat$data_clstlvl$X1k==1])

  # cdist <- c(cdist[1:4],sum(cdist[5:length(cdist)])) ## collapse 5+
  mdist <- c(mdist[1:4],sum(mdist[5:length(mdist)]));mdist0 <- c(mdist0[1:4],sum(mdist0[5:length(mdist0)]));mdist1 <- c(mdist1[1:4],sum(mdist1[5:length(mdist1)]))

  # names(cdist) <-  <- c("1","2","3","4","5+")  ## rename
  names(mdist) <- c("1","2","3","4","5+");names(mdist0) <- paste("X=0,",c("1","2","3","4","5+")); names(mdist1) <- paste("X=1,",c("1","2","3","4","5+"))    ## rename


  # cdist <- round(100*cdist/KK) ## get percentages
  mdist <- round(100*mdist/KK);mdist0 <- round(100*mdist0/sum(mdat$data_clstlvl$X1k==0));mdist1 <- round(100*mdist1/sum(mdat$data_clstlvl$X1k==1))

  # dists <- data.frame(cbind(mdist,cdist)) ## combine

  ## report means
  meanmdist <- mean(mdat$data_clstlvl$Nk);meanmdist0 <- mean(mdat$data_clstlvl$Nk[mdat$data_clstlvl$X1k==0]);meanmdist1 <- mean(mdat$data_clstlvl$Nk[mdat$data_clstlvl$X1k==1])

  return(list(dist=mdist,dist0=mdist0,dist1=mdist1,
              mean=meanmdist,mean0=meanmdist0,mean1=meanmdist1))

}



##############################################################
## generate distributions for each simulation setting

dist_gam00 <- get_sizeDist(gam=0.00) ## for sim 1+2
dist_gam25 <- get_sizeDist(gam=-0.25) ## for sim 3+4
dist_gamP25 <- get_sizeDist(gam=+0.25)  ## for sim 5+6


##############################################################
## Make distribution tables

tab_Ndist <- rbind(dist_gam00$dist,dist_gam25$dist,dist_gamP25$dist)
tab_Ndist0 <- rbind(dist_gam00$dist0,dist_gam25$dist0,dist_gamP25$dist0)
tab_Ndist1 <- rbind(dist_gam00$dist1,dist_gam25$dist1,dist_gamP25$dist1)
rownames(tab_Ndist) <- c("NonICS","ICS Gam=-0.25","ICS Gam=0.25")
rownames(tab_Ndist0) <- c("NonICS","ICS Gam=-0.25","ICS Gam=0.25")
rownames(tab_Ndist1) <- c("NonICS","ICS Gam=-0.25","ICS Gam=0.25")

tabNdist <- cbind(tab_Ndist,tab_Ndist0,tab_Ndist1)
xtable(tabNdist,digits = 0)

## mean table
tab_means <- rbind(c(dist_gam00$mean,dist_gam00$mean0,dist_gam00$mean1),
                   c(dist_gam25$mean,dist_gam25$mean0,dist_gam25$mean1),
                   c(dist_gamP25$mean,dist_gamP25$mean0,dist_gamP25$mean1))
rownames(tab_means) <- c("NonICS","ICS Gam=-0.25","ICS Gam=0.25")
colnames(tab_means) <- c("Overall","X=0","X=1")
xtable(tab_means,digits = 2)





##############################################################################
## function to generate conditional and marginal distributions of outcome rate
get_outDist <- function(gam=0,m=1,ZIP=FALSE,slopes=TRUE,dig=1){
  # cdat <- genICS(gamma_bk=gam,gamma1_bk=gam,min_Nk=m,slopes=slopes,cond=TRUE,ZIP=ZIP,   K=KK,Xrate=Xprev,alpha0=a0,alpha1=a1,beta0=b0,beta1=b1,beta2=b2,eps0=e0,eps1=e1,sig.bk=sigma0,sig1.bk=sigma1) ## conditional
  mdat <- genICS(gamma_bk=gam,gamma1_bk=gam,min_Nk=m,slopes=slopes,cond=FALSE,ZIP=ZIP,  K=KK,Xrate=Xprev,alpha0=a0,alpha1=a1,beta0=b0,beta1=b1,beta2=b2,eps0=e0,eps1=e1,sig.bk=sigma0,sig1.bk=sigma1) ## marginal

  # cprev <- round(100*mean(cdat$data$Yki),digits=dig) ## get mean outcome
  mprev <- round(100*mean(mdat$data$Yki),digits=dig)
  # c0prev <- round(100*mean(cdat$data$Yki[cdat$data$X1ki==0]),digits=dig) ## get mean outcome among unexposed
  m0prev <- round(100*mean(mdat$data$Yki[mdat$data$X1ki==0]),digits=dig)
  # c1prev <- round(100*mean(cdat$data$Yki[cdat$data$X1ki==1]),digits=dig) ## get mean outcome among exposed
  m1prev <- round(100*mean(mdat$data$Yki[mdat$data$X1ki==1]),digits=dig)



  # cN1prev <- round(100*mean(cdat$data$Yki[cdat$data$Nki==1]),digits=dig) ## get mean outcome among exposed
  mN1prev <- round(100*mean(mdat$data$Yki[mdat$data$Nki==1]),digits=dig)
  # cN2prev <- round(100*mean(cdat$data$Yki[cdat$data$Nki==2]),digits=dig) ## get mean outcome among exposed
  mN2prev <- round(100*mean(mdat$data$Yki[mdat$data$Nki==2]),digits=dig)
  # cN3prev <- round(100*mean(cdat$data$Yki[cdat$data$Nki==3]),digits=dig) ## get mean outcome among exposed
  mN3prev <- round(100*mean(mdat$data$Yki[mdat$data$Nki==3]),digits=dig)
  # cN4prev <- round(100*mean(cdat$data$Yki[cdat$data$Nki==4]),digits=dig) ## get mean outcome among exposed
  mN4prev <- round(100*mean(mdat$data$Yki[mdat$data$Nki==4]),digits=dig)
  # cN5prev <- round(100*mean(cdat$data$Yki[cdat$data$Nki>=5]),digits=dig) ## get mean outcome among exposed
  mN5prev <- round(100*mean(mdat$data$Yki[mdat$data$Nki>=5]),digits=dig)

  return(list(#cprev=cprev,
              mprev=mprev,
              #c01prev=c(c0prev,c1prev),
              m01prev=c(m0prev,m1prev),
              #cNprev=c(cN1prev,cN2prev,cN3prev,cN4prev,cN5prev),
              mNprev=c(mN1prev,mN2prev,mN3prev,mN4prev,mN5prev)))

}

##############################################################
## generate outcome prevalences for each simulation setting

prev_gam00 <- get_outDist(gam=0.00)
prev_gam25 <- get_outDist(gam=-0.25)
prev_gamP25 <- get_outDist(gam=+0.25)


##############################################################
## Make prevalence tables for exposed and unexposed

## marginal
mprev <- rbind(prev_gam00$mprev,prev_gam25$mprev,prev_gamP25$mprev)
m01prev <- rbind(prev_gam00$m01prev,prev_gam25$m01prev,prev_gamP25$m01prev)
prev_tab <- cbind(mprev,m01prev)
rownames(prev_tab) <- c("NonICS","ICS Gam=-0.25","ICS Gam=0.25")
colnames(prev_tab) <- c("Overall","X=0","X=1")
xtable(prev_tab,digits = 1)




##############################################################
## Make prevalence vs SIZE tables

## marginal
mNprev <- rbind(prev_gam00$mNprev,prev_gam25$mNprev,prev_gamP25$mNprev)
rownames(mNprev) <- c("NonICS","ICS Gam=-0.25","ICS Gam=0.25")
colnames(mNprev) <- c("Nk=1","2","3","4","5+")
xtable(mNprev,digits=1)


