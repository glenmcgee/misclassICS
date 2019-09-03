###############################################################
##              Produce Results Tables and Plots             ##
###############################################################
## Edited: 05/29/2018
library(data.table)
library(xtable)
library(wesanderson)
library(tidyverse)
library(gridExtra)
library(kableExtra)
library(dplyr)

## True parameter values
a0 <- 0.6; a1 <- -0.2 ## volume model coefficients
b0 <- -4; b1 <- 0.5; b2 <- 0.2 ## outcome model coefficients
e0 <- -2; e1 <- 1 ## zero-inflation model coefficients
gamma <- 0#,-0.1,-0.25 ## scaling factor for random intercept in volume model
## intercepts model
sigma <- 2 ## SD of rand int
## slopes model
sigma0 <- 2 ## SD of rand eff for X=0
sigma1 <- 1.5 ## SD of rand eff for X=1


## set path
path <- "../Simulations/Results/"
outpath <- "../Simulations/Plots/"

## colors
wes_red <- wes_palette(n=5, name="Darjeeling")[1]
wes_green <- wes_palette(n=5, name="Darjeeling")[2]
wes_gold <- wes_palette(n=5, name="Darjeeling")[3]
wes_blue <- wes_palette(n=5, name="Darjeeling")[5]



################################################
##                  Load Data                 ##
################################################
results_simG0_Fixed <- read.table(paste0(path,"results_simG0_Fixed.txt"),header = TRUE);
results_simG0_BothDec <- read.table(paste0(path,"results_simG0_BothDec.txt"),header = TRUE);
results_simG0_BothInc <- read.table(paste0(path,"results_simG0_BothInc.txt"),header = TRUE);
results_simG0_SensDec <- read.table(paste0(path,"results_simG0_SensDec.txt"),header = TRUE);
results_simGN25_Fixed <- read.table(paste0(path,"results_simGN25_Fixed.txt"),header = TRUE);
results_simGN25_BothDec <- read.table(paste0(path,"results_simGN25_BothDec.txt"),header = TRUE);
results_simGN25_BothInc <- read.table(paste0(path,"results_simGN25_BothInc.txt"),header = TRUE);
results_simGN25_SensDec <- read.table(paste0(path,"results_simGN25_SensDec.txt"),header = TRUE);
results_simG25_Fixed <- read.table(paste0(path,"results_simG25_Fixed.txt"),header = TRUE);
results_simG25_BothDec <- read.table(paste0(path,"results_simG25_BothDec.txt"),header = TRUE);
results_simG25_BothInc <- read.table(paste0(path,"results_simG25_BothInc.txt"),header = TRUE);
results_simG25_SensDec <- read.table(paste0(path,"results_simG25_SensDec.txt"),header = TRUE);
#
results_simGN25_BothDecCOND <- read.table(paste0(path,"results_simGN25_BothDecCOND.txt"),header = TRUE); ## CONDITIONAL MODELS ONLY ## gamma=-0.25, sens/spec decrease in N
#
results_simG0_FixedINDUCED <- read.table(paste0(path,"results_simG0_FixedINDUCED.txt"),header = TRUE); ## INDUCED INFORMATIVENESS: a1=-1, b1=1, gamma=0, simple misclass
results_simG0_BothDecINDUCED <- read.table(paste0(path,"results_simG0_BothDecINDUCED.txt"),header = TRUE); ## INDUCED INFORMATIVENESS: a1=-1, b1=1, gamma=0, simple misclass
results_simG0_BothIncINDUCED <- read.table(paste0(path,"results_simG0_BothIncINDUCED.txt"),header = TRUE); ## INDUCED INFORMATIVENESS: a1=-1, b1=1, gamma=0, simple misclass

################################################
##        Functions for Data Processing       ##
################################################
## function to identify outliers
get_outliers <- function(df,no_sd=4){
  ## get data
  df <- df[,order(colnames(df))] ## order data (just in case)
  ests <- df[,grep(pattern="est_",names(df))] ## get estimates
  meanmat <- matrix(apply(ests,2,mean,na.rm=T),byrow=T,nrow=nrow(ests),ncol=ncol(ests)) ## matrix of columnwise SDs
  sdmat <- matrix(apply(ests,2,sd,na.rm=T),byrow=T,nrow=nrow(ests),ncol=ncol(ests)) ## matrix of columnwise SDs
  # ## get parameters
  # true_vals <- get_truevals(gamma=gamma,ZIP=ZIP,slopes=slopes)
  # true_mat <- matrix(true_vals,nrow=nrow(ests),ncol=length(true_vals),byrow=TRUE) ## matrix of true parameters

  res <- abs(ests-meanmat)>no_sd*sdmat ## outliers
  which_outliers <- which(apply(res,1,sum)>0) ## which rows contain outliers
  return(which_outliers)
}

## function to get true parameters
get_truevals <- function(gamma=0,ZIP=FALSE,slopes=FALSE,suppsim=FALSE,COND=FALSE){

  ## add alphas and betas
  if(suppsim==FALSE & COND==FALSE){
    true_vals <- c(rep(c(a0,a1),each=4),  ## alpha only in joint models
                   rep(c(b0,b1,b2),each=(3+7+7+4+4))) ## beta in all models
  }else if(COND==TRUE){
    true_vals <- c(rep(c(a0,a1),each=4),  ## alpha only in joint models
                   rep(c(b0,b1,b2),each=(4+4))) ## beta in all models
  }else{## stronger a1/b1 for supplementary induced informativeness simulation
    true_vals <- c(rep(c(a0,-1),each=4),  ## alpha only in joint models
                   rep(c(b0,1,b2),each=(3+7+7+4+4))) ## beta in all models
  }


  ## add epsilons if ZIP model
  if(ZIP==TRUE){
    true_vals <- c(true_vals,
                   rep(c(e0,e1),each=4)) ## joint models EXCLUDING no0s (since it cant have zeros)
  }

  ## add gammas and sigmas as required
  if(slopes==FALSE){ ## random intercepts
    true_vals <- c(true_vals,
                   rep(gamma,each=4), ## only joint models
                   rep(sigma,each=8)) ## joint models AND naive (likelihood based methods)

  }else{ ## random slopes
    true_vals <- c(true_vals,
                   rep(c(gamma,gamma),each=4), ## only joint models ## currently gamma0=gamma1
                   rep(c(sigma0,sigma1),each=8)) ## joint models AND naive (likelihood based methods)
  }
  return(true_vals)
}




## function to compute pctbias and coverage (to later be averaged)
process_data <- function(df,gamma=0,ZIP=FALSE,slopes=TRUE,suppsim=FALSE,COND=FALSE){

  ## get data
  df <- df[,order(colnames(df))] ## order data (just in case)
  ests <- df[,grep(pattern="est_",names(df))] ## get estimates
  SEs <- df[,grep(pattern="SE_",names(df))] ## get standard errors

  ## get parameters
  true_vals <- get_truevals(gamma=gamma,ZIP=ZIP,slopes=slopes,suppsim=suppsim,COND=COND)
  true_mat <- matrix(true_vals,nrow=nrow(ests),ncol=length(true_vals),byrow=TRUE) ## matrix of true parameters

  ## get pctbias
  bias <- 100*(ests-true_mat)/true_mat   ## is true value covered by 95% Wald interval?
  colnames(bias) <- gsub(pattern="est",replacement="pctbias",colnames(bias))   ## rename columns

  ## get MSE
  mse <- (ests-true_mat)^2   ## squared difference
  colnames(mse) <- gsub(pattern="est",replacement="mse",colnames(mse))   ## rename columns

  ## get coveraage
  cvg <- 100*(true_mat>=ests-1.96*SEs & true_mat<=ests+1.96*SEs)   ## is true value covered by 95% Wald interval?
  colnames(cvg) <- gsub(pattern="est",replacement="cvg",colnames(cvg))   ## rename columns

  ## comine data
  df <- cbind(df,bias,mse,cvg)

  return(df)
}


################################################
##              Data Processing               ##
################################################
## get and exclude outliers
# outliers_pois_gam25_m1 <- get_outliers(res_pois_gam25_m1); res_pois_gam25_m1 <- res_pois_gam25_m1[-outliers_pois_gam25_m1,]


res_simG0_Fixed <- process_data(results_simG0_Fixed,gamma=-0.00,ZIP=FALSE,slopes=TRUE)
res_simG0_BothDec <- process_data(results_simG0_BothDec,gamma=-0.00,ZIP=FALSE,slopes=TRUE)
res_simG0_BothInc <- process_data(results_simG0_BothInc,gamma=-0.00,ZIP=FALSE,slopes=TRUE)
res_simG0_SensDec <- process_data(results_simG0_SensDec,gamma=-0.00,ZIP=FALSE,slopes=TRUE)
res_simGN25_Fixed <- process_data(results_simGN25_Fixed,gamma=-0.25,ZIP=FALSE,slopes=TRUE)
res_simGN25_BothDec <- process_data(results_simGN25_BothDec,gamma=-0.25,ZIP=FALSE,slopes=TRUE)
res_simGN25_BothInc <- process_data(results_simGN25_BothInc,gamma=-0.25,ZIP=FALSE,slopes=TRUE)
res_simGN25_SensDec <- process_data(results_simGN25_SensDec,gamma=-0.25,ZIP=FALSE,slopes=TRUE)
res_simG25_Fixed <- process_data(results_simG25_Fixed,gamma=0.25,ZIP=FALSE,slopes=TRUE)
res_simG25_BothDec <- process_data(results_simG25_BothDec,gamma=0.25,ZIP=FALSE,slopes=TRUE)
res_simG25_BothInc <- process_data(results_simG25_BothInc,gamma=0.25,ZIP=FALSE,slopes=TRUE)
res_simG25_SensDec <- process_data(results_simG25_SensDec,gamma=0.25,ZIP=FALSE,slopes=TRUE)
#
res_simGN25_BothDecCOND <- process_data(results_simGN25_BothDecCOND,gamma=-0.25,ZIP=FALSE,slopes=TRUE,COND=TRUE)
#
res_simG0_FixedINDUCED <- process_data(results_simG0_FixedINDUCED,gamma=-0.00,ZIP=FALSE,slopes=TRUE,suppsim=TRUE)
res_simG0_BothDecINDUCED <- process_data(results_simG0_BothDecINDUCED,gamma=-0.00,ZIP=FALSE,slopes=TRUE,suppsim=TRUE)
res_simG0_BothIncINDUCED <- process_data(results_simG0_BothIncINDUCED,gamma=-0.00,ZIP=FALSE,slopes=TRUE,suppsim=TRUE)

# ## list of data files
# pois_ls <- list(res_pois_gam00_m1,res_pois_gam10_m1,res_pois_gam25_m1,
#                 res_pois_gam00_m0,res_pois_gam10_m0,res_pois_gam25_m0)


################################################
##                   Tables                   ##
################################################


## handle NaNs
is.nan.data.frame <- function(x){
  do.call(cbind, lapply(x, is.nan))
}


make_tab <- function(df,stat="est",func=mean,ZIP=FALSE,slopes=TRUE,print=TRUE,dig=2,scale=1,COND=FALSE){#,root=FALSE){

  if(COND==FALSE){
    ## gee
    B0 <- apply(cbind(df[,paste0("B0",stat,"_gee","_complete")],df[,paste0("B0",stat,"_gee","_mis")],df[,paste0("B0",stat,"_gee","_valid")]),2,func,na.rm=TRUE)
    B1 <- apply(cbind(df[,paste0("B1",stat,"_gee","_complete")],df[,paste0("B1",stat,"_gee","_mis")],df[,paste0("B1",stat,"_gee","_valid")]),2,func,na.rm=TRUE)
    B2 <- apply(cbind(df[,paste0("B2",stat,"_gee","_complete")],df[,paste0("B2",stat,"_gee","_mis")],df[,paste0("B2",stat,"_gee","_valid")]),2,func,na.rm=TRUE)
    gee <- rbind(B0,B1,B2)
    colnames(gee) <- c("Comp","Mis","Valid")
    gee[is.nan(gee)] <- NA

    ## iee
    B0 <- apply(cbind(df[,paste0("B0",stat,"_iee","_complete")],df[,paste0("B0",stat,"_iee","_mis")],df[,paste0("B0",stat,"_iee","_valid")],df[,paste0("B0",stat,"_iee","_simple_EEEnoN")],df[,paste0("B0",stat,"_iee","_simple_EEE")],df[,paste0("B0",stat,"_iee","_EEEnoN")],df[,paste0("B0",stat,"_iee","_EEE")]),2,func,na.rm=TRUE)
    B1 <- apply(cbind(df[,paste0("B1",stat,"_iee","_complete")],df[,paste0("B1",stat,"_iee","_mis")],df[,paste0("B1",stat,"_iee","_valid")],df[,paste0("B1",stat,"_iee","_simple_EEEnoN")],df[,paste0("B1",stat,"_iee","_simple_EEE")],df[,paste0("B1",stat,"_iee","_EEEnoN")],df[,paste0("B1",stat,"_iee","_EEE")]),2,func,na.rm=TRUE)
    B2 <- apply(cbind(df[,paste0("B2",stat,"_iee","_complete")],df[,paste0("B2",stat,"_iee","_mis")],df[,paste0("B2",stat,"_iee","_valid")],df[,paste0("B2",stat,"_iee","_simple_EEEnoN")],df[,paste0("B2",stat,"_iee","_simple_EEE")],df[,paste0("B2",stat,"_iee","_EEEnoN")],df[,paste0("B2",stat,"_iee","_EEE")]),2,func,na.rm=TRUE)
    iee <- rbind(B0,B1,B2)
    colnames(iee) <- c("Comp","Mis","Valid","Simple NoN","Simple","EEE NoN","EEE")
    iee[is.nan(iee)] <- NA

    ## wee
    B0 <- apply(cbind(df[,paste0("B0",stat,"_wee","_complete")],df[,paste0("B0",stat,"_wee","_mis")],df[,paste0("B0",stat,"_wee","_valid")],df[,paste0("B0",stat,"_wee","_simple_EEEnoN")],df[,paste0("B0",stat,"_wee","_simple_EEE")],df[,paste0("B0",stat,"_wee","_EEEnoN")],df[,paste0("B0",stat,"_wee","_EEE")]),2,func,na.rm=TRUE)
    B1 <- apply(cbind(df[,paste0("B1",stat,"_wee","_complete")],df[,paste0("B1",stat,"_wee","_mis")],df[,paste0("B1",stat,"_wee","_valid")],df[,paste0("B1",stat,"_wee","_simple_EEEnoN")],df[,paste0("B1",stat,"_wee","_simple_EEE")],df[,paste0("B1",stat,"_wee","_EEEnoN")],df[,paste0("B1",stat,"_wee","_EEE")]),2,func,na.rm=TRUE)
    B2 <- apply(cbind(df[,paste0("B2",stat,"_wee","_complete")],df[,paste0("B2",stat,"_wee","_mis")],df[,paste0("B2",stat,"_wee","_valid")],df[,paste0("B2",stat,"_wee","_simple_EEEnoN")],df[,paste0("B2",stat,"_wee","_simple_EEE")],df[,paste0("B2",stat,"_wee","_EEEnoN")],df[,paste0("B2",stat,"_wee","_EEE")]),2,func,na.rm=TRUE)
    wee <- rbind(B0,B1,B2)
    colnames(wee) <- c("Comp","Mis","Valid","Simple NoN","Simple","EEE NoN","EEE")
    wee[is.nan(wee)] <- NA
  }

  ## jmm
  # E0 <- E1 <- NULL
  A0 <- apply(cbind(df[,paste0("A0",stat,"_jmm","_complete")],df[,paste0("A0",stat,"_jmm","_mis")],df[,paste0("A0",stat,"_jmm","_valid")],df[,paste0("A0",stat,"_jmm","_obsLik")]),2,func,na.rm=TRUE)
  A1 <- apply(cbind(df[,paste0("A1",stat,"_jmm","_complete")],df[,paste0("A1",stat,"_jmm","_mis")],df[,paste0("A1",stat,"_jmm","_valid")],df[,paste0("A1",stat,"_jmm","_obsLik")]),2,func,na.rm=TRUE)
  B0 <- apply(cbind(df[,paste0("B0",stat,"_jmm","_complete")],df[,paste0("B0",stat,"_jmm","_mis")],df[,paste0("B0",stat,"_jmm","_valid")],df[,paste0("B0",stat,"_jmm","_obsLik")]),2,func,na.rm=TRUE)
  B1 <- apply(cbind(df[,paste0("B1",stat,"_jmm","_complete")],df[,paste0("B1",stat,"_jmm","_mis")],df[,paste0("B1",stat,"_jmm","_valid")],df[,paste0("B1",stat,"_jmm","_obsLik")]),2,func,na.rm=TRUE)
  B2 <- apply(cbind(df[,paste0("B2",stat,"_jmm","_complete")],df[,paste0("B2",stat,"_jmm","_mis")],df[,paste0("B2",stat,"_jmm","_valid")],df[,paste0("B2",stat,"_jmm","_obsLik")]),2,func,na.rm=TRUE)
  if(slopes==FALSE){
    Gam <- apply(cbind(df[,paste0("Gam",stat,"_jmm","_complete")],df[,paste0("Gam",stat,"_jmm","_mis")],df[,paste0("Gam",stat,"_jmm","_valid")],df[,paste0("Gam",stat,"_jmm","_obsLik")]),2,func,na.rm=TRUE)
    Sig <- apply(cbind(df[,paste0("Sig",stat,"_jmm","_complete")],df[,paste0("Sig",stat,"_jmm","_mis")],df[,paste0("Sig",stat,"_jmm","_valid")],df[,paste0("Sig",stat,"_jmm","_obsLik")]),2,func,na.rm=TRUE)
    jmm <- rbind(A0,A1,B0,B1,B2,Sig,Gam)
  }else{
    Gam0 <- apply(cbind(df[,paste0("Gam0",stat,"_jmm","_complete")],df[,paste0("Gam0",stat,"_jmm","_mis")],df[,paste0("Gam0",stat,"_jmm","_valid")],df[,paste0("Gam0",stat,"_jmm","_obsLik")]),2,func,na.rm=TRUE)
    Gam1 <- apply(cbind(df[,paste0("Gam1",stat,"_jmm","_complete")],df[,paste0("Gam1",stat,"_jmm","_mis")],df[,paste0("Gam1",stat,"_jmm","_valid")],df[,paste0("Gam1",stat,"_jmm","_obsLik")]),2,func,na.rm=TRUE)
    Sig0 <- apply(cbind(df[,paste0("Sig0",stat,"_jmm","_complete")],df[,paste0("Sig0",stat,"_jmm","_mis")],df[,paste0("Sig0",stat,"_jmm","_valid")],df[,paste0("Sig0",stat,"_jmm","_obsLik")]),2,func,na.rm=TRUE)
    Sig1 <- apply(cbind(df[,paste0("Sig1",stat,"_jmm","_complete")],df[,paste0("Sig1",stat,"_jmm","_mis")],df[,paste0("Sig1",stat,"_jmm","_valid")],df[,paste0("Sig1",stat,"_jmm","_obsLik")]),2,func,na.rm=TRUE)
    jmm <- rbind(A0,A1,B0,B1,B2,Sig0,Sig1,Gam0,Gam1)
  }
  # if(ZIP==TRUE){
  #   E0 <- apply(cbind(df[,paste0("E0",stat,"_jmm","_complete")],df[,paste0("E0",stat,"_jmm","_mis")],df[,paste0("E0",stat,"_jmm","_valid")],df[,paste0("E0",stat,"_jmm","_obsLik")]),2,func,na.rm=TRUE)
  #   E1 <- apply(cbind(df[,paste0("E1",stat,"_jmm","_complete")],df[,paste0("E1",stat,"_jmm","_mis")],df[,paste0("E1",stat,"_jmm","_valid")],df[,paste0("E1",stat,"_jmm","_obsLik")]),2,func,na.rm=TRUE)
  # }
  colnames(jmm) <- c("Comp","Mis","Valid","ObsLik")
  jmm[is.nan(jmm)] <- NA

  ## mm
  B0 <- apply(cbind(df[,paste0("B0",stat,"_mm","_complete")],df[,paste0("B0",stat,"_mm","_mis")],df[,paste0("B0",stat,"_mm","_valid")],df[,paste0("B0",stat,"_mm","_obsLik")]),2,func,na.rm=TRUE)
  B1 <- apply(cbind(df[,paste0("B1",stat,"_mm","_complete")],df[,paste0("B1",stat,"_mm","_mis")],df[,paste0("B1",stat,"_mm","_valid")],df[,paste0("B1",stat,"_mm","_obsLik")]),2,func,na.rm=TRUE)
  B2 <- apply(cbind(df[,paste0("B2",stat,"_mm","_complete")],df[,paste0("B2",stat,"_mm","_mis")],df[,paste0("B2",stat,"_mm","_valid")],df[,paste0("B2",stat,"_mm","_obsLik")]),2,func,na.rm=TRUE)
  if(slopes==FALSE){
    Sig <- apply(cbind(df[,paste0("Sig",stat,"_mm","_complete")],df[,paste0("Sig",stat,"_mm","_mis")],df[,paste0("Sig",stat,"_mm","_valid")],df[,paste0("Sig",stat,"_mm","_obsLik")]),2,func,na.rm=TRUE)
    mm <- rbind(B0,B1,B2,Sig,Gam)
  }else{
    Sig0 <- apply(cbind(df[,paste0("Sig0",stat,"_mm","_complete")],df[,paste0("Sig0",stat,"_mm","_mis")],df[,paste0("Sig0",stat,"_mm","_valid")],df[,paste0("Sig0",stat,"_mm","_obsLik")]),2,func,na.rm=TRUE)
    Sig1 <- apply(cbind(df[,paste0("Sig1",stat,"_mm","_complete")],df[,paste0("Sig1",stat,"_mm","_mis")],df[,paste0("Sig1",stat,"_mm","_valid")],df[,paste0("Sig1",stat,"_mm","_obsLik")]),2,func,na.rm=TRUE)
    mm <- rbind(B0,B1,B2,Sig0,Sig1,Gam0,Gam1)
  }
  colnames(mm) <- c("Comp","Mis","Valid","ObsLik")
  mm[is.nan(mm)] <- NA

  ## convert to string to preserve rounding
  if(COND==FALSE){
    ls <- list(gee=data.frame(gee),iee=data.frame(iee),wee=data.frame(wee),mm=data.frame(mm),jmm=data.frame(jmm))
  }else{
    ls <- list(mm=data.frame(mm),jmm=data.frame(jmm))
  }

  # newfun
  for(ll in 1:length(ls)){
    # if(root==TRUE){ls[[ll]] <- 1000*data.frame(ls[[ll]]^2)}  ## for RMSE
    ls[[ll]] <- data.frame(scale*ls[[ll]])
    test <- ls[[ll]]%>% mutate_all(function(x) format(round(x,dig),nsmall=dig))
    rownames(test) <- rownames(ls[[ll]]) ## preserve row names
    ls[[ll]] <- test
  }

  return(ls)
}

master_tab <- function(ls,nams=c("Mean Est","SD","100xMSE","Cvg")){

  for(ll in 1:length(ls)){
    # ls[[ll]] <- round(ls[[ll]],dig) ## round
    ls[[ll]] <- data.frame(cbind(index=seq(1,nrow(ls[[ll]])),param=rownames(ls[[ll]]),stat=nams[ll],ls[[ll]])) ## add index for reordering
  }

  res <- data.frame(rbindlist(ls)[order(index),-1]) ## combine
  # rownames(res) <- paste(res$Param,res$Gamma) ## unique row name
  res <- data.frame(res) ## remove redundant column

  return(res)
}




mean_simG0_Fixed <- make_tab(res_simG0_Fixed,stat="est",func=mean,ZIP=FALSE,slopes=TRUE,print=TRUE,dig=2)
median_simG0_Fixed <- make_tab(res_simG0_Fixed,stat="est",func=median,ZIP=FALSE,slopes=TRUE,print=TRUE,dig=2)
pctbias_simG0_Fixed <- make_tab(res_simG0_Fixed,stat="pctbias",func=mean,ZIP=FALSE,slopes=TRUE,print=TRUE,dig=0)
median_simG0_Fixed <- make_tab(res_simG0_Fixed,stat="pctbias",func=median,ZIP=FALSE,slopes=TRUE,print=TRUE,dig=0)
cvg_simG0_Fixed <- make_tab(res_simG0_Fixed,stat="cvg",func=mean,ZIP=FALSE,slopes=TRUE,print=TRUE,dig=0)
sd_simG0_Fixed <- make_tab(res_simG0_Fixed,stat="est",func=sd,ZIP=FALSE,slopes=TRUE,print=TRUE,dig=2)
mse_simG0_Fixed <- make_tab(res_simG0_Fixed,stat="mse",func=mean,ZIP=FALSE,slopes=TRUE,print=TRUE,dig=2,scale=100)

mean_simG0_BothDec <- make_tab(res_simG0_BothDec,stat="est",func=mean,ZIP=FALSE,slopes=TRUE,print=TRUE,dig=2)
median_simG0_BothDec <- make_tab(res_simG0_BothDec,stat="est",func=median,ZIP=FALSE,slopes=TRUE,print=TRUE,dig=2)
pctbias_simG0_BothDec <- make_tab(res_simG0_BothDec,stat="pctbias",func=mean,ZIP=FALSE,slopes=TRUE,print=TRUE,dig=0)
median_simG0_BothDec <- make_tab(res_simG0_BothDec,stat="pctbias",func=median,ZIP=FALSE,slopes=TRUE,print=TRUE,dig=0)
cvg_simG0_BothDec <- make_tab(res_simG0_BothDec,stat="cvg",func=mean,ZIP=FALSE,slopes=TRUE,print=TRUE,dig=0)
sd_simG0_BothDec <- make_tab(res_simG0_BothDec,stat="est",func=sd,ZIP=FALSE,slopes=TRUE,print=TRUE,dig=2)
mse_simG0_BothDec <- make_tab(res_simG0_BothDec,stat="mse",func=mean,ZIP=FALSE,slopes=TRUE,print=TRUE,dig=2,scale=100)

mean_simG0_BothInc <- make_tab(res_simG0_BothInc,stat="est",func=mean,ZIP=FALSE,slopes=TRUE,print=TRUE,dig=2)
median_simG0_BothInc <- make_tab(res_simG0_BothInc,stat="est",func=median,ZIP=FALSE,slopes=TRUE,print=TRUE,dig=2)
pctbias_simG0_BothInc <- make_tab(res_simG0_BothInc,stat="pctbias",func=mean,ZIP=FALSE,slopes=TRUE,print=TRUE,dig=0)
median_simG0_BothInc <- make_tab(res_simG0_BothInc,stat="pctbias",func=median,ZIP=FALSE,slopes=TRUE,print=TRUE,dig=0)
cvg_simG0_BothInc <- make_tab(res_simG0_BothInc,stat="cvg",func=mean,ZIP=FALSE,slopes=TRUE,print=TRUE,dig=0)
sd_simG0_BothInc <- make_tab(res_simG0_BothInc,stat="est",func=sd,ZIP=FALSE,slopes=TRUE,print=TRUE,dig=2)
mse_simG0_BothInc <- make_tab(res_simG0_BothInc,stat="mse",func=mean,ZIP=FALSE,slopes=TRUE,print=TRUE,dig=2,scale=100)

mean_simG0_SensDec <- make_tab(res_simG0_SensDec,stat="est",func=mean,ZIP=FALSE,slopes=TRUE,print=TRUE,dig=2)
median_simG0_SensDec <- make_tab(res_simG0_SensDec,stat="est",func=median,ZIP=FALSE,slopes=TRUE,print=TRUE,dig=2)
pctbias_simG0_SensDec <- make_tab(res_simG0_SensDec,stat="pctbias",func=mean,ZIP=FALSE,slopes=TRUE,print=TRUE,dig=0)
median_simG0_SensDec <- make_tab(res_simG0_SensDec,stat="pctbias",func=median,ZIP=FALSE,slopes=TRUE,print=TRUE,dig=0)
cvg_simG0_SensDec <- make_tab(res_simG0_SensDec,stat="cvg",func=mean,ZIP=FALSE,slopes=TRUE,print=TRUE,dig=0)
sd_simG0_SensDec <- make_tab(res_simG0_SensDec,stat="est",func=sd,ZIP=FALSE,slopes=TRUE,print=TRUE,dig=2)
mse_simG0_SensDec <- make_tab(res_simG0_SensDec,stat="mse",func=mean,ZIP=FALSE,slopes=TRUE,print=TRUE,dig=2,scale=100)



mean_simGN25_Fixed <- make_tab(res_simGN25_Fixed,stat="est",func=mean,ZIP=FALSE,slopes=TRUE,print=TRUE,dig=2)
median_simGN25_Fixed <- make_tab(res_simGN25_Fixed,stat="est",func=median,ZIP=FALSE,slopes=TRUE,print=TRUE,dig=2)
pctbias_simGN25_Fixed <- make_tab(res_simGN25_Fixed,stat="pctbias",func=mean,ZIP=FALSE,slopes=TRUE,print=TRUE,dig=0)
median_simGN25_Fixed <- make_tab(res_simGN25_Fixed,stat="pctbias",func=median,ZIP=FALSE,slopes=TRUE,print=TRUE,dig=0)
cvg_simGN25_Fixed <- make_tab(res_simGN25_Fixed,stat="cvg",func=mean,ZIP=FALSE,slopes=TRUE,print=TRUE,dig=0)
sd_simGN25_Fixed <- make_tab(res_simGN25_Fixed,stat="est",func=sd,ZIP=FALSE,slopes=TRUE,print=TRUE,dig=2)
mse_simGN25_Fixed <- make_tab(res_simGN25_Fixed,stat="mse",func=mean,ZIP=FALSE,slopes=TRUE,print=TRUE,dig=2,scale=100)

mean_simGN25_BothDec <- make_tab(res_simGN25_BothDec,stat="est",func=mean,ZIP=FALSE,slopes=TRUE,print=TRUE,dig=2)
median_simGN25_BothDec <- make_tab(res_simGN25_BothDec,stat="est",func=median,ZIP=FALSE,slopes=TRUE,print=TRUE,dig=2)
pctbias_simGN25_BothDec <- make_tab(res_simGN25_BothDec,stat="pctbias",func=mean,ZIP=FALSE,slopes=TRUE,print=TRUE,dig=0)
median_simGN25_BothDec <- make_tab(res_simGN25_BothDec,stat="pctbias",func=median,ZIP=FALSE,slopes=TRUE,print=TRUE,dig=0)
cvg_simGN25_BothDec <- make_tab(res_simGN25_BothDec,stat="cvg",func=mean,ZIP=FALSE,slopes=TRUE,print=TRUE,dig=0)
sd_simGN25_BothDec <- make_tab(res_simGN25_BothDec,stat="est",func=sd,ZIP=FALSE,slopes=TRUE,print=TRUE,dig=2)
mse_simGN25_BothDec <- make_tab(res_simGN25_BothDec,stat="mse",func=mean,ZIP=FALSE,slopes=TRUE,print=TRUE,dig=2,scale=100)

mean_simGN25_BothInc <- make_tab(res_simGN25_BothInc,stat="est",func=mean,ZIP=FALSE,slopes=TRUE,print=TRUE,dig=2)
median_simGN25_BothInc <- make_tab(res_simGN25_BothInc,stat="est",func=median,ZIP=FALSE,slopes=TRUE,print=TRUE,dig=2)
pctbias_simGN25_BothInc <- make_tab(res_simGN25_BothInc,stat="pctbias",func=mean,ZIP=FALSE,slopes=TRUE,print=TRUE,dig=0)
median_simGN25_BothInc <- make_tab(res_simGN25_BothInc,stat="pctbias",func=median,ZIP=FALSE,slopes=TRUE,print=TRUE,dig=0)
cvg_simGN25_BothInc <- make_tab(res_simGN25_BothInc,stat="cvg",func=mean,ZIP=FALSE,slopes=TRUE,print=TRUE,dig=0)
sd_simGN25_BothInc <- make_tab(res_simGN25_BothInc,stat="est",func=sd,ZIP=FALSE,slopes=TRUE,print=TRUE,dig=2)
mse_simGN25_BothInc <- make_tab(res_simGN25_BothInc,stat="mse",func=mean,ZIP=FALSE,slopes=TRUE,print=TRUE,dig=2,scale=100)

mean_simGN25_SensDec <- make_tab(res_simGN25_SensDec,stat="est",func=mean,ZIP=FALSE,slopes=TRUE,print=TRUE,dig=2)
median_simGN25_SensDec <- make_tab(res_simGN25_SensDec,stat="est",func=median,ZIP=FALSE,slopes=TRUE,print=TRUE,dig=2)
pctbias_simGN25_SensDec <- make_tab(res_simGN25_SensDec,stat="pctbias",func=mean,ZIP=FALSE,slopes=TRUE,print=TRUE,dig=0)
median_simGN25_SensDec <- make_tab(res_simGN25_SensDec,stat="pctbias",func=median,ZIP=FALSE,slopes=TRUE,print=TRUE,dig=0)
cvg_simGN25_SensDec <- make_tab(res_simGN25_SensDec,stat="cvg",func=mean,ZIP=FALSE,slopes=TRUE,print=TRUE,dig=0)
sd_simGN25_SensDec <- make_tab(res_simGN25_SensDec,stat="est",func=sd,ZIP=FALSE,slopes=TRUE,print=TRUE,dig=2)
mse_simGN25_SensDec <- make_tab(res_simGN25_SensDec,stat="mse",func=mean,ZIP=FALSE,slopes=TRUE,print=TRUE,dig=2,scale=100)



mean_simG25_Fixed <- make_tab(res_simG25_Fixed,stat="est",func=mean,ZIP=FALSE,slopes=TRUE,print=TRUE,dig=2)
median_simG25_Fixed <- make_tab(res_simG25_Fixed,stat="est",func=median,ZIP=FALSE,slopes=TRUE,print=TRUE,dig=2)
pctbias_simG25_Fixed <- make_tab(res_simG25_Fixed,stat="pctbias",func=mean,ZIP=FALSE,slopes=TRUE,print=TRUE,dig=0)
median_simG25_Fixed <- make_tab(res_simG25_Fixed,stat="pctbias",func=median,ZIP=FALSE,slopes=TRUE,print=TRUE,dig=0)
cvg_simG25_Fixed <- make_tab(res_simG25_Fixed,stat="cvg",func=mean,ZIP=FALSE,slopes=TRUE,print=TRUE,dig=0)
sd_simG25_Fixed <- make_tab(res_simG25_Fixed,stat="est",func=sd,ZIP=FALSE,slopes=TRUE,print=TRUE,dig=2)
mse_simG25_Fixed <- make_tab(res_simG25_Fixed,stat="mse",func=mean,ZIP=FALSE,slopes=TRUE,print=TRUE,dig=2,scale=100)

mean_simG25_BothDec <- make_tab(res_simG25_BothDec,stat="est",func=mean,ZIP=FALSE,slopes=TRUE,print=TRUE,dig=2)
median_simG25_BothDec <- make_tab(res_simG25_BothDec,stat="est",func=median,ZIP=FALSE,slopes=TRUE,print=TRUE,dig=2)
pctbias_simG25_BothDec <- make_tab(res_simG25_BothDec,stat="pctbias",func=mean,ZIP=FALSE,slopes=TRUE,print=TRUE,dig=0)
median_simG25_BothDec <- make_tab(res_simG25_BothDec,stat="pctbias",func=median,ZIP=FALSE,slopes=TRUE,print=TRUE,dig=0)
cvg_simG25_BothDec <- make_tab(res_simG25_BothDec,stat="cvg",func=mean,ZIP=FALSE,slopes=TRUE,print=TRUE,dig=0)
sd_simG25_BothDec <- make_tab(res_simG25_BothDec,stat="est",func=sd,ZIP=FALSE,slopes=TRUE,print=TRUE,dig=2)
mse_simG25_BothDec <- make_tab(res_simG25_BothDec,stat="mse",func=mean,ZIP=FALSE,slopes=TRUE,print=TRUE,dig=2,scale=100)

mean_simG25_BothInc <- make_tab(res_simG25_BothInc,stat="est",func=mean,ZIP=FALSE,slopes=TRUE,print=TRUE,dig=2)
median_simG25_BothInc <- make_tab(res_simG25_BothInc,stat="est",func=median,ZIP=FALSE,slopes=TRUE,print=TRUE,dig=2)
pctbias_simG25_BothInc <- make_tab(res_simG25_BothInc,stat="pctbias",func=mean,ZIP=FALSE,slopes=TRUE,print=TRUE,dig=0)
median_simG25_BothInc <- make_tab(res_simG25_BothInc,stat="pctbias",func=median,ZIP=FALSE,slopes=TRUE,print=TRUE,dig=0)
cvg_simG25_BothInc <- make_tab(res_simG25_BothInc,stat="cvg",func=mean,ZIP=FALSE,slopes=TRUE,print=TRUE,dig=0)
sd_simG25_BothInc <- make_tab(res_simG25_BothInc,stat="est",func=sd,ZIP=FALSE,slopes=TRUE,print=TRUE,dig=2)
mse_simG25_BothInc <- make_tab(res_simG25_BothInc,stat="mse",func=mean,ZIP=FALSE,slopes=TRUE,print=TRUE,dig=2,scale=100)

mean_simG25_SensDec <- make_tab(res_simG25_SensDec,stat="est",func=mean,ZIP=FALSE,slopes=TRUE,print=TRUE,dig=2)
median_simG25_SensDec <- make_tab(res_simG25_SensDec,stat="est",func=median,ZIP=FALSE,slopes=TRUE,print=TRUE,dig=2)
pctbias_simG25_SensDec <- make_tab(res_simG25_SensDec,stat="pctbias",func=mean,ZIP=FALSE,slopes=TRUE,print=TRUE,dig=0)
median_simG25_SensDec <- make_tab(res_simG25_SensDec,stat="pctbias",func=median,ZIP=FALSE,slopes=TRUE,print=TRUE,dig=0)
cvg_simG25_SensDec <- make_tab(res_simG25_SensDec,stat="cvg",func=mean,ZIP=FALSE,slopes=TRUE,print=TRUE,dig=0)
sd_simG25_SensDec <- make_tab(res_simG25_SensDec,stat="est",func=sd,ZIP=FALSE,slopes=TRUE,print=TRUE,dig=2)
mse_simG25_SensDec <- make_tab(res_simG25_SensDec,stat="mse",func=mean,ZIP=FALSE,slopes=TRUE,print=TRUE,dig=2,scale=100)

## conditional
mean_simGN25_BothDecCOND <- make_tab(res_simGN25_BothDecCOND,stat="est",func=mean,ZIP=FALSE,slopes=TRUE,print=TRUE,COND=TRUE,dig=2)
median_simGN25_BothDecCOND <- make_tab(res_simGN25_BothDecCOND,stat="est",func=median,ZIP=FALSE,slopes=TRUE,print=TRUE,COND=TRUE,dig=2)
pctbias_simGN25_BothDecCOND <- make_tab(res_simGN25_BothDecCOND,stat="pctbias",func=mean,ZIP=FALSE,slopes=TRUE,print=TRUE,COND=TRUE,dig=0)
median_simGN25_BothDecCOND <- make_tab(res_simGN25_BothDecCOND,stat="pctbias",func=median,ZIP=FALSE,slopes=TRUE,print=TRUE,COND=TRUE,dig=0)
cvg_simGN25_BothDecCOND <- make_tab(res_simGN25_BothDecCOND,stat="cvg",func=mean,ZIP=FALSE,slopes=TRUE,print=TRUE,COND=TRUE,dig=0)
sd_simGN25_BothDecCOND <- make_tab(res_simGN25_BothDecCOND,stat="est",func=sd,ZIP=FALSE,slopes=TRUE,print=TRUE,COND=TRUE,dig=2)
mse_simGN25_BothDecCOND <- make_tab(res_simGN25_BothDecCOND,stat="mse",func=mean,ZIP=FALSE,slopes=TRUE,print=TRUE,COND=TRUE,dig=2,scale=100)
## induced
mean_simG0_FixedINDUCED <- make_tab(res_simG0_FixedINDUCED,stat="est",func=mean,ZIP=FALSE,slopes=TRUE,print=TRUE,dig=2)
median_simG0_FixedINDUCED <- make_tab(res_simG0_FixedINDUCED,stat="est",func=median,ZIP=FALSE,slopes=TRUE,print=TRUE,dig=2)
pctbias_simG0_FixedINDUCED <- make_tab(res_simG0_FixedINDUCED,stat="pctbias",func=mean,ZIP=FALSE,slopes=TRUE,print=TRUE,dig=0)
median_simG0_FixedINDUCED <- make_tab(res_simG0_FixedINDUCED,stat="pctbias",func=median,ZIP=FALSE,slopes=TRUE,print=TRUE,dig=0)
cvg_simG0_FixedINDUCED <- make_tab(res_simG0_FixedINDUCED,stat="cvg",func=mean,ZIP=FALSE,slopes=TRUE,print=TRUE,dig=0)
sd_simG0_FixedINDUCED <- make_tab(res_simG0_FixedINDUCED,stat="est",func=sd,ZIP=FALSE,slopes=TRUE,print=TRUE,dig=2)
mse_simG0_FixedINDUCED <- make_tab(res_simG0_FixedINDUCED,stat="mse",func=mean,ZIP=FALSE,slopes=TRUE,print=TRUE,dig=2,scale=100)

mean_simG0_BothDecINDUCED <- make_tab(res_simG0_BothDecINDUCED,stat="est",func=mean,ZIP=FALSE,slopes=TRUE,print=TRUE,dig=2)
median_simG0_BothDecINDUCED <- make_tab(res_simG0_BothDecINDUCED,stat="est",func=median,ZIP=FALSE,slopes=TRUE,print=TRUE,dig=2)
pctbias_simG0_BothDecINDUCED <- make_tab(res_simG0_BothDecINDUCED,stat="pctbias",func=mean,ZIP=FALSE,slopes=TRUE,print=TRUE,dig=0)
median_simG0_BothDecINDUCED <- make_tab(res_simG0_BothDecINDUCED,stat="pctbias",func=median,ZIP=FALSE,slopes=TRUE,print=TRUE,dig=0)
cvg_simG0_BothDecINDUCED <- make_tab(res_simG0_BothDecINDUCED,stat="cvg",func=mean,ZIP=FALSE,slopes=TRUE,print=TRUE,dig=0)
sd_simG0_BothDecINDUCED <- make_tab(res_simG0_BothDecINDUCED,stat="est",func=sd,ZIP=FALSE,slopes=TRUE,print=TRUE,dig=2)
mse_simG0_BothDecINDUCED <- make_tab(res_simG0_BothDecINDUCED,stat="mse",func=mean,ZIP=FALSE,slopes=TRUE,print=TRUE,dig=2,scale=100)

mean_simG0_BothIncINDUCED <- make_tab(res_simG0_BothIncINDUCED,stat="est",func=mean,ZIP=FALSE,slopes=TRUE,print=TRUE,dig=2)
median_simG0_BothIncINDUCED <- make_tab(res_simG0_BothIncINDUCED,stat="est",func=median,ZIP=FALSE,slopes=TRUE,print=TRUE,dig=2)
pctbias_simG0_BothIncINDUCED <- make_tab(res_simG0_BothIncINDUCED,stat="pctbias",func=mean,ZIP=FALSE,slopes=TRUE,print=TRUE,dig=0)
median_simG0_BothIncINDUCED <- make_tab(res_simG0_BothIncINDUCED,stat="pctbias",func=median,ZIP=FALSE,slopes=TRUE,print=TRUE,dig=0)
cvg_simG0_BothIncINDUCED <- make_tab(res_simG0_BothIncINDUCED,stat="cvg",func=mean,ZIP=FALSE,slopes=TRUE,print=TRUE,dig=0)
sd_simG0_BothIncINDUCED <- make_tab(res_simG0_BothIncINDUCED,stat="est",func=sd,ZIP=FALSE,slopes=TRUE,print=TRUE,dig=2)
mse_simG0_BothIncINDUCED <- make_tab(res_simG0_BothIncINDUCED,stat="mse",func=mean,ZIP=FALSE,slopes=TRUE,print=TRUE,dig=2,scale=100)


## Master tables
simG0_Fixedtab_gee <- master_tab(ls=list(mean_simG0_Fixed$gee,sd_simG0_Fixed$gee,mse_simG0_Fixed$gee,cvg_simG0_Fixed$gee))
simG0_Fixedtab_iee <- master_tab(ls=list(mean_simG0_Fixed$iee,sd_simG0_Fixed$iee,mse_simG0_Fixed$iee,cvg_simG0_Fixed$iee))
simG0_Fixedtab_wee <- master_tab(ls=list(mean_simG0_Fixed$wee,sd_simG0_Fixed$wee,mse_simG0_Fixed$wee,cvg_simG0_Fixed$wee))
simG0_Fixedtab_mm <- master_tab(ls=list(mean_simG0_Fixed$mm,sd_simG0_Fixed$mm,mse_simG0_Fixed$mm,cvg_simG0_Fixed$mm))
simG0_Fixedtab_jmm <- master_tab(ls=list(mean_simG0_Fixed$jmm,sd_simG0_Fixed$jmm,mse_simG0_Fixed$jmm,cvg_simG0_Fixed$jmm))

simG0_BothDectab_gee <- master_tab(ls=list(mean_simG0_BothDec$gee,sd_simG0_BothDec$gee,mse_simG0_BothDec$gee,cvg_simG0_BothDec$gee))
simG0_BothDectab_iee <- master_tab(ls=list(mean_simG0_BothDec$iee,sd_simG0_BothDec$iee,mse_simG0_BothDec$iee,cvg_simG0_BothDec$iee))
simG0_BothDectab_wee <- master_tab(ls=list(mean_simG0_BothDec$wee,sd_simG0_BothDec$wee,mse_simG0_BothDec$wee,cvg_simG0_BothDec$wee))
simG0_BothDectab_mm <- master_tab(ls=list(mean_simG0_BothDec$mm,sd_simG0_BothDec$mm,mse_simG0_BothDec$mm,cvg_simG0_BothDec$mm))
simG0_BothDectab_jmm <- master_tab(ls=list(mean_simG0_BothDec$jmm,sd_simG0_BothDec$jmm,mse_simG0_BothDec$jmm,cvg_simG0_BothDec$jmm))

simG0_BothInctab_gee <- master_tab(ls=list(mean_simG0_BothInc$gee,sd_simG0_BothInc$gee,mse_simG0_BothInc$gee,cvg_simG0_BothInc$gee))
simG0_BothInctab_iee <- master_tab(ls=list(mean_simG0_BothInc$iee,sd_simG0_BothInc$iee,mse_simG0_BothInc$iee,cvg_simG0_BothInc$iee))
simG0_BothInctab_wee <- master_tab(ls=list(mean_simG0_BothInc$wee,sd_simG0_BothInc$wee,mse_simG0_BothInc$wee,cvg_simG0_BothInc$wee))
simG0_BothInctab_mm <- master_tab(ls=list(mean_simG0_BothInc$mm,sd_simG0_BothInc$mm,mse_simG0_BothInc$mm,cvg_simG0_BothInc$mm))
simG0_BothInctab_jmm <- master_tab(ls=list(mean_simG0_BothInc$jmm,sd_simG0_BothInc$jmm,mse_simG0_BothInc$jmm,cvg_simG0_BothInc$jmm))

simG0_SensDectab_gee <- master_tab(ls=list(mean_simG0_SensDec$gee,sd_simG0_SensDec$gee,mse_simG0_SensDec$gee,cvg_simG0_SensDec$gee))
simG0_SensDectab_iee <- master_tab(ls=list(mean_simG0_SensDec$iee,sd_simG0_SensDec$iee,mse_simG0_SensDec$iee,cvg_simG0_SensDec$iee))
simG0_SensDectab_wee <- master_tab(ls=list(mean_simG0_SensDec$wee,sd_simG0_SensDec$wee,mse_simG0_SensDec$wee,cvg_simG0_SensDec$wee))
simG0_SensDectab_mm <- master_tab(ls=list(mean_simG0_SensDec$mm,sd_simG0_SensDec$mm,mse_simG0_SensDec$mm,cvg_simG0_SensDec$mm))
simG0_SensDectab_jmm <- master_tab(ls=list(mean_simG0_SensDec$jmm,sd_simG0_SensDec$jmm,mse_simG0_SensDec$jmm,cvg_simG0_SensDec$jmm))


simGN25_Fixedtab_gee <- master_tab(ls=list(mean_simGN25_Fixed$gee,sd_simGN25_Fixed$gee,mse_simGN25_Fixed$gee,cvg_simGN25_Fixed$gee))
simGN25_Fixedtab_iee <- master_tab(ls=list(mean_simGN25_Fixed$iee,sd_simGN25_Fixed$iee,mse_simGN25_Fixed$iee,cvg_simGN25_Fixed$iee))
simGN25_Fixedtab_wee <- master_tab(ls=list(mean_simGN25_Fixed$wee,sd_simGN25_Fixed$wee,mse_simGN25_Fixed$wee,cvg_simGN25_Fixed$wee))
simGN25_Fixedtab_mm <- master_tab(ls=list(mean_simGN25_Fixed$mm,sd_simGN25_Fixed$mm,mse_simGN25_Fixed$mm,cvg_simGN25_Fixed$mm))
simGN25_Fixedtab_jmm <- master_tab(ls=list(mean_simGN25_Fixed$jmm,sd_simGN25_Fixed$jmm,mse_simGN25_Fixed$jmm,cvg_simGN25_Fixed$jmm))

simGN25_BothDectab_gee <- master_tab(ls=list(mean_simGN25_BothDec$gee,sd_simGN25_BothDec$gee,mse_simGN25_BothDec$gee,cvg_simGN25_BothDec$gee))
simGN25_BothDectab_iee <- master_tab(ls=list(mean_simGN25_BothDec$iee,sd_simGN25_BothDec$iee,mse_simGN25_BothDec$iee,cvg_simGN25_BothDec$iee))
simGN25_BothDectab_wee <- master_tab(ls=list(mean_simGN25_BothDec$wee,sd_simGN25_BothDec$wee,mse_simGN25_BothDec$wee,cvg_simGN25_BothDec$wee))
simGN25_BothDectab_mm <- master_tab(ls=list(mean_simGN25_BothDec$mm,sd_simGN25_BothDec$mm,mse_simGN25_BothDec$mm,cvg_simGN25_BothDec$mm))
simGN25_BothDectab_jmm <- master_tab(ls=list(mean_simGN25_BothDec$jmm,sd_simGN25_BothDec$jmm,mse_simGN25_BothDec$jmm,cvg_simGN25_BothDec$jmm))

simGN25_BothInctab_gee <- master_tab(ls=list(mean_simGN25_BothInc$gee,sd_simGN25_BothInc$gee,mse_simGN25_BothInc$gee,cvg_simGN25_BothInc$gee))
simGN25_BothInctab_iee <- master_tab(ls=list(mean_simGN25_BothInc$iee,sd_simGN25_BothInc$iee,mse_simGN25_BothInc$iee,cvg_simGN25_BothInc$iee))
simGN25_BothInctab_wee <- master_tab(ls=list(mean_simGN25_BothInc$wee,sd_simGN25_BothInc$wee,mse_simGN25_BothInc$wee,cvg_simGN25_BothInc$wee))
simGN25_BothInctab_mm <- master_tab(ls=list(mean_simGN25_BothInc$mm,sd_simGN25_BothInc$mm,mse_simGN25_BothInc$mm,cvg_simGN25_BothInc$mm))
simGN25_BothInctab_jmm <- master_tab(ls=list(mean_simGN25_BothInc$jmm,sd_simGN25_BothInc$jmm,mse_simGN25_BothInc$jmm,cvg_simGN25_BothInc$jmm))

simGN25_SensDectab_gee <- master_tab(ls=list(mean_simGN25_SensDec$gee,sd_simGN25_SensDec$gee,mse_simGN25_SensDec$gee,cvg_simGN25_SensDec$gee))
simGN25_SensDectab_iee <- master_tab(ls=list(mean_simGN25_SensDec$iee,sd_simGN25_SensDec$iee,mse_simGN25_SensDec$iee,cvg_simGN25_SensDec$iee))
simGN25_SensDectab_wee <- master_tab(ls=list(mean_simGN25_SensDec$wee,sd_simGN25_SensDec$wee,mse_simGN25_SensDec$wee,cvg_simGN25_SensDec$wee))
simGN25_SensDectab_mm <- master_tab(ls=list(mean_simGN25_SensDec$mm,sd_simGN25_SensDec$mm,mse_simGN25_SensDec$mm,cvg_simGN25_SensDec$mm))
simGN25_SensDectab_jmm <- master_tab(ls=list(mean_simGN25_SensDec$jmm,sd_simGN25_SensDec$jmm,mse_simGN25_SensDec$jmm,cvg_simGN25_SensDec$jmm))


simG25_Fixedtab_gee <- master_tab(ls=list(mean_simG25_Fixed$gee,sd_simG25_Fixed$gee,mse_simG25_Fixed$gee,cvg_simG25_Fixed$gee))
simG25_Fixedtab_iee <- master_tab(ls=list(mean_simG25_Fixed$iee,sd_simG25_Fixed$iee,mse_simG25_Fixed$iee,cvg_simG25_Fixed$iee))
simG25_Fixedtab_wee <- master_tab(ls=list(mean_simG25_Fixed$wee,sd_simG25_Fixed$wee,mse_simG25_Fixed$wee,cvg_simG25_Fixed$wee))
simG25_Fixedtab_mm <- master_tab(ls=list(mean_simG25_Fixed$mm,sd_simG25_Fixed$mm,mse_simG25_Fixed$mm,cvg_simG25_Fixed$mm))
simG25_Fixedtab_jmm <- master_tab(ls=list(mean_simG25_Fixed$jmm,sd_simG25_Fixed$jmm,mse_simG25_Fixed$jmm,cvg_simG25_Fixed$jmm))

simG25_BothDectab_gee <- master_tab(ls=list(mean_simG25_BothDec$gee,sd_simG25_BothDec$gee,mse_simG25_BothDec$gee,cvg_simG25_BothDec$gee))
simG25_BothDectab_iee <- master_tab(ls=list(mean_simG25_BothDec$iee,sd_simG25_BothDec$iee,mse_simG25_BothDec$iee,cvg_simG25_BothDec$iee))
simG25_BothDectab_wee <- master_tab(ls=list(mean_simG25_BothDec$wee,sd_simG25_BothDec$wee,mse_simG25_BothDec$wee,cvg_simG25_BothDec$wee))
simG25_BothDectab_mm <- master_tab(ls=list(mean_simG25_BothDec$mm,sd_simG25_BothDec$mm,mse_simG25_BothDec$mm,cvg_simG25_BothDec$mm))
simG25_BothDectab_jmm <- master_tab(ls=list(mean_simG25_BothDec$jmm,sd_simG25_BothDec$jmm,mse_simG25_BothDec$jmm,cvg_simG25_BothDec$jmm))

simG25_BothInctab_gee <- master_tab(ls=list(mean_simG25_BothInc$gee,sd_simG25_BothInc$gee,mse_simG25_BothInc$gee,cvg_simG25_BothInc$gee))
simG25_BothInctab_iee <- master_tab(ls=list(mean_simG25_BothInc$iee,sd_simG25_BothInc$iee,mse_simG25_BothInc$iee,cvg_simG25_BothInc$iee))
simG25_BothInctab_wee <- master_tab(ls=list(mean_simG25_BothInc$wee,sd_simG25_BothInc$wee,mse_simG25_BothInc$wee,cvg_simG25_BothInc$wee))
simG25_BothInctab_mm <- master_tab(ls=list(mean_simG25_BothInc$mm,sd_simG25_BothInc$mm,mse_simG25_BothInc$mm,cvg_simG25_BothInc$mm))
simG25_BothInctab_jmm <- master_tab(ls=list(mean_simG25_BothInc$jmm,sd_simG25_BothInc$jmm,mse_simG25_BothInc$jmm,cvg_simG25_BothInc$jmm))

simG25_SensDectab_gee <- master_tab(ls=list(mean_simG25_SensDec$gee,sd_simG25_SensDec$gee,mse_simG25_SensDec$gee,cvg_simG25_SensDec$gee))
simG25_SensDectab_iee <- master_tab(ls=list(mean_simG25_SensDec$iee,sd_simG25_SensDec$iee,mse_simG25_SensDec$iee,cvg_simG25_SensDec$iee))
simG25_SensDectab_wee <- master_tab(ls=list(mean_simG25_SensDec$wee,sd_simG25_SensDec$wee,mse_simG25_SensDec$wee,cvg_simG25_SensDec$wee))
simG25_SensDectab_mm <- master_tab(ls=list(mean_simG25_SensDec$mm,sd_simG25_SensDec$mm,mse_simG25_SensDec$mm,cvg_simG25_SensDec$mm))
simG25_SensDectab_jmm <- master_tab(ls=list(mean_simG25_SensDec$jmm,sd_simG25_SensDec$jmm,mse_simG25_SensDec$jmm,cvg_simG25_SensDec$jmm))

## cond
simGN25_BothDecCONDtab_mm <- master_tab(ls=list(mean_simGN25_BothDecCOND$mm,sd_simGN25_BothDecCOND$mm,mse_simGN25_BothDecCOND$mm,cvg_simGN25_BothDecCOND$mm))
simGN25_BothDecCONDtab_jmm <- master_tab(ls=list(mean_simGN25_BothDecCOND$jmm,sd_simGN25_BothDecCOND$jmm,mse_simGN25_BothDecCOND$jmm,cvg_simGN25_BothDecCOND$jmm))
## induced
simG0_FixedINDUCEDtab_gee <- master_tab(ls=list(mean_simG0_FixedINDUCED$gee,sd_simG0_FixedINDUCED$gee,mse_simG0_FixedINDUCED$gee,cvg_simG0_FixedINDUCED$gee))
simG0_FixedINDUCEDtab_iee <- master_tab(ls=list(mean_simG0_FixedINDUCED$iee,sd_simG0_FixedINDUCED$iee,mse_simG0_FixedINDUCED$iee,cvg_simG0_FixedINDUCED$iee))
simG0_FixedINDUCEDtab_wee <- master_tab(ls=list(mean_simG0_FixedINDUCED$wee,sd_simG0_FixedINDUCED$wee,mse_simG0_FixedINDUCED$wee,cvg_simG0_FixedINDUCED$wee))
simG0_FixedINDUCEDtab_mm <- master_tab(ls=list(mean_simG0_FixedINDUCED$mm,sd_simG0_FixedINDUCED$mm,mse_simG0_FixedINDUCED$mm,cvg_simG0_FixedINDUCED$mm))
simG0_FixedINDUCEDtab_jmm <- master_tab(ls=list(mean_simG0_FixedINDUCED$jmm,sd_simG0_FixedINDUCED$jmm,mse_simG0_FixedINDUCED$jmm,cvg_simG0_FixedINDUCED$jmm))

simG0_BothDecINDUCEDtab_gee <- master_tab(ls=list(mean_simG0_BothDecINDUCED$gee,sd_simG0_BothDecINDUCED$gee,mse_simG0_BothDecINDUCED$gee,cvg_simG0_BothDecINDUCED$gee))
simG0_BothDecINDUCEDtab_iee <- master_tab(ls=list(mean_simG0_BothDecINDUCED$iee,sd_simG0_BothDecINDUCED$iee,mse_simG0_BothDecINDUCED$iee,cvg_simG0_BothDecINDUCED$iee))
simG0_BothDecINDUCEDtab_wee <- master_tab(ls=list(mean_simG0_BothDecINDUCED$wee,sd_simG0_BothDecINDUCED$wee,mse_simG0_BothDecINDUCED$wee,cvg_simG0_BothDecINDUCED$wee))
simG0_BothDecINDUCEDtab_mm <- master_tab(ls=list(mean_simG0_BothDecINDUCED$mm,sd_simG0_BothDecINDUCED$mm,mse_simG0_BothDecINDUCED$mm,cvg_simG0_BothDecINDUCED$mm))
simG0_BothDecINDUCEDtab_jmm <- master_tab(ls=list(mean_simG0_BothDecINDUCED$jmm,sd_simG0_BothDecINDUCED$jmm,mse_simG0_BothDecINDUCED$jmm,cvg_simG0_BothDecINDUCED$jmm))

simG0_BothIncINDUCEDtab_gee <- master_tab(ls=list(mean_simG0_BothIncINDUCED$gee,sd_simG0_BothIncINDUCED$gee,mse_simG0_BothIncINDUCED$gee,cvg_simG0_BothIncINDUCED$gee))
simG0_BothIncINDUCEDtab_iee <- master_tab(ls=list(mean_simG0_BothIncINDUCED$iee,sd_simG0_BothIncINDUCED$iee,mse_simG0_BothIncINDUCED$iee,cvg_simG0_BothIncINDUCED$iee))
simG0_BothIncINDUCEDtab_wee <- master_tab(ls=list(mean_simG0_BothIncINDUCED$wee,sd_simG0_BothIncINDUCED$wee,mse_simG0_BothIncINDUCED$wee,cvg_simG0_BothIncINDUCED$wee))
simG0_BothIncINDUCEDtab_mm <- master_tab(ls=list(mean_simG0_BothIncINDUCED$mm,sd_simG0_BothIncINDUCED$mm,mse_simG0_BothIncINDUCED$mm,cvg_simG0_BothIncINDUCED$mm))
simG0_BothIncINDUCEDtab_jmm <- master_tab(ls=list(mean_simG0_BothIncINDUCED$jmm,sd_simG0_BothIncINDUCED$jmm,mse_simG0_BothIncINDUCED$jmm,cvg_simG0_BothIncINDUCED$jmm))



################################################
##              Print tables                  ##
################################################
##  PRINT ALL COMPONENTS
# kable(simG0_Fixedtab_gee,format="latex",booktabs=TRUE,linesep="",caption="simG0_Fixed -- GEE",align=rep(c("l","r"),c(2,3)))%>%collapse_rows(columns = 1,latex_hline="major")
# kable(simG0_Fixedtab_iee,format="latex",booktabs=TRUE,linesep="",caption="simG0_Fixed -- IEE",align=rep(c("l","r"),c(2,7)))%>%collapse_rows(columns = 1,latex_hline="major")
# kable(simG0_Fixedtab_wee,format="latex",booktabs=TRUE,linesep="",caption="simG0_Fixed -- WEE",align=rep(c("l","r"),c(2,7)))%>%collapse_rows(columns = 1,latex_hline="major")
# kable(simG0_Fixedtab_mm,format="latex",booktabs=TRUE,linesep="",caption="simG0_Fixed -- MM",align=rep(c("l","r"),c(2,4)))%>%collapse_rows(columns = 1,latex_hline="major")
# kable(simG0_Fixedtab_jmm,format="latex",booktabs=TRUE,linesep="",caption="simG0_Fixed -- JMM",align=rep(c("l","r"),c(2,4)))%>%collapse_rows(columns = 1,latex_hline="major")
#
# kable(simG0_BothDectab_gee,format="latex",booktabs=TRUE,linesep="",caption="simG0_BothDec -- GEE",align=rep(c("l","r"),c(2,3)))%>%collapse_rows(columns = 1,latex_hline="major")
# kable(simG0_BothDectab_iee,format="latex",booktabs=TRUE,linesep="",caption="simG0_BothDec -- IEE",align=rep(c("l","r"),c(2,7)))%>%collapse_rows(columns = 1,latex_hline="major")
# kable(simG0_BothDectab_wee,format="latex",booktabs=TRUE,linesep="",caption="simG0_BothDec -- WEE",align=rep(c("l","r"),c(2,7)))%>%collapse_rows(columns = 1,latex_hline="major")
# kable(simG0_BothDectab_mm,format="latex",booktabs=TRUE,linesep="",caption="simG0_BothDec -- MM",align=rep(c("l","r"),c(2,4)))%>%collapse_rows(columns = 1,latex_hline="major")
# kable(simG0_BothDectab_jmm,format="latex",booktabs=TRUE,linesep="",caption="simG0_BothDec -- JMM",align=rep(c("l","r"),c(2,4)))%>%collapse_rows(columns = 1,latex_hline="major")
#
# kable(simG0_BothInctab_gee,format="latex",booktabs=TRUE,linesep="",caption="simG0_BothInc -- GEE",align=rep(c("l","r"),c(2,3)))%>%collapse_rows(columns = 1,latex_hline="major")
# kable(simG0_BothInctab_iee,format="latex",booktabs=TRUE,linesep="",caption="simG0_BothInc -- IEE",align=rep(c("l","r"),c(2,7)))%>%collapse_rows(columns = 1,latex_hline="major")
# kable(simG0_BothInctab_wee,format="latex",booktabs=TRUE,linesep="",caption="simG0_BothInc -- WEE",align=rep(c("l","r"),c(2,7)))%>%collapse_rows(columns = 1,latex_hline="major")
# kable(simG0_BothInctab_mm,format="latex",booktabs=TRUE,linesep="",caption="simG0_BothInc -- MM",align=rep(c("l","r"),c(2,4)))%>%collapse_rows(columns = 1,latex_hline="major")
# kable(simG0_BothInctab_jmm,format="latex",booktabs=TRUE,linesep="",caption="simG0_BothInc -- JMM",align=rep(c("l","r"),c(2,4)))%>%collapse_rows(columns = 1,latex_hline="major")
#
# kable(simG0_SensDectab_gee,format="latex",booktabs=TRUE,linesep="",caption="simG0_SensDec -- GEE",align=rep(c("l","r"),c(2,3)))%>%collapse_rows(columns = 1,latex_hline="major")
# kable(simG0_SensDectab_iee,format="latex",booktabs=TRUE,linesep="",caption="simG0_SensDec -- IEE",align=rep(c("l","r"),c(2,7)))%>%collapse_rows(columns = 1,latex_hline="major")
# kable(simG0_SensDectab_wee,format="latex",booktabs=TRUE,linesep="",caption="simG0_SensDec -- WEE",align=rep(c("l","r"),c(2,7)))%>%collapse_rows(columns = 1,latex_hline="major")
# kable(simG0_SensDectab_mm,format="latex",booktabs=TRUE,linesep="",caption="simG0_SensDec -- MM",align=rep(c("l","r"),c(2,4)))%>%collapse_rows(columns = 1,latex_hline="major")
# kable(simG0_SensDectab_jmm,format="latex",booktabs=TRUE,linesep="",caption="simG0_SensDec -- JMM",align=rep(c("l","r"),c(2,4)))%>%collapse_rows(columns = 1,latex_hline="major")
#
#
# kable(simGN25_Fixedtab_gee,format="latex",booktabs=TRUE,linesep="",caption="simGN25_Fixed -- GEE",align=rep(c("l","r"),c(2,3)))%>%collapse_rows(columns = 1,latex_hline="major")
# kable(simGN25_Fixedtab_iee,format="latex",booktabs=TRUE,linesep="",caption="simGN25_Fixed -- IEE",align=rep(c("l","r"),c(2,7)))%>%collapse_rows(columns = 1,latex_hline="major")
# kable(simGN25_Fixedtab_wee,format="latex",booktabs=TRUE,linesep="",caption="simGN25_Fixed -- WEE",align=rep(c("l","r"),c(2,7)))%>%collapse_rows(columns = 1,latex_hline="major")
# kable(simGN25_Fixedtab_mm,format="latex",booktabs=TRUE,linesep="",caption="simGN25_Fixed -- MM",align=rep(c("l","r"),c(2,4)))%>%collapse_rows(columns = 1,latex_hline="major")
# kable(simGN25_Fixedtab_jmm,format="latex",booktabs=TRUE,linesep="",caption="simGN25_Fixed -- JMM",align=rep(c("l","r"),c(2,4)))%>%collapse_rows(columns = 1,latex_hline="major")
#
# kable(simGN25_BothDectab_gee,format="latex",booktabs=TRUE,linesep="",caption="simGN25_BothDec -- GEE",align=rep(c("l","r"),c(2,3)))%>%collapse_rows(columns = 1,latex_hline="major")
# kable(simGN25_BothDectab_iee,format="latex",booktabs=TRUE,linesep="",caption="simGN25_BothDec -- IEE",align=rep(c("l","r"),c(2,7)))%>%collapse_rows(columns = 1,latex_hline="major")
# kable(simGN25_BothDectab_wee,format="latex",booktabs=TRUE,linesep="",caption="simGN25_BothDec -- WEE",align=rep(c("l","r"),c(2,7)))%>%collapse_rows(columns = 1,latex_hline="major")
# kable(simGN25_BothDectab_mm,format="latex",booktabs=TRUE,linesep="",caption="simGN25_BothDec -- MM",align=rep(c("l","r"),c(2,4)))%>%collapse_rows(columns = 1,latex_hline="major")
# kable(simGN25_BothDectab_jmm,format="latex",booktabs=TRUE,linesep="",caption="simGN25_BothDec -- JMM",align=rep(c("l","r"),c(2,4)))%>%collapse_rows(columns = 1,latex_hline="major")
#
# kable(simGN25_BothInctab_gee,format="latex",booktabs=TRUE,linesep="",caption="simGN25_BothInc -- GEE",align=rep(c("l","r"),c(2,3)))%>%collapse_rows(columns = 1,latex_hline="major")
# kable(simGN25_BothInctab_iee,format="latex",booktabs=TRUE,linesep="",caption="simGN25_BothInc -- IEE",align=rep(c("l","r"),c(2,7)))%>%collapse_rows(columns = 1,latex_hline="major")
# kable(simGN25_BothInctab_wee,format="latex",booktabs=TRUE,linesep="",caption="simGN25_BothInc -- WEE",align=rep(c("l","r"),c(2,7)))%>%collapse_rows(columns = 1,latex_hline="major")
# kable(simGN25_BothInctab_mm,format="latex",booktabs=TRUE,linesep="",caption="simGN25_BothInc -- MM",align=rep(c("l","r"),c(2,4)))%>%collapse_rows(columns = 1,latex_hline="major")
# kable(simGN25_BothInctab_jmm,format="latex",booktabs=TRUE,linesep="",caption="simGN25_BothInc -- JMM",align=rep(c("l","r"),c(2,4)))%>%collapse_rows(columns = 1,latex_hline="major")
#
# kable(simGN25_SensDectab_gee,format="latex",booktabs=TRUE,linesep="",caption="simGN25_SensDec -- GEE",align=rep(c("l","r"),c(2,3)))%>%collapse_rows(columns = 1,latex_hline="major")
# kable(simGN25_SensDectab_iee,format="latex",booktabs=TRUE,linesep="",caption="simGN25_SensDec -- IEE",align=rep(c("l","r"),c(2,7)))%>%collapse_rows(columns = 1,latex_hline="major")
# kable(simGN25_SensDectab_wee,format="latex",booktabs=TRUE,linesep="",caption="simGN25_SensDec -- WEE",align=rep(c("l","r"),c(2,7)))%>%collapse_rows(columns = 1,latex_hline="major")
# kable(simGN25_SensDectab_mm,format="latex",booktabs=TRUE,linesep="",caption="simGN25_SensDec -- MM",align=rep(c("l","r"),c(2,4)))%>%collapse_rows(columns = 1,latex_hline="major")
# kable(simGN25_SensDectab_jmm,format="latex",booktabs=TRUE,linesep="",caption="simGN25_SensDec -- JMM",align=rep(c("l","r"),c(2,4)))%>%collapse_rows(columns = 1,latex_hline="major")
#
#
# kable(simG25_Fixedtab_gee,format="latex",booktabs=TRUE,linesep="",caption="simG25_Fixed -- GEE",align=rep(c("l","r"),c(2,3)))%>%collapse_rows(columns = 1,latex_hline="major")
# kable(simG25_Fixedtab_iee,format="latex",booktabs=TRUE,linesep="",caption="simG25_Fixed -- IEE",align=rep(c("l","r"),c(2,7)))%>%collapse_rows(columns = 1,latex_hline="major")
# kable(simG25_Fixedtab_wee,format="latex",booktabs=TRUE,linesep="",caption="simG25_Fixed -- WEE",align=rep(c("l","r"),c(2,7)))%>%collapse_rows(columns = 1,latex_hline="major")
# kable(simG25_Fixedtab_mm,format="latex",booktabs=TRUE,linesep="",caption="simG25_Fixed -- MM",align=rep(c("l","r"),c(2,4)))%>%collapse_rows(columns = 1,latex_hline="major")
# kable(simG25_Fixedtab_jmm,format="latex",booktabs=TRUE,linesep="",caption="simG25_Fixed -- JMM",align=rep(c("l","r"),c(2,4)))%>%collapse_rows(columns = 1,latex_hline="major")
#
# kable(simG25_BothDectab_gee,format="latex",booktabs=TRUE,linesep="",caption="simG25_BothDec -- GEE",align=rep(c("l","r"),c(2,3)))%>%collapse_rows(columns = 1,latex_hline="major")
# kable(simG25_BothDectab_iee,format="latex",booktabs=TRUE,linesep="",caption="simG25_BothDec -- IEE",align=rep(c("l","r"),c(2,7)))%>%collapse_rows(columns = 1,latex_hline="major")
# kable(simG25_BothDectab_wee,format="latex",booktabs=TRUE,linesep="",caption="simG25_BothDec -- WEE",align=rep(c("l","r"),c(2,7)))%>%collapse_rows(columns = 1,latex_hline="major")
# kable(simG25_BothDectab_mm,format="latex",booktabs=TRUE,linesep="",caption="simG25_BothDec -- MM",align=rep(c("l","r"),c(2,4)))%>%collapse_rows(columns = 1,latex_hline="major")
# kable(simG25_BothDectab_jmm,format="latex",booktabs=TRUE,linesep="",caption="simG25_BothDec -- JMM",align=rep(c("l","r"),c(2,4)))%>%collapse_rows(columns = 1,latex_hline="major")
#
# kable(simG25_BothInctab_gee,format="latex",booktabs=TRUE,linesep="",caption="simG25_BothInc -- GEE",align=rep(c("l","r"),c(2,3)))%>%collapse_rows(columns = 1,latex_hline="major")
# kable(simG25_BothInctab_iee,format="latex",booktabs=TRUE,linesep="",caption="simG25_BothInc -- IEE",align=rep(c("l","r"),c(2,7)))%>%collapse_rows(columns = 1,latex_hline="major")
# kable(simG25_BothInctab_wee,format="latex",booktabs=TRUE,linesep="",caption="simG25_BothInc -- WEE",align=rep(c("l","r"),c(2,7)))%>%collapse_rows(columns = 1,latex_hline="major")
# kable(simG25_BothInctab_mm,format="latex",booktabs=TRUE,linesep="",caption="simG25_BothInc -- MM",align=rep(c("l","r"),c(2,4)))%>%collapse_rows(columns = 1,latex_hline="major")
# kable(simG25_BothInctab_jmm,format="latex",booktabs=TRUE,linesep="",caption="simG25_BothInc -- JMM",align=rep(c("l","r"),c(2,4)))%>%collapse_rows(columns = 1,latex_hline="major")
#
# kable(simG25_SensDectab_gee,format="latex",booktabs=TRUE,linesep="",caption="simG25_SensDec -- GEE",align=rep(c("l","r"),c(2,3)))%>%collapse_rows(columns = 1,latex_hline="major")
# kable(simG25_SensDectab_iee,format="latex",booktabs=TRUE,linesep="",caption="simG25_SensDec -- IEE",align=rep(c("l","r"),c(2,7)))%>%collapse_rows(columns = 1,latex_hline="major")
# kable(simG25_SensDectab_wee,format="latex",booktabs=TRUE,linesep="",caption="simG25_SensDec -- WEE",align=rep(c("l","r"),c(2,7)))%>%collapse_rows(columns = 1,latex_hline="major")
# kable(simG25_SensDectab_mm,format="latex",booktabs=TRUE,linesep="",caption="simG25_SensDec -- MM",align=rep(c("l","r"),c(2,4)))%>%collapse_rows(columns = 1,latex_hline="major")
# kable(simG25_SensDectab_jmm,format="latex",booktabs=TRUE,linesep="",caption="simG25_SensDec -- JMM",align=rep(c("l","r"),c(2,4)))%>%collapse_rows(columns = 1,latex_hline="major")
#
## cond
kable(simGN25_BothDecCONDtab_mm,format="latex",booktabs=TRUE,linesep="",caption="simGN25_BothDecCOND -- MM",align=rep(c("l","r"),c(2,4)))%>%collapse_rows(columns = 1,latex_hline="major")
kable(simGN25_BothDecCONDtab_jmm,format="latex",booktabs=TRUE,linesep="",caption="simGN25_BothDecCOND -- JMM",align=rep(c("l","r"),c(2,4)))%>%collapse_rows(columns = 1,latex_hline="major")
# ## induced
# kable(simG0_FixedINDUCEDtab_gee,format="latex",booktabs=TRUE,linesep="",caption="simG0_FixedINDUCED -- GEE",align=rep(c("l","r"),c(2,3)))%>%collapse_rows(columns = 1,latex_hline="major")
# kable(simG0_FixedINDUCEDtab_iee,format="latex",booktabs=TRUE,linesep="",caption="simG0_FixedINDUCED -- IEE",align=rep(c("l","r"),c(2,7)))%>%collapse_rows(columns = 1,latex_hline="major")
# kable(simG0_FixedINDUCEDtab_wee,format="latex",booktabs=TRUE,linesep="",caption="simG0_FixedINDUCED -- WEE",align=rep(c("l","r"),c(2,7)))%>%collapse_rows(columns = 1,latex_hline="major")
# kable(simG0_FixedINDUCEDtab_mm,format="latex",booktabs=TRUE,linesep="",caption="simG0_FixedINDUCED -- MM",align=rep(c("l","r"),c(2,4)))%>%collapse_rows(columns = 1,latex_hline="major")
# kable(simG0_FixedINDUCEDtab_jmm,format="latex",booktabs=TRUE,linesep="",caption="simG0_FixedINDUCED -- JMM",align=rep(c("l","r"),c(2,4)))%>%collapse_rows(columns = 1,latex_hline="major")
#
# kable(simG0_BothDecINDUCEDtab_gee,format="latex",booktabs=TRUE,linesep="",caption="simG0_BothDecINDUCED -- GEE",align=rep(c("l","r"),c(2,3)))%>%collapse_rows(columns = 1,latex_hline="major")
# kable(simG0_BothDecINDUCEDtab_iee,format="latex",booktabs=TRUE,linesep="",caption="simG0_BothDecINDUCED -- IEE",align=rep(c("l","r"),c(2,7)))%>%collapse_rows(columns = 1,latex_hline="major")
# kable(simG0_BothDecINDUCEDtab_wee,format="latex",booktabs=TRUE,linesep="",caption="simG0_BothDecINDUCED -- WEE",align=rep(c("l","r"),c(2,7)))%>%collapse_rows(columns = 1,latex_hline="major")
# kable(simG0_BothDecINDUCEDtab_mm,format="latex",booktabs=TRUE,linesep="",caption="simG0_BothDecINDUCED -- MM",align=rep(c("l","r"),c(2,4)))%>%collapse_rows(columns = 1,latex_hline="major")
# kable(simG0_BothDecINDUCEDtab_jmm,format="latex",booktabs=TRUE,linesep="",caption="simG0_BothDecINDUCED -- JMM",align=rep(c("l","r"),c(2,4)))%>%collapse_rows(columns = 1,latex_hline="major")
#
# kable(simG0_BothIncINDUCEDtab_gee,format="latex",booktabs=TRUE,linesep="",caption="simG0_BothIncINDUCED -- GEE",align=rep(c("l","r"),c(2,3)))%>%collapse_rows(columns = 1,latex_hline="major")
# kable(simG0_BothIncINDUCEDtab_iee,format="latex",booktabs=TRUE,linesep="",caption="simG0_BothIncINDUCED -- IEE",align=rep(c("l","r"),c(2,7)))%>%collapse_rows(columns = 1,latex_hline="major")
# kable(simG0_BothIncINDUCEDtab_wee,format="latex",booktabs=TRUE,linesep="",caption="simG0_BothIncINDUCED -- WEE",align=rep(c("l","r"),c(2,7)))%>%collapse_rows(columns = 1,latex_hline="major")
# kable(simG0_BothIncINDUCEDtab_mm,format="latex",booktabs=TRUE,linesep="",caption="simG0_BothIncINDUCED -- MM",align=rep(c("l","r"),c(2,4)))%>%collapse_rows(columns = 1,latex_hline="major")
# kable(simG0_BothIncINDUCEDtab_jmm,format="latex",booktabs=TRUE,linesep="",caption="simG0_BothIncINDUCED -- JMM",align=rep(c("l","r"),c(2,4)))%>%collapse_rows(columns = 1,latex_hline="major")
#



##  BETA TABLES
beta_tab <- function(ls){
  B0 <- B1 <- B2 <- NULL
  for(ll in 1:length(ls)){
    B0 <- rbind(B0,cbind(scenario=ll,ls[[ll]][ls[[ll]]$param=="B0",-1]))
    B1 <- rbind(B1,cbind(scenario=ll,ls[[ll]][ls[[ll]]$param=="B1",-1]))
    B2 <- rbind(B2,cbind(scenario=ll,ls[[ll]][ls[[ll]]$param=="B2",-1]))
  }
  return(list(B0=B0,B1=B1,B2=B2))
}

beta_tab_gee <- beta_tab(list(simG0_Fixedtab_gee,simG0_BothDectab_gee,simG0_BothInctab_gee,simG0_SensDectab_gee,simGN25_Fixedtab_gee,simGN25_BothDectab_gee,simGN25_BothInctab_gee,simGN25_SensDectab_gee,simG25_Fixedtab_gee,simG25_BothDectab_gee,simG25_BothInctab_gee,simG25_SensDectab_gee))
beta_tab_iee <- beta_tab(list(simG0_Fixedtab_iee,simG0_BothDectab_iee,simG0_BothInctab_iee,simG0_SensDectab_iee,simGN25_Fixedtab_iee,simGN25_BothDectab_iee,simGN25_BothInctab_iee,simGN25_SensDectab_iee,simG25_Fixedtab_iee,simG25_BothDectab_iee,simG25_BothInctab_iee,simG25_SensDectab_iee))
beta_tab_wee <- beta_tab(list(simG0_Fixedtab_wee,simG0_BothDectab_wee,simG0_BothInctab_wee,simG0_SensDectab_wee,simGN25_Fixedtab_wee,simGN25_BothDectab_wee,simGN25_BothInctab_wee,simGN25_SensDectab_wee,simG25_Fixedtab_wee,simG25_BothDectab_wee,simG25_BothInctab_wee,simG25_SensDectab_wee))
beta_tab_mm <- beta_tab(list(simG0_Fixedtab_mm,simG0_BothDectab_mm,simG0_BothInctab_mm,simG0_SensDectab_mm,simGN25_Fixedtab_mm,simGN25_BothDectab_mm,simGN25_BothInctab_mm,simGN25_SensDectab_mm,simG25_Fixedtab_mm,simG25_BothDectab_mm,simG25_BothInctab_mm,simG25_SensDectab_mm))
beta_tab_jmm <- beta_tab(list(simG0_Fixedtab_jmm,simG0_BothDectab_jmm,simG0_BothInctab_jmm,simG0_SensDectab_jmm,simGN25_Fixedtab_jmm,simGN25_BothDectab_jmm,simGN25_BothInctab_jmm,simGN25_SensDectab_jmm,simG25_Fixedtab_jmm,simG25_BothDectab_jmm,simG25_BothInctab_jmm,simG25_SensDectab_jmm))

B0_tab_gee <- beta_tab_gee$B0; B1_tab_gee <- beta_tab_gee$B1; B2_tab_gee <- beta_tab_gee$B2
B0_tab_iee <- beta_tab_iee$B0; B1_tab_iee <- beta_tab_iee$B1; B2_tab_iee <- beta_tab_iee$B2
B0_tab_wee <- beta_tab_wee$B0; B1_tab_wee <- beta_tab_wee$B1; B2_tab_wee <- beta_tab_wee$B2
B0_tab_mm <- beta_tab_mm$B0; B1_tab_mm <- beta_tab_mm$B1; B2_tab_mm <- beta_tab_mm$B2
B0_tab_jmm <- beta_tab_jmm$B0; B1_tab_jmm <- beta_tab_jmm$B1; B2_tab_jmm <- beta_tab_jmm$B2

## Final tables
B0tab <- (cbind(B0_tab_jmm,B0_tab_wee[,-(1:2)]));rownames(B0tab) <- NULL
B1tab <- (cbind(B1_tab_jmm,B1_tab_wee[,-(1:2)]));rownames(B1tab) <- NULL
B2tab <- (cbind(B2_tab_jmm,B2_tab_wee[,-(1:2)]));rownames(B2tab) <- NULL

## Final
kable(B0tab,format="latex",booktabs=TRUE,linesep="",caption="B0 -- JMM/WEE",align=rep(c("l","r"),c(1,1+7+4)))%>%collapse_rows(columns = 1,latex_hline="major")
kable(B1tab,format="latex",booktabs=TRUE,linesep="",caption="B1 -- JMM/WEE",align=rep(c("l","r"),c(1,1+7+4)))%>%collapse_rows(columns = 1,latex_hline="major")
kable(B2tab,format="latex",booktabs=TRUE,linesep="",caption="B2 -- JMM/WEE",align=rep(c("l","r"),c(1,1+7+4)))%>%collapse_rows(columns = 1,latex_hline="major")


## IEE/MMM table
B0tab_noninf <- data.frame(cbind(B0_tab_mm,B0_tab_iee[,-(1:2)]));rownames(B0tab_noninf) <- NULL
B1tab_noninf <- data.frame(cbind(B1_tab_mm,B1_tab_iee[,-(1:2)]));rownames(B1tab_noninf) <- NULL
B2tab_noninf <- data.frame(cbind(B2_tab_mm,B2_tab_iee[,-(1:2)]));rownames(B2tab_noninf) <- NULL
kable(B0tab_noninf,format="latex",booktabs=TRUE,linesep="",caption="B0 -- MM/IEE",align=rep(c("l","r","l","r"),c(2,7,1,4)))%>%collapse_rows(columns = 1,latex_hline="major")
kable(B1tab_noninf,format="latex",booktabs=TRUE,linesep="",caption="B1 -- MM/IEE",align=rep(c("l","r","l","r"),c(2,7,1,4)))%>%collapse_rows(columns = 1,latex_hline="major")
kable(B2tab_noninf,format="latex",booktabs=TRUE,linesep="",caption="B2 -- MM/IEE",align=rep(c("l","r","l","r"),c(2,7,1,4)))%>%collapse_rows(columns = 1,latex_hline="major")

### Induced Informativeness
beta_inducedtab_iee <- beta_tab(list(simG0_Fixedtab_iee,simG0_BothDectab_iee,simG0_BothInctab_iee,simG0_FixedINDUCEDtab_iee,simG0_BothDecINDUCEDtab_iee,simG0_BothIncINDUCEDtab_iee))
beta_inducedtab_wee <- beta_tab(list(simG0_Fixedtab_wee,simG0_BothDectab_wee,simG0_BothInctab_wee,simG0_FixedINDUCEDtab_wee,simG0_BothDecINDUCEDtab_wee,simG0_BothIncINDUCEDtab_wee))
beta_inducedtab_mm <- beta_tab(list(simG0_Fixedtab_mm,simG0_BothDectab_mm,simG0_BothInctab_mm,simG0_FixedINDUCEDtab_mm,simG0_BothDecINDUCEDtab_mm,simG0_BothIncINDUCEDtab_mm))
beta_inducedtab_jmm <- beta_tab(list(simG0_Fixedtab_jmm,simG0_BothDectab_jmm,simG0_BothInctab_jmm,simG0_FixedINDUCEDtab_jmm,simG0_BothDecINDUCEDtab_jmm,simG0_BothIncINDUCEDtab_jmm))
B0_inducedtab_iee <- beta_inducedtab_iee$B0; B1_inducedtab_iee <- beta_inducedtab_iee$B1; B2_inducedtab_iee <- beta_inducedtab_iee$B2
B0_inducedtab_wee <- beta_inducedtab_wee$B0; B1_inducedtab_wee <- beta_inducedtab_wee$B1; B2_inducedtab_wee <- beta_inducedtab_wee$B2
B0_inducedtab_mm <- beta_inducedtab_mm$B0; B1_inducedtab_mm <- beta_inducedtab_mm$B1; B2_inducedtab_mm <- beta_inducedtab_mm$B2
B0_inducedtab_jmm <- beta_inducedtab_jmm$B0; B1_inducedtab_jmm <- beta_inducedtab_jmm$B1; B2_inducedtab_jmm <- beta_inducedtab_jmm$B2
B0tab_induced <- data.frame(cbind(B0_inducedtab_mm,B0_inducedtab_iee[,-(1:2)]));rownames(B0tab_induced) <- NULL
B1tab_induced <- data.frame(cbind(B1_inducedtab_mm,B1_inducedtab_iee[,-(1:2)]));rownames(B1tab_induced) <- NULL
B2tab_induced <- data.frame(cbind(B2_inducedtab_mm,B2_inducedtab_iee[,-(1:2)]));rownames(B2tab_induced) <- NULL
kable(B0tab_induced,format="latex",booktabs=TRUE,linesep="",caption="B0 -- IEE/MM",align=rep(c("l","r","l","r"),c(2,7,1,4)))%>%collapse_rows(columns = 1,latex_hline="major")
kable(B1tab_induced,format="latex",booktabs=TRUE,linesep="",caption="B1 -- IEE/MM",align=rep(c("l","r","l","r"),c(2,7,1,4)))%>%collapse_rows(columns = 1,latex_hline="major")
kable(B2tab_induced,format="latex",booktabs=TRUE,linesep="",caption="B2 -- IEE/MM",align=rep(c("l","r","l","r"),c(2,7,1,4)))%>%collapse_rows(columns = 1,latex_hline="major")
B0tab_induced2 <- data.frame(cbind(B0_inducedtab_jmm,B0_inducedtab_wee[,-(1:2)]));rownames(B0tab_induced2) <- NULL
B1tab_induced2 <- data.frame(cbind(B1_inducedtab_jmm,B1_inducedtab_wee[,-(1:2)]));rownames(B1tab_induced2) <- NULL
B2tab_induced2 <- data.frame(cbind(B2_inducedtab_jmm,B2_inducedtab_wee[,-(1:2)]));rownames(B2tab_induced2) <- NULL
kable(B0tab_induced2,format="latex",booktabs=TRUE,linesep="",caption="B0 -- WEE/JMM",align=rep(c("l","r","l","r"),c(2,7,1,4)))%>%collapse_rows(columns = 1,latex_hline="major")
kable(B1tab_induced2,format="latex",booktabs=TRUE,linesep="",caption="B1 -- WEE/JMM",align=rep(c("l","r","l","r"),c(2,7,1,4)))%>%collapse_rows(columns = 1,latex_hline="major")
kable(B2tab_induced2,format="latex",booktabs=TRUE,linesep="",caption="B2 -- WEE/JMM",align=rep(c("l","r","l","r"),c(2,7,1,4)))%>%collapse_rows(columns = 1,latex_hline="major")



################################################
##                  Boxplots                  ##
################################################
make_box <- function(df,method="jmm",caption=""){#,stat="est",func=mean,ZIP=FALSE,slopes=FALSE){

  if(method=="gee"){
    nams <- c("3Comp","1Mis","2Valid")
    nams2 <- c("Comp","Mis","Valid")
  }else if(method=="iee"|method=="wee"){
    nams <- c("7Comp","6EEE-N","5EEE","1Mis","4SEEE-N","3SEEE","2Valid")
    nams2 <- c("Comp","E4","E3","Mis","E2","E1","Valid")
  }else if(method=="mm"|method=="jmm"){
    nams <- c("4Comp","1Mis","3ObsLik","2Valid")
    nams2 <- c("Comp","Mis","ObsLik","Valid")
  }
  B0 <- df[,grep(pattern=paste0("B0est_",method),names(df))]
  colnames(B0) <- nams
  B0 <- B0%>%gather("Group","Est",1:ncol(B0))
  B0$color <- 4; B0$color[grepl("Mis", B0$Group)] <- 1; B0$color[grepl("Val", B0$Group)] <- 2; B0$color[grepl("Comp", B0$Group)] <- 3
  B0$opacity <- 0.95; #B0$opacity[B0$Group %in% c("4JMM","7Joint")] <- 1

  B1 <- df[,grep(pattern=paste0("B1est_",method),names(df))]
  colnames(B1) <- nams
  B1 <- B1%>%gather("Group","Est",1:ncol(B1))
  B1$color <- 4; B1$color[grepl("Mis", B1$Group)] <- 1; B1$color[grepl("Val", B1$Group)] <- 2; B1$color[grepl("Comp", B1$Group)] <- 3
  B1$opacity <- 0.95; #B1$opacity[B1$Group %in% c("4JMM","7Joint")] <- 1

  B2 <- df[,grep(pattern=paste0("B2est_",method),names(df))]
  colnames(B2) <- nams
  B2 <- B2%>%gather("Group","Est",1:ncol(B2))
  B2$color <- 4; B2$color[grepl("Mis", B2$Group)] <- 1; B2$color[grepl("Val", B2$Group)] <- 2; B2$color[grepl("Comp", B2$Group)] <- 3
  B2$opacity <- 0.95; #B2$opacity[B2$Group %in% c("4JMM","7Joint")] <- 1



  box_B0 <- ggplot(B0,aes(factor(Group),Est,fill=as.factor(color),alpha=as.numeric(opacity)))+
    geom_boxplot(aes(),outlier.size=0.5)+
    geom_hline(yintercept=b0)+
    scale_y_continuous(limits=c(b0-0.55,b0+0.55))+
    scale_fill_manual(values = c(wes_red,wes_gold,wes_green,wes_blue),guide=F)+
    scale_alpha_continuous(guide=F,range=c(0.5,1))+
    ylab("Estimates")+xlab("")+
    scale_x_discrete(breaks=nams,labels=nams2)+
    theme_classic() +
    ggtitle(caption) +
    theme(plot.title=element_text(hjust=0.99,margin=margin(t=5,b=-20),size=10))


  box_B1 <- ggplot(B1,aes(factor(Group),Est,fill=as.factor(color),alpha=as.numeric(opacity)))+
    geom_boxplot(aes(),outlier.size=0.5)+
    geom_hline(yintercept=b1)+
    scale_y_continuous(limits=c(b1-0.9,b1+0.99))+
    scale_fill_manual(values = c(wes_red,wes_gold,wes_green,wes_blue),guide=F)+
    scale_alpha_continuous(guide=F,range=c(0.5,1))+
    ylab("Estimates")+xlab("")+
    scale_x_discrete(breaks=nams,labels=nams2)+
    theme_classic()+
    ggtitle(caption) +
    theme(plot.title=element_text(hjust=0.99,margin=margin(t=5,b=-20),size=10))

  box_B2 <- ggplot(B2,aes(factor(Group),Est,fill=as.factor(color),alpha=as.numeric(opacity)))+
    geom_boxplot(aes(),outlier.size=0.5)+
    geom_hline(yintercept=b2)+
    scale_y_continuous(limits=c(b2-0.4,b2+0.4))+
    scale_fill_manual(values = c(wes_red,wes_gold,wes_green,wes_blue),guide=F)+
    scale_alpha_continuous(guide=F,range=c(0.5,1))+
    ylab("Estimates")+xlab("")+
    scale_x_discrete(breaks=nams,labels=nams2)+
    theme_classic()+
    ggtitle(caption) +
    theme(plot.title=element_text(hjust=0.99,margin=margin(t=5,b=-20),size=10))

  return(list(B0=box_B0,B1=box_B1,B2=box_B2))
}


box_gee_simG0_Fixed <- make_box(res_simG0_Fixed,method="gee",caption="Non-ICS \n Simple Mis")
box_iee_simG0_Fixed <- make_box(res_simG0_Fixed,method="iee",caption="Non-ICS \n Simple Mis")
box_wee_simG0_Fixed <- make_box(res_simG0_Fixed,method="wee",caption="Non-ICS \n Simple Mis")
box_mm_simG0_Fixed <- make_box(res_simG0_Fixed,method="mm",caption="Non-ICS \n Simple Mis")
box_jmm_simG0_Fixed <- make_box(res_simG0_Fixed,method="jmm",caption="Non-ICS \n Simple Mis")

box_gee_simG0_BothDec <- make_box(res_simG0_BothDec,method="gee",caption="Non-ICS \n Sens/Spec Dec")
box_iee_simG0_BothDec <- make_box(res_simG0_BothDec,method="iee",caption="Non-ICS \n Sens/Spec Dec")
box_wee_simG0_BothDec <- make_box(res_simG0_BothDec,method="wee",caption="Non-ICS \n Sens/Spec Dec")
box_mm_simG0_BothDec <- make_box(res_simG0_BothDec,method="mm",caption="Non-ICS \n Sens/Spec Dec")
box_jmm_simG0_BothDec <- make_box(res_simG0_BothDec,method="jmm",caption="Non-ICS \n Sens/Spec Dec")

box_gee_simG0_BothInc <- make_box(res_simG0_BothInc,method="gee",caption="Non-ICS \n Sens/Spec Inc")
box_iee_simG0_BothInc <- make_box(res_simG0_BothInc,method="iee",caption="Non-ICS \n Sens/Spec Inc")
box_wee_simG0_BothInc <- make_box(res_simG0_BothInc,method="wee",caption="Non-ICS \n Sens/Spec Inc")
box_mm_simG0_BothInc <- make_box(res_simG0_BothInc,method="mm",caption="Non-ICS \n Sens/Spec Inc")
box_jmm_simG0_BothInc <- make_box(res_simG0_BothInc,method="jmm",caption="Non-ICS \n Sens/Spec Inc")

box_gee_simG0_SensDec <- make_box(res_simG0_SensDec,method="gee",caption="Non-ICS \n Sens Dec")
box_iee_simG0_SensDec <- make_box(res_simG0_SensDec,method="iee",caption="Non-ICS \n Sens Dec")
box_wee_simG0_SensDec <- make_box(res_simG0_SensDec,method="wee",caption="Non-ICS \n Sens Dec")
box_mm_simG0_SensDec <- make_box(res_simG0_SensDec,method="mm",caption="Non-ICS \n Sens Dec")
box_jmm_simG0_SensDec <- make_box(res_simG0_SensDec,method="jmm",caption="Non-ICS \n Sens Dec")

box_gee_simGN25_Fixed <- make_box(res_simGN25_Fixed,method="gee",caption=expression(atop(paste("ICS (",gamma,"<0)"),"Simple Mis")))
box_iee_simGN25_Fixed <- make_box(res_simGN25_Fixed,method="iee",caption=expression(atop(paste("ICS (",gamma,"<0)"),"Simple Mis")))
box_wee_simGN25_Fixed <- make_box(res_simGN25_Fixed,method="wee",caption=expression(atop(paste("ICS (",gamma,"<0)"),"Simple Mis")))
box_mm_simGN25_Fixed <- make_box(res_simGN25_Fixed,method="mm",caption=expression(atop(paste("ICS (",gamma,"<0)"),"Simple Mis")))
box_jmm_simGN25_Fixed <- make_box(res_simGN25_Fixed,method="jmm",caption=expression(atop(paste("ICS (",gamma,"<0)"),"Simple Mis")))

box_gee_simGN25_BothDec <- make_box(res_simGN25_BothDec,method="gee",caption=expression(atop(paste("ICS (",gamma,"<0)"),"Sens/Spec Dec")))
box_iee_simGN25_BothDec <- make_box(res_simGN25_BothDec,method="iee",caption=expression(atop(paste("ICS (",gamma,"<0)"),"Sens/Spec Dec")))
box_wee_simGN25_BothDec <- make_box(res_simGN25_BothDec,method="wee",caption=expression(atop(paste("ICS (",gamma,"<0)"),"Sens/Spec Dec")))
box_mm_simGN25_BothDec <- make_box(res_simGN25_BothDec,method="mm",caption=expression(atop(paste("ICS (",gamma,"<0)"),"Sens/Spec Dec")))
box_jmm_simGN25_BothDec <- make_box(res_simGN25_BothDec,method="jmm",caption=expression(atop(paste("ICS (",gamma,"<0)"),"Sens/Spec Dec")))

box_gee_simGN25_BothInc <- make_box(res_simGN25_BothInc,method="gee",caption=expression(atop(paste("ICS (",gamma,"<0)"),"Sens/Spec Inc")))
box_iee_simGN25_BothInc <- make_box(res_simGN25_BothInc,method="iee",caption=expression(atop(paste("ICS (",gamma,"<0)"),"Sens/Spec Inc")))
box_wee_simGN25_BothInc <- make_box(res_simGN25_BothInc,method="wee",caption=expression(atop(paste("ICS (",gamma,"<0)"),"Sens/Spec Inc")))
box_mm_simGN25_BothInc <- make_box(res_simGN25_BothInc,method="mm",caption=expression(atop(paste("ICS (",gamma,"<0)"),"Sens/Spec Inc")))
box_jmm_simGN25_BothInc <- make_box(res_simGN25_BothInc,method="jmm",caption=expression(atop(paste("ICS (",gamma,"<0)"),"Sens/Spec Inc")))

box_gee_simGN25_SensDec <- make_box(res_simGN25_SensDec,method="gee",caption=expression(atop(paste("ICS (",gamma,"<0)"),"Sens Dec")))
box_iee_simGN25_SensDec <- make_box(res_simGN25_SensDec,method="iee",caption=expression(atop(paste("ICS (",gamma,"<0)"),"Sens Dec")))
box_wee_simGN25_SensDec <- make_box(res_simGN25_SensDec,method="wee",caption=expression(atop(paste("ICS (",gamma,"<0)"),"Sens Dec")))
box_mm_simGN25_SensDec <- make_box(res_simGN25_SensDec,method="mm",caption=expression(atop(paste("ICS (",gamma,"<0)"),"Sens Dec")))
box_jmm_simGN25_SensDec <- make_box(res_simGN25_SensDec,method="jmm",caption=expression(atop(paste("ICS (",gamma,"<0)"),"Sens Dec")))

box_gee_simG25_Fixed <- make_box(res_simG25_Fixed,method="gee",caption=expression(atop(paste("ICS (",gamma,">0)"),"Simple Mis")))
box_iee_simG25_Fixed <- make_box(res_simG25_Fixed,method="iee",caption=expression(atop(paste("ICS (",gamma,">0)"),"Simple Mis")))
box_wee_simG25_Fixed <- make_box(res_simG25_Fixed,method="wee",caption=expression(atop(paste("ICS (",gamma,">0)"),"Simple Mis")))
box_mm_simG25_Fixed <- make_box(res_simG25_Fixed,method="mm",caption=expression(atop(paste("ICS (",gamma,">0)"),"Simple Mis")))
box_jmm_simG25_Fixed <- make_box(res_simG25_Fixed,method="jmm",caption=expression(atop(paste("ICS (",gamma,">0)"),"Simple Mis")))

box_gee_simG25_BothDec <- make_box(res_simG25_BothDec,method="gee",caption=expression(atop(paste("ICS (",gamma,">0)"),"Sens/Spec Dec")))
box_iee_simG25_BothDec <- make_box(res_simG25_BothDec,method="iee",caption=expression(atop(paste("ICS (",gamma,">0)"),"Sens/Spec Dec")))
box_wee_simG25_BothDec <- make_box(res_simG25_BothDec,method="wee",caption=expression(atop(paste("ICS (",gamma,">0)"),"Sens/Spec Dec")))
box_mm_simG25_BothDec <- make_box(res_simG25_BothDec,method="mm",caption=expression(atop(paste("ICS (",gamma,">0)"),"Sens/Spec Dec")))
box_jmm_simG25_BothDec <- make_box(res_simG25_BothDec,method="jmm",caption=expression(atop(paste("ICS (",gamma,">0)"),"Sens/Spec Dec")))

box_gee_simG25_BothInc <- make_box(res_simG25_BothInc,method="gee",caption=expression(atop(paste("ICS (",gamma,">0)"),"Sens/Spec Inc")))
box_iee_simG25_BothInc <- make_box(res_simG25_BothInc,method="iee",caption=expression(atop(paste("ICS (",gamma,">0)"),"Sens/Spec Inc")))
box_wee_simG25_BothInc <- make_box(res_simG25_BothInc,method="wee",caption=expression(atop(paste("ICS (",gamma,">0)"),"Sens/Spec Inc")))
box_mm_simG25_BothInc <- make_box(res_simG25_BothInc,method="mm",caption=expression(atop(paste("ICS (",gamma,">0)"),"Sens/Spec Inc")))
box_jmm_simG25_BothInc <- make_box(res_simG25_BothInc,method="jmm",caption=expression(atop(paste("ICS (",gamma,">0)"),"Sens/Spec Inc")))

box_gee_simG25_SensDec <- make_box(res_simG25_SensDec,method="gee",caption=expression(atop(paste("ICS (",gamma,">0)"),"Sens Dec")))
box_iee_simG25_SensDec <- make_box(res_simG25_SensDec,method="iee",caption=expression(atop(paste("ICS (",gamma,">0)"),"Sens Dec")))
box_wee_simG25_SensDec <- make_box(res_simG25_SensDec,method="wee",caption=expression(atop(paste("ICS (",gamma,">0)"),"Sens Dec")))
box_mm_simG25_SensDec <- make_box(res_simG25_SensDec,method="mm",caption=expression(atop(paste("ICS (",gamma,">0)"),"Sens Dec")))
box_jmm_simG25_SensDec <- make_box(res_simG25_SensDec,method="jmm",caption=expression(atop(paste("ICS (",gamma,">0)"),"Sens Dec")))

##cond
box_mm_simGN25_BothDecCOND <- make_box(res_simGN25_BothDecCOND,method="mm",caption=expression(atop(paste("ICS (",gamma,"<0)"),"Sens/Spec Dec")))
box_jmm_simGN25_BothDecCOND <- make_box(res_simGN25_BothDecCOND,method="jmm",caption=expression(atop(paste("ICS (",gamma,"<0)"),"Sens/Spec Dec")))
## induced
box_gee_simG0_FixedINDUCED <- make_box(res_simG0_FixedINDUCED,method="gee",caption="Non-ICS \n Simple Mis")
box_iee_simG0_FixedINDUCED <- make_box(res_simG0_FixedINDUCED,method="iee",caption="Non-ICS \n Simple Mis")
box_wee_simG0_FixedINDUCED <- make_box(res_simG0_FixedINDUCED,method="wee",caption="Non-ICS \n Simple Mis")
box_mm_simG0_FixedINDUCED <- make_box(res_simG0_FixedINDUCED,method="mm",caption="Non-ICS \n Simple Mis")
box_jmm_simG0_FixedINDUCED <- make_box(res_simG0_FixedINDUCED,method="jmm",caption="Non-ICS \n Simple Mis")

box_gee_simG0_BothDecINDUCED <- make_box(res_simG0_BothDecINDUCED,method="gee",caption="Non-ICS \n Sens/Spec Dec")
box_iee_simG0_BothDecINDUCED <- make_box(res_simG0_BothDecINDUCED,method="iee",caption="Non-ICS \n Sens/Spec Dec")
box_wee_simG0_BothDecINDUCED <- make_box(res_simG0_BothDecINDUCED,method="wee",caption="Non-ICS \n Sens/Spec Dec")
box_mm_simG0_BothDecINDUCED <- make_box(res_simG0_BothDecINDUCED,method="mm",caption="Non-ICS \n Sens/Spec Dec")
box_jmm_simG0_BothDecINDUCED <- make_box(res_simG0_BothDecINDUCED,method="jmm",caption="Non-ICS \n Sens/Spec Dec")

box_gee_simG0_BothIncINDUCED <- make_box(res_simG0_BothIncINDUCED,method="gee",caption="Non-ICS \n Sens/Spec Inc")
box_iee_simG0_BothIncINDUCED <- make_box(res_simG0_BothIncINDUCED,method="iee",caption="Non-ICS \n Sens/Spec Inc")
box_wee_simG0_BothIncINDUCED <- make_box(res_simG0_BothIncINDUCED,method="wee",caption="Non-ICS \n Sens/Spec Inc")
box_mm_simG0_BothIncINDUCED <- make_box(res_simG0_BothIncINDUCED,method="mm",caption="Non-ICS \n Sens/Spec Inc")
box_jmm_simG0_BothIncINDUCED <- make_box(res_simG0_BothIncINDUCED,method="jmm",caption="Non-ICS \n Sens/Spec Inc")



ggsave(filename=paste0(outpath,"B0box_gee_simG0_Fixed.pdf"),plot=box_gee_simG0_Fixed$B0,width=4,height=4);ggsave(filename=paste0(outpath,"B1box_gee_simG0_Fixed.pdf"),plot=box_gee_simG0_Fixed$B1,width=4,height=4);ggsave(filename=paste0(outpath,"B2box_gee_simG0_Fixed.pdf"),plot=box_gee_simG0_Fixed$B2,width=4,height=4)
ggsave(filename=paste0(outpath,"B0box_iee_simG0_Fixed.pdf"),plot=box_iee_simG0_Fixed$B0,width=4,height=4);ggsave(filename=paste0(outpath,"B1box_iee_simG0_Fixed.pdf"),plot=box_iee_simG0_Fixed$B1,width=4,height=4);ggsave(filename=paste0(outpath,"B2box_iee_simG0_Fixed.pdf"),plot=box_iee_simG0_Fixed$B2,width=4,height=4)
ggsave(filename=paste0(outpath,"B0box_wee_simG0_Fixed.pdf"),plot=box_wee_simG0_Fixed$B0,width=4,height=4);ggsave(filename=paste0(outpath,"B1box_wee_simG0_Fixed.pdf"),plot=box_wee_simG0_Fixed$B1,width=4,height=4);ggsave(filename=paste0(outpath,"B2box_wee_simG0_Fixed.pdf"),plot=box_wee_simG0_Fixed$B2,width=4,height=4)
ggsave(filename=paste0(outpath,"B0box_mm_simG0_Fixed.pdf"),plot=box_mm_simG0_Fixed$B0,width=4,height=4);ggsave(filename=paste0(outpath,"B1box_mm_simG0_Fixed.pdf"),plot=box_mm_simG0_Fixed$B1,width=4,height=4);ggsave(filename=paste0(outpath,"B2box_mm_simG0_Fixed.pdf"),plot=box_mm_simG0_Fixed$B2,width=4,height=4)
ggsave(filename=paste0(outpath,"B0box_jmm_simG0_Fixed.pdf"),plot=box_jmm_simG0_Fixed$B0,width=4,height=4);ggsave(filename=paste0(outpath,"B1box_jmm_simG0_Fixed.pdf"),plot=box_jmm_simG0_Fixed$B1,width=4,height=4);ggsave(filename=paste0(outpath,"B2box_jmm_simG0_Fixed.pdf"),plot=box_jmm_simG0_Fixed$B2,width=4,height=4)

ggsave(filename=paste0(outpath,"B0box_gee_simG0_BothDec.pdf"),plot=box_gee_simG0_BothDec$B0,width=4,height=4);ggsave(filename=paste0(outpath,"B1box_gee_simG0_BothDec.pdf"),plot=box_gee_simG0_BothDec$B1,width=4,height=4);ggsave(filename=paste0(outpath,"B2box_gee_simG0_BothDec.pdf"),plot=box_gee_simG0_BothDec$B2,width=4,height=4)
ggsave(filename=paste0(outpath,"B0box_iee_simG0_BothDec.pdf"),plot=box_iee_simG0_BothDec$B0,width=4,height=4);ggsave(filename=paste0(outpath,"B1box_iee_simG0_BothDec.pdf"),plot=box_iee_simG0_BothDec$B1,width=4,height=4);ggsave(filename=paste0(outpath,"B2box_iee_simG0_BothDec.pdf"),plot=box_iee_simG0_BothDec$B2,width=4,height=4)
ggsave(filename=paste0(outpath,"B0box_wee_simG0_BothDec.pdf"),plot=box_wee_simG0_BothDec$B0,width=4,height=4);ggsave(filename=paste0(outpath,"B1box_wee_simG0_BothDec.pdf"),plot=box_wee_simG0_BothDec$B1,width=4,height=4);ggsave(filename=paste0(outpath,"B2box_wee_simG0_BothDec.pdf"),plot=box_wee_simG0_BothDec$B2,width=4,height=4)
ggsave(filename=paste0(outpath,"B0box_mm_simG0_BothDec.pdf"),plot=box_mm_simG0_BothDec$B0,width=4,height=4);ggsave(filename=paste0(outpath,"B1box_mm_simG0_BothDec.pdf"),plot=box_mm_simG0_BothDec$B1,width=4,height=4);ggsave(filename=paste0(outpath,"B2box_mm_simG0_BothDec.pdf"),plot=box_mm_simG0_BothDec$B2,width=4,height=4)
ggsave(filename=paste0(outpath,"B0box_jmm_simG0_BothDec.pdf"),plot=box_jmm_simG0_BothDec$B0,width=4,height=4);ggsave(filename=paste0(outpath,"B1box_jmm_simG0_BothDec.pdf"),plot=box_jmm_simG0_BothDec$B1,width=4,height=4);ggsave(filename=paste0(outpath,"B2box_jmm_simG0_BothDec.pdf"),plot=box_jmm_simG0_BothDec$B2,width=4,height=4)

ggsave(filename=paste0(outpath,"B0box_gee_simG0_BothInc.pdf"),plot=box_gee_simG0_BothInc$B0,width=4,height=4);ggsave(filename=paste0(outpath,"B1box_gee_simG0_BothInc.pdf"),plot=box_gee_simG0_BothInc$B1,width=4,height=4);ggsave(filename=paste0(outpath,"B2box_gee_simG0_BothInc.pdf"),plot=box_gee_simG0_BothInc$B2,width=4,height=4)
ggsave(filename=paste0(outpath,"B0box_iee_simG0_BothInc.pdf"),plot=box_iee_simG0_BothInc$B0,width=4,height=4);ggsave(filename=paste0(outpath,"B1box_iee_simG0_BothInc.pdf"),plot=box_iee_simG0_BothInc$B1,width=4,height=4);ggsave(filename=paste0(outpath,"B2box_iee_simG0_BothInc.pdf"),plot=box_iee_simG0_BothInc$B2,width=4,height=4)
ggsave(filename=paste0(outpath,"B0box_wee_simG0_BothInc.pdf"),plot=box_wee_simG0_BothInc$B0,width=4,height=4);ggsave(filename=paste0(outpath,"B1box_wee_simG0_BothInc.pdf"),plot=box_wee_simG0_BothInc$B1,width=4,height=4);ggsave(filename=paste0(outpath,"B2box_wee_simG0_BothInc.pdf"),plot=box_wee_simG0_BothInc$B2,width=4,height=4)
ggsave(filename=paste0(outpath,"B0box_mm_simG0_BothInc.pdf"),plot=box_mm_simG0_BothInc$B0,width=4,height=4);ggsave(filename=paste0(outpath,"B1box_mm_simG0_BothInc.pdf"),plot=box_mm_simG0_BothInc$B1,width=4,height=4);ggsave(filename=paste0(outpath,"B2box_mm_simG0_BothInc.pdf"),plot=box_mm_simG0_BothInc$B2,width=4,height=4)
ggsave(filename=paste0(outpath,"B0box_jmm_simG0_BothInc.pdf"),plot=box_jmm_simG0_BothInc$B0,width=4,height=4);ggsave(filename=paste0(outpath,"B1box_jmm_simG0_BothInc.pdf"),plot=box_jmm_simG0_BothInc$B1,width=4,height=4);ggsave(filename=paste0(outpath,"B2box_jmm_simG0_BothInc.pdf"),plot=box_jmm_simG0_BothInc$B2,width=4,height=4)

ggsave(filename=paste0(outpath,"B0box_gee_simG0_SensDec.pdf"),plot=box_gee_simG0_SensDec$B0,width=4,height=4);ggsave(filename=paste0(outpath,"B1box_gee_simG0_SensDec.pdf"),plot=box_gee_simG0_SensDec$B1,width=4,height=4);ggsave(filename=paste0(outpath,"B2box_gee_simG0_SensDec.pdf"),plot=box_gee_simG0_SensDec$B2,width=4,height=4)
ggsave(filename=paste0(outpath,"B0box_iee_simG0_SensDec.pdf"),plot=box_iee_simG0_SensDec$B0,width=4,height=4);ggsave(filename=paste0(outpath,"B1box_iee_simG0_SensDec.pdf"),plot=box_iee_simG0_SensDec$B1,width=4,height=4);ggsave(filename=paste0(outpath,"B2box_iee_simG0_SensDec.pdf"),plot=box_iee_simG0_SensDec$B2,width=4,height=4)
ggsave(filename=paste0(outpath,"B0box_wee_simG0_SensDec.pdf"),plot=box_wee_simG0_SensDec$B0,width=4,height=4);ggsave(filename=paste0(outpath,"B1box_wee_simG0_SensDec.pdf"),plot=box_wee_simG0_SensDec$B1,width=4,height=4);ggsave(filename=paste0(outpath,"B2box_wee_simG0_SensDec.pdf"),plot=box_wee_simG0_SensDec$B2,width=4,height=4)
ggsave(filename=paste0(outpath,"B0box_mm_simG0_SensDec.pdf"),plot=box_mm_simG0_SensDec$B0,width=4,height=4);ggsave(filename=paste0(outpath,"B1box_mm_simG0_SensDec.pdf"),plot=box_mm_simG0_SensDec$B1,width=4,height=4);ggsave(filename=paste0(outpath,"B2box_mm_simG0_SensDec.pdf"),plot=box_mm_simG0_SensDec$B2,width=4,height=4)
ggsave(filename=paste0(outpath,"B0box_jmm_simG0_SensDec.pdf"),plot=box_jmm_simG0_SensDec$B0,width=4,height=4);ggsave(filename=paste0(outpath,"B1box_jmm_simG0_SensDec.pdf"),plot=box_jmm_simG0_SensDec$B1,width=4,height=4);ggsave(filename=paste0(outpath,"B2box_jmm_simG0_SensDec.pdf"),plot=box_jmm_simG0_SensDec$B2,width=4,height=4)

ggsave(filename=paste0(outpath,"B0box_gee_simGN25_Fixed.pdf"),plot=box_gee_simGN25_Fixed$B0,width=4,height=4);ggsave(filename=paste0(outpath,"B1box_gee_simGN25_Fixed.pdf"),plot=box_gee_simGN25_Fixed$B1,width=4,height=4);ggsave(filename=paste0(outpath,"B2box_gee_simGN25_Fixed.pdf"),plot=box_gee_simGN25_Fixed$B2,width=4,height=4)
ggsave(filename=paste0(outpath,"B0box_iee_simGN25_Fixed.pdf"),plot=box_iee_simGN25_Fixed$B0,width=4,height=4);ggsave(filename=paste0(outpath,"B1box_iee_simGN25_Fixed.pdf"),plot=box_iee_simGN25_Fixed$B1,width=4,height=4);ggsave(filename=paste0(outpath,"B2box_iee_simGN25_Fixed.pdf"),plot=box_iee_simGN25_Fixed$B2,width=4,height=4)
ggsave(filename=paste0(outpath,"B0box_wee_simGN25_Fixed.pdf"),plot=box_wee_simGN25_Fixed$B0,width=4,height=4);ggsave(filename=paste0(outpath,"B1box_wee_simGN25_Fixed.pdf"),plot=box_wee_simGN25_Fixed$B1,width=4,height=4);ggsave(filename=paste0(outpath,"B2box_wee_simGN25_Fixed.pdf"),plot=box_wee_simGN25_Fixed$B2,width=4,height=4)
ggsave(filename=paste0(outpath,"B0box_mm_simGN25_Fixed.pdf"),plot=box_mm_simGN25_Fixed$B0,width=4,height=4);ggsave(filename=paste0(outpath,"B1box_mm_simGN25_Fixed.pdf"),plot=box_mm_simGN25_Fixed$B1,width=4,height=4);ggsave(filename=paste0(outpath,"B2box_mm_simGN25_Fixed.pdf"),plot=box_mm_simGN25_Fixed$B2,width=4,height=4)
ggsave(filename=paste0(outpath,"B0box_jmm_simGN25_Fixed.pdf"),plot=box_jmm_simGN25_Fixed$B0,width=4,height=4);ggsave(filename=paste0(outpath,"B1box_jmm_simGN25_Fixed.pdf"),plot=box_jmm_simGN25_Fixed$B1,width=4,height=4);ggsave(filename=paste0(outpath,"B2box_jmm_simGN25_Fixed.pdf"),plot=box_jmm_simGN25_Fixed$B2,width=4,height=4)

ggsave(filename=paste0(outpath,"B0box_gee_simGN25_BothDec.pdf"),plot=box_gee_simGN25_BothDec$B0,width=4,height=4);ggsave(filename=paste0(outpath,"B1box_gee_simGN25_BothDec.pdf"),plot=box_gee_simGN25_BothDec$B1,width=4,height=4);ggsave(filename=paste0(outpath,"B2box_gee_simGN25_BothDec.pdf"),plot=box_gee_simGN25_BothDec$B2,width=4,height=4)
ggsave(filename=paste0(outpath,"B0box_iee_simGN25_BothDec.pdf"),plot=box_iee_simGN25_BothDec$B0,width=4,height=4);ggsave(filename=paste0(outpath,"B1box_iee_simGN25_BothDec.pdf"),plot=box_iee_simGN25_BothDec$B1,width=4,height=4);ggsave(filename=paste0(outpath,"B2box_iee_simGN25_BothDec.pdf"),plot=box_iee_simGN25_BothDec$B2,width=4,height=4)
ggsave(filename=paste0(outpath,"B0box_wee_simGN25_BothDec.pdf"),plot=box_wee_simGN25_BothDec$B0,width=4,height=4);ggsave(filename=paste0(outpath,"B1box_wee_simGN25_BothDec.pdf"),plot=box_wee_simGN25_BothDec$B1,width=4,height=4);ggsave(filename=paste0(outpath,"B2box_wee_simGN25_BothDec.pdf"),plot=box_wee_simGN25_BothDec$B2,width=4,height=4)
ggsave(filename=paste0(outpath,"B0box_mm_simGN25_BothDec.pdf"),plot=box_mm_simGN25_BothDec$B0,width=4,height=4);ggsave(filename=paste0(outpath,"B1box_mm_simGN25_BothDec.pdf"),plot=box_mm_simGN25_BothDec$B1,width=4,height=4);ggsave(filename=paste0(outpath,"B2box_mm_simGN25_BothDec.pdf"),plot=box_mm_simGN25_BothDec$B2,width=4,height=4)
ggsave(filename=paste0(outpath,"B0box_jmm_simGN25_BothDec.pdf"),plot=box_jmm_simGN25_BothDec$B0,width=4,height=4);ggsave(filename=paste0(outpath,"B1box_jmm_simGN25_BothDec.pdf"),plot=box_jmm_simGN25_BothDec$B1,width=4,height=4);ggsave(filename=paste0(outpath,"B2box_jmm_simGN25_BothDec.pdf"),plot=box_jmm_simGN25_BothDec$B2,width=4,height=4)

ggsave(filename=paste0(outpath,"B0box_gee_simGN25_BothInc.pdf"),plot=box_gee_simGN25_BothInc$B0,width=4,height=4);ggsave(filename=paste0(outpath,"B1box_gee_simGN25_BothInc.pdf"),plot=box_gee_simGN25_BothInc$B1,width=4,height=4);ggsave(filename=paste0(outpath,"B2box_gee_simGN25_BothInc.pdf"),plot=box_gee_simGN25_BothInc$B2,width=4,height=4)
ggsave(filename=paste0(outpath,"B0box_iee_simGN25_BothInc.pdf"),plot=box_iee_simGN25_BothInc$B0,width=4,height=4);ggsave(filename=paste0(outpath,"B1box_iee_simGN25_BothInc.pdf"),plot=box_iee_simGN25_BothInc$B1,width=4,height=4);ggsave(filename=paste0(outpath,"B2box_iee_simGN25_BothInc.pdf"),plot=box_iee_simGN25_BothInc$B2,width=4,height=4)
ggsave(filename=paste0(outpath,"B0box_wee_simGN25_BothInc.pdf"),plot=box_wee_simGN25_BothInc$B0,width=4,height=4);ggsave(filename=paste0(outpath,"B1box_wee_simGN25_BothInc.pdf"),plot=box_wee_simGN25_BothInc$B1,width=4,height=4);ggsave(filename=paste0(outpath,"B2box_wee_simGN25_BothInc.pdf"),plot=box_wee_simGN25_BothInc$B2,width=4,height=4)
ggsave(filename=paste0(outpath,"B0box_mm_simGN25_BothInc.pdf"),plot=box_mm_simGN25_BothInc$B0,width=4,height=4);ggsave(filename=paste0(outpath,"B1box_mm_simGN25_BothInc.pdf"),plot=box_mm_simGN25_BothInc$B1,width=4,height=4);ggsave(filename=paste0(outpath,"B2box_mm_simGN25_BothInc.pdf"),plot=box_mm_simGN25_BothInc$B2,width=4,height=4)
ggsave(filename=paste0(outpath,"B0box_jmm_simGN25_BothInc.pdf"),plot=box_jmm_simGN25_BothInc$B0,width=4,height=4);ggsave(filename=paste0(outpath,"B1box_jmm_simGN25_BothInc.pdf"),plot=box_jmm_simGN25_BothInc$B1,width=4,height=4);ggsave(filename=paste0(outpath,"B2box_jmm_simGN25_BothInc.pdf"),plot=box_jmm_simGN25_BothInc$B2,width=4,height=4)

ggsave(filename=paste0(outpath,"B0box_gee_simGN25_SensDec.pdf"),plot=box_gee_simGN25_SensDec$B0,width=4,height=4);ggsave(filename=paste0(outpath,"B1box_gee_simGN25_SensDec.pdf"),plot=box_gee_simGN25_SensDec$B1,width=4,height=4);ggsave(filename=paste0(outpath,"B2box_gee_simGN25_SensDec.pdf"),plot=box_gee_simGN25_SensDec$B2,width=4,height=4)
ggsave(filename=paste0(outpath,"B0box_iee_simGN25_SensDec.pdf"),plot=box_iee_simGN25_SensDec$B0,width=4,height=4);ggsave(filename=paste0(outpath,"B1box_iee_simGN25_SensDec.pdf"),plot=box_iee_simGN25_SensDec$B1,width=4,height=4);ggsave(filename=paste0(outpath,"B2box_iee_simGN25_SensDec.pdf"),plot=box_iee_simGN25_SensDec$B2,width=4,height=4)
ggsave(filename=paste0(outpath,"B0box_wee_simGN25_SensDec.pdf"),plot=box_wee_simGN25_SensDec$B0,width=4,height=4);ggsave(filename=paste0(outpath,"B1box_wee_simGN25_SensDec.pdf"),plot=box_wee_simGN25_SensDec$B1,width=4,height=4);ggsave(filename=paste0(outpath,"B2box_wee_simGN25_SensDec.pdf"),plot=box_wee_simGN25_SensDec$B2,width=4,height=4)
ggsave(filename=paste0(outpath,"B0box_mm_simGN25_SensDec.pdf"),plot=box_mm_simGN25_SensDec$B0,width=4,height=4);ggsave(filename=paste0(outpath,"B1box_mm_simGN25_SensDec.pdf"),plot=box_mm_simGN25_SensDec$B1,width=4,height=4);ggsave(filename=paste0(outpath,"B2box_mm_simGN25_SensDec.pdf"),plot=box_mm_simGN25_SensDec$B2,width=4,height=4)
ggsave(filename=paste0(outpath,"B0box_jmm_simGN25_SensDec.pdf"),plot=box_jmm_simGN25_SensDec$B0,width=4,height=4);ggsave(filename=paste0(outpath,"B1box_jmm_simGN25_SensDec.pdf"),plot=box_jmm_simGN25_SensDec$B1,width=4,height=4);ggsave(filename=paste0(outpath,"B2box_jmm_simGN25_SensDec.pdf"),plot=box_jmm_simGN25_SensDec$B2,width=4,height=4)

ggsave(filename=paste0(outpath,"B0box_gee_simG25_Fixed.pdf"),plot=box_gee_simG25_Fixed$B0,width=4,height=4);ggsave(filename=paste0(outpath,"B1box_gee_simG25_Fixed.pdf"),plot=box_gee_simG25_Fixed$B1,width=4,height=4);ggsave(filename=paste0(outpath,"B2box_gee_simG25_Fixed.pdf"),plot=box_gee_simG25_Fixed$B2,width=4,height=4)
ggsave(filename=paste0(outpath,"B0box_iee_simG25_Fixed.pdf"),plot=box_iee_simG25_Fixed$B0,width=4,height=4);ggsave(filename=paste0(outpath,"B1box_iee_simG25_Fixed.pdf"),plot=box_iee_simG25_Fixed$B1,width=4,height=4);ggsave(filename=paste0(outpath,"B2box_iee_simG25_Fixed.pdf"),plot=box_iee_simG25_Fixed$B2,width=4,height=4)
ggsave(filename=paste0(outpath,"B0box_wee_simG25_Fixed.pdf"),plot=box_wee_simG25_Fixed$B0,width=4,height=4);ggsave(filename=paste0(outpath,"B1box_wee_simG25_Fixed.pdf"),plot=box_wee_simG25_Fixed$B1,width=4,height=4);ggsave(filename=paste0(outpath,"B2box_wee_simG25_Fixed.pdf"),plot=box_wee_simG25_Fixed$B2,width=4,height=4)
ggsave(filename=paste0(outpath,"B0box_mm_simG25_Fixed.pdf"),plot=box_mm_simG25_Fixed$B0,width=4,height=4);ggsave(filename=paste0(outpath,"B1box_mm_simG25_Fixed.pdf"),plot=box_mm_simG25_Fixed$B1,width=4,height=4);ggsave(filename=paste0(outpath,"B2box_mm_simG25_Fixed.pdf"),plot=box_mm_simG25_Fixed$B2,width=4,height=4)
ggsave(filename=paste0(outpath,"B0box_jmm_simG25_Fixed.pdf"),plot=box_jmm_simG25_Fixed$B0,width=4,height=4);ggsave(filename=paste0(outpath,"B1box_jmm_simG25_Fixed.pdf"),plot=box_jmm_simG25_Fixed$B1,width=4,height=4);ggsave(filename=paste0(outpath,"B2box_jmm_simG25_Fixed.pdf"),plot=box_jmm_simG25_Fixed$B2,width=4,height=4)

ggsave(filename=paste0(outpath,"B0box_gee_simG25_BothDec.pdf"),plot=box_gee_simG25_BothDec$B0,width=4,height=4);ggsave(filename=paste0(outpath,"B1box_gee_simG25_BothDec.pdf"),plot=box_gee_simG25_BothDec$B1,width=4,height=4);ggsave(filename=paste0(outpath,"B2box_gee_simG25_BothDec.pdf"),plot=box_gee_simG25_BothDec$B2,width=4,height=4)
ggsave(filename=paste0(outpath,"B0box_iee_simG25_BothDec.pdf"),plot=box_iee_simG25_BothDec$B0,width=4,height=4);ggsave(filename=paste0(outpath,"B1box_iee_simG25_BothDec.pdf"),plot=box_iee_simG25_BothDec$B1,width=4,height=4);ggsave(filename=paste0(outpath,"B2box_iee_simG25_BothDec.pdf"),plot=box_iee_simG25_BothDec$B2,width=4,height=4)
ggsave(filename=paste0(outpath,"B0box_wee_simG25_BothDec.pdf"),plot=box_wee_simG25_BothDec$B0,width=4,height=4);ggsave(filename=paste0(outpath,"B1box_wee_simG25_BothDec.pdf"),plot=box_wee_simG25_BothDec$B1,width=4,height=4);ggsave(filename=paste0(outpath,"B2box_wee_simG25_BothDec.pdf"),plot=box_wee_simG25_BothDec$B2,width=4,height=4)
ggsave(filename=paste0(outpath,"B0box_mm_simG25_BothDec.pdf"),plot=box_mm_simG25_BothDec$B0,width=4,height=4);ggsave(filename=paste0(outpath,"B1box_mm_simG25_BothDec.pdf"),plot=box_mm_simG25_BothDec$B1,width=4,height=4);ggsave(filename=paste0(outpath,"B2box_mm_simG25_BothDec.pdf"),plot=box_mm_simG25_BothDec$B2,width=4,height=4)
ggsave(filename=paste0(outpath,"B0box_jmm_simG25_BothDec.pdf"),plot=box_jmm_simG25_BothDec$B0,width=4,height=4);ggsave(filename=paste0(outpath,"B1box_jmm_simG25_BothDec.pdf"),plot=box_jmm_simG25_BothDec$B1,width=4,height=4);ggsave(filename=paste0(outpath,"B2box_jmm_simG25_BothDec.pdf"),plot=box_jmm_simG25_BothDec$B2,width=4,height=4)

ggsave(filename=paste0(outpath,"B0box_gee_simG25_BothInc.pdf"),plot=box_gee_simG25_BothInc$B0,width=4,height=4);ggsave(filename=paste0(outpath,"B1box_gee_simG25_BothInc.pdf"),plot=box_gee_simG25_BothInc$B1,width=4,height=4);ggsave(filename=paste0(outpath,"B2box_gee_simG25_BothInc.pdf"),plot=box_gee_simG25_BothInc$B2,width=4,height=4)
ggsave(filename=paste0(outpath,"B0box_iee_simG25_BothInc.pdf"),plot=box_iee_simG25_BothInc$B0,width=4,height=4);ggsave(filename=paste0(outpath,"B1box_iee_simG25_BothInc.pdf"),plot=box_iee_simG25_BothInc$B1,width=4,height=4);ggsave(filename=paste0(outpath,"B2box_iee_simG25_BothInc.pdf"),plot=box_iee_simG25_BothInc$B2,width=4,height=4)
ggsave(filename=paste0(outpath,"B0box_wee_simG25_BothInc.pdf"),plot=box_wee_simG25_BothInc$B0,width=4,height=4);ggsave(filename=paste0(outpath,"B1box_wee_simG25_BothInc.pdf"),plot=box_wee_simG25_BothInc$B1,width=4,height=4);ggsave(filename=paste0(outpath,"B2box_wee_simG25_BothInc.pdf"),plot=box_wee_simG25_BothInc$B2,width=4,height=4)
ggsave(filename=paste0(outpath,"B0box_mm_simG25_BothInc.pdf"),plot=box_mm_simG25_BothInc$B0,width=4,height=4);ggsave(filename=paste0(outpath,"B1box_mm_simG25_BothInc.pdf"),plot=box_mm_simG25_BothInc$B1,width=4,height=4);ggsave(filename=paste0(outpath,"B2box_mm_simG25_BothInc.pdf"),plot=box_mm_simG25_BothInc$B2,width=4,height=4)
ggsave(filename=paste0(outpath,"B0box_jmm_simG25_BothInc.pdf"),plot=box_jmm_simG25_BothInc$B0,width=4,height=4);ggsave(filename=paste0(outpath,"B1box_jmm_simG25_BothInc.pdf"),plot=box_jmm_simG25_BothInc$B1,width=4,height=4);ggsave(filename=paste0(outpath,"B2box_jmm_simG25_BothInc.pdf"),plot=box_jmm_simG25_BothInc$B2,width=4,height=4)

ggsave(filename=paste0(outpath,"B0box_gee_simG25_SensDec.pdf"),plot=box_gee_simG25_SensDec$B0,width=4,height=4);ggsave(filename=paste0(outpath,"B1box_gee_simG25_SensDec.pdf"),plot=box_gee_simG25_SensDec$B1,width=4,height=4);ggsave(filename=paste0(outpath,"B2box_gee_simG25_SensDec.pdf"),plot=box_gee_simG25_SensDec$B2,width=4,height=4)
ggsave(filename=paste0(outpath,"B0box_iee_simG25_SensDec.pdf"),plot=box_iee_simG25_SensDec$B0,width=4,height=4);ggsave(filename=paste0(outpath,"B1box_iee_simG25_SensDec.pdf"),plot=box_iee_simG25_SensDec$B1,width=4,height=4);ggsave(filename=paste0(outpath,"B2box_iee_simG25_SensDec.pdf"),plot=box_iee_simG25_SensDec$B2,width=4,height=4)
ggsave(filename=paste0(outpath,"B0box_wee_simG25_SensDec.pdf"),plot=box_wee_simG25_SensDec$B0,width=4,height=4);ggsave(filename=paste0(outpath,"B1box_wee_simG25_SensDec.pdf"),plot=box_wee_simG25_SensDec$B1,width=4,height=4);ggsave(filename=paste0(outpath,"B2box_wee_simG25_SensDec.pdf"),plot=box_wee_simG25_SensDec$B2,width=4,height=4)
ggsave(filename=paste0(outpath,"B0box_mm_simG25_SensDec.pdf"),plot=box_mm_simG25_SensDec$B0,width=4,height=4);ggsave(filename=paste0(outpath,"B1box_mm_simG25_SensDec.pdf"),plot=box_mm_simG25_SensDec$B1,width=4,height=4);ggsave(filename=paste0(outpath,"B2box_mm_simG25_SensDec.pdf"),plot=box_mm_simG25_SensDec$B2,width=4,height=4)
ggsave(filename=paste0(outpath,"B0box_jmm_simG25_SensDec.pdf"),plot=box_jmm_simG25_SensDec$B0,width=4,height=4);ggsave(filename=paste0(outpath,"B1box_jmm_simG25_SensDec.pdf"),plot=box_jmm_simG25_SensDec$B1,width=4,height=4);ggsave(filename=paste0(outpath,"B2box_jmm_simG25_SensDec.pdf"),plot=box_jmm_simG25_SensDec$B2,width=4,height=4)

## cond
ggsave(filename=paste0(outpath,"B0box_gee_simGN25_BothDecCOND.pdf"),plot=box_gee_simGN25_BothDecCOND$B0,width=4,height=4);ggsave(filename=paste0(outpath,"B1box_gee_simGN25_BothDecCOND.pdf"),plot=box_gee_simGN25_BothDecCOND$B1,width=4,height=4);ggsave(filename=paste0(outpath,"B2box_gee_simGN25_BothDecCOND.pdf"),plot=box_gee_simGN25_BothDecCOND$B2,width=4,height=4)
ggsave(filename=paste0(outpath,"B0box_iee_simGN25_BothDecCOND.pdf"),plot=box_iee_simGN25_BothDecCOND$B0,width=4,height=4);ggsave(filename=paste0(outpath,"B1box_iee_simGN25_BothDecCOND.pdf"),plot=box_iee_simGN25_BothDecCOND$B1,width=4,height=4);ggsave(filename=paste0(outpath,"B2box_iee_simGN25_BothDecCOND.pdf"),plot=box_iee_simGN25_BothDecCOND$B2,width=4,height=4)
ggsave(filename=paste0(outpath,"B0box_wee_simGN25_BothDecCOND.pdf"),plot=box_wee_simGN25_BothDecCOND$B0,width=4,height=4);ggsave(filename=paste0(outpath,"B1box_wee_simGN25_BothDecCOND.pdf"),plot=box_wee_simGN25_BothDecCOND$B1,width=4,height=4);ggsave(filename=paste0(outpath,"B2box_wee_simGN25_BothDecCOND.pdf"),plot=box_wee_simGN25_BothDecCOND$B2,width=4,height=4)
ggsave(filename=paste0(outpath,"B0box_mm_simGN25_BothDecCOND.pdf"),plot=box_mm_simGN25_BothDecCOND$B0,width=4,height=4);ggsave(filename=paste0(outpath,"B1box_mm_simGN25_BothDecCOND.pdf"),plot=box_mm_simGN25_BothDecCOND$B1,width=4,height=4);ggsave(filename=paste0(outpath,"B2box_mm_simGN25_BothDecCOND.pdf"),plot=box_mm_simGN25_BothDecCOND$B2,width=4,height=4)
ggsave(filename=paste0(outpath,"B0box_jmm_simGN25_BothDecCOND.pdf"),plot=box_jmm_simGN25_BothDecCOND$B0,width=4,height=4);ggsave(filename=paste0(outpath,"B1box_jmm_simGN25_BothDecCOND.pdf"),plot=box_jmm_simGN25_BothDecCOND$B1,width=4,height=4);ggsave(filename=paste0(outpath,"B2box_jmm_simGN25_BothDecCOND.pdf"),plot=box_jmm_simGN25_BothDecCOND$B2,width=4,height=4)
## induced
ggsave(filename=paste0(outpath,"B0box_gee_simG0_FixedINDUCED.pdf"),plot=box_gee_simG0_FixedINDUCED$B0,width=4,height=4);ggsave(filename=paste0(outpath,"B1box_gee_simG0_FixedINDUCED.pdf"),plot=box_gee_simG0_FixedINDUCED$B1,width=4,height=4);ggsave(filename=paste0(outpath,"B2box_gee_simG0_FixedINDUCED.pdf"),plot=box_gee_simG0_FixedINDUCED$B2,width=4,height=4)
ggsave(filename=paste0(outpath,"B0box_iee_simG0_FixedINDUCED.pdf"),plot=box_iee_simG0_FixedINDUCED$B0,width=4,height=4);ggsave(filename=paste0(outpath,"B1box_iee_simG0_FixedINDUCED.pdf"),plot=box_iee_simG0_FixedINDUCED$B1,width=4,height=4);ggsave(filename=paste0(outpath,"B2box_iee_simG0_FixedINDUCED.pdf"),plot=box_iee_simG0_FixedINDUCED$B2,width=4,height=4)
ggsave(filename=paste0(outpath,"B0box_wee_simG0_FixedINDUCED.pdf"),plot=box_wee_simG0_FixedINDUCED$B0,width=4,height=4);ggsave(filename=paste0(outpath,"B1box_wee_simG0_FixedINDUCED.pdf"),plot=box_wee_simG0_FixedINDUCED$B1,width=4,height=4);ggsave(filename=paste0(outpath,"B2box_wee_simG0_FixedINDUCED.pdf"),plot=box_wee_simG0_FixedINDUCED$B2,width=4,height=4)
ggsave(filename=paste0(outpath,"B0box_mm_simG0_FixedINDUCED.pdf"),plot=box_mm_simG0_FixedINDUCED$B0,width=4,height=4);ggsave(filename=paste0(outpath,"B1box_mm_simG0_FixedINDUCED.pdf"),plot=box_mm_simG0_FixedINDUCED$B1,width=4,height=4);ggsave(filename=paste0(outpath,"B2box_mm_simG0_FixedINDUCED.pdf"),plot=box_mm_simG0_FixedINDUCED$B2,width=4,height=4)
ggsave(filename=paste0(outpath,"B0box_jmm_simG0_FixedINDUCED.pdf"),plot=box_jmm_simG0_FixedINDUCED$B0,width=4,height=4);ggsave(filename=paste0(outpath,"B1box_jmm_simG0_FixedINDUCED.pdf"),plot=box_jmm_simG0_FixedINDUCED$B1,width=4,height=4);ggsave(filename=paste0(outpath,"B2box_jmm_simG0_FixedINDUCED.pdf"),plot=box_jmm_simG0_FixedINDUCED$B2,width=4,height=4)

ggsave(filename=paste0(outpath,"B0box_gee_simG0_BothDecINDUCED.pdf"),plot=box_gee_simG0_BothDecINDUCED$B0,width=4,height=4);ggsave(filename=paste0(outpath,"B1box_gee_simG0_BothDecINDUCED.pdf"),plot=box_gee_simG0_BothDecINDUCED$B1,width=4,height=4);ggsave(filename=paste0(outpath,"B2box_gee_simG0_BothDecINDUCED.pdf"),plot=box_gee_simG0_BothDecINDUCED$B2,width=4,height=4)
ggsave(filename=paste0(outpath,"B0box_iee_simG0_BothDecINDUCED.pdf"),plot=box_iee_simG0_BothDecINDUCED$B0,width=4,height=4);ggsave(filename=paste0(outpath,"B1box_iee_simG0_BothDecINDUCED.pdf"),plot=box_iee_simG0_BothDecINDUCED$B1,width=4,height=4);ggsave(filename=paste0(outpath,"B2box_iee_simG0_BothDecINDUCED.pdf"),plot=box_iee_simG0_BothDecINDUCED$B2,width=4,height=4)
ggsave(filename=paste0(outpath,"B0box_wee_simG0_BothDecINDUCED.pdf"),plot=box_wee_simG0_BothDecINDUCED$B0,width=4,height=4);ggsave(filename=paste0(outpath,"B1box_wee_simG0_BothDecINDUCED.pdf"),plot=box_wee_simG0_BothDecINDUCED$B1,width=4,height=4);ggsave(filename=paste0(outpath,"B2box_wee_simG0_BothDecINDUCED.pdf"),plot=box_wee_simG0_BothDecINDUCED$B2,width=4,height=4)
ggsave(filename=paste0(outpath,"B0box_mm_simG0_BothDecINDUCED.pdf"),plot=box_mm_simG0_BothDecINDUCED$B0,width=4,height=4);ggsave(filename=paste0(outpath,"B1box_mm_simG0_BothDecINDUCED.pdf"),plot=box_mm_simG0_BothDecINDUCED$B1,width=4,height=4);ggsave(filename=paste0(outpath,"B2box_mm_simG0_BothDecINDUCED.pdf"),plot=box_mm_simG0_BothDecINDUCED$B2,width=4,height=4)
ggsave(filename=paste0(outpath,"B0box_jmm_simG0_BothDecINDUCED.pdf"),plot=box_jmm_simG0_BothDecINDUCED$B0,width=4,height=4);ggsave(filename=paste0(outpath,"B1box_jmm_simG0_BothDecINDUCED.pdf"),plot=box_jmm_simG0_BothDecINDUCED$B1,width=4,height=4);ggsave(filename=paste0(outpath,"B2box_jmm_simG0_BothDecINDUCED.pdf"),plot=box_jmm_simG0_BothDecINDUCED$B2,width=4,height=4)

ggsave(filename=paste0(outpath,"B0box_gee_simG0_BothIncINDUCED.pdf"),plot=box_gee_simG0_BothIncINDUCED$B0,width=4,height=4);ggsave(filename=paste0(outpath,"B1box_gee_simG0_BothIncINDUCED.pdf"),plot=box_gee_simG0_BothIncINDUCED$B1,width=4,height=4);ggsave(filename=paste0(outpath,"B2box_gee_simG0_BothIncINDUCED.pdf"),plot=box_gee_simG0_BothIncINDUCED$B2,width=4,height=4)
ggsave(filename=paste0(outpath,"B0box_iee_simG0_BothIncINDUCED.pdf"),plot=box_iee_simG0_BothIncINDUCED$B0,width=4,height=4);ggsave(filename=paste0(outpath,"B1box_iee_simG0_BothIncINDUCED.pdf"),plot=box_iee_simG0_BothIncINDUCED$B1,width=4,height=4);ggsave(filename=paste0(outpath,"B2box_iee_simG0_BothIncINDUCED.pdf"),plot=box_iee_simG0_BothIncINDUCED$B2,width=4,height=4)
ggsave(filename=paste0(outpath,"B0box_wee_simG0_BothIncINDUCED.pdf"),plot=box_wee_simG0_BothIncINDUCED$B0,width=4,height=4);ggsave(filename=paste0(outpath,"B1box_wee_simG0_BothIncINDUCED.pdf"),plot=box_wee_simG0_BothIncINDUCED$B1,width=4,height=4);ggsave(filename=paste0(outpath,"B2box_wee_simG0_BothIncINDUCED.pdf"),plot=box_wee_simG0_BothIncINDUCED$B2,width=4,height=4)
ggsave(filename=paste0(outpath,"B0box_mm_simG0_BothIncINDUCED.pdf"),plot=box_mm_simG0_BothIncINDUCED$B0,width=4,height=4);ggsave(filename=paste0(outpath,"B1box_mm_simG0_BothIncINDUCED.pdf"),plot=box_mm_simG0_BothIncINDUCED$B1,width=4,height=4);ggsave(filename=paste0(outpath,"B2box_mm_simG0_BothIncINDUCED.pdf"),plot=box_mm_simG0_BothIncINDUCED$B2,width=4,height=4)
ggsave(filename=paste0(outpath,"B0box_jmm_simG0_BothIncINDUCED.pdf"),plot=box_jmm_simG0_BothIncINDUCED$B0,width=4,height=4);ggsave(filename=paste0(outpath,"B1box_jmm_simG0_BothIncINDUCED.pdf"),plot=box_jmm_simG0_BothIncINDUCED$B1,width=4,height=4);ggsave(filename=paste0(outpath,"B2box_jmm_simG0_BothIncINDUCED.pdf"),plot=box_jmm_simG0_BothIncINDUCED$B2,width=4,height=4)


## 1x6 boxplots
B0box_iee <- arrangeGrob(box_iee_simG0_Fixed$B0, box_iee_simGN25_Fixed$B0, box_iee_simG25_Fixed$B0, box_iee_simG0_BothDec$B0, box_iee_simGN25_BothDec$B0, box_iee_simG25_BothDec$B0,nrow=2,ncol=3)
B1box_iee <- arrangeGrob(box_iee_simG0_Fixed$B1, box_iee_simGN25_Fixed$B1, box_iee_simG25_Fixed$B1, box_iee_simG0_BothDec$B1, box_iee_simGN25_BothDec$B1, box_iee_simG25_BothDec$B1,nrow=2,ncol=3)
B2box_iee <- arrangeGrob(box_iee_simG0_Fixed$B2, box_iee_simGN25_Fixed$B2, box_iee_simG25_Fixed$B2, box_iee_simG0_BothDec$B2, box_iee_simGN25_BothDec$B2, box_iee_simG25_BothDec$B2,nrow=2,ncol=3)
ggsave(filename=paste0(outpath,"B0box_iee.pdf"),plot=B0box_iee,width=8,height=6)
ggsave(filename=paste0(outpath,"B1box_iee.pdf"),plot=B1box_iee,width=8,height=6)
ggsave(filename=paste0(outpath,"B2box_iee.pdf"),plot=B2box_iee,width=8,height=6)

B0box_wee <- arrangeGrob(box_wee_simG0_Fixed$B0, box_wee_simGN25_Fixed$B0, box_wee_simG25_Fixed$B0, box_wee_simG0_BothDec$B0, box_wee_simGN25_BothDec$B0, box_wee_simG25_BothDec$B0,nrow=2,ncol=3)
B1box_wee <- arrangeGrob(box_wee_simG0_Fixed$B1, box_wee_simGN25_Fixed$B1, box_wee_simG25_Fixed$B1, box_wee_simG0_BothDec$B1, box_wee_simGN25_BothDec$B1, box_wee_simG25_BothDec$B1,nrow=2,ncol=3)
B2box_wee <- arrangeGrob(box_wee_simG0_Fixed$B2, box_wee_simGN25_Fixed$B2, box_wee_simG25_Fixed$B2, box_wee_simG0_BothDec$B2, box_wee_simGN25_BothDec$B2, box_wee_simG25_BothDec$B2,nrow=2,ncol=3)
ggsave(filename=paste0(outpath,"B0box_wee.pdf"),plot=B0box_wee,width=8,height=6)
ggsave(filename=paste0(outpath,"B1box_wee.pdf"),plot=B1box_wee,width=8,height=6)
ggsave(filename=paste0(outpath,"B2box_wee.pdf"),plot=B2box_wee,width=8,height=6)

B0box_mm <- arrangeGrob(box_mm_simG0_Fixed$B0, box_mm_simGN25_Fixed$B0, box_mm_simG25_Fixed$B0, box_mm_simG0_BothDec$B0, box_mm_simGN25_BothDec$B0, box_mm_simG25_BothDec$B0,nrow=2,ncol=3)
B1box_mm <- arrangeGrob(box_mm_simG0_Fixed$B1, box_mm_simGN25_Fixed$B1, box_mm_simG25_Fixed$B1, box_mm_simG0_BothDec$B1, box_mm_simGN25_BothDec$B1, box_mm_simG25_BothDec$B1,nrow=2,ncol=3)
B2box_mm <- arrangeGrob(box_mm_simG0_Fixed$B2, box_mm_simGN25_Fixed$B2, box_mm_simG25_Fixed$B2, box_mm_simG0_BothDec$B2, box_mm_simGN25_BothDec$B2, box_mm_simG25_BothDec$B2,nrow=2,ncol=3)
ggsave(filename=paste0(outpath,"B0box_mm.pdf"),plot=B0box_mm,width=8,height=6)
ggsave(filename=paste0(outpath,"B1box_mm.pdf"),plot=B1box_mm,width=8,height=6)
ggsave(filename=paste0(outpath,"B2box_mm.pdf"),plot=B2box_mm,width=8,height=6)

B0box_jmm <- arrangeGrob(box_jmm_simG0_Fixed$B0, box_jmm_simGN25_Fixed$B0, box_jmm_simG25_Fixed$B0, box_jmm_simG0_BothDec$B0, box_jmm_simGN25_BothDec$B0, box_jmm_simG25_BothDec$B0,nrow=2,ncol=3)
B1box_jmm <- arrangeGrob(box_jmm_simG0_Fixed$B1, box_jmm_simGN25_Fixed$B1, box_jmm_simG25_Fixed$B1, box_jmm_simG0_BothDec$B1, box_jmm_simGN25_BothDec$B1, box_jmm_simG25_BothDec$B1,nrow=2,ncol=3)
B2box_jmm <- arrangeGrob(box_jmm_simG0_Fixed$B2, box_jmm_simGN25_Fixed$B2, box_jmm_simG25_Fixed$B2, box_jmm_simG0_BothDec$B2, box_jmm_simGN25_BothDec$B2, box_jmm_simG25_BothDec$B2,nrow=2,ncol=3)
ggsave(filename=paste0(outpath,"B0box_jmm.pdf"),plot=B0box_jmm,width=8,height=6)
ggsave(filename=paste0(outpath,"B1box_jmm.pdf"),plot=B1box_jmm,width=8,height=6)
ggsave(filename=paste0(outpath,"B2box_jmm.pdf"),plot=B2box_jmm,width=8,height=6)


## 3x3 plot
B0box_wee <- arrangeGrob(box_wee_simG0_Fixed$B0, box_wee_simGN25_Fixed$B0, box_wee_simG25_Fixed$B0, box_wee_simG0_BothDec$B0, box_wee_simGN25_BothDec$B0, box_wee_simG25_BothDec$B0, box_wee_simG0_BothInc$B0, box_wee_simGN25_BothInc$B0, box_wee_simG25_BothInc$B0,nrow=3,ncol=3)
B1box_wee <- arrangeGrob(box_wee_simG0_Fixed$B1, box_wee_simGN25_Fixed$B1, box_wee_simG25_Fixed$B1, box_wee_simG0_BothDec$B1, box_wee_simGN25_BothDec$B1, box_wee_simG25_BothDec$B1, box_wee_simG0_BothInc$B1, box_wee_simGN25_BothInc$B1, box_wee_simG25_BothInc$B1,nrow=3,ncol=3)
B2box_wee <- arrangeGrob(box_wee_simG0_Fixed$B2, box_wee_simGN25_Fixed$B2, box_wee_simG25_Fixed$B2, box_wee_simG0_BothDec$B2, box_wee_simGN25_BothDec$B2, box_wee_simG25_BothDec$B2, box_wee_simG0_BothInc$B2, box_wee_simGN25_BothInc$B2, box_wee_simG25_BothInc$B2,nrow=3,ncol=3)
ggsave(filename=paste0(outpath,"B0box_wee_3by3.pdf"),plot=B0box_wee,width=8,height=9)
ggsave(filename=paste0(outpath,"B1box_wee_3by3.pdf"),plot=B1box_wee,width=8,height=9)
ggsave(filename=paste0(outpath,"B2box_wee_3by3.pdf"),plot=B2box_wee,width=8,height=9)

## 3x3 plot
B0box_jmm <- arrangeGrob(box_jmm_simG0_Fixed$B0, box_jmm_simGN25_Fixed$B0, box_jmm_simG25_Fixed$B0, box_jmm_simG0_BothDec$B0, box_jmm_simGN25_BothDec$B0, box_jmm_simG25_BothDec$B0, box_jmm_simG0_BothInc$B0, box_jmm_simGN25_BothInc$B0, box_jmm_simG25_BothInc$B0,nrow=3,ncol=3)
B1box_jmm <- arrangeGrob(box_jmm_simG0_Fixed$B1, box_jmm_simGN25_Fixed$B1, box_jmm_simG25_Fixed$B1, box_jmm_simG0_BothDec$B1, box_jmm_simGN25_BothDec$B1, box_jmm_simG25_BothDec$B1, box_jmm_simG0_BothInc$B1, box_jmm_simGN25_BothInc$B1, box_jmm_simG25_BothInc$B1,nrow=3,ncol=3)
B2box_jmm <- arrangeGrob(box_jmm_simG0_Fixed$B2, box_jmm_simGN25_Fixed$B2, box_jmm_simG25_Fixed$B2, box_jmm_simG0_BothDec$B2, box_jmm_simGN25_BothDec$B2, box_jmm_simG25_BothDec$B2, box_jmm_simG0_BothInc$B2, box_jmm_simGN25_BothInc$B2, box_jmm_simG25_BothInc$B2,nrow=3,ncol=3)
ggsave(filename=paste0(outpath,"B0box_jmm_3by3.pdf"),plot=B0box_jmm,width=8,height=9)
ggsave(filename=paste0(outpath,"B1box_jmm_3by3.pdf"),plot=B1box_jmm,width=8,height=9)
ggsave(filename=paste0(outpath,"B2box_jmm_3by3.pdf"),plot=B2box_jmm,width=8,height=9)


bias_box <- function(df,caption=""){
  df <- cbind(df[,grepl("comp", names(df))],df[,grepl("mis", names(df))])
  nams <- c("6GEEComp","7IEEComp","99JMMComp","9MMComp","8WEEComp",
            "1GEEMis","2IEEMis","5JMMMis","4MMMis","3WEEMis")
  # nams2 <- c("GEE-Comp","IEE-Comp","JMM-Comp","MM-Comp","WEE-Comp",
  #           "GEE-Mis","IEE-Mis","JMM-Mis","MM-Mis","WEE-Mis")
  nams2 <- c("GEE","IEE","JMM","MM","WEE",
             "GEE","IEE","JMM","MM","WEE")

  B0 <- df[,grep(pattern=paste0("B0est_"),names(df))]
  colnames(B0) <- nams
  B0 <- B0%>%gather("Group","Est",1:ncol(B0))
  B0$color <- 2; B0$color[grepl("Mis", B0$Group)] <- 1
  B0$opacity <- 0.95; #B0$opacity[B0$Group %in% c("4JMM","7Joint")] <- 1

  B1 <- df[,grep(pattern=paste0("B1est_"),names(df))]
  colnames(B1) <- nams
  B1 <- B1%>%gather("Group","Est",1:ncol(B1))
  B1$color <- 2; B1$color[grepl("Mis", B1$Group)] <- 1
  B1$opacity <- 0.95; #B1$opacity[B1$Group %in% c("4JMM","7Joint")] <- 1

  B2 <- df[,grep(pattern=paste0("B2est_"),names(df))]
  colnames(B2) <- nams
  B2 <- B2%>%gather("Group","Est",1:ncol(B2))
  B2$color <- 2; B2$color[grepl("Mis", B2$Group)] <- 1
  B2$opacity <- 0.95; #B2$opacity[B2$Group %in% c("4JMM","7Joint")] <- 1



  box_B0 <- ggplot(B0,aes(factor(Group),Est,fill=as.factor(color),alpha=as.numeric(opacity)))+
    geom_boxplot(aes())+
    geom_hline(yintercept=b0)+
    scale_y_continuous(limits=c(b0-0.55,b0+0.55))+
    scale_fill_manual(values = c(wes_red,wes_green),guide=F)+
    scale_alpha_continuous(guide=F,range=c(0.5,1))+
    ylab("Estimates")+xlab("")+
    scale_x_discrete(breaks=nams,labels=nams2)+
    theme_classic() +
    ggtitle(paste(caption)) +
    theme(plot.title = element_text(hjust = 0.5,margin = margin(t = 10, b = -20)))


  box_B1 <- ggplot(B1,aes(factor(Group),Est,fill=as.factor(color),alpha=as.numeric(opacity)))+
    geom_boxplot(aes())+
    geom_hline(yintercept=b1)+
    scale_y_continuous(limits=c(b1-0.9,b1+0.9))+
    scale_fill_manual(values = c(wes_red,wes_green),guide=F)+
    scale_alpha_continuous(guide=F,range=c(0.5,1))+
    ylab("Estimates")+xlab("")+
    scale_x_discrete(breaks=nams,labels=nams2)+
    theme_classic()+
    ggtitle(paste(caption)) +
    theme(plot.title = element_text(hjust = 0.5,margin = margin(t = 10, b = -20)))

  box_B2 <- ggplot(B2,aes(factor(Group),Est,fill=as.factor(color),alpha=as.numeric(opacity)))+
    geom_boxplot(aes())+
    geom_hline(yintercept=b2)+
    scale_y_continuous(limits=c(b2-0.4,b2+0.4))+
    scale_fill_manual(values = c(wes_red,wes_green),guide=F)+
    scale_alpha_continuous(guide=F,range=c(0.5,1))+
    ylab("Estimates")+xlab("")+
    scale_x_discrete(breaks=nams,labels=nams2)+
    theme_classic()+
    ggtitle(paste(caption)) +
    theme(plot.title = element_text(hjust = 0.5,margin = margin(t = 10, b = -20)))

  return(list(B0=box_B0,B1=box_B1,B2=box_B2))
}


biasbox_simG0_Fixed <- bias_box(res_simG0_Fixed,caption="")
biasbox_simG0_BothDec <- bias_box(res_simG0_BothDec,caption="")
biasbox_simGN25_Fixed <- bias_box(res_simGN25_Fixed,caption="")
biasbox_simGN25_BothDec <- bias_box(res_simGN25_BothDec,caption="")
biasbox_simG25_Fixed <- bias_box(res_simG25_Fixed,caption="")
biasbox_simG25_BothDec <- bias_box(res_simG25_BothDec,caption="")








##########################
## RESULTS FOR ALPHA


##  BETA TABLES
alpha_tab <- function(ls){
  A0 <- A1 <- NULL
  for(ll in 1:length(ls)){
    A0 <- rbind(A0,cbind(scenario=ll,ls[[ll]][ls[[ll]]$param=="A0",-1]))
    A1 <- rbind(A1,cbind(scenario=ll,ls[[ll]][ls[[ll]]$param=="A1",-1]))
  }
  return(list(A0=A0,A1=A1))
}

alpha_tab_jmm <- alpha_tab(list(simG0_Fixedtab_jmm,simG0_BothDectab_jmm,simG0_BothInctab_jmm,simG0_SensDectab_jmm,simGN25_Fixedtab_jmm,simGN25_BothDectab_jmm,simGN25_BothInctab_jmm,simGN25_SensDectab_jmm,simG25_Fixedtab_jmm,simG25_BothDectab_jmm,simG25_BothInctab_jmm,simG25_SensDectab_jmm))

A0_tab_jmm <- alpha_tab_jmm$A0; A1_tab_jmm <- alpha_tab_jmm$A1; 

## Final tables
A0tab <- (cbind(A0_tab_jmm));rownames(A0tab) <- NULL
A1tab <- (cbind(A1_tab_jmm));rownames(A1tab) <- NULL


## Final
kable(A0tab,format="latex",booktabs=TRUE,linesep="",caption="A0 -- JMM",align=rep(c("l","r"),c(1,1+4)))%>%collapse_rows(columns = 1,latex_hline="major")
kable(A1tab,format="latex",booktabs=TRUE,linesep="",caption="A1 -- JMM",align=rep(c("l","r"),c(1,1+4)))%>%collapse_rows(columns = 1,latex_hline="major")



