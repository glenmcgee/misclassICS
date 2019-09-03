##############################
##        Simulation 1      ##
##############################
## Updated: 03/31/2019
## generating data:
## with random slopes (ICS and not)
## and MISCLASSIFICATION

## misclassification scenario (see manuscript)
scenario <-  "GN25_BothInc"     #"G0_Fixed" "G0_BothDec" "G0_BothInc" "G0_SensDec"
                            #"GN25_Fixed" "GN25_BothDec" "GN25_BothInc" "GN25_SensDec"
                            #"G25_Fixed" "G25_BothDec" "G25_BothInc" "G25_SensDec"

runLOCAL=TRUE ## set to TRUE to run locally ## FALSE is on cluster

#####################
## setting parameters
RR <- 2000 ## No. of iterations
KK <- 8000 ## No. of clusters
Vpercent <- 0.20 ## Validation sampling fraction
# VV <- 1000  ## No. of clusters in (internal) validation sample
JJ <- 500 ## Number of jobs in array (only matters on cluster)
Xprev <- 0.25 ## Exposure prevalence
a0 <- 0.6; a1 <- -0.2 ## volume model coefficients
b0 <- -4; b1 <- 0.5; b2 <- 0.2 ## outcome model coefficients
#gamma <- -0.1 ## scaling factor for random intercept in volume model
sigma0 <- 2 ## SD of rand eff for X=0
sigma1 <- 1.5 ## SD of rand eff for X=1
nb <- 3 ## number of betas
na <- 2 ## number of alphas
ng <- 2 ## number of gammas (2 for random slopes)
ns <- 2 ## number of sigmas (2 for random slopes)

# varnames <- c("B0est_naive","B1est_naive","B2est_naive","Sig0est_naive","Sig1est_naive","B0SE_naive","B1SE_naive","B2SE_naive","Sig0SE_naive","Sig1SE_naive","MSEP_naive",
#               # "A0est_jno0s","A1est_jno0s","G0est_jno0s","G1est_jno0s","B0est_jno0s","B1est_jno0s","B2est_jno0s","Sig0est_jno0s","Sig1est_jno0s","A0SE_jno0s","A1SE_jno0s","G0SE_jno0s","G1SE_jno0s","B0SE_jno0s","B1SE_jno0s","B2SE_jno0s","Sig0SE_jno0s","Sig1SE_jno0s","MSEP_jno0s",
#               "A0est_joint","A1est_joint","G0est_joint","G1est_joint","B0est_joint","B1est_joint","B2est_joint","Sig0est_joint","Sig1est_joint","A0SE_joint","A1SE_joint","G0SE_joint","G1SE_joint","B0SE_joint","B1SE_joint","B2SE_joint","Sig0SE_joint","Sig1SE_joint","MSEP_joint",
#               "B0est_gee","B1est_gee","B2est_gee","B0SE_gee","B1SE_gee","B2SE_gee",
#               "B0est_iee","B1est_iee","B2est_iee","B0SE_iee","B1SE_iee","B2SE_iee",
#               "B0est_wgee","B1est_wgee","B2est_wgee","B0SE_wgee","B1SE_wgee","B2SE_wgee",
#               "A0est_jmm","A1est_jmm","G0est_jmm","G1est_jmm","B0est_jmm","B1est_jmm","B2est_jmm","Sig0est_jmm","Sig1est_jmm","A0SE_jmm","A1SE_jmm","G0SE_jmm","G1SE_jmm","B0SE_jmm","B1SE_jmm","B2SE_jmm","Sig0SE_jmm","Sig1SE_jmm","MSEP_jmm")
# outputsize <- length(varnames) #67 ## width of output (ie number of variables reported) DO NOT CHANGE UNLESS MODELS/FITS ARE CHANGED
##
expit <- function(x){exp(x)/(1+exp(x))}
logit <- function(x){log(x/(1-x))}
# ## Misclassification parameters
# mis_varies <- c(100,logit(seq(0.90,0.70,by=-0.10))) ## varying sens/spec
# mis_100 <- rep(100,length(mis_varies)) ## constant/perfect classification
# NU01 <- rep(0,length(mis_varies)) ## Nk vs specificity
# NU11 <- rep(0,length(mis_varies)) ## Nk vs sensitivity
# ## Misclassification parameters
# NU00 <- logit(0.90)
# NU10 <- logit(0.90)
# NU00 <- logit(0.90)
# NU10 <- logit(0.90)



###############################################
## set folder to save data and load functions
if(runLOCAL==TRUE){
  path <- "../Simulations/Results/" ## path for results
  funpath <- "../Functions/" ## path for functions
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
# source(paste0(funpath,"fitEEE.R"))
require(lme4)
require(geepack)
require(JMMICSpack)

####################
## set up parallel
if(runLOCAL==TRUE){
  loopvec <- 1:RR ## run entire loop
  suffix <- "" ## output file doesnt have a suffix
}else{
  args <- commandArgs(trailingOnly = TRUE) ## collect arguments from shell script
  iter_no <- as.integer(args[1])  ## here there is only one argument, the iteration number
  loopvec <- ((RR/JJ)*(iter_no-1) +1):((RR/JJ)*iter_no) ## portion of the loop to run
  suffix <- paste0("_iter",iter_no) ## append output file names with the iteration number, to be combined later
}


#######################################
## Generate data with misclassification
gendat <- function(nu_00=logit(0.90), ## baseline specificities
                   nu_01=0, ## Nk vs specificities
                   nu_10=logit(0.90), ## baseline sensitivities
                   nu_11=0, ## Nk vs sensitivity
                   gam=0, ## scaling factor in volume model
                   m=1, ## minimum cluster size (0 or 1)
                   VVpercent=0.2,## iteration number used to set seed for reproducibility
                   iter=0){
  ## set seed
  if(iter!=0){ set.seed(1234+10*iter) }

  ## generate marginal data
  mdat <- genICS(gamma_bk=gam,gamma1_bk=gam,min_Nk=m,slopes=TRUE,cond=FALSE,  K=KK,Xrate=Xprev,alpha0=a0,alpha1=a1,beta0=b0,beta1=b1,beta2=b2,sig.bk=sigma0,sig1.bk=sigma1) ## marginal
  mdatMIS <- genMIS(mdat,nu00=nu_00,nu01=nu_01,nu10=nu_10,nu11=nu_11) ## marginal

  ## select validation sample
  VALIDk <- rep(c(1,0),c(VVpercent*KK,(1-VVpercent)*KK)); mdatMIS$data_clstlvl$VALIDk <- VALIDk
  VALIDki <- rep(VALIDk,mdatMIS$data_clstlvl$Nk); mdatMIS$data$VALIDki <- VALIDki

  return(mdatMIS)
}

#######################################
## fit complete data models
fit_complete <- function(datMIS,valid=FALSE,startpos=FALSE){

  ## collect data
  df <- datMIS$data; df2 <- datMIS$data_clstlvl
  if(valid==TRUE){ ## for validation sample only
    df <- df[df$VALIDki==1,]; df2 <- df2[df2$VALIDk==1,]
  }

  ## GEE-EXCH
  gee <- try(geeglm(Yki ~ X1ki + X2ki,id=IDki,data=df,family=binomial,corstr="exchangeable"),silent=T)
  if(length(gee)>1){
    est_gee <- gee$coef; SE_gee <- summary(gee)$coef[,2]
  }else{est_gee <- SE_gee <- rep(NA,nb)}
  if(length(SE_gee)<length(est_gee)){SE_gee <- rep(NA,length(est_gee))} ## error handling
  names(est_gee) <- paste0(c("B0","B1","B2"),"est_","gee",sep="")
  names(SE_gee) <- paste0(c("B0","B1","B2"),"SE_","gee",sep="")

  ## IEE
  iee <- try(geeglm(Yki ~ X1ki + X2ki,id=IDki,data=df,family=binomial,corstr="independence"),silent=T)
  if(length(iee)>1){
    est_iee <- iee$coef; SE_iee <- summary(iee)$coef[,2]
  }else{est_iee <- SE_iee <- rep(NA,nb)}
  if(length(SE_iee)<length(est_iee)){SE_iee <- rep(NA,length(est_iee))} ## error handling
  names(est_iee) <- paste0(c("B0","B1","B2"),"est_","iee",sep="")
  names(SE_iee) <- paste0(c("B0","B1","B2"),"SE_","iee",sep="")

  ## WGEE
  wee <- try(suppressWarnings(geeglm(Yki ~ X1ki + X2ki,id=IDki,data=df,family=binomial,corstr="independence",weights=ICWki)),silent=T) ## without suppress warnings, get warning about non integer successes in binomial, due to weights
  if(length(wee)>1){
    est_wee <- wee$coef; SE_wee <- summary(wee)$coef[,2]
  }else{est_wee <- SE_wee <- rep(NA,nb)}
  if(length(SE_wee)<length(est_wee)){SE_wee <- rep(NA,length(est_wee))} ## error handling
  names(est_wee) <- paste0(c("B0","B1","B2"),"est_","wee",sep="")
  names(SE_wee) <- paste0(c("B0","B1","B2"),"SE_","wee",sep="")

  ## MM
  mm <- try(JMMICS_fit(Nk=df2$Nk,Zk=cbind(1,df2$X1k),Yki=df$Yki,Xki=cbind(1,df$X1ki,df$X2ki),IDk=df2$IDk,Z0k=NULL,minNk=1,NegBin=FALSE,ZIP=FALSE,slopes=TRUE,slope_col=2,condSize=TRUE,condOut=FALSE,joint=FALSE,startpos=startpos,nquad=60),silent=T)
  if(length(mm)>1){
    est_mm <- mm$ests; SE_mm <- mm$SEs
  }else{est_mm <- SE_mm <- rep(NA,nb+ns)}
  if(length(SE_mm)<length(est_mm)){SE_mm <- rep(NA,length(est_mm))} ## error handling
  names(est_mm) <- paste0(c("B0","B1","B2","Sig0","Sig1"),"est_","mm",sep="")
  names(SE_mm) <- paste0(c("B0","B1","B2","Sig0","Sig1"),"SE_","mm",sep="")

  ## JMM
  jmm <- try(JMMICS_fit(Nk=df2$Nk,Zk=cbind(1,df2$X1k),Yki=df$Yki,Xki=cbind(1,df$X1ki,df$X2ki),IDk=df2$IDk,Z0k=NULL,minNk=1,NegBin=FALSE,ZIP=FALSE,slopes=TRUE,slope_col=2,condSize=TRUE,condOut=FALSE,joint=TRUE,startpos=startpos,nquad=60),silent=T)
  if(length(jmm)>1){
    est_jmm <- jmm$ests; SE_jmm <- jmm$SEs
  }else{est_jmm <- SE_jmm <- rep(NA,na+ng+nb+ns)}
  if(length(SE_jmm)<length(est_jmm)){SE_jmm <- rep(NA,length(est_jmm))} ## error handling
  names(est_jmm) <- paste0(c("A0","A1","Gam0","Gam1","B0","B1","B2","Sig0","Sig1"),"est_","jmm",sep="")
  names(SE_jmm) <- paste0(c("A0","A1","Gam0","Gam1","B0","B1","B2","Sig0","Sig1"),"SE_","jmm",sep="")

  ##
  ests <- c(est_gee,est_iee,est_wee,est_mm,est_jmm)
  SEs <- c(SE_gee,SE_iee,SE_wee,SE_mm,SE_jmm)
  if(valid==TRUE){
    names(ests) <- paste0(names(ests),"_valid",sep="")
    names(SEs) <- paste0(names(SEs),"_valid",sep="")
  }else{
    names(ests) <- paste0(names(ests),"_complete",sep="")
    names(SEs) <- paste0(names(SEs),"_complete",sep="")
  }

  # return(list(gee=c(est_gee,SE_gee),iee=c(est_iee,SE_iee),wee=c(est_wee,SE_wee),mm=c(est_mm,SE_mm),jmm=c(est_jmm,SE_jmm)))
  return(list(ests=ests,SEs=SEs))
}

#######################################
## fit misclassified models
fit_mis <- function(datMIS,startpos=FALSE){
  df <- datMIS$data; df2 <- datMIS$data_clstlvl

  ## GEE-EXCH
  gee <- try(geeglm(Yki ~ W1ki + X2ki,id=IDki,data=df,family=binomial,corstr="exchangeable"),silent=T)
  if(length(gee)>1){
    est_gee <- gee$coef; SE_gee <- summary(gee)$coef[,2]
  }else{est_gee <- SE_gee <- rep(NA,nb)}
  if(length(SE_gee)<length(est_gee)){SE_gee <- rep(NA,length(est_gee))} ## error handling
  names(est_gee) <- paste0(c("B0","B1","B2"),"est_","gee",sep="")
  names(SE_gee) <- paste0(c("B0","B1","B2"),"SE_","gee",sep="")

  ## IEE
  iee <- try(geeglm(Yki ~ W1ki + X2ki,id=IDki,data=df,family=binomial,corstr="independence"),silent=T)
  if(length(iee)>1){
    est_iee <- iee$coef; SE_iee <- summary(iee)$coef[,2]
  }else{est_iee <- SE_iee <- rep(NA,nb)}
  if(length(SE_iee)<length(est_iee)){SE_iee <- rep(NA,length(est_iee))} ## error handling
  names(est_iee) <- paste0(c("B0","B1","B2"),"est_","iee",sep="")
  names(SE_iee) <- paste0(c("B0","B1","B2"),"SE_","iee",sep="")

  ## WGEE
  wee <- try(suppressWarnings(geeglm(Yki ~ W1ki + X2ki,id=IDki,data=df,family=binomial,corstr="independence",weights=ICWki)),silent=T) ## without suppress warnings, get warning about non integer successes in binomial, due to weights
  if(length(wee)>1){
    est_wee <- wee$coef; SE_wee <- summary(wee)$coef[,2]
  }else{est_wee <- SE_wee <- rep(NA,nb)}
  if(length(SE_wee)<length(est_wee)){SE_wee <- rep(NA,length(est_wee))} ## error handling
  names(est_wee) <- paste0(c("B0","B1","B2"),"est_","wee",sep="")
  names(SE_wee) <- paste0(c("B0","B1","B2"),"SE_","wee",sep="")

  ## MM
  mm <- try(JMMICS_fit(Nk=df2$Nk,Zk=cbind(1,df2$W1k),Yki=df$Yki,Xki=cbind(1,df$W1ki,df$X2ki),IDk=df2$IDk,Z0k=NULL,minNk=1,NegBin=FALSE,ZIP=FALSE,slopes=TRUE,slope_col=2,condSize=TRUE,condOut=FALSE,joint=FALSE,startpos=startpos,nquad=60),silent=T)
  if(length(mm)>1){
    est_mm <- mm$ests; SE_mm <- mm$SEs
  }else{est_mm <- SE_mm <- rep(NA,nb+ns)}
  if(length(SE_mm)<length(est_mm)){SE_mm <- rep(NA,length(est_mm))} ## error handling
  names(est_mm) <- paste0(c("B0","B1","B2","Sig0","Sig1"),"est_","mm",sep="")
  names(SE_mm) <- paste0(c("B0","B1","B2","Sig0","Sig1"),"SE_","mm",sep="")

  ## JMM
  jmm <- try(JMMICS_fit(Nk=df2$Nk,Zk=cbind(1,df2$W1k),Yki=df$Yki,Xki=cbind(1,df$W1ki,df$X2ki),IDk=df2$IDk,Z0k=NULL,minNk=1,NegBin=FALSE,ZIP=FALSE,slopes=TRUE,slope_col=2,condSize=TRUE,condOut=FALSE,joint=TRUE,startpos=startpos,nquad=60),silent=T)
  if(length(jmm)>1){
    est_jmm <- jmm$ests; SE_jmm <- jmm$SEs
  }else{est_jmm <- SE_jmm <- rep(NA,na+ng+nb+ns)}
  if(length(SE_jmm)<length(est_jmm)){SE_jmm <- rep(NA,length(est_jmm))} ## error handling
  names(est_jmm) <- paste0(c("A0","A1","Gam0","Gam1","B0","B1","B2","Sig0","Sig1"),"est_","jmm",sep="")
  names(SE_jmm) <- paste0(c("A0","A1","Gam0","Gam1","B0","B1","B2","Sig0","Sig1"),"SE_","jmm",sep="")

  ##
  ests <- c(est_gee,est_iee,est_wee,est_mm,est_jmm)
  SEs <- c(SE_gee,SE_iee,SE_wee,SE_mm,SE_jmm)
  names(ests) <- paste0(names(ests),"_mis",sep="")
  names(SEs) <- paste0(names(SEs),"_mis",sep="")


    return(list(ests=ests,SEs=SEs))
}

################################################
## fit expected estimating equations approaches
fit_EEE <- function(datMIS,dependN=TRUE){
  ## collect data
  df <- datMIS$data; df2 <- datMIS$data_clstlvl
  ## placeholder exposure in main sample to ensure we're not using true exposure
  df$X1ki[df$VALIDki==0] <- 0
  df2$X1k[df2$VALIDk==0] <- 0

  if(dependN==FALSE){ ## exposure model not including N
    ## IEE -- simplified
    iee_simple <- try(fit_naive2stage_R(Xki=cbind(1,df$X1ki,df$X2ki),Zxki=model.matrix(X1ki~Yki+W1ki+X2ki,data=df),Nk_actual=df2$Nk,Yki=df$Yki,Xk=df2$X1k,IDk=df2$IDk,WEIGHTk=NA,VALIDk=df2$VALIDk,nboot=50),silent=T)
    if(length(iee_simple)>1){
      est_iee_simple <- iee_simple$ests; SE_iee_simple <- iee_simple$SEs
    }else{est_iee_simple <- SE_iee_simple <- rep(NA,nb)}
    if(length(SE_iee_simple)<length(est_iee_simple)){SE_iee_simple <- rep(NA,length(est_iee_simple))} ## error handling
    names(est_iee_simple) <- paste0(c("B0","B1","B2"),"est_","iee_simple",sep="")
    names(SE_iee_simple) <- paste0(c("B0","B1","B2"),"SE_","iee_simple",sep="")

    ## WEE -- simplified
    wee_simple <- try(fit_naive2stage_R(Xki=cbind(1,df$X1ki,df$X2ki),Zxki=model.matrix(X1ki~Yki+W1ki+X2ki,data=df),Nk_actual=df2$Nk,Yki=df$Yki,Xk=df2$X1k,IDk=df2$IDk,WEIGHTk=(1/df2$Nk),VALIDk=df2$VALIDk,nboot=50),silent=T)
    if(length(wee_simple)>1){
      est_wee_simple <- wee_simple$ests; SE_wee_simple <- wee_simple$SEs
    }else{est_wee_simple <- SE_wee_simple <- rep(NA,nb)}
    if(length(SE_wee_simple)<length(est_wee_simple)){SE_wee_simple <- rep(NA,length(est_wee_simple))} ## error handling
    names(est_wee_simple) <- paste0(c("B0","B1","B2"),"est_","wee_simple",sep="")
    names(SE_wee_simple) <- paste0(c("B0","B1","B2"),"SE_","wee_simple",sep="")

    ## IEE
    iee <- try(fit_naive2stage_R(Xki=cbind(1,df$X1ki,df$X2ki),Zxki=model.matrix(X1ki~Yki*W1ki+X2ki,data=df),Nk_actual=df2$Nk,Yki=df$Yki,Xk=df2$X1k,IDk=df2$IDk,WEIGHTk=NA,VALIDk=df2$VALIDk,nboot=50),silent=T)
    if(length(iee)>1){
      est_iee <- iee$ests; SE_iee <- iee$SEs
    }else{est_iee <- SE_iee <- rep(NA,nb)}
    if(length(SE_iee)<length(est_iee)){SE_iee <- rep(NA,length(est_iee))} ## error handling
    names(est_iee) <- paste0(c("B0","B1","B2"),"est_","iee",sep="")
    names(SE_iee) <- paste0(c("B0","B1","B2"),"SE_","iee",sep="")

    ## WEE
    wee <- try(fit_naive2stage_R(Xki=cbind(1,df$X1ki,df$X2ki),Zxki=model.matrix(X1ki~Yki*W1ki+X2ki,data=df),Nk_actual=df2$Nk,Yki=df$Yki,Xk=df2$X1k,IDk=df2$IDk,WEIGHTk=(1/df2$Nk),VALIDk=df2$VALIDk,nboot=50),silent=T)
    if(length(wee)>1){
      est_wee <- wee$ests; SE_wee <- wee$SEs
    }else{est_wee <- SE_wee <- rep(NA,nb)}
    if(length(SE_wee)<length(est_wee)){SE_wee <- rep(NA,length(est_wee))} ## error handling
    names(est_wee) <- paste0(c("B0","B1","B2"),"est_","wee",sep="")
    names(SE_wee) <- paste0(c("B0","B1","B2"),"SE_","wee",sep="")

    ##
    ests <- c(est_iee_simple,est_wee_simple,est_iee,est_wee)
    SEs <- c(SE_iee_simple,SE_wee_simple,SE_iee,SE_wee)
    names(ests) <- paste0(names(ests),"_EEEnoN",sep="")
    names(SEs) <- paste0(names(SEs),"_EEEnoN",sep="")

  }else{ ## exposure model includes N
    ## IEE -- simplified
    iee_simple <- try(fit_naive2stage_R(Xki=cbind(1,df$X1ki,df$X2ki),Zxki=model.matrix(X1ki~Yki+W1ki+Nki+X2ki,data=df),Nk_actual=df2$Nk,Yki=df$Yki,Xk=df2$X1k,IDk=df2$IDk,WEIGHTk=NA,VALIDk=df2$VALIDk,nboot=50),silent=T)
    if(length(iee_simple)>1){
      est_iee_simple <- iee_simple$ests; SE_iee_simple <- iee_simple$SEs
    }else{est_iee_simple <- SE_iee_simple <- rep(NA,nb)}
    if(length(SE_iee_simple)<length(est_iee_simple)){SE_iee_simple <- rep(NA,length(est_iee_simple))} ## error handling
    names(est_iee_simple) <- paste0(c("B0","B1","B2"),"est_","iee_simple",sep="")
    names(SE_iee_simple) <- paste0(c("B0","B1","B2"),"SE_","iee_simple",sep="")

    ## WEE -- simplified
    wee_simple <- try(fit_naive2stage_R(Xki=cbind(1,df$X1ki,df$X2ki),Zxki=model.matrix(X1ki~Yki+W1ki+Nki+X2ki,data=df),Nk_actual=df2$Nk,Yki=df$Yki,Xk=df2$X1k,IDk=df2$IDk,WEIGHTk=(1/df2$Nk),VALIDk=df2$VALIDk,nboot=50),silent=T)
    if(length(wee_simple)>1){
      est_wee_simple <- wee_simple$ests; SE_wee_simple <- wee_simple$SEs
    }else{est_wee_simple <- SE_wee_simple <- rep(NA,nb)}
    if(length(SE_wee_simple)<length(est_wee_simple)){SE_wee_simple <- rep(NA,length(est_wee_simple))} ## error handling
    names(est_wee_simple) <- paste0(c("B0","B1","B2"),"est_","wee_simple",sep="")
    names(SE_wee_simple) <- paste0(c("B0","B1","B2"),"SE_","wee_simple",sep="")

    ## IEE
    iee <- try(fit_naive2stage_R(Xki=cbind(1,df$X1ki,df$X2ki),Zxki=model.matrix(X1ki~Yki*W1ki*Nki+X2ki,data=df),Nk_actual=df2$Nk,Yki=df$Yki,Xk=df2$X1k,IDk=df2$IDk,WEIGHTk=NA,VALIDk=df2$VALIDk,nboot=50),silent=T)
    if(length(iee)>1){
      est_iee <- iee$ests; SE_iee <- iee$SEs
    }else{est_iee <- SE_iee <- rep(NA,nb)}
    if(length(SE_iee)<length(est_iee)){SE_iee <- rep(NA,length(est_iee))} ## error handling
    names(est_iee) <- paste0(c("B0","B1","B2"),"est_","iee",sep="")
    names(SE_iee) <- paste0(c("B0","B1","B2"),"SE_","iee",sep="")

    ## WEE
    wee <- try(fit_naive2stage_R(Xki=cbind(1,df$X1ki,df$X2ki),Zxki=model.matrix(X1ki~Yki*W1ki*Nki+X2ki,data=df),Nk_actual=df2$Nk,Yki=df$Yki,Xk=df2$X1k,IDk=df2$IDk,WEIGHTk=(1/df2$Nk),VALIDk=df2$VALIDk,nboot=50),silent=T)
    if(length(wee)>1){
      est_wee <- wee$ests; SE_wee <- wee$SEs
    }else{est_wee <- SE_wee <- rep(NA,nb)}
    if(length(SE_wee)<length(est_wee)){SE_wee <- rep(NA,length(est_wee))} ## error handling
    names(est_wee) <- paste0(c("B0","B1","B2"),"est_","wee",sep="")
    names(SE_wee) <- paste0(c("B0","B1","B2"),"SE_","wee",sep="")

    ##
    ests <- c(est_iee_simple,est_wee_simple,est_iee,est_wee)
    SEs <- c(SE_iee_simple,SE_wee_simple,SE_iee,SE_wee)
    names(ests) <- paste0(names(ests),"_EEE",sep="")
    names(SEs) <- paste0(names(SEs),"_EEE",sep="")
  }

  return(list(ests=ests,SEs=SEs))
}


#######################################
## fit observed likelihood approaches
fit_obsLik <- function(datMIS,dependN=TRUE,startpos=FALSE){
  ## collect data
  df <- datMIS$data; df2 <- datMIS$data_clstlvl
  ## placeholder exposure in main sample to ensure we're not using true exposure
  df$X1ki[df$VALIDki==0] <- 0
  df2$X1k[df2$VALIDk==0] <- 0

  if(dependN==FALSE){ ## misclassification model does not include N
    ## MM
    mm <- try(MISclassICS_fit(Nk=df2$Nk,Zk=cbind(1,df2$X1k),Yki=df$Yki,Xki=cbind(1,df$X1ki,df$X2ki),Wk=df2$W1k,Zwk=cbind(1,df2$X1k),Xk=df2$X1k,Zxk=matrix(1,nrow=nrow(df2),ncol=1),VALIDk=df2$VALIDk,IDk=df2$IDk,Z0k=NULL,minNk=1,NegBin=FALSE,ZIP=FALSE,slopes=TRUE,slope_col=2,condSize=TRUE,condOut=FALSE,joint=FALSE,startpos=startpos,nquad=60),silent=T)
    if(length(mm)>1){
      est_mm <- mm$ests[1:5]; SE_mm <- mm$SEs[1:5]
    }else{est_mm <- SE_mm <- rep(NA,nb+ns)}
    if(length(SE_mm)<length(est_mm)){SE_mm <- rep(NA,length(est_mm))} ## error handling
    names(est_mm) <- paste0(c("B0","B1","B2","Sig0","Sig1"),"est_","mm",sep="")
    names(SE_mm) <- paste0(c("B0","B1","B2","Sig0","Sig1"),"SE_","mm",sep="")


    ## JMM
    jmm <- try(MISclassICS_fit(Nk=df2$Nk,Zk=cbind(1,df2$X1k),Yki=df$Yki,Xki=cbind(1,df$X1ki,df$X2ki),Wk=df2$W1k,Zwk=cbind(1,df2$X1k),Xk=df2$X1k,Zxk=matrix(1,nrow=nrow(df2),ncol=1),VALIDk=df2$VALIDk,IDk=df2$IDk,Z0k=NULL,minNk=1,NegBin=FALSE,ZIP=FALSE,slopes=TRUE,slope_col=2,condSize=TRUE,condOut=FALSE,joint=TRUE,startpos=startpos,nquad=60),silent=T)
    if(length(jmm)>1){
      est_jmm <- jmm$ests[1:9]; SE_jmm <- jmm$SEs[1:9]
    }else{est_jmm <- SE_jmm <- rep(NA,na+ng+nb+ns)}
    if(length(SE_jmm)<length(est_jmm)){SE_jmm <- rep(NA,length(est_jmm))} ## error handling
    names(est_jmm) <- paste0(c("A0","A1","Gam0","Gam1","B0","B1","B2","Sig0","Sig1"),"est_","jmm",sep="")
    names(SE_jmm) <- paste0(c("A0","A1","Gam0","Gam1","B0","B1","B2","Sig0","Sig1"),"SE_","jmm",sep="")

    ##
    ests <- c(est_mm,est_jmm)
    SEs <- c(SE_mm,SE_jmm)
    names(ests) <- paste0(names(ests),"_obsLik",sep="")
    names(SEs) <- paste0(names(SEs),"_obsLik",sep="")

  }else{ ## misclassification model includes N
    ## MM
    mm <- try(MISclassICS_fit(Nk=df2$Nk,Zk=cbind(1,df2$X1k),Yki=df$Yki,Xki=cbind(1,df$X1ki,df$X2ki),Wk=df2$W1k,Zwk=cbind(1,df2$X1k,df2$Nk,df2$X1k*df2$Nk),Xk=df2$X1k,Zxk=matrix(1,nrow=nrow(df2),ncol=1),VALIDk=df2$VALIDk,IDk=df2$IDk,Z0k=NULL,minNk=1,NegBin=FALSE,ZIP=FALSE,slopes=TRUE,slope_col=2,condSize=TRUE,condOut=FALSE,joint=FALSE,interac=TRUE,inter_id=4,main_id=3,startpos=startpos,nquad=60),silent=T)
    if(length(mm)>1){
      est_mm <- mm$ests[1:5]; SE_mm <- mm$SEs[1:5]
    }else{est_mm <- SE_mm <- rep(NA,nb+ns)}
    if(length(SE_mm)<length(est_mm)){SE_mm <- rep(NA,length(est_mm))} ## error handling
    names(est_mm) <- paste0(c("B0","B1","B2","Sig0","Sig1"),"est_","mm",sep="")
    names(SE_mm) <- paste0(c("B0","B1","B2","Sig0","Sig1"),"SE_","mm",sep="")

    ## JMM
    jmm <- try(MISclassICS_fit(Nk=df2$Nk,Zk=cbind(1,df2$X1k),Yki=df$Yki,Xki=cbind(1,df$X1ki,df$X2ki),Wk=df2$W1k,Zwk=cbind(1,df2$X1k,df2$Nk,df2$X1k*df2$Nk),Xk=df2$X1k,Zxk=matrix(1,nrow=nrow(df2),ncol=1),VALIDk=df2$VALIDk,IDk=df2$IDk,Z0k=NULL,minNk=1,NegBin=FALSE,ZIP=FALSE,slopes=TRUE,slope_col=2,condSize=TRUE,condOut=FALSE,joint=TRUE,interac=TRUE,inter_id=4,main_id=3,startpos=startpos,nquad=60),silent=T)
    if(length(jmm)>1){
      est_jmm <- jmm$ests[1:9]; SE_jmm <- jmm$SEs[1:9]
    }else{est_jmm <- SE_jmm <- rep(NA,na+ng+nb+ns)}
    if(length(SE_jmm)<length(est_jmm)){SE_jmm <- rep(NA,length(est_jmm))} ## error handling
    names(est_jmm) <- paste0(c("A0","A1","Gam0","Gam1","B0","B1","B2","Sig0","Sig1"),"est_","jmm",sep="")
    names(SE_jmm) <- paste0(c("A0","A1","Gam0","Gam1","B0","B1","B2","Sig0","Sig1"),"SE_","jmm",sep="")

    ##
    ests <- c(est_mm,est_jmm)
    SEs <- c(SE_mm,SE_jmm)
    names(ests) <- paste0(names(ests),"_obsLik",sep="")
    names(SEs) <- paste0(names(SEs),"_obsLik",sep="")

  }

  return(list(ests=ests,SEs=SEs))
}


#################
## fit all models
fitmods <- function(dat,depN=FALSE,startpos=FALSE){
  g_complete <- fit_complete(dat,startpos=startpos)                ## complete data
  g_mis <- fit_mis(dat,startpos=startpos)                          ## misclassified data
  g_valid <- fit_complete(dat,valid=TRUE,startpos=startpos)        ## validation set only
  g_EEE_noN <- fit_EEE(dat,dependN=FALSE)                          ## EEE (exposure model doesnt depend on N)
  g_EEE <- fit_EEE(dat,dependN=TRUE)                               ## EEE
  g_obsLik <- fit_obsLik(dat,dependN=depN,startpos=startpos) ## obsLik

  return(list(ests=c(g_complete$ests,g_mis$ests,g_valid$ests,g_EEE_noN$ests,g_EEE$ests,g_obsLik$ests),
              SEs=c(g_complete$SEs,g_mis$SEs,g_valid$SEs,g_EEE_noN$SEs,g_EEE$SEs,g_obsLik$SEs)))

}


################
##    Loop    ##
################
print("starting loop")
res_ests <- res_SEs <- c()
## loop over datasets
for (rr in loopvec){

  if(scenario=="G0_Fixed"){
    dat <- gendat(nu_00=logit(0.85),nu_01=0,nu_10=logit(0.75),nu_11=0,gam=0,m=1,iter=rr)
    fit <- fitmods(dat,depN=FALSE,startpos=FALSE)
  }else if(scenario=="G0_BothDec"){
    dat <- gendat(nu_00=logit(0.9),nu_01=-0.7,nu_10=logit(0.9),nu_11=-0.7,gam=0,m=1,iter=rr)
    fit <- fitmods(dat,depN=TRUE,startpos=FALSE)
  }else if(scenario=="G0_BothInc"){
    dat <- gendat(nu_00=logit(0.7),nu_01=0.7,nu_10=logit(0.7),nu_11=0.7,gam=0,m=1,iter=rr)
    fit <- fitmods(dat,depN=TRUE,startpos=FALSE)
  }else if(scenario=="G0_SensDec"){
    dat <- gendat(nu_00=logit(0.99),nu_01=0,nu_10=logit(0.65),nu_11=-0.3,gam=0,m=1,iter=rr)
    fit <- fitmods(dat,depN=TRUE,startpos=FALSE)
  }else if(scenario=="GN25_Fixed"){
    dat <- gendat(nu_00=logit(0.85),nu_01=0,nu_10=logit(0.75),nu_11=0,gam=-0.25,m=1,iter=rr)
    fit <- fitmods(dat,depN=FALSE,startpos=FALSE)
  }else if(scenario=="GN25_BothDec"){
    dat <- gendat(nu_00=logit(0.9),nu_01=-0.7,nu_10=logit(0.9),nu_11=-0.7,gam=-0.25,m=1,iter=rr)
    fit <- fitmods(dat,depN=TRUE,startpos=FALSE)
  }else if(scenario=="GN25_BothInc"){
    dat <- gendat(nu_00=logit(0.7),nu_01=0.7,nu_10=logit(0.7),nu_11=0.7,gam=-0.25,m=1,iter=rr)
    fit <- fitmods(dat,depN=TRUE,startpos=FALSE)
  }else if(scenario=="GN25_SensDec"){
    dat <- gendat(nu_00=logit(0.99),nu_01=0,nu_10=logit(0.65),nu_11=-0.3,gam=-0.25,m=1,iter=rr)
    fit <- fitmods(dat,depN=TRUE,startpos=FALSE)
  }else if(scenario=="G25_Fixed"){
    dat <- gendat(nu_00=logit(0.85),nu_01=0,nu_10=logit(0.75),nu_11=0,gam=0.25,m=1,iter=rr)
    fit <- fitmods(dat,depN=FALSE,startpos=TRUE)
  }else if(scenario=="G25_BothDec"){
    dat <- gendat(nu_00=logit(0.9),nu_01=-0.7,nu_10=logit(0.9),nu_11=-0.7,gam=0.25,m=1,iter=rr)
    fit <- fitmods(dat,depN=TRUE,startpos=TRUE)
  }else if(scenario=="G25_BothInc"){
    dat <- gendat(nu_00=logit(0.7),nu_01=0.7,nu_10=logit(0.7),nu_11=0.7,gam=0.25,m=1,iter=rr)
    fit <- fitmods(dat,depN=TRUE,startpos=TRUE)
  }else if(scenario=="G25_SensDec"){
    dat <- gendat(nu_00=logit(0.99),nu_01=0,nu_10=logit(0.65),nu_11=-0.3,gam=0.25,m=1,iter=rr)
    fit <- fitmods(dat,depN=TRUE,startpos=TRUE)
  }



  res_ests <- rbind(res_ests,fit$ests)
  res_SEs <- rbind(res_SEs,fit$SEs)


  print(rr)

}

results <- cbind(res_ests,res_SEs)

## save
write.table(results,file=paste0(path,"results_sim",scenario,suffix,".txt"),col.names=TRUE)









