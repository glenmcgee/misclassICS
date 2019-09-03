#######################
##    NHS Analysis   ##
#######################
## 05/28/2019
## Fit models to NHS data
## prepare dataframes for JMMICS pack
## then fit:
 ## Conditional:
  ## joint conditional (no 0s)
  ## joint conditional
  ## outcome-only glmm
 ## Marginal:
  ## GEE-Exch
  ## IEE
  ## wee
  ## JMM
  ## marginal size model (for comparing joint models)

## load packages
#library(tidyr)
#library(dplyr)
#library(xtable)
library(glmmML)
library(lme4)
library(geepack)
library(ggplot2)
library(JMMICSpack)
library(ResourceSelection) ## for hosmer lemeshow
expit <- function(x)exp(x)/(1+exp(x))
#library(MISclassICSpack) ## called by "MISclassICS_C.R"
source("../Functions/JMMICS.R")
source("../Functions/MISclassICS_C.R") ## for observed likelihood
source("../Functions/MISclassICS.R") ## for EEE
## load and clean data
#source("~/code/1Mis_clean_data.R") ## called by "2Mis_EDA.R"
## prelim tables
source("../Data Analysis/2Mis_EDA.R")


##################################
## ignore non-mothers
datG1 <- datG1_no0
datG2 <- datG2_no0

##################################
## avoid NAs
## just a placeholder; 0 not used.
datG1$mdes1[is.na(datG1$mdes1)] <- 0
datG2$mdes1[is.na(datG2$mdes1)] <- 0

################################
### prepare data for JMMICSpack
YY <- datG2$adhd ## outcome (ki)
XX <- cbind(1,datG2$desqx1,datG2$msmk2,datG2$yob89_5155,datG2$yob89_5660,datG2$yob89_61plus) ## outcome-related covariates (ki)    #,datG2$raceWhite,datG2$momed2,datG2$momed3,datG2$momed4)
NN <- datG1$totalkids ## cluster size (k)
ZZ <- cbind(1,datG1$desqx1,datG1$msmk2,datG1$yob89_5155,datG1$yob89_5660,datG1$yob89_61plus) ## size-related covaraites (k)
IDD <- datG1$id2 ## cluster id (k)
IDV <- datG1$validID ## member of validation set (k)

### prepare validation data
vYY <- datG2$adhd[datG2$validID==1]
vXX <- cbind(1,datG2$mdes1[datG2$validID==1],datG2$msmk2[datG2$validID==1],datG2$yob89_5155[datG2$validID==1],datG2$yob89_5660[datG2$validID==1],datG2$yob89_61plus[datG2$validID==1])#,datG2$raceWhite,datG2$momed2,datG2$momed3,datG2$momed4)
vNN <- datG1$totalkids[datG1$validID==1]
vZZ <- cbind(1,datG1$mdes1[datG1$validID==1],datG1$msmk2[datG1$validID==1],datG1$yob89_5155[datG1$validID==1],datG1$yob89_5660[datG1$validID==1],datG1$yob89_61plus[datG1$validID==1])
# ZZ0 <- ZZ[,1:2]
vIDD <- as.numeric(factor(datG1$id2[datG1$validID==1])) #datG2$id2[datG2$validID==1]

## prepare for EEE
XX1ki <- cbind(1,datG2$mdes1,XX[,-(1:2)]) ## outcome model with true exposure (ki)
ZZxki <- cbind(1,datG2$desqx1,datG2$totalkids,YY,datG2$desqx1*datG2$totalkids,datG2$desqx1*YY,YY*datG2$totalkids,datG2$desqx1*YY*datG2$totalkids,XX[,-(1:2)]) ## conditional exposure model
XX1k <- datG1$mdes1 ## true exposure (k) (where observed)
## version 2 for simpler EEE with no interactions
ZZxki_2 <- cbind(1,datG2$desqx1,datG2$totalkids,YY,XX[,-(1:2)]) ## conditional exposure model


## prepare for obsLik
ZZ1 <- cbind(1,XX1k,ZZ[,-(1:2)]) ## size-related covariates with true exposure (k)
#XX1ki as above
WW1k <- datG1$desqx1 ## misclassified exposure (k)
ZZw <- cbind(1,XX1k,NN,XX1k*NN,ZZ[,-(1:2)]) ## misclassification model


########################################
##            Model Fits              ##
########################################

########################################
## Proposed corrections
## EEE
iee_corr <- suppressWarnings(fit_naive2stage_R(Xki=XX1ki,Zxki=ZZxki,Nk_actual=NN,Yki=YY,Xk=XX1k,IDk=IDD,WEIGHTk=NA,VALIDk=IDV,nboot=200))
iee_corr_est <- iee_corr$ests
iee_corr_SE <- iee_corr$SEs
print("iee_corr")
wee_corr <- suppressWarnings(fit_naive2stage_R(Xki=XX1ki,Zxki=ZZxki,Nk_actual=NN,Yki=YY,Xk=XX1k,IDk=IDD,WEIGHTk=(1/NN),VALIDk=IDV,nboot=200))
wee_corr_est <- wee_corr$ests
wee_corr_SE <- wee_corr$SEs
print("wee_corr")
## EEE ## simplified exposure model
iee_corr_2 <- suppressWarnings(fit_naive2stage_R(Xki=XX1ki,Zxki=ZZxki_2,Nk_actual=NN,Yki=YY,Xk=XX1k,IDk=IDD,WEIGHTk=NA,VALIDk=IDV,nboot=200))
iee_corr_2_est <- iee_corr_2$ests
iee_corr_2_SE <- iee_corr_2$SEs
print("iee_corr_2")
wee_corr_2 <- suppressWarnings(fit_naive2stage_R(Xki=XX1ki,Zxki=ZZxki_2,Nk_actual=NN,Yki=YY,Xk=XX1k,IDk=IDD,WEIGHTk=(1/NN),VALIDk=IDV,nboot=200))
wee_corr_2_est <- wee_corr_2$ests
wee_corr_2_SE <- wee_corr_2$SEs
print("wee_corr_2")
## ObsLik
jmm_corr <- suppressWarnings(MISclassICS_fit(Nk=NN,Zk=ZZ1,Yki=YY,Xki=XX1ki,Wk=WW1k,Zwk=ZZw,Xk=XX1k,Zxk=ZZ[,-2],VALIDk=IDV,IDk=IDD,Z0k=NULL,minNk=1,NegBin=FALSE,ZIP=FALSE,slopes=TRUE,slope_col=2,condSize=TRUE,condOut=FALSE,joint=TRUE,interac=FALSE,inter_id=4,main_id=3,startpos=FALSE,nquad=50))
jmm_corr_est <- jmm_corr$ests; print(jmm_corr_est)
jmm_corr_SE <- jmm_corr$SEs; print(jmm_corr_SE)


########################################
## Validation Set Only
##
iee_valid <- suppressWarnings(geeglm(adhd~mdes1+msmk2+yob89_5155+yob89_5660+yob89_61plus,id=id2,data=datG2[datG2$validID==1,],family=binomial,corstr="independence"))
iee_valid_est <- iee_valid$coef
iee_valid_SE <- summary(iee_valid)$coef[,2]
print("iee_valid")
##
wee_valid <- suppressWarnings(geeglm(adhd~mdes1+msmk2+yob89_5155+yob89_5660+yob89_61plus,id=id2,data=datG2[datG2$validID==1,],family=binomial,corstr="independence",weights=(1/totalkids)))
wee_valid_est <- wee_valid$coef
wee_valid_SE <- summary(wee_valid)$coef[,2]
print("wee_valid")
##
jmm_valid <- suppressWarnings(JMMICS_fit(Nk=vNN,Zk=vZZ,Yki=vYY,Xki=vXX,IDk=vIDD,Z0k=vZZ,weights=NA,minNk=1,NegBin=FALSE,ZIP=FALSE,slopes=TRUE,slope_col=2,condSize=TRUE,condOut=FALSE,joint=TRUE,nquad=50))
jmm_valid_est <- jmm_valid[[1]]
jmm_valid_SE <- jmm_valid[[2]]
print("jmm_valid")

########################################
## Naive Models
##
iee_naive <- suppressWarnings(geeglm(adhd~desqx1+msmk2+yob89_5155+yob89_5660+yob89_61plus,id=id2,data=datG2,family=binomial,corstr="independence"))
iee_naive_est <- iee_naive$coef
iee_naive_SE <- summary(iee_naive)$coef[,2]
print("iee_naive")
##
wee_naive <- suppressWarnings(geeglm(adhd~desqx1+msmk2+yob89_5155+yob89_5660+yob89_61plus,id=id2,data=datG2,family=binomial,corstr="independence",weights=(1/totalkids)))
wee_naive_est <- wee_naive$coef
wee_naive_SE <- summary(wee_naive)$coef[,2]
print("wee_naive")
##
jmm_naive <- suppressWarnings(JMMICS_fit(Nk=NN,Zk=ZZ,Yki=YY,Xki=XX,IDk=IDD,Z0k=ZZ,weights=NA,minNk=1,NegBin=FALSE,ZIP=FALSE,slopes=TRUE,slope_col=2,condSize=TRUE,condOut=FALSE,joint=TRUE,nquad=50))
jmm_naive_est <- jmm_naive[[1]]
jmm_naive_SE <- jmm_naive[[2]]
print("jmm_naive")




########################################
##          Collect Results           ##
########################################
n_a <- ncol(ZZ)        ## no. of alphas
n_b <- ncol(XX)        ## no. of betas

## note: not including nuisance parameters epsW and epsX, but they get printed above

## order of parameters:    alphas,        betas,    sigma0, sigma1,   gamma0, gamma1
                    ##  rep(NA,n_a),  rep(NA,n_b),   rep(NA,2),        rep(NA,2))

## order that parameters are output from JMMICSpack: alpha gamma0 gamma1 beta sigma0 sigma1
## order that parameters are output from MISclassICSpack: alpha gamma0 gamma1 beta sigma0 sigma1 epsW epsX

## combine results
ests <- cbind(c(rep(NA,n_a),iee_naive_est,rep(NA,2),rep(NA,2)),
              c(rep(NA,n_a),wee_naive_est,rep(NA,2),rep(NA,2)),
              c(jmm_naive_est[1:n_a],jmm_naive_est[(n_a+2)+(1:n_b)],jmm_naive_est[(n_a+2+n_b)+(1:2)],jmm_naive_est[(n_a)+(1:2)]),
              c(rep(NA,n_a),iee_valid_est,rep(NA,2),rep(NA,2)),
              c(rep(NA,n_a),wee_valid_est,rep(NA,2),rep(NA,2)),
              c(jmm_valid_est[1:n_a],jmm_valid_est[(n_a+2)+(1:n_b)],jmm_valid_est[(n_a+2+n_b)+(1:2)],jmm_valid_est[(n_a)+(1:2)]),
              c(rep(NA,n_a),iee_corr_est,rep(NA,2),rep(NA,2)),
              c(rep(NA,n_a),wee_corr_est,rep(NA,2),rep(NA,2)),
              c(rep(NA,n_a),iee_corr_2_est,rep(NA,2),rep(NA,2)),
              c(rep(NA,n_a),wee_corr_2_est,rep(NA,2),rep(NA,2)),
              c(jmm_corr_est[1:n_a],jmm_corr_est[(n_a+2)+(1:n_b)],jmm_corr_est[(n_a+2+n_b)+(1:2)],jmm_corr_est[(n_a)+(1:2)]) )

SEs <- cbind(c(rep(NA,n_a),iee_naive_SE,rep(NA,2),rep(NA,2)),
              c(rep(NA,n_a),wee_naive_SE,rep(NA,2),rep(NA,2)),
              c(jmm_naive_SE[1:n_a],jmm_naive_SE[(n_a+2)+(1:n_b)],jmm_naive_SE[(n_a+2+n_b)+(1:2)],jmm_naive_SE[(n_a)+(1:2)]),
              c(rep(NA,n_a),iee_valid_SE,rep(NA,2),rep(NA,2)),
              c(rep(NA,n_a),wee_valid_SE,rep(NA,2),rep(NA,2)),
              c(jmm_valid_SE[1:n_a],jmm_valid_SE[(n_a+2)+(1:n_b)],jmm_valid_SE[(n_a+2+n_b)+(1:2)],jmm_valid_SE[(n_a)+(1:2)]),
              c(rep(NA,n_a),iee_corr_SE,rep(NA,2),rep(NA,2)),
              c(rep(NA,n_a),wee_corr_SE,rep(NA,2),rep(NA,2)),
              c(rep(NA,n_a),iee_corr_2_SE,rep(NA,2),rep(NA,2)),
              c(rep(NA,n_a),wee_corr_2_SE,rep(NA,2),rep(NA,2)),
              c(jmm_corr_SE[1:n_a],jmm_corr_SE[(n_a+2)+(1:n_b)],jmm_corr_SE[(n_a+2+n_b)+(1:2)],jmm_corr_SE[(n_a)+(1:2)]) )

## label results
rownames(ests) <- rownames(SEs) <- c("a0","a1 DES","a2 msmk","a3 yob5155","a4 yob5660","a5 yob61plus",
                                     "b0","b1 DES","b2 msmk","b3 yob5155","b4 yob5660","b5 yob61plus",
                                     "Sigma0","Sigma1",
                                     "Gamma0","Gamma1" )
colnames(ests) <- colnames(SEs) <- c("IEE-Naive","WEE-Naive","JMM-Naive",
                                     "IEE-Valid","WEE-Valid","JMM-Valid",
                                     "IEE-EEE","WEE-EEE",
                                     "IEE-EEE-Simple","WEE-EEE-Simple","JMM-ObsLik" )

## save output
write.table(ests,file="Mis_ests_2.txt")
write.table(SEs,file="Mis_SEs_2.txt")



## make tables for latex
xtable(ests)
xtable(SEs)


############################################
## Extension: hosmer lemeshow
## goodness of fit check for EEE

## fit conditional exposure model with interactions
gX_EEE <- glm(datG2$mdes1[datG2$validID==1]~ZZxki[datG2$validID==1,]-1,family="binomial")
summary(gX_EEE)
## fit conditional exposure model with main effects only
gX_EEE_2 <- glm(datG2$mdes1[datG2$validID==1]~ZZxki_2[datG2$validID==1,]-1,family="binomial")
summary(gX_EEE_2)

## LRT (simple model nested within more complex)
anova(gX_EEE,gX_EEE_2,test="LRT")

### EEE hosmer lemeshow test
hoslem.test(x=datG2$mdes1[datG2$validID==1],y=fitted(gX_EEE),g=10)
##EEE simple
hoslem.test(x=datG2$mdes1[datG2$validID==1],y=fitted(gX_EEE_2),g=10)



