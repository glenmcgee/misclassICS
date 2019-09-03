##############################
##        Simulation 1      ##
##############################
## Updated: 04/07/2019
## generating data:
## with random slopes (ICS and not)
## and MISCLASSIFICATION

## misclassification scenario (see manuscript)
# scenario <- 1 # 2 3 4 5 6 ()

runLOCAL=TRUE ## set to TRUE to run locally ## FALSE is on cluster
GRAY=TRUE ## for black and white plots

#####################
## setting parameters
RR <- 2000 ## No. of iterations
KK <- 8000 ## No. of clusters
Vpercent <- 0.20 ## Validation sampling fraction
# VV <- 1000  ## No. of clusters in (internal) validation sample
JJ <- 500 ## Number of jobs in array (only matters on cluster)
Xprev <- 0.25 ## Exposure prevalence
a0 <- 0.6; ## volume model coefficients
# a1 <- -0.2
b0 <- -4;  b2 <- 0.2 ## outcome model coefficients
# b1 <- 0.5;
#gamma <- -0.1 ## scaling factor for random intercept in volume model
sigma0 <- 2 ## SD of rand eff for X=0
sigma1 <- 1.5 ## SD of rand eff for X=1
nb <- 3 ## number of betas
na <- 2 ## number of alphas
ng <- 2 ## number of gammas (2 for random slopes)
ns <- 2 ## number of sigmas (2 for random slopes)

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
# source(paste0(funpath,"fitEEE.R"))
require(lme4)
require(geepack)
require(JMMICSpack)
require(wesanderson)

## colors
wes_red <- wes_palette(n=5, name="Darjeeling")[1]
wes_green <- wes_palette(n=5, name="Darjeeling")[2]
wes_gold <- wes_palette(n=5, name="Darjeeling")[3]
wes_orange <- wes_palette(n=5, name="Darjeeling")[4]
wes_blue <- wes_palette(n=5, name="Darjeeling")[5]
if(GRAY==TRUE){
  wes_orange <- "gray45"
  wes_red <- "gray98"
  wes_gold <- "gray75"
}



#######################################
## Generate data with misclassification
gendat <- function(a1=-0.2,
                   b1=0.5,
                   nu_00=logit(0.90), ## baseline specificities
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

  # ## MM
  # mm <- try(JMMICS_fit(Nk=df2$Nk,Zk=cbind(1,df2$X1k),Yki=df$Yki,Xki=cbind(1,df$X1ki,df$X2ki),IDk=df2$IDk,Z0k=NULL,minNk=1,NegBin=FALSE,ZIP=FALSE,slopes=TRUE,slope_col=2,condSize=TRUE,condOut=FALSE,joint=FALSE,startpos=startpos,nquad=60),silent=T)
  # if(length(mm)>1){
  #   est_mm <- mm$ests; SE_mm <- mm$SEs
  # }else{est_mm <- SE_mm <- rep(NA,nb+ns)}
  # if(length(SE_mm)<length(est_mm)){SE_mm <- rep(NA,length(est_mm))} ## error handling
  # names(est_mm) <- paste0(c("B0","B1","B2","Sig0","Sig1"),"est_","mm",sep="")
  # names(SE_mm) <- paste0(c("B0","B1","B2","Sig0","Sig1"),"SE_","mm",sep="")
  #
  # ## JMM
  # jmm <- try(JMMICS_fit(Nk=df2$Nk,Zk=cbind(1,df2$X1k),Yki=df$Yki,Xki=cbind(1,df$X1ki,df$X2ki),IDk=df2$IDk,Z0k=NULL,minNk=1,NegBin=FALSE,ZIP=FALSE,slopes=TRUE,slope_col=2,condSize=TRUE,condOut=FALSE,joint=TRUE,startpos=startpos,nquad=60),silent=T)
  # if(length(jmm)>1){
  #   est_jmm <- jmm$ests; SE_jmm <- jmm$SEs
  # }else{est_jmm <- SE_jmm <- rep(NA,na+ng+nb+ns)}
  # if(length(SE_jmm)<length(est_jmm)){SE_jmm <- rep(NA,length(est_jmm))} ## error handling
  # names(est_jmm) <- paste0(c("A0","A1","Gam0","Gam1","B0","B1","B2","Sig0","Sig1"),"est_","jmm",sep="")
  # names(SE_jmm) <- paste0(c("A0","A1","Gam0","Gam1","B0","B1","B2","Sig0","Sig1"),"SE_","jmm",sep="")

  ##
  ests <- c(est_gee,est_iee,est_wee)#,est_mm,est_jmm)
  SEs <- c(SE_gee,SE_iee,SE_wee)#,SE_mm,SE_jmm)
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

  # ## MM
  # mm <- try(JMMICS_fit(Nk=df2$Nk,Zk=cbind(1,df2$W1k),Yki=df$Yki,Xki=cbind(1,df$W1ki,df$X2ki),IDk=df2$IDk,Z0k=NULL,minNk=1,NegBin=FALSE,ZIP=FALSE,slopes=TRUE,slope_col=2,condSize=TRUE,condOut=FALSE,joint=FALSE,startpos=startpos,nquad=60),silent=T)
  # if(length(mm)>1){
  #   est_mm <- mm$ests; SE_mm <- mm$SEs
  # }else{est_mm <- SE_mm <- rep(NA,nb+ns)}
  # if(length(SE_mm)<length(est_mm)){SE_mm <- rep(NA,length(est_mm))} ## error handling
  # names(est_mm) <- paste0(c("B0","B1","B2","Sig0","Sig1"),"est_","mm",sep="")
  # names(SE_mm) <- paste0(c("B0","B1","B2","Sig0","Sig1"),"SE_","mm",sep="")
  #
  # ## JMM
  # jmm <- try(JMMICS_fit(Nk=df2$Nk,Zk=cbind(1,df2$W1k),Yki=df$Yki,Xki=cbind(1,df$W1ki,df$X2ki),IDk=df2$IDk,Z0k=NULL,minNk=1,NegBin=FALSE,ZIP=FALSE,slopes=TRUE,slope_col=2,condSize=TRUE,condOut=FALSE,joint=TRUE,startpos=startpos,nquad=60),silent=T)
  # if(length(jmm)>1){
  #   est_jmm <- jmm$ests; SE_jmm <- jmm$SEs
  # }else{est_jmm <- SE_jmm <- rep(NA,na+ng+nb+ns)}
  # if(length(SE_jmm)<length(est_jmm)){SE_jmm <- rep(NA,length(est_jmm))} ## error handling
  # names(est_jmm) <- paste0(c("A0","A1","Gam0","Gam1","B0","B1","B2","Sig0","Sig1"),"est_","jmm",sep="")
  # names(SE_jmm) <- paste0(c("A0","A1","Gam0","Gam1","B0","B1","B2","Sig0","Sig1"),"SE_","jmm",sep="")
  #
  ##
  ests <- c(est_gee,est_iee,est_wee)#,est_mm,est_jmm)
  SEs <- c(SE_gee,SE_iee,SE_wee)#,SE_mm,SE_jmm)
  names(ests) <- paste0(names(ests),"_mis",sep="")
  names(SEs) <- paste0(names(SEs),"_mis",sep="")


  return(list(ests=ests,SEs=SEs))
}


#################
## fit all models
fitmods <- function(dat){ #},depN=FALSE,startpos=FALSE){
  g_complete <- fit_complete(dat,startpos=startpos)                ## complete data
  g_mis <- fit_mis(dat,startpos=startpos)                          ## misclassified data
  # g_valid <- fit_complete(dat,valid=TRUE,startpos=startpos)        ## validation set only
  # g_EEE_noN <- fit_EEE(dat,dependN=FALSE)                          ## EEE (exposure model doesnt depend on N)
  # g_EEE <- fit_EEE(dat,dependN=TRUE)                               ## EEE
  # g_obsLik <- fit_obsLik(dat,dependN=depN,startpos=startpos) ## obsLik

  return(list(ests=c(g_complete$ests,g_mis$ests),#g_valid$ests,g_EEE_noN$ests,g_EEE$ests,g_obsLik$ests),
              SEs=c(g_complete$SEs,g_mis$SEs)))#,g_valid$SEs,g_EEE_noN$SEs,g_EEE$SEs,g_obsLik$SEs)))

}


################
##    Loop    ##
################
print("starting loop")
# res_ests <- res_SEs <- c()
res_00 <- res_02 <- res_06 <- res_1 <- c()
resdep_00 <- resdep_02 <- resdep_06 <- resdep_1 <- c()
## loop over datasets
for (rr in loopvec){

  ## simple misclassification
  dat <- gendat(a1=-0,b1=1,nu_00=logit(0.85),nu_01=0,nu_10=logit(0.75),nu_11=0,gam=0,m=1,iter=rr)
  res_00 <- rbind(res_00,fitmods(dat)$ests) #,depN=FALSE,startpos=FALSE)

  dat <- gendat(a1=-0.2,b1=1,nu_00=logit(0.85),nu_01=0,nu_10=logit(0.75),nu_11=0,gam=0,m=1,iter=rr)
  res_02 <- rbind(res_02,fitmods(dat)$ests) #,depN=FALSE,startpos=FALSE)

  dat <- gendat(a1=-0.6,b1=1,nu_00=logit(0.85),nu_01=0,nu_10=logit(0.75),nu_11=0,gam=0,m=1,iter=rr)
  res_06 <- rbind(res_06,fitmods(dat)$ests) #,depN=FALSE,startpos=FALSE)

  dat <- gendat(a1=-1,b1=1,nu_00=logit(0.85),nu_01=0,nu_10=logit(0.75),nu_11=0,gam=0,m=1,iter=rr)
  res_1 <- rbind(res_1,fitmods(dat)$ests) #,depN=FALSE,startpos=FALSE)



  print(rr)

}

colnames(res_00) <- paste0(colnames(res_00),"_00")
colnames(res_02) <- paste0(colnames(res_02),"_02")
colnames(res_06) <- paste0(colnames(res_06),"_06")
colnames(res_1) <- paste0(colnames(res_1),"_1")

res <- cbind(res_00,res_02,res_06,res_1)
write.table(res,file=paste0(path,"inducedICS_results",suffix,".txt"),col.names=TRUE)
# res <- read.table(paste0(path,"inducedICS_results.txt"),header = TRUE);





## make boxplots
geenams <- paste0("B1est_gee_mis_",c("00","02","06","1"))
ieenams <- paste0("B1est_iee_mis_",c("00","02","06","1"))
weenams <- paste0("B1est_wee_mis_",c("00","02","06","1"))
B1 <- data.frame(res[,colnames(res)%in%c(geenams,ieenams,weenams)])
B1 <- B1%>%gather("Group","Est",1:ncol(B1))
B1$color <- 1; B1$color[grepl("iee", B1$Group)] <- 2; B1$color[grepl("wee", B1$Group)] <- 3
B1$opacity <- 0.95; #B1$opacity[B1$Group %in% c("4JMM","7Joint")] <- 1
B1$lvl <- 0; B1$lvl[grepl("_02", B1$Group)] <- 1; B1$lvl[grepl("_06", B1$Group)] <- 2; B1$lvl[grepl("_1", B1$Group)] <- 3;


box_B1 <- ggplot(B1,aes(x=factor(lvl),Est,fill=as.factor(color),alpha=as.numeric(opacity)))+
  geom_boxplot(aes(),outlier.size=0.5)+
  geom_hline(yintercept=1)+
  scale_y_continuous(limits=c(0.15,1.35),expand=c(0.01,0.01))+
  scale_fill_manual(values = c(wes_orange,wes_red,wes_gold),guide=F)+
  scale_alpha_continuous(guide=F,range=c(0.8,1))+
  ylab("Naive Estimates")+xlab(expression(paste(X[1][k],"~",N[k]," (log RR)")))+
  scale_x_discrete(breaks=c(0,1,2,3),labels=c(0,-0.2,-0.6,-1.0))+ #,
  theme_classic() +
  ggtitle("") +
  theme(plot.title=element_text(hjust=0.99,margin=margin(t=5,b=-20),size=10))
ggsave(filename=paste0(outpath,"inducedICS_plot_gee.pdf"),plot=box_B1,width=3,height=6)


box_B1_nogee <- ggplot(B1[B1$color!=1,],aes(x=factor(lvl),Est,fill=as.factor(color),alpha=as.numeric(opacity)))+
  geom_boxplot(aes(),outlier.size=0.5)+
  geom_hline(yintercept=1)+
  scale_y_continuous(limits=c(0.15,1.35),expand=c(0.01,0.01))+
  scale_fill_manual(values = c(wes_red,wes_gold),guide=F)+
  scale_alpha_continuous(guide=F,range=c(0.8,1))+
  ylab("Naive Estimates")+xlab(expression(paste(X[1][k],"~",N[k]," (log RR)")))+
  scale_x_discrete(breaks=c(0,1,2,3),labels=c(0,-0.2,-0.6,-1.0))+ #,
  theme_classic() +
  ggtitle("") +
  theme(plot.title=element_text(hjust=0.99,margin=margin(t=5,b=-20),size=10))
ggsave(filename=paste0(outpath,"inducedICS_plot.pdf"),plot=box_B1_nogee,width=3,height=6)




## make boxplots
geenams <- paste0("B1est_gee_complete_",c("00","02","06","1"))
ieenams <- paste0("B1est_iee_complete_",c("00","02","06","1"))
weenams <- paste0("B1est_wee_complete_",c("00","02","06","1"))
B1 <- data.frame(res[,colnames(res)%in%c(geenams,ieenams,weenams)])
B1 <- B1%>%gather("Group","Est",1:ncol(B1))
B1$color <- 1; B1$color[grepl("iee", B1$Group)] <- 2; B1$color[grepl("wee", B1$Group)] <- 3
B1$opacity <- 0.95; #B1$opacity[B1$Group %in% c("4JMM","7Joint")] <- 1
B1$lvl <- 0; B1$lvl[grepl("_02", B1$Group)] <- 1; B1$lvl[grepl("_06", B1$Group)] <- 2; B1$lvl[grepl("_1", B1$Group)] <- 3;


box_B1 <- ggplot(B1,aes(x=factor(lvl),Est,fill=as.factor(color),alpha=as.numeric(opacity)))+
  geom_boxplot(aes(),outlier.size=0.5)+
  geom_hline(yintercept=1)+
  scale_y_continuous(limits=c(0.15,1.35),expand=c(0.01,0.01))+
  scale_fill_manual(values = c(wes_orange,wes_red,wes_gold),guide=F)+
  scale_alpha_continuous(guide=F,range=c(0.8,1))+
  ylab("Estimates")+xlab(expression(paste(X[1][k],"~",N[k]," (log RR)")))+
  scale_x_discrete(breaks=c(0,1,2,3),labels=c(0,-0.2,-0.6,-1.0))+ #,
  theme_classic() +
  ggtitle("") +
  theme(plot.title=element_text(hjust=0.99,margin=margin(t=5,b=-20),size=10))
ggsave(filename=paste0(outpath,"inducedICS_plot_gee_TRUE.pdf"),plot=box_B1,width=3,height=6)


box_B1_nogee <- ggplot(B1[B1$color!=1,],aes(x=factor(lvl),Est,fill=as.factor(color),alpha=as.numeric(opacity)))+
  geom_boxplot(aes(),outlier.size=0.5)+
  geom_hline(yintercept=1)+
  scale_y_continuous(limits=c(0.15,1.35),expand=c(0.01,0.01))+
  scale_fill_manual(values = c(wes_red,wes_gold),guide=F)+
  scale_alpha_continuous(guide=F,range=c(0.8,1))+
  ylab("Estimates")+xlab(expression(paste(X[1][k],"~",N[k]," (log RR)")))+
  scale_x_discrete(breaks=c(0,1,2,3),labels=c(0,-0.2,-0.6,-1.0))+ #,
  theme_classic() +
  ggtitle("") +
  theme(plot.title=element_text(hjust=0.99,margin=margin(t=5,b=-20),size=10))
ggsave(filename=paste0(outpath,"inducedICS_plot_TRUE.pdf"),plot=box_B1_nogee,width=3,height=6)

