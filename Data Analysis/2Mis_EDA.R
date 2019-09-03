##################
##    NHS EDA   ##
##################
## 05/07/2019
## Generate tables and plots for EDA of NHS data
## 1 - generate table of ADHD and DES by cluster size
## 2 - generate table of emptiness rates in each covariate category
## 3 - some prelim tables for ADHD/DES rates by covariates
## 4 - run  preliminary regressions on size, emptiness, and outcome to determine presence of ICS, etc.

## load packages
#library(tidyr)
#library(dplyr)
library(xtable)
## load clean data
source("~/code/1Mis_clean_data.R")


######################
## some summary stats
## no. nurses
length(unique(dat$id))
## rate of adhd diagnosis
mean(datG2$adhd)
## rate of "having at least one child with adhd given that you had at least one child"
mean((datG1_no0$adhd>0))
## rate of multiple diagnoses, given at least one diagnosis
mean((datG1$adhd*datG1$totalkids>1)[datG1_no0$adhd>0])
## rate of DES nurse-reported exposure
mean((datG1$desqx1))
## rate of DES mother-reported exposure
mean((datG1$mdes1))
## misclassification and cluster size
table(datG1$mdes1,datG1$desqx1,datG1$totalkids)


###################
## 1 Nk table
#### collapse large cluster sizes
dat$Nk_collapse4plus <- dat$totalkids; dat$Nk_collapse4plus[dat$Nk_collapse4plus>=4] <- "4+"
datG1$Nk_collapse4plus <- datG1$totalkids; datG1$Nk_collapse4plus[datG1$Nk_collapse4plus>=4] <- "4+"

tabNk_1 <- dat %>% group_by(Nk_collapse4plus) %>% summarise(K=length(unique(id)),
                                                     Out_pct=round(100*mean(adhd),2))
tabNk_2 <- datG1 %>% group_by(Nk_collapse4plus) %>% summarise(K=length(unique(id)),
                                                     Exp_pct=round(100*mean(desqx1),2))

## characteristics about validation sample
VALIDtabNk_1 <- dat %>% filter(!is.na(mdes1)) %>% group_by(Nk_collapse4plus) %>% summarise(K=length(unique(id)),
                                                                                            ValidOut_pct=round(100*mean(adhd),2))
VALIDtabNk_2 <- datG1 %>% filter(!is.na(mdes1)) %>% group_by(Nk_collapse4plus) %>% summarise(K=length(unique(id)),
                                                                                             ValidExp_pct=round(100*mean(desqx1),2),
                                                                                             ValidMExp_pct=round(100*mean(mdes1),2))
## summary of concordance
CONCtabNk_1 <- datG1 %>% filter(mdes1==1) %>% group_by(Nk_collapse4plus) %>% summarise(K=length(unique(id)),
                                                                                             SENS_pct=round(100*mean(desqx1),2))
CONCtabNk_2 <- datG1 %>% filter(mdes1==0) %>% group_by(Nk_collapse4plus) %>% summarise(K=length(unique(id)),
                                                                                       SPEC_pct=round(100*mean(desqx0),2))


tabNk <- cbind(tabNk_1,tabNk_2[,3],VALIDtabNk_1[,-1],VALIDtabNk_2[,3:4],CONCtabNk_1[,3],CONCtabNk_2[,3])
colnames(tabNk) <- c("N_k","No. G1s","% ADHD","% DES","Valid % ADHD","Valid % DES","Valid % MDES","Sensitivity","Specificity")
tabNk[is.na(tabNk)] <- "--"
xtable(tabNk,include.rownames = FALSE)
xtable(tabNk[,-(6:8)],include.rownames = FALSE)
write.table(tabNk,"~/results/Mis_tabNk.txt")


###################
##  Same table but this time only collapse clusters of size 5 or larger
dat$Nk_collapse5plus <- dat$totalkids; dat$Nk_collapse5plus[dat$Nk_collapse5plus>=5] <- "5+"
datG1$Nk_collapse5plus <- datG1$totalkids; datG1$Nk_collapse5plus[datG1$Nk_collapse5plus>=5] <- "5+"

tabNk_1 <- dat %>% group_by(Nk_collapse5plus) %>% summarise(K=length(unique(id)),
                                                            Out_pct=round(100*mean(adhd),2))
tabNk_2 <- datG1 %>% group_by(Nk_collapse5plus) %>% summarise(K=length(unique(id)),
                                                              Exp_pct=round(100*mean(desqx1),2))

## characteristics about validation sample
VALIDtabNk_1 <- dat %>% filter(!is.na(mdes1)) %>% group_by(Nk_collapse5plus) %>% summarise(K=length(unique(id)),
                                                                                           ValidOut_pct=round(100*mean(adhd),2))
VALIDtabNk_2 <- datG1 %>% filter(!is.na(mdes1)) %>% group_by(Nk_collapse5plus) %>% summarise(K=length(unique(id)),
                                                                                             ValidExp_pct=round(100*mean(desqx1),2),
                                                                                             ValidMExp_pct=round(100*mean(mdes1),2))
## summary of concordance
CONCtabNk_1 <- datG1 %>% filter(mdes1==1) %>% group_by(Nk_collapse5plus) %>% summarise(K=length(unique(id)),
                                                                                       SENS_pct=round(100*mean(desqx1),2))
CONCtabNk_2 <- datG1 %>% filter(mdes1==0) %>% group_by(Nk_collapse5plus) %>% summarise(K=length(unique(id)),
                                                                                       SPEC_pct=round(100*mean(desqx0),2))


tabNk <- cbind(tabNk_1,tabNk_2[,3],VALIDtabNk_1[,-1],VALIDtabNk_2[,3:4],CONCtabNk_1[,3],CONCtabNk_2[,3])
colnames(tabNk) <- c("N_k","No. G1s","% ADHD","% DES","Valid % ADHD","Valid % DES","Valid % MDES","Sensitivity","Specificity")
tabNk[is.na(tabNk)] <- "--"
xtable(tabNk,include.rownames = FALSE)
xtable(tabNk[,-(6:8)],include.rownames = FALSE)
write.table(tabNk,"~/results/Mis_tabNk5.txt")


################################################
## 2 Table of overall, mothers, and non-mothers
## setting rate=TRUE gives the rates of emptiness in each covariate category ## This one probably more useful
## setting rate=FALSE just gives overall cell percentages
nm_table <- function(EDAvar="N0",cats,names=NA,df=datG1,G1s=TRUE,digits=0,rate=FALSE){
  if(length(names)<length(cats)){
    rnames <- c("Overall",cats)
  }else{
    rnames <- c("Overall",names)
  }

  ## column of variable of interest
  EDAvec <- df[,EDAvar]

  ## first get overall number of G1s
  total <- c(nrow(df))
  EDA_1 <- c(sum(EDAvec!=1)) ## N>=1
  EDA_0 <- c(sum(EDAvec==1)) ##N=0


  ## get category specific numbers
  for(var in cats){
    varvec <- df[,var]
    total <- c(total,sum(varvec==1))
    EDA_1 <- c(EDA_1,sum((EDAvec!=1) & (varvec==1)))
    EDA_0 <- c(EDA_0,sum((EDAvec==1) & (varvec==1)))
  }


  if(rate==FALSE){
    # percentage
    total_pct <- round(100*total/nrow(df),digits)
    EDA_1_pct <- round(100*EDA_1/nrow(df),digits)
    EDA_0_pct <- round(100*EDA_0/nrow(df),digits)
    ## make table
    tab <- cbind(total,total_pct,EDA_1,EDA_1_pct,EDA_0,EDA_0_pct)
    colnames(tab) <- c("Overall No.","Overall %","N>=1 No.","N>=1 %","N=0 No.","N=0 %")

  }else {
    # percentage
    EDA_1_pct <- round(100*EDA_1/sum(EDAvec!=1),digits)
    EDA_0_pct <- round(100*EDA_0/sum(EDAvec==1),digits)
    ## make table
    tab <- cbind(total,EDA_1,EDA_1_pct,EDA_0,EDA_0_pct)
    colnames(tab) <- c("Overall No.","N>=1 No.","N>=1 Rate %","N=0 No.","N=0 Rate %")
  }


  #
  rownames(tab) <- rnames
  return(tab)

}

nmcats <- c("desqx1","desqx0",
            "mdes1","mdes0",
           "msmk2","msmk0",
           "yob89_4650","yob89_5155","yob89_5660","yob89_61plus",
           "raceWhite","raceAA","raceAsian","raceOther",
           "hisp89","nonhisp",
           "momed2","momed3","momed4","momed5","momed6",
           "daded2","daded3","daded4","daded5","daded6",
           "mprof05","mexec05","msale05","mmech05","mmach05","mserv05","mlabo05","mfarm05","mmili05","mdk05",
           "fprof05","fexec05","fsale05","fmech05","fmach05","fserv05","flabo05","ffarm05","fmili05","fdk05",
           "phome05","nophome")

nmnames <- c("DES Yes","DES No",
             "MDES Yes","MDES No",
              "Smoke Yes","Smoke No",
              "YOB 45-60","51-55","56-60","61+",
              "White","African American","Asian","Other",
              "Hispanic","Nonhisp",
              "MomEd 2","MomEd 3","MomEd 4 ","MomEd 5 ","MomEd 6 ",
              "DadEd 2","DadEd 3","DadEd 4 ","DadEd 5 ","DadEd 6 ",
              "mprof05","mexec05","msale05","mmech05","mmach05","mserv05","mlabo05","mfarm05","mmili05","mdk05",
              "fprof05","fexec05","fsale05","fmech05","fmach05","fserv05","flabo05","ffarm05","fmili05","fdk05",
              "PHOME","PHOME No")

## table of overall/mothers/nonmothers (cell percentages)
nonmothers_tab <- nm_table("N0",cats=nmcats,names=nmnames,df=datG1,digits=1)
xtable(nonmothers_tab,digits=c(0,0,1,0,1,0,1))
write.table(nonmothers_tab,file="~/results/Mis_nonmothers_tab.txt")
## table of overall/mothers/nonmothers (rate of emptiness) ## this one is probably more meaninful if emptiness is the focus
nonmothers_rate_tab <- nm_table("N0",cats=nmcats,names=nmnames,df=datG1,digits=1,rate=TRUE)
xtable(nonmothers_rate_tab,digits=c(0,0,0,1,0,1))
write.table(nonmothers_rate_tab,file="~/results/Mis_nonmothers_rate_tab.txt")


###########################
#### 3 EDA Table
## (less complete but still interesting plots of DES and ADHD rates by covariate categories)
EDA_table <- function(EDAvar,cats,names=NA,df=datG1,G1s=TRUE,digits=0){
  if(length(names)<length(cats)){
    rnames <- c("Overall",cats)
  }else{
    rnames <- c("Overall",names)
  }

  ## column of variable of interest
  EDAvec <- df[,EDAvar]

  ## first get overall number of G1s
  EDA_K <- c(sum(EDAvec))
  total_K <- c(nrow(df))

  ## get category specific numbers
  for(var in cats){
    total_K <- c(total_K,sum(df[,var]==1))
    EDA_K <- c(EDA_K,sum(EDAvec[df[,var]==1]))
  }

  # percentage
  EDA_pct <- round(100*EDA_K/total_K,digits)

  ## make table
  tab <- cbind(total_K,EDA_K,EDA_pct)
  rownames(tab) <- rnames
  if(G1s==TRUE){
    colnames(tab) <- c("No. of G1s","No. Exposed","%")
  }else{
    colnames(tab) <- c("No. of G2s","No. Exposed","%")
  }


  return(tab)

}

DEScats<-c("N0","N1","N2","N3","N4plus",
           "yob89_4650","yob89_5155","yob89_5660","yob89_61plus",
           "scand89","ocauc89","afric89","hisp89","asian89","oanc89",
           "mprof05","mexec05","msale05","mmech05","mmach05","mserv05","mlabo05","mfarm05","mmili05","mdk05",
           "fprof05","fexec05","fsale05","fmech05","fmach05","fserv05","flabo05","ffarm05","fmili05","fdk05",
           "momed2","momed3","momed4","momed5","momed6",
           "daded2","daded3","daded4","daded5","daded6",
           "msmk2","msmk3","phome05")
DESnames <- c("Nk=0","1","2","3","4+",
              "YOB 45-60","51-55","56-60","61+",
              "Scand","Ocauc","Afric","Hisp","Asian","Oanc",
              "mprof05","mexec05","msale05","mmech05","mmach05","mserv05","mlabo05","mfarm05","mmili05","mdk05",
              "fprof05","fexec05","fsale05","fmech05","fmach05","fserv05","flabo05","ffarm05","fmili05","fdk05",
              "MomEd 2","MomEd 3","MomEd 4 ","MomEd 5 ","MomEd 6 ",
              "DadEd 2","DadEd 3","DadEd 4 ","DadEd 5 ","DadEd 6 ",
              "MSMK 2","MSMK 3","PHOME")

## DES table (G1 level)
DES_tab <- EDA_table("desqx1",cats=DEScats,names=DESnames,df=datG1,G1s=TRUE,digits=2)
xtable(DES_tab)
write.table(DES_tab,file="~/results/Mis_DES_tab.txt")
MDES_tab <- EDA_table("mdes1",cats=DEScats,names=DESnames,df=datG1,G1s=TRUE,digits=2)
xtable(MDES_tab)
write.table(MDES_tab,file="~/results/Mis_MDES_tab.txt")

## ADHD table (G2 level)
ADHDcats <- c("desqx1","mdes1",DEScats)
ADHDnames <- c("DES","MDES",DESnames)
ADHD_tab <- EDA_table("adhd",cats=ADHDcats,names=ADHDnames,df=dat[dat$N0==0,],G1s=FALSE,digits=2)
xtable(ADHD_tab)
write.table(ADHD_tab,file="~/results/Mis_ADHD_tab.txt")


####################################
## 4 Prelim Regressions

###################################
## function to clean up reg results
clean_reg <- function(g,glmm=FALSE,digits=2){
  if(glmm==FALSE){
    coef <- exp(g$coef)
    # conf <- exp(confint(g))
    SE <- sqrt(diag(summary(g)$cov.scaled))
    conf <- exp(cbind(g$coef-1.96*SE,g$coef+1.96*SE))
  }else{
    coef <- exp(summary(g)$coefficients[,1])
    SE <- summary(g)$coefficients[,2]
    conf <- exp(cbind(summary(g)$coefficients[,1]-1.96*SE,summary(g)$coefficients[,1]+1.96*SE))
  }
  coef <- round(coef,digits)
  conf <- round(conf,digits)
  conf <- paste("(",conf[,1],", ",conf[,2],")",sep="")
  return(cbind(coef,conf))
}


#############################################
## preliminary size regressions
g_pois <- glm(totalkids~desqx1+msmk2+adhd+yob89_5155+yob89_5660+yob89_61plus,data=datG1,family=poisson)
g_nonzero <- glm((1-N0)~desqx1+msmk2+adhd+yob89_5155+yob89_5660+yob89_61plus,data=datG1,family=poisson)

results_Nkreg <- cbind(clean_reg(g_pois),clean_reg(g_nonzero))
colnames(results_Nkreg) <- c("RR","CI","OR","CI")
#rowname(results_Nkreg) <-
xtable(results_Nkreg)


#
# ######## G1 level plot
# plot_ADHD_Nk <- ggplot(datG1,aes(x=totalkids,y=adhd,color=DES))+
#   geom_point()+
#   xlab('Nk') +ylab ('Effect Estimate')
#
#
#









# ### Table 1
# vars_binary <- colnames(datG1)[!colnames(datG1)%in%c("id","adhd","yob89","yob89sq","totalkids")]
# tab1_G1 <- NULL
# for(var in vars_binary){
#   tab1_G1 <- rbind(tab1_G1,c(sum(datG1[,var],na.rm=TRUE),100*mean(datG1[,var],na.rm=TRUE)))
# }
# tab1_G1 <- rbind(c(nrow(datG1),100),tab1_G1)
# rownames(tab1_G1) <- c("Total",vars_binary)
# colnames(tab1_G1) <- c("K","%")
#
# write.table(tab1_G1,"/udd/n2glm/results/tab1_G1.txt")


# ##########################################################
# ## G2 level data table 1
# tab1_G2 <- dat %>% group_by(totalkids) %>% summarise(K=length(unique(id)),
#                                                   Out_pct=100*mean(adhd),
#                                                   Exp_pct=100*mean(desqx1))
# colnames(tab1_G2) <- c("Nk","No. G1s","% ADHD","% DES")
#
# write.table(tab1_G2,"/udd/n2glm/results/tab1_G2.txt")

#
#
#
# ##### EDA for DES
# DES_K <- c(sum(datG1$desqx1),
#             sum(datG1$desqx1[datG1$N0==1]),
#             sum(datG1$desqx1[datG1$N1==1]),
#             sum(datG1$desqx1[datG1$N2==1]),
#             sum(datG1$desqx1[datG1$N3==1]),
#             sum(datG1$desqx1[datG1$N4plus==1]),
#             sum(datG1$desqx1[datG1$yob89_4650==1]),
#             sum(datG1$desqx1[datG1$yob89_5155==1]),
#             sum(datG1$desqx1[datG1$yob89_5660==1]),
#             sum(datG1$desqx1[datG1$yob89_61plus==1]),
#            sum(datG1$desqx1[datG1$scand89==1]),
#            sum(datG1$desqx1[datG1$ocauc89==1]),
#            sum(datG1$desqx1[datG1$afric89==1]),
#            sum(datG1$desqx1[datG1$hisp89==1]),
#            sum(datG1$desqx1[datG1$asian89==1]),
#            sum(datG1$desqx1[datG1$oanc89==1]),
#            sum(datG1$desqx1[datG1$momed2==1]),
#            sum(datG1$desqx1[datG1$momed3==1]),
#            sum(datG1$desqx1[datG1$momed4==1]),
#            sum(datG1$desqx1[datG1$momed5==1]),
#            sum(datG1$desqx1[datG1$momed6==1]),
#            sum(datG1$desqx1[datG1$daded2==1]),
#            sum(datG1$desqx1[datG1$daded3==1]),
#            sum(datG1$desqx1[datG1$daded4==1]),
#            sum(datG1$desqx1[datG1$daded5==1]),
#            sum(datG1$desqx1[datG1$daded6==1]),
#            sum(datG1$desqx1[datG1$msmk2==1]),
#            sum(datG1$desqx1[datG1$msmk3==1]),
#            sum(datG1$desqx1[datG1$phome05==1]),
#
#
#           )
# tab_DES <- cbind(DES_K,DES_K/nrow(datG1))
# colnames(tab_DES) <- c("No. G1s","%")
# rownames(tab_DES) <- c("Overall",
#                        "Nk=0","1","2","3","4+",
#                        "YOB 45-60","51-55","56-60","61+",
#                        "Scand","Ocauc","Afric","Hisp","Asian","Oanc",
#                        "mprof05","mexec05","msale05","mmech05","mmach05","mserv05","mlabo05","mfarm05","mmili05","mdk05",
#                        "fprof05","fexec05","fsale05","fmech05","fmach05","fserv05","flabo05","ffarm05","fmili05","fdk05",
#                        "MomEd 2","MomEd 3","MomEd 4 ","MomEd 5 ","MomEd 6 ",
#                        "DadEd 2","DadEd 3","DadEd 4 ","DadEd 5 ","DadEd 6 ",
#                        "MSMK 2","MSMK 3","PHOME")
#
#
#
#
#
#















## tibble method
# datG1 %>% select(-c(id,adhd)) %>%
#   gather(desqx1:msmk3,key="var",value="value") %>%
#   group_by(var) %>%
#   summarise_all(funs(K = sum,pctK = 100*mean))
#                                                                                                                          pctK = 100*mean))
#
# datG1 %>% select(-c(id,adhd)) %>% summarise_each(funs(K = sum,
#                                                     pctK = 100*mean))



# dim(dat)
# [1] 106198     45

# names(dat)
# [1] "id"        "desqx1"    "yob89"     "scand89"   "ocauc89"   "afric89"
# [7] "hisp89"    "asian89"   "oanc89"    "mprof05"   "mexec05"   "msale05"
# [13] "mmech05"   "mmach05"   "mserv05"   "mlabo05"   "mfarm05"   "mmili05"
# [19] "mdk05"     "fprof05"   "fexec05"   "fsale05"   "fmech05"   "fmach05"
# [25] "fserv05"   "flabo05"   "ffarm05"   "fmili05"   "fdk05"     "phome05"
# [31] "adhd"      "totalkids" "yob89sq"   "momed2"    "momed3"    "momed4"
# [37] "momed5"    "momed6"    "daded2"    "daded3"    "daded4"    "daded5"
# [43] "daded6"    "msmk2"     "msmk3"
