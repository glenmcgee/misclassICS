#########################
##    LOAD NHS DATA    ##
#########################
## 05/07/2019
## Needs to be updated after consulting with Marianthi about data coding issues

## load packages
library(tidyr)
library(dplyr)

######################
## Load and clean data

dat <- read.csv("~/data/GlenData.csv")

## clean coding
dat$validID <- 0
dat$validID[dat$validated==1 & !is.na(dat$mdes)] <- 1

## clean coding (1,2) --> (0,1)
dat$mdes1[dat$mdes==1] <- 0 ## VALIDATED DES EXPOSURE
dat$mdes1[dat$mdes==2] <- 1
dat <- dat[,!(colnames(dat) %in% c("mdes") )] ## discard old coding

## clean data coding ("0","yes") --> (0,1)
vars_clean <- c("scand89","ocauc89","afric89","hisp89","asian89","oanc89","phome05")
for (var in vars_clean){
  dat[,var] <- as.numeric(dat[,var]=="yes")
}

## make indicators for strata of non-binary vars
dat$N0 <- as.numeric(dat$totalkids==0)
dat$N1 <- as.numeric(dat$totalkids==1)
dat$N2 <- as.numeric(dat$totalkids==2)
dat$N3 <- as.numeric(dat$totalkids==3)
dat$N4plus <- as.numeric(dat$totalkids>=4)

## categorize year of birth
dat$yob89_4650 <- as.numeric(dat$yob89<=50)
dat$yob89_5155 <- as.numeric(dat$yob89>50 & dat$yob89<=55)
dat$yob89_5660 <- as.numeric(dat$yob89>55 & dat$yob89<=60)
dat$yob89_61plus <- as.numeric(dat$yob89>60)

# complements of categories for building tables
dat$desqx0 <- 1-dat$desqx1
dat$mdes0 <- 1-dat$mdes1
dat$msmk0 <- 1-dat$msmk2 ## assuming msmk2 is smoking... not sure what msmk3 is? --> ask MARIANTHI
dat$nophome <- 1-dat$phome05
dat$raceWhite <- as.numeric((dat$scand89==1) | (dat$ocauc89==1)) ## assuming these are white --> what does it mean to have NO race??
dat$raceAA <- dat$afric89
dat$raceAsian <- dat$asian89
dat$raceOther <- dat$oanc89
dat$nonhisp <- 1-dat$hisp89


#### FIX DATA CODING ERROR -- ASK MARIANTHI. is it adoptions?
dat$adhd[!is.na(dat$adhd) & dat$totalkids==0] <- NA

### renumber clusters (for use in JMMMICSpack)
dat$id2 <- as.numeric(factor(dat$id))


##########################################################
## G1 level dataset
datG1 <- dat %>% group_by(id) %>% summarise_all(mean)
datG1 <- data.frame(datG1)

## make variable for at least one adhd diagnosis
datG1$adhd_1 <- as.numeric(datG1$adhd>0)


##########################################################
## G2 level dataset
datG2 <- dat[dat$totalkids>0,]


##########################################################
## G1 and G2 datasets -- excluding empty clusters (no 0s)
dat_no0 <- dat[dat$totalkids>0,]
dat_no0$id2 <- as.numeric(factor(dat_no0$id))
datG1_no0 <- dat_no0 %>% group_by(id) %>% summarise_all(mean)
datG1_no0 <- data.frame(datG1_no0)
datG2_no0 <- dat_no0


