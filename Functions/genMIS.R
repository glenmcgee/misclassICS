########################################
## Function to misclassify exposure   ##
########################################
## Updated: 02/25/2019
  ## Function takes result of genICS() as variable
  ## Misclassifies exposure (X1k-->W1k) based on the following model:
    ## logit(E[1(W1k=X1k)])=(1-X1k)*(nu00+nu01*(Nk-2))+
    ##                           X1k*(nu10+nu11*(Nk-2))
  ## Adds W1k and W1ki to datasets as in genICS()

## set functions
expit <- function(y){  return(exp(y)/(1+exp(y)))}
logit <- function(y){  return(log(y/(1-y)))}


genMIS <- function(datICS, ## object produced by genICS()
                   nu00=logit(0.90), ## baseline specificity (for Nk=1)
                   nu10=logit(0.90), ## baseline sensitivity
                   nu01=0, ## logOR for Nk on specificity
                   nu11=0 ## logOR for Nk on sensitivity
                   ){

  ## get data
  data <- datICS$data ## unit-level data
  data_clstlvl <- datICS$data_clstlvl ## cluster-level data
    Nk <- data_clstlvl$Nk ## cluster sizes
    X1k <- data_clstlvl$X1k ## cluster-level exposure
  bk <- datICS$true_randeff ## random effects

  #########################################
  ## Simulate Exposure Misclassification ##
  #########################################

  ## concordance between W1k and X1k
  IW1k_X1k <- rbinom(length(X1k),1,expit((1-X1k)*(nu00+nu01*(Nk-2))+X1k*(nu10+nu11*(Nk-2))))

  ## cluster-level misclassified exposure
  W1k <- IW1k_X1k*X1k+(1-IW1k_X1k)*(1-X1k)

  ## obs-level
  W1ki <- rep(W1k,times=Nk)


  #################
  ## Output Data ##
  #################
  data_clstlvl$W1k <- W1k
  data$W1ki <- W1ki

  return(list(data=data,data_clstlvl=data_clstlvl,true_randeff=bk))
}

