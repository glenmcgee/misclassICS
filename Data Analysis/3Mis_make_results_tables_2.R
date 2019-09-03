

###########################
##     Report Results    ##
###########################
## 05/29/2019

library(xtable)

## set number of digits for final table
dig <- 2

## load data
## Server:
ests <- data.frame(read.table(file="Mis_ests_2.txt",header=T,row.names=1))
SEs <- data.frame(read.table(file="Mis_SEs_2.txt",header=T,row.names=1))


## create intervals
CIs_lo <- ests-1.96*SEs
CIs_hi <- ests+1.96*SEs

## exponentiate
ests[1:(nrow(ests)-4),] <- exp(ests[1:(nrow(ests)-4),])
CIs_lo[1:(nrow(CIs_lo)-4),] <- exp(CIs_lo[1:(nrow(CIs_lo)-4),])
CIs_hi[1:(nrow(CIs_hi)-4),] <- exp(CIs_hi[1:(nrow(CIs_hi)-4),])

## make full results table
tab_res <- c()
for(cc in 1:ncol(ests)){ ## loop over analyses
  tab_res <- cbind(tab_res,
                   round(ests[,cc],dig), ## estimates
                   paste0("(",round(CIs_lo[,cc],dig),",",round(CIs_hi[,cc],dig),")") ## CIs
                   )
}

colnames(tab_res) <- paste0(rep(colnames(ests),each=2),
                            rep(c("_Est","_CI"),times=ncol(ests)) )
rownames(tab_res) <- c("a0","a1 DES","a2 msmk","a3 yob5155","a4 yob5660","a5 yob61plus",
                       "b0","b1 DES","b2 msmk","b3 yob5155","b4 yob5660","b5 yob61plus",
                       "Sigma0","Sigma1",
                       "Gamma0","Gamma1" )


## NAs
tab_res[tab_res=="(NA,NA)"] <- ""
tab_res[is.na(tab_res)] <- ""



## marginal slopes results
xtable(tab_res)
write.table(tab_res,"Mis_results_2.txt")


