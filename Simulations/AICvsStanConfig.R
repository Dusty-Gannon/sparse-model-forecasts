# AICvsStan config file prep
library(here)

#1. short TS, without seasons, trend, or correlation
trial1=c(100,100,50,5,0,0,1,0.5,0,0,0.0)

#2. with seasons and trend, but no correlation
trial2=c(100,100,50,5,0.5,0.5,1,0.5,0,0,0.0)

#3. with season and trend and correlation
trial3=c(100,100,50,5,0.5,0.5,1,0.5,0.2,0.4,0.8)

# number of trials per category
n=50

configx=as.data.frame(matrix(data=c(rep(trial1,n),rep(trial2,n),rep(trial3,n)),nrow=3*n,ncol=11,byrow=T))

configx=cbind(1:(3*n),configx)

configx[,13]=paste0("ACPtrial",1:(3*n))

colnames(configx)=c("ArrayTaskID","numTrials","n","K","num_strong","prob_cycle","trend_fraction","freq","sigma","probWeakCorr","probStrongCorr","corrLevel","nameID")
rownames(configx)=NULL

write.table(configx,file=here("Simulations/AICvsStanConfig.txt"),row.names = F,quote=F)
