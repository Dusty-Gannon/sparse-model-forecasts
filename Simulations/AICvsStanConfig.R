# AICvsStan config file prep
library(here)

#1. short TS, without seasons or trend
trial1=c(100,100,50,5,0,0,1,0.5,0,0,0.0)

#2. with seasons and trend
trial2=c(100,100,50,5,0.5,0.5,1,0.5,0,0,0.0)

#3. with seasons and trend and medium correlation
trial3=c(100,100,50,5,0.5,0.5,1,0.5,0.2,0.4,0.4)

#4. with season and trend and larger correlation
trial4=c(100,100,50,5,0.5,0.5,1,0.5,0.2,0.4,0.8)

configx=as.data.frame(matrix(data=c(rep(trial1,10),rep(trial2,10),rep(trial3,10),rep(trial4,10)),nrow=40,ncol=11,byrow=T))

configx=cbind(1:4,configx)

configx[,13]=paste0("ACPtrial",1:4)

colnames(configx)=c("ArrayTaskID","numTrials","n","K","num_strong","prob_cycle","trend_fraction","freq","sigma","probWeakCorr","probStrongCorr","corrLevel","nameID")
rownames(configx)=NULL

write.table(configx,file=here("Simulations/AICvsStanConfig.txt"),row.names = F)
