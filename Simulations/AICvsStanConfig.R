# AICvsStan config file prep
library(here)

#1. short TS, without seasons or trend
trial1=c(1000,100,50,5,0,0,1,0.5,0,0,0.0,"ACPtrial1")

#2. with seasons and trend
trial2=c(1000,100,50,5,0.5,0.5,1,0.5,0,0,0.0,"ACPtrial2")

#3. with seasons and trend and medium correlation
trial3=c(1000,100,50,5,0.5,0.5,1,0.5,0.2,0.4,0.4,"ACPtrial3")

#4. with season and trend and larger correlation
trial4=c(1000,100,50,5,0.5,0.5,1,0.5,0.2,0.4,0.8,"ACPtrial4")

configx=as.data.frame(matrix(data=c(trial1,trial2,trial3,trial4),nrow=4,ncol=12,byrow=T))

configx=cbind(1:4,configx)

colnames(configx)=c("ArrayTaskID","numTrials","n","K","num_strong","prob_cycle","trend_fraction","freq","sigma","probWeakCorr","probStrongCorr","corrLevel","nameID")
rownames(configx)=NULL

write.table(configx,file=here("Simulations/AICvsStanConfig.txt"))
