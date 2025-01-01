# AICvsStan config file prep
library(here)

# Dec 31- adding in what would happen for different values of correlation (besides 0.9) to match with challenge 1
#numTrials,n=100,K=50, num_strong=5, prob_cycle=0, trend_fraction=0,freq = 1, sigma = 0.5,probWeakCorr=0,numStrongCorr=0,strongSelf=F,corrLevel=0
total_trials=120

# first columns
numTrials=rep(50,total_trials)
n=rep(100,total_trials)
K=rep(50,total_trials)
numstrong=rep(5,total_trials)
prob_cycle=rep(0.5,total_trials)
trend_fraction=rep(0.5,total_trials)
freq=rep(1,total_trials)
sigma=rep(0.5,total_trials)

# last columns
probWeakCorr=rep(0.8,total_trials)
numStrongCorr=rep(3,total_trials)
strongSelf=rep(T,total_trials)
corrLevel=c(rep(0.9,total_trials/3),rep(0.5,total_trials/3),rep(0.1,total_trials/3))

# de-correlation
corrChange=rep(c(rep(T,total_trials/6),rep(F,total_trials/6)),3)
propChange=rep(0.75,total_trials)
changeSize=rep(0,total_trials)
changeTimeVar=rep(0,total_trials)

# combine them
configx=as.data.frame(matrix(data=c(numTrials,n,K,numstrong,prob_cycle,trend_fraction,freq,sigma,probWeakCorr,numStrongCorr,strongSelf,corrLevel,corrChange,propChange,changeSize,changeTimeVar),nrow=total_trials,ncol=16,byrow=F))
configx[,11]<-as.logical(configx[,11])
configx[,13]<-as.logical(configx[,13])



configx=cbind(1:total_trials,configx)

configx[,18]=paste0("decorr",1:total_trials)

colnames(configx)=c("ArrayTaskID","numTrials","n","K","num_strong","prob_cycle","trend_fraction","freq","sigma","probWeakCorr","numStrongCorr","strongSelf","corrLevel","corrChange","propChange","changeSize","changeTimeVar","nameID")
rownames(configx)=NULL

write.table(configx,file=here("Simulations/DecorrConfig.txt"),row.names = F,quote=F)
