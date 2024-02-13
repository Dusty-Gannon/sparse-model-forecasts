# AICvsStan config file prep
library(here)


#numTrials,n=100,K=50, num_strong=5, prob_cycle=0, trend_fraction=0,freq = 1, sigma = 0.5,probWeakCorr=0,numStrongCorr=0,strongSelf=F,corrLevel=0

# first columns
numTrials=rep(50,36)
n=rep(100,36)
K=rep(50,36)
numstrong=rep(5,36)
prob_cycle=rep(0.5,36)
trend_fraction=rep(0.5,36)
freq=rep(1,36)
sigma=rep(0.5,36)

# last columns
probWeakCorr=rep(0.2,36)
numStrongCorr=rep(3,36)
strongSelf=rep(T,36)
corrLevel=rep(0.7,36)

# de-correlation
corrChange=rep(T,36)
propChange=c(rep(0.25,12),rep(0.5,12),rep(0.75,12))
changeSize=rep(c(rep(0.5,6),rep(0.2,6)),3)
changeTimeVar=rep(c(rep(1,3),rep(5,3)),6)

# combine them
configx=as.data.frame(matrix(data=c(numTrials,n,K,numstrong,prob_cycle,trend_fraction,freq,sigma,probWeakCorr,numStrongCorr,strongSelf,corrLevel,corrChange,propChange,changeSize,changeTimeVar),nrow=36,ncol=16,byrow=F))
configx[,11]<-as.logical(configx[,11])
configx[,13]<-as.logical(configx[,13])



configx=cbind(1:36,configx)

configx[,18]=paste0("ACPtrial",1:36)

colnames(configx)=c("ArrayTaskID","numTrials","n","K","num_strong","prob_cycle","trend_fraction","freq","sigma","probWeakCorr","numStrongCorr","strongSelf","corrLevel","corrChange","propChange","changeSize","changeTimeVar","nameID")
rownames(configx)=NULL

write.table(configx,file=here("Simulations/AICvsStanConfig.txt"),row.names = F,quote=F)
