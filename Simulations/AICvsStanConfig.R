# AICvsStan config file prep
library(here)


#numTrials,n=100,K=50, num_strong=5, prob_cycle=0, trend_fraction=0,freq = 1, sigma = 0.5,probWeakCorr=0,numStrongCorr=0,strongSelf=F,corrLevel=0

# first columns
numTrials=rep(50,72)
n=rep(100,72)
K=rep(50,72)
numstrong=rep(5,72)
prob_cycle=rep(0.5,72)
trend_fraction=rep(0.5,72)
freq=rep(1,72)
sigma=rep(0.5,72)

# last columns
probWeakCorr=rep(0.8,72)
numStrongCorr=rep(3,72)
strongSelf=rep(T,72)
corrLevel=rep(0.9,72)

# de-correlation
corrChange=c(rep(T,36),rep(F,36))
propChange=rep(0.75,72)
changeSize=rep(0,72)
changeTimeVar=rep(0,72)

# combine them
configx=as.data.frame(matrix(data=c(numTrials,n,K,numstrong,prob_cycle,trend_fraction,freq,sigma,probWeakCorr,numStrongCorr,strongSelf,corrLevel,corrChange,propChange,changeSize,changeTimeVar),nrow=72,ncol=16,byrow=F))
configx[,11]<-as.logical(configx[,11])
configx[,13]<-as.logical(configx[,13])



configx=cbind(1:72,configx)

configx[,18]=paste0("ACPtrial",1:72)

colnames(configx)=c("ArrayTaskID","numTrials","n","K","num_strong","prob_cycle","trend_fraction","freq","sigma","probWeakCorr","numStrongCorr","strongSelf","corrLevel","corrChange","propChange","changeSize","changeTimeVar","nameID")
rownames(configx)=NULL

write.table(configx,file=here("Simulations/AICvsStanConfig.txt"),row.names = F,quote=F)
