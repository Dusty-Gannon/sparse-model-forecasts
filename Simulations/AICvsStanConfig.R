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
probWeakCorr=c(rep(0.2,12),rep(0.5,12),rep(0.8,12))
numStrongCorr=rep(c(rep(1,6),rep(3,6)),3)
strongSelf=rep(c(rep(T,3),rep(F,3)),6)
corrLevel=rep(c(0.5,0.7,0.9),12)

# no de-correlation
corrChange=rep(F,36)
propChange=rep(0,36)
changeSize=rep(0,36)
changeTimeVar=rep(0,36)

# combine them
configx=as.data.frame(matrix(data=c(numTrials,n,K,numstrong,prob_cycle,trend_fraction,freq,sigma,probWeakCorr,numStrongCorr,strongSelf,corrLevel,corrChange,propChange,changeSize,changeTimeVar),nrow=36,ncol=16,byrow=F))
configx[,11]<-as.logical(configx[,11])
configx[,13]<-as.logical(configx[,13])
# number of trials per category
nx=3

configx=rbind(configx,configx,configx)


configx=cbind(1:(36*nx),configx)

configx[,18]=paste0("AICSTANredo",1:(36*nx))

colnames(configx)=c("ArrayTaskID","numTrials","n","K","num_strong","prob_cycle","trend_fraction","freq","sigma","probWeakCorr","numStrongCorr","strongSelf","corrLevel","corrChange","propChange","changeSize","changeTimeVar","nameID")
rownames(configx)=NULL

write.table(configx,file=here("Simulations/AICvsStanConfig.txt"),row.names = F,quote=F)
