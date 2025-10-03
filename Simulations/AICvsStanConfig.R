# AICvsStan config file prep
library(here)


#numTrials,n=100,K=50, num_strong=5, prob_cycle=0, trend_fraction=0,freq = 1, sigma = 0.5,probWeakCorr=0,numStrongCorr=0,strongSelf=F,corrLevel=0

# for 300 runs per case, we have 3 cases, so need 24 runs of 50 trails each
totalruns=24
# first columns everything held constant
numTrials=rep(50,totalruns) # number of trials per run (split up for parallel processing)
n=rep(100,totalruns) # length of time series
K=rep(50,totalruns) # number of betas/covariates
numstrong=rep(5,totalruns) # number of strong predictors
prob_cycle=rep(0.5,totalruns) # probability of a covariate having a yearly cycle
trend_fraction=rep(0.5,totalruns) # probability of a covariate having a trend
freq=rep(1,totalruns) # frequency of observations per time step
sigma=rep(0.5,totalruns) #Standard deviation of the noise component

# last columns
probWeakCorr=rep(0.5,totalruns) # probability of a weak covariate correlating with a strong covariate
numStrongCorr=rep(3,totalruns) # number of strong covariates that correlate with weak covariates
strongSelf=rep(F,totalruns) # whether strong covariates correlate with other strong covariates
corrLevel=rep(c(0.1,0.5,0.9),8) # correlation levels

# no de-correlation
corrChange=rep(F,totalruns)
propChange=rep(0,totalruns)
changeSize=rep(0,totalruns)
changeTimeVar=rep(0,totalruns)

# combine them
configx=as.data.frame(matrix(data=c(numTrials,n,K,numstrong,prob_cycle,trend_fraction,freq,sigma,probWeakCorr,numStrongCorr,strongSelf,corrLevel,corrChange,propChange,changeSize,changeTimeVar),nrow=totalruns,ncol=16,byrow=F))
configx[,11]<-as.logical(configx[,11])
configx[,13]<-as.logical(configx[,13])
# number of trials per category
configx=cbind(1:totalruns,configx)

configx[,18]=paste0("AICSTAN300_",1:totalruns)

colnames(configx)=c("ArrayTaskID","numTrials","n","K","num_strong","prob_cycle","trend_fraction","freq","sigma","probWeakCorr","numStrongCorr","strongSelf","corrLevel","corrChange","propChange","changeSize","changeTimeVar","nameID")
rownames(configx)=NULL

write.table(configx,file=here("Simulations/AICvsStanConfig.txt"),row.names = F,quote=F)
