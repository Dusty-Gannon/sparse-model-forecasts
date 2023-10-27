# get all the functions and stuff that we need
source("Simulations/AICvsStanRMSE.R")

# Arguments we need for the get time series function...
# numTrials,n=100,K=50, num_strong=5, prob_cycle=0, trend_fraction=0,freq = 1, sigma = 0.5,correlated=F,rateCorr=2

# this gets the arguments from the shell script
args<-commandArgs(TRUE)

# args= c(20, 100, 60, 5, 0.2, 0.2, 1, 0.5, F, 2)
# arguments come in as strings
numTrials=as.numeric(args[1])
n=as.numeric(args[2])
K=as.numeric(args[3])
num_strong=as.numeric(args[4])
prob_cycle=as.numeric(args[5])
trend_fraction=as.numeric(args[6])
freq=as.numeric(args[7])
sigma=as.numeric(args[8])
correlated=as.logical(args[9])
rateCorr=as.numeric(args[10])
nameID=args[11]

ts1=getTS(numTrials = numTrials, n = n, K = K, num_strong = num_strong,prob_cycle = prob_cycle, trend_fraction = trend_fraction,freq=freq,sigma=sigma,correlated=correlated,rateCorr=rateCorr) # make timeseries


#trials
ts1clean=lapply(ts1,FUN=function(x) cleanTS(x)) # clean them
ts1test=lapply(ts1clean,FUN=function(x) splitTS(x,set="test",n=n,nfit=0.6*n)) #testing set
ts1train=lapply(ts1clean,FUN=function(x) splitTS(x,set="train",n=n,nfit=0.6*n)) #training set
AICmodlist=lapply(ts1train,FUN=function(x) AICselect(x)) # do model selection AIC
RMSEAIClist=mapply(function(x,y) RMSE_AIC(x,y), AICmodlist,ts1test) # get RMSE AIC
STANmodlist=mapply(function(x,y) STANselect(x,y,nfit=60,n=n,K=K),ts1train,ts1test) # do model selection STAN
STANpredlist=lapply(STANmodlist,FUN=function(x) STANgetpredict(x)) # get STAN predictions for y
STANy=lapply(ts1test,"[",,1)


RMSESTANlistraw=mapply(function(x,y) RMSE_bayes(x,y),STANy,STANpredlist) # get RMSE STAN
RMSESTANlist=colMeans(RMSESTANlistraw)


AICconfusion=mapply(function(x,y) AICconfusionRates(x,y),ts1,AICmodlist)
mean(AICconfusion[1,])# TPR is 0.992
mean(AICconfusion[2,])# TNR is 0.498

STANbetalist=lapply(STANmodlist, FUN=function(x) STANbetapost(x)) # get summary of stan predictions for beta
# ts_betas=lapply(ts1....)  # How do I extract the beta values from the timeseries list?
STANconfusion=mapply(function(x,y) STANconfusionRates(x,y), STANbetalist, ts_betas) #returns a matrix of false pos, false neg,
## NB: this is VERY close to Alice C's function, not 100% sure it's the best fit for what we want here.

#hist(RMSESTANlist,xlim=c(0,2.5),ylim=c(0,50),breaks=seq(0,2.5,by=0.1),xlab="RMSE",main="RMSE using horseshoe priors")
#hist(RMSEAIClist,xlim=c(0,2.5),ylim=c(0,50),breaks=seq(0,2.5,by=0.1),xlab="RMSE",main="RMSE using stepwise AIC")

# Output- want to save the confusion values, and save the RMSE values
# Put each into dataframes and save as .csv?

# Save RMSESTANlist
# Save RMSEAIClist
saveRDS(RMSESTANlist,file=paste0("RMSESTAN_",nameID,".rds")) # customize these, allow args input?
saveRDS(RMSEAIClist,file="RMSEAIC_",nameID,".rds")

# Save the confusion matrices
saveRDS(AICconfusion,file="AICconfusion_",nameID,".rds")
saveRDS(STANconfusion,file="STANconfusion_",nameID,".rds")

# Save the lists of models for possible future analysis
saveRDS(AICmodlist,file="AICmodlist_",nameID,".rds")
saveRDS(STANmodlist,file="STANmodlist_",nameID,".rds")


