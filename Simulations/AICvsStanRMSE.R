# For RMSE of predictions, we need to output a longer time series than we actually use to fit the model
#
# Fit the step-wise AIC method
# Extract RMSE
# Fit the STAN model
# Extract RMSE
#
# Plot and compare

library(MASS)
library(rstan)
library(here)
devtools::load_all()

# compile the stan model
stanFHS<-stan_model(here("Stan/sparse_reg_FHS.stan"))


#' Title AIC Stepwise Model Selection
#'
#' @param dataSet The dataset used to fit the model, response variable is called y
#'
#' @return The model fit by stepAIC forward selection
#' @export
#'
#' @examples
AICselect=function(dataSet){
  mod_full=lm(y~.,data=dataSet)
  mod_null=lm(y~1,data=dataSet)
  stepAIC<-MASS::stepAIC(mod_null, direction = "forward", scope = list(lower = mod_null,                                                                     upper = mod_full))
  return(stepAIC)
}


#' Title Calculate prediction RMSE for standard lm model
#'
#' @param model A model to be evaluated
#' @param testData The test dataset to be used for prediction and model evaluation, response variable is called y
#'
#' @return The root mean squared error for the model prediction
#' @export
#'
#' @examples
RMSE_AIC=function(model,testData){
  predValues<-predict(model,testData)
  RMSE_AIC<-sqrt(mean((predValues-testData$y)^2))
  return(RMSE_AIC)
}


# function to collect true/false positives/negatives from the AIC based model
#' Title
#'
#' @param timeseries timeseries that model is based on, needed to extract true betas
#' @param model step AIC lm result model to be evaluated
#' @param strongCutoff cutoff for true positive
#'
#' @return true positive and negative rates
#' @export
#'
#' @examples
AICconfusionRates<-function(timeseries,model,strongCutoff=0.3){

  strong1<-which(abs(timeseries$beta)>strongCutoff)
  strong2<-paste0("driver_",strong1)
  strongaic<-names(model$coefficients)

  positives=length(strong1)
  negatives=(length(timeseries$beta)-length(strong1)-1)

  truePos=length(intersect(strongaic,strong2))
  trueNeg=negatives-length(strongaic)+length(intersect(strongaic,strong2))

  TPR=truePos/positives
  TNR=trueNeg/negatives

  return(c(TPR,TNR))
}


# function to get a list of time series
#' Title
#'
#' @param numTrials number of timeseries that you want
#' @param n length of timeseries
#' @param K number of betas/covariates
#'
#' @return list of numTrials timeseries of length n
#' @export
#'
#' @examples
getTS<-function(numTrials,n=100,K=50, num_strong=5, prob_cycle=0, trend_fraction=0,freq = 1, sigma = 0.5,correlated=F,rateCorr=2){
  ts_AICselection <- lapply(as.list(1:numTrials),
                            FUN = function(x)
                              basic_timeseries(K = K, num_strong = num_strong, n = n, freq = freq, prob_cycle = prob_cycle, trend_fraction = trend_fraction, sigma = sigma, correlated=correlated, rateCorr = rateCorr)
  )
  return(ts_AICselection)
}

# function to clean up time series
#' Title
#'
#' @param ts raw timeseries to be cleaned
#'
#' @return data frame with the y and x values for the timeseries- ready for modeling
#' @export
#'
#' @examples
cleanTS<-function(ts){
  K=length(ts$beta)-1
  simDat_df <- data.frame("y" = ts$y, ts$X[,2:(K+1)])
  names(simDat_df) <- c("y",paste0("driver_",2:(K+1)))
  return(simDat_df)
}

# function to split the cleaned time series into training and test sets
#' Title
#'
#' @param tsclean cleaned time series with y in column 1 and x in the other columns
#' @param set either "test" or "train"
#' @param nfit the number of observations to be used to fit the model (to train the data)
#'
#' @return the testing or training set you requested
#' @export
#'
#' @examples
splitTS<-function(tsclean,set,nfit=60,n=100){
  if(set=="test"){
    ans=tsclean[(nfit+1):n,]
  }
  if(set=="train"){
    ans=tsclean[1:nfit,]
  }
  return(ans)
}

# function to give stan based models
#' Title
#'
#' @param dataSet a training data set, cleaned, and formatted with y in column 1 and x's in other columns
#' @param testSet a testing data set to predict from in order to test the model
#' @param m0 predicted number of strong predictors
#' @param K number of covariates
#' @param n sum of lengths of training plus testing data set
#' @param nfit length of training data set
#' @param slab_scl estimate for value of strong predictors
#' @param slab_df degrees of freedom, determining amount of regularization, with higher df meaning large effects are more regularized but model converges better , and lower meaning large effects closer to true value but model converges worse
#'
#' @return the model fit
#' @export
#'
#' @examples
STANselect<-function(dataSet,testSet,m0=5,K=50,n=100,nfit=60,slab_scl=1,slab_df=10){
  tau_0=tau0(
    y=dataSet$y,
    m0=m0,
    M=K,
    N=nfit,
    fam="gaussian"
  )

  # compile data for stan
  datlist<-list(
    N=nfit, # length of time series
    P0=1, # for intercept
    P=K, # number of covariates
    y=dataSet$y,
    X=dataSet[,2:51],
    tau0=tau_0,
    slab_scl=slab_scl,
    slab_df=slab_df,
    N_new=n-nfit,
    X_new=testSet[,2:51]
  )

  # sample the posterior
  mfit_FHScompare<-sampling(
    stanFHS,
    data=datlist,
    chains=3,
    cores=3
  )

  return(mfit_FHScompare)
}


# function to extract posterior predictions for beta and place them in a summary data frame
#' Title
#'
#' @param modelfit a stan model fit
#'
#' @return the summary data frame for the posterior beta predictions
#' @export
#'
#' @examples

STANbetapost <- function(modelfit) {

  # extract the posterior betas
  beta_pred=rstan::extract(modelfit, pars = 'beta')$beta

  # create a data frame w summary stats for these posterior betas
  beta_post <- data.frame(
    mean = apply(beta_pred, 2, mean),
    median = apply(beta_pred, 2, median),
    min = apply(beta_pred, 2, min),
    max = apply(beta_pred, 2, max),
    low = apply(beta_pred, 2, quantile, probs = 0.025),
    high = apply(beta_pred, 2, quantile, probs = 0.975),
    q0.01 = apply(beta_pred, 2, quantile, probs = 0.01),
    q0.05 = apply(beta_pred, 2, quantile, probs = 0.05),
    q0.1 = apply(beta_pred, 2, quantile, probs = 0.1),
    q0.9 = apply(beta_pred, 2, quantile, probs = 0.9),
    q0.95 = apply(beta_pred, 2, quantile, probs = 0.95),
    q0.99 = apply(beta_pred, 2, quantile, probs = 0.99)

  )

  return(beta_post)
}


## The STANconfusionRates() code leans HEAVILY on Alice C's 'calculate_true_pos_rate()', but some things don't make sense to me"
##    - the 'zero' values for beta
##    - there is one definition of 'true positive' in the description, but it looks like the code is doing something else...

#' function to provide confusion metrics (rates for true pos, false pos, false neg) for stan based models
#'
#' For a model fit to simulated AR-p data, a true positive is a correctly
#' identified beta or phi parameter as non-zero. **We define a true positive as
#' any parameter estimate for which 90% of the posterior mass lies above or
#' below zero, given that the true parameter is non-zero.**  CHECK THIS IN AM
#'
#' Title
#'
#' @param par_post a summary of the posterior estimates of the parameter values (beta is the default)
#' @param par_vals the data set ('true') values for the parameter
#' @param threshold the posterior mass threshold that defines a 'true positive' result;
#' 0.9 (90%) is the default value
#'
#' @return a matrix of the true positive, false positive and false negative rates for the model
#' @export
#'
#' @examples
STANconfusionRates<-function(par_post, par_vals, par='beta', threshold=0.90) {
  min_pos = 0
  if(par == 'beta') min_pos = 1.96 * 0.05 #this is from Alice C's code, not sure why this was done/ why these values for 'zero' - CT

  # calculate the boundaries of the posterior mass based on the threshold
  bottom <- paste0('q', (1-threshold))
  top <- paste0('q', (threshold))

  ests <- par_post[,c('mean', 'median', bottom, top)]
  colnames(ests) <- c('mean', 'median', 'bottom', 'top')

  ests <- ests %>%
    mutate(value = par_vals,
           effect = case_when(abs(value) <= min_pos ~ 'zero',
                              value > min_pos ~ 'positive',
                              value < -min_pos ~ 'negative'),
           post_mass = case_when(bottom > 0 ~ 'positive',
                                 top < 0 ~ 'negative',
                                 TRUE ~ 'zero'))
  # True Positive Rate
  true_pos = sum((ests$effect == 'positive' & ests$post_mass == 'positive')|
                   (ests$effect == 'negative' & ests$post_mass == 'negative'))/
    sum(ests$effect != 'zero')

  # False Positive Rate
  false_pos = sum((ests$effect != 'negative' & ests$post_mass == 'negative')|
                    (ests$effect != 'positive' & ests$post_mass == 'positive'))/
    sum(ests$effect == 'zero')

  # True Negative Rate
  true_neg = sum(ests$effect == 'zero' & ests$post_mass == 'zero')/
    sum(ests$effect == 'zero')

  # False Negative Rate
  false_neg = sum((ests$effect == 'negative' & ests$post_mass != 'negative')|
                    (ests$effect == 'positive' & ests$post_mass != 'positive'))/
    sum(ests$effect != 'zero')


  TFP <- data.frame(TPR = true_pos,
                    TNR = true_neg,
                    FPR = false_pos,
                    FNR = false_neg) %>%
    rename_with(function(x) paste0(par, '_', x))

  return( TFP )

}


#trials
start=Sys.time()
ts1=getTS(100) # make timeseries
ts1clean=lapply(ts1,FUN=function(x) cleanTS(x)) # clean them
ts1test=lapply(ts1clean,FUN=function(x) splitTS(x,set="test")) #testing set
ts1train=lapply(ts1clean,FUN=function(x) splitTS(x,set="train")) #training set
AICmodlist=lapply(ts1train,FUN=function(x) AICselect(x)) # do model selection AIC
RMSEAIClist=mapply(function(x,y) RMSE_AIC(x,y), AICmodlist,ts1test) # get RMSE AIC
STANmodlist=mapply(function(x,y) STANselect(x,y,nfit=60),ts1train,ts1test) # do model selection STAN
STANpredlist=lapply(STANmodlist,FUN=function(x) STANgetpredict(x)) # get STAN predictions for y
STANy=lapply(ts1test,"[",,1)

STANbetalist=lapply(STANmodlist, FUN=function(x) STANbetapost(x)) # get summary of stan predictions for beta
ts_betas=lapply(ts1, FUN=function(x) x$beta[2:51])  # removing the beta[1], 0 in all cases, so there are 50 ts_betas to match the 50 modeled betas

RMSESTANlistraw=mapply(function(x,y) RMSE_bayes(x,y),STANy,STANpredlist) # get RMSE STAN
RMSESTANlist=colMeans(RMSESTANlistraw)
Sys.time()-start

AICconfusion=mapply(function(x,y) AICconfusionRates(x,y),ts1,AICmodlist)
mean(AICconfusion[1,])# TPR is 0.992
mean(AICconfusion[2,])# TNR is 0.498

STANconfusion=mapply(function(x,y) STANconfusionRates(x,y), STANbetalist, ts_betas) #returns a matrix of TPR, TNR, FPR, FNR (redundant but added so I could check myself)


TPR_STANbeta <- mean(unlist(STANconfusion[1,]))   # 0.759
TNR_STANbeta <- mean(unlist(STANconfusion[2,]))   # 0.946
FPR_STANbeta <- mean(unlist(STANconfusion[3,]))   # 0.054
FNR_STANbeta <- mean(unlist(STANconfusion[4,]))   # 0.241


hist(RMSESTANlist,xlim=c(0,2.5),ylim=c(0,50),breaks=seq(0,2.5,by=0.1),xlab="RMSE",main="RMSE using horseshoe priors")
hist(RMSEAIClist,xlim=c(0,2.5),ylim=c(0,50),breaks=seq(0,2.5,by=0.1),xlab="RMSE",main="RMSE using stepwise AIC")

mean(RMSESTANlist)
mean(RMSEAIClist)

