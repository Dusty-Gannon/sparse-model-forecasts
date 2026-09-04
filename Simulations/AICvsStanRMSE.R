# For RMSE of predictions, we need to output a longer time series than we actually use to fit the model
#
# Fit the step-wise AIC method
# Extract RMSE
# Fit the STAN model
# Extract RMSE
#
# Plot and compare

library(dplyr)
library(MASS)
library(rstan)
library(here)
library(purrr)
source(here("R/sim1_temporal_drivers.R"))
source(here("R/tau0_from_data.R"))

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
  all_vars <- names(dataSet)[-1]
  K <- length(all_vars)
  mod_full=lm(y~.,data=dataSet)
  mod_null=lm(y~1,data=dataSet)
  result<-MASS::stepAIC(mod_null, direction = "forward",
    scope = list(lower = mod_null, upper = mod_full),
    keep = function(model, aic) list(model = model, aic = aic))

  chain_models <- result$keep["model", ]
  chain_aic    <- unlist(result$keep["aic", ])
  delta        <- chain_aic - min(chain_aic)
  weights      <- exp(-0.5 * delta) / sum(exp(-0.5 * delta))

  avg_coef <- setNames(numeric(K + 1), c("(Intercept)", all_vars))
  for (i in seq_along(chain_models)) {
    cf <- coef(chain_models[[i]])
    avg_coef[names(cf)] <- avg_coef[names(cf)] + weights[i] * cf
  }

  result$avg_coef    <- avg_coef
  result$avg_weights <- weights
  return(result)
}



#' Title Calculate prediction RMSE for standard linear model with AIC selection
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


#' Title Calculate prediction RMSE for full GLM model
#'
#' @param model The model to be evaluated
#' @param testData The test dataset to be used for prediction and model evaluation, response variable is called y
#'
#' @returns The root mean squared error for the model prediction
#' @export
#'
#' @examples
RMSE_GLM=function(model,testData){
  predValues<-predict(model,testData)
  RMSE_GLM<-sqrt(mean((predValues-testData$y)^2))
  return(RMSE_GLM)
}


#' Title Calculate prediction RMSE using model-averaged coefficients from AICselect
#'
#' @param aicModel AICselect result (must have $avg_coef attached)
#' @param testData The test dataset, response in column y
#'
#' @return scalar RMSE
#' @export
RMSE_modelAvg=function(aicModel, testData){
  avg_coef <- aicModel$avg_coef
  X <- model.matrix(~ ., data = testData[, -1, drop = FALSE])
  preds <- as.vector(X[, names(avg_coef)] %*% avg_coef)
  sqrt(mean((preds - testData$y)^2))
}


#' Title Function to collect true/false positives/negatives from the full GLM model
#'
#' @param timeseries timeseries that model is based on; must contain strong_ids
#' @param model glm result model to be evaluated
#'
#' @returns true positive and negative rates based on p<0.05
#' @export
#'
#' @examples
GLMconfusionRates<-function(timeseries, model){

  strong2 <- paste0("driver_", timeseries$strong_ids)
  strongaic <- names(which(summary(model)$coefficients[,4]<0.05))

  positives <- length(timeseries$strong_ids)
  negatives <- length(timeseries$beta) - 1 - positives

  truePos <- length(intersect(strongaic, strong2))
  trueNeg <- negatives - length(strongaic) + length(intersect(strongaic, strong2))

  TPR <- truePos/positives
  TNR <- trueNeg/negatives

  return(c(TPR,TNR))
}

# function to collect true/false positives/negatives from the AIC based model
#' Title
#'
#' @param timeseries timeseries that model is based on; must contain strong_ids
#' @param model step AIC lm result model to be evaluated
#'
#' @return true positive and negative rates
#' @export
#'
#' @examples
AICconfusionRates<-function(timeseries, model){

  strong2 <- paste0("driver_", timeseries$strong_ids)
  strongaic <- names(model$coefficients)

  positives <- length(timeseries$strong_ids)
  negatives <- length(timeseries$beta) - 1 - positives

  truePos <- length(intersect(strongaic, strong2))
  trueNeg <- negatives - length(strongaic) + length(intersect(strongaic, strong2))

  TPR <- truePos/positives
  TNR <- trueNeg/negatives

  return(c(TPR,TNR))
}


#' Title Confusion rates from model-averaged AIC variable importances
#'
#' @param timeseries Simulated timeseries with strong_ids
#' @param aicModel AICselect result (must have $keep, $avg_weights)
#'
#' @return c(TPR, TNR); detection defined as relative importance > 0.5
#' @export
modelAvg_confusionRates<-function(timeseries, aicModel){
  strong2 <- paste0("driver_", timeseries$strong_ids)
  K <- length(timeseries$beta) - 1
  all_vars <- paste0("driver_", 1:K)
  positives <- length(timeseries$strong_ids)
  negatives <- K - positives

  chain_models <- aicModel$keep["model", ]
  weights <- aicModel$avg_weights

  ri <- setNames(numeric(K), all_vars)
  for (i in seq_along(chain_models)) {
    in_model <- setdiff(names(coef(chain_models[[i]])), "(Intercept)")
    present <- intersect(in_model, all_vars)
    ri[present] <- ri[present] + weights[i]
  }
  detected <- names(ri)[ri > 0.5]

  truePos <- length(intersect(detected, strong2))
  trueNeg <- negatives - length(detected) + length(intersect(detected, strong2))

  return(c(TPR = truePos/positives, TNR = trueNeg/negatives))
}


# function to get a list of time series
#' Title
#'
#' @param numTrials number of timeseries that you want
#' @param n length of timeseries
#' @param K number of betas/covariates
#' @param num_strong number of strong predictors
#' @param prob_cycle probability of a covariate having a yearly cycle
#' @param trend_fraction probability of a covariate having a trend
#' @param freq frequency of observations per time step
#' @param sigma Standard deviation of the noise component
#' @param correlated Logical value giving whether correlated covariates should be used
#' @param rateCorr Amount of correlation to use, making the rate parameter smaller increases the possible cov values
#'
#' @return list of numTrials timeseries of length n
#' @export
#'
#' @examples
getTS<-function(numTrials,n=100,K=50, num_strong=5, prob_cycle=0, trend_fraction=0,freq = 1, sigma = 0.5,probWeakCorr=0,numStrongCorr=0,strongSelf=F,corrLevel=0,corrChange=F, propChange=0.5, changeSize=0.5, changeTimeVar=5){
  ts_AICselection <- lapply(as.list(1:numTrials),
                            FUN = function(x)
                              basic_timeseries(K = K, num_strong = num_strong, n = n, freq = freq, prob_cycle = prob_cycle, trend_fraction = trend_fraction, sigma = sigma, probWeakCorr=probWeakCorr, numStrongCorr = numStrongCorr, strongSelf=strongSelf, corrLevel=corrLevel, corrChange=corrChange, propChange=propChange, changeSize=changeSize, changeTimeVar=changeTimeVar)
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
  names(simDat_df) <- c("y",paste0("driver_",1:(K)))
  return(simDat_df)
}

# function to split the cleaned time series into training and test sets
#' Title
#'
#' @param tsclean cleaned time series with y in column 1 and x in the other columns
#' @param set either "test" or "train"
#' @param nfit the number of observations to be used to fit the model (to train the data)
#' @param n full length of time series
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

  # standardize predictor columns using training-set statistics
  X_raw     <- as.matrix(dataSet[, 2:(K+1)])
  X_new_raw <- as.matrix(testSet[, 2:(K+1)])
  col_means <- colMeans(X_raw)
  col_sds   <- apply(X_raw, 2, sd)
  X_std     <- scale(X_raw, center = col_means, scale = col_sds)
  X_new_std <- scale(X_new_raw, center = col_means, scale = col_sds)

  # compile data for stan
  datlist<-list(
    N=nfit, # length of time series
    P0=1, # for intercept
    P=K+1, # number of covariates
    y=dataSet$y,
    X=cbind(rep(1,nfit), X_std),
    tau0=tau_0,
    slab_scl=slab_scl,
    slab_df=slab_df,
    N_new=n - nfit,
    X_new=cbind(rep(1, n - nfit), X_new_std)
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

#' Title Function to get predictions from stan model
#'
#' @param model model fit from stan
#'
#' @return the predictions the model made for the testing data
#' @export
#'
#' @examples
STANgetpredict=function(model){
  predictionInfo=rstan::extract(model, pars = "y_rep")$y_rep
  return(predictionInfo)
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

STANbetapost <- function(modelfit, scaled = T, sd_x = NULL) {

  if(scaled & is.null(sd_x)){
    stop("Were columns of X_train standardized?
         If so, we need to rescale beta with sd_x for them to be compared to the simulation values.")
  }
  # extract the posterior betas
  beta_post=rstan::extract(modelfit, pars = 'beta')$beta
  # rescale all columns except the intercept
  beta_post[, 2:ncol(beta_post)] <- t(t(beta_post[, 2:ncol(beta_post)]) / sd_x)

  # create a data frame w summary stats for these posterior betas
  beta_df <- data.frame(
    mean = apply(beta_post, 2, mean),
    median = apply(beta_post, 2, median),
    min = apply(beta_post, 2, min),
    max = apply(beta_post, 2, max),
    low = apply(beta_post, 2, quantile, probs = 0.025),
    high = apply(beta_post, 2, quantile, probs = 0.975),
    q0.01 = apply(beta_post, 2, quantile, probs = 0.01),
    q0.05 = apply(beta_post, 2, quantile, probs = 0.05),
    q0.1 = apply(beta_post, 2, quantile, probs = 0.1),
    q0.9 = apply(beta_post, 2, quantile, probs = 0.9),
    q0.95 = apply(beta_post, 2, quantile, probs = 0.95),
    q0.99 = apply(beta_post, 2, quantile, probs = 0.99)

  )

  return(beta_df)
}



#' Compute unconditional standard errors for full model-averaged coefficients
#'
#' Uses the Burnham & Anderson (2004) formula: SE_j = sum_i w_i * sqrt(var_ij + (beta_ij - beta_bar_j)^2)
#' where beta_ij = 0 and var_ij = 0 for models that do not contain variable j (full averaging).
#'
#' @param aicModel AICselect result with $keep, $avg_weights, and $avg_coef attached
#' @param all_vars Character vector of predictor names (e.g. "driver_1" ... "driver_K")
#'
#' @return Named numeric vector of unconditional SEs, one per variable in all_vars
modelAvg_uncond_SE <- function(aicModel, all_vars) {
  chain_models <- aicModel$keep["model", ]
  weights      <- aicModel$avg_weights
  avg_coef     <- aicModel$avg_coef

  se <- setNames(numeric(length(all_vars)), all_vars)
  for (j in all_vars) {
    beta_bar <- avg_coef[j]
    contrib <- vapply(seq_along(chain_models), function(i) {
      cf <- coef(chain_models[[i]])
      if (j %in% names(cf)) {
        beta_ij <- cf[j]
        var_ij  <- vcov(chain_models[[i]])[j, j]
      } else {
        beta_ij <- 0
        var_ij  <- 0
      }
      weights[i] * sqrt(var_ij + (beta_ij - beta_bar)^2)
    }, numeric(1))
    se[j] <- sum(contrib)
  }
  se
}


#' Compute interval coverage rates for AIC, Stan, and full-GLM models
#'
#' For each model, checks what fraction of the true beta values (from the
#' simulation) fall inside the estimated 95% confidence or credible interval.
#'
#' AIC: coverage is computed only for the variables that were selected.
#' Stan: coverage is split by whether the 95% CrI excludes zero ("escapes") or
#'   not.  Stan posteriors are back-transformed from the standardised predictor
#'   scale to the original scale before comparison.
#' GLM: coverage is computed for all K predictor variables.
#'
#' @param timeseries Simulated timeseries object containing true beta values
#'   (intercept first, then K driver betas).
#' @param trainData Training data frame used to fit the models (y in column 1,
#'   driver_1 … driver_K in the remaining columns).  Used to recover the
#'   column standard deviations needed to back-transform Stan posteriors.
#' @param aicModel Stepwise-AIC lm model fit returned by AICselect().
#' @param stanModel Stan model fit returned by STANselect().
#' @param glmModel Full lm model fit (all K predictors).
#' @param K Number of predictor variables.
#' @param conf Nominal coverage level (default 0.95).
#'
#' @return A one-row data frame with four coverage rates:
#'   \describe{
#'     \item{AIC_coverage}{Proportion of *selected* variables whose true beta
#'       falls inside the 95\% CI.}
#'     \item{GLM_coverage}{Proportion of *all* K variables whose true beta
#'       falls inside the 95\% CI.}
#'     \item{Stan_coverage_escape}{Coverage among posteriors whose 95\% CrI
#'       excludes zero.}
#'     \item{Stan_coverage_noescape}{Coverage among posteriors whose 95\% CrI
#'       includes zero.}
#'   }
#' @export
coverageRates <- function(timeseries, trainData, aicModel, stanModel, glmModel,
                          K = 50, conf = 0.95) {

  true_betas <- timeseries$beta[2:(K + 1)]  # driver betas, original scale
  driver_names <- paste0("driver_", 1:K)

  ## ---- AIC: CI for selected variables only ----
  aic_selected <- setdiff(names(aicModel$coefficients), "(Intercept)")

  if (length(aic_selected) == 0) {
    aic_coverage <- NA_real_
  } else {
    aic_ci <- confint(aicModel, parm = aic_selected, level = conf)
    aic_idx <- as.integer(sub("driver_", "", aic_selected))
    aic_true <- true_betas[aic_idx]
    aic_coverage <- mean(aic_true >= aic_ci[, 1] & aic_true <= aic_ci[, 2])
  }

  ## ---- GLM: CI for all K variables ----
  glm_ci <- confint(glmModel, level = conf)[driver_names, ]
  glm_coverage <- mean(true_betas >= glm_ci[, 1] & true_betas <= glm_ci[, 2])

  ## ---- Stan: CrI back-transformed to original predictor scale ----
  # Stan fits on standardised X; dividing by col_sds recovers the original scale.
  col_sds <- apply(as.matrix(trainData[, 2:(K + 1)]), 2, sd)
  beta_post <- STANbetapost(stanModel, sd_x = col_sds)
  bp <- beta_post[2:(K + 1), ]  # drop intercept row

  escape_zero <- bp$low > 0 | bp$high < 0

  if (sum(escape_zero) == 0) {
    stan_coverage_escape <- NA_real_
  } else {
    stan_coverage_escape <- mean(
      true_betas[escape_zero] >= bp$low[escape_zero] &
      true_betas[escape_zero] <= bp$high[escape_zero]
    )
  }

  if (sum(!escape_zero) == 0) {
    stan_coverage_noescape <- NA_real_
  } else {
    stan_coverage_noescape <- mean(
      true_betas[!escape_zero] >= bp$low[!escape_zero] &
      true_betas[!escape_zero] <= bp$high[!escape_zero]
    )
  }

  ## ---- Model averaging: unconditional CIs (Burnham & Anderson 2004) ----
  z <- qnorm((1 + conf) / 2)
  mavg_se <- modelAvg_uncond_SE(aicModel, driver_names)
  mavg_lo <- aicModel$avg_coef[driver_names] - z * mavg_se
  mavg_hi <- aicModel$avg_coef[driver_names] + z * mavg_se
  escape_mavg <- mavg_lo > 0 | mavg_hi < 0
  in_ci_mavg  <- true_betas >= mavg_lo & true_betas <= mavg_hi

  mavg_coverage_escape <- if (sum(escape_mavg) == 0) NA_real_ else
    mean(in_ci_mavg[escape_mavg])
  mavg_coverage_noescape <- if (sum(!escape_mavg) == 0) NA_real_ else
    mean(in_ci_mavg[!escape_mavg])

  data.frame(
    AIC_coverage           = aic_coverage,
    GLM_coverage           = glm_coverage,
    Stan_coverage_escape   = stan_coverage_escape,
    Stan_coverage_noescape = stan_coverage_noescape,
    MAvg_coverage_escape   = mavg_coverage_escape,
    MAvg_coverage_noescape = mavg_coverage_noescape
  )
}



## The STANconfusionRates() code leans HEAVILY on Alice C's 'calculate_true_pos_rate()', but some things don't make sense to me"
##    - the 'zero' values for beta
##    - there is one definition of 'true positive' in the description, but it looks like the code is doing something else...


#' function to provide confusion metrics
#'
#' A true positive is a driver whose posterior mass (defined by threshold)
#' lies entirely above or below zero, given that it is a true large-effect
#' predictor (i.e. its index is in timeseries$strong_ids).  The intercept row
#' of par_post is excluded; only the K driver rows are evaluated.
#'
#' @param par_post a summary of the posterior estimates from STANbetapost()
#' @param timeseries the simulated timeseries object; must contain strong_ids
#' @param par label prefix for the returned column names (default 'beta')
#' @param threshold posterior mass threshold defining a 'detected' effect (default 0.90)
#'
#' @return a data frame of TPR, TNR, FPR, FNR prefixed by par
#' @export
#'
#' @examples
STANconfusionRates<-function(par_post, timeseries, par='beta', threshold=0.90) {

  K <- length(timeseries$beta) - 1
  strong_ids <- timeseries$strong_ids
  true_betas <- timeseries$beta[2:(K + 1)]

  bottom <- paste0('q', (1-threshold))
  top <- paste0('q', (threshold))

  # operate on driver rows only (skip intercept at row 1)
  bp <- par_post[2:(K + 1), c('mean', 'median', bottom, top)]
  colnames(bp) <- c('mean', 'median', 'bottom', 'top')

  is_strong <- seq_len(K) %in% strong_ids

  ests <- bp %>%
    mutate(effect = case_when(!is_strong      ~ 'zero',
                              true_betas > 0  ~ 'positive',
                              true_betas < 0  ~ 'negative'),
           post_mass = case_when(bottom > 0 ~ 'positive',
                                 top < 0    ~ 'negative',
                                 TRUE       ~ 'zero'))

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

  return(TFP)

}

