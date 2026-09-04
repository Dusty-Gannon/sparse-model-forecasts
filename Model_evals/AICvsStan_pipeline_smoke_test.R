# Smoke test for the full AICvsStanBatch pipeline.
# Runs a single simulated trial end-to-end so each step can be verified by hand.
# Run from the project root: Rscript Model_evals/pipeline_smoke_test.R

devtools::load_all()
source(here::here("Simulations/AICvsStanRMSE.R"))

# ---- Parameters ----
set.seed(260904)
K <- 10
n <- 50
nfit <- 30   # train: rows 1-30, test: rows 31-50
num_strong <- 3

# ---- Simulate one time series ----
ts_raw <- getTS(
  numTrials=1, n=n, K=K, num_strong=num_strong,
  prob_cycle=0, trend_fraction=0, freq=1, sigma=0.5,
  probWeakCorr=0, numStrongCorr=0, strongSelf=FALSE,
  corrLevel=0, corrChange=FALSE, propChange=0,
  changeSize=0, changeTimeVar=0
)[[1]]

cat("\n=== TRUE SIMULATION VALUES ===\n")
cat("strong_ids (1-based driver columns):", ts_raw$strong_ids, "\n")
cat("beta (intercept + drivers 1-K):\n")
bdf <- data.frame(
  index = 0:K,
  beta  = round(ts_raw$beta, 4),
  strong = c(FALSE, seq_len(K) %in% ts_raw$strong_ids)
)
print(bdf)

# ---- Clean and split ----
ts_clean <- cleanTS(ts_raw)
train    <- splitTS(ts_clean, set="train", n=n, nfit=nfit)
test     <- splitTS(ts_clean, set="test",  n=n, nfit=nfit)

# ---- AIC model ----
cat("\n=== AIC STEPWISE ===\n")
aic_mod  <- AICselect(train)
rmse_aic <- RMSE_AIC(aic_mod, test)
aic_conf <- AICconfusionRates(ts_raw, aic_mod)
cat("Selected variables:", paste(setdiff(names(aic_mod$coefficients), "(Intercept)"), collapse=", "), "\n")
cat("RMSE:", round(rmse_aic, 4), "\n")
cat("TPR:", round(aic_conf[1], 4), "  TNR:", round(aic_conf[2], 4), "\n")

# ---- AIC model averaging ----
cat("\n=== AIC MODEL AVERAGING ===\n")
rmse_mavg <- RMSE_modelAvg(aic_mod, test)
mavg_conf <- modelAvg_confusionRates(ts_raw, aic_mod)
cat("Chain length (models in average):", length(aic_mod$avg_weights), "\n")
cat("AIC weights:", round(aic_mod$avg_weights, 3), "\n")
cat("Model-averaged coefficients:\n")
print(round(aic_mod$avg_coef, 4))
cat("RMSE:", round(rmse_mavg, 4), "\n")
cat("TPR:", round(mavg_conf["TPR"], 4), "  TNR:", round(mavg_conf["TNR"], 4), "\n")

# ---- Full GLM ----
cat("\n=== FULL GLM ===\n")
glm_mod  <- glm(y~., data=train)
rmse_glm <- RMSE_GLM(glm_mod, test)
glm_conf <- GLMconfusionRates(ts_raw, glm_mod)
cat("RMSE:", round(rmse_glm, 4), "\n")
cat("TPR:", round(glm_conf[1], 4), "  TNR:", round(glm_conf[2], 4), "\n")

# ---- Stan (RHS) ----
cat("\n=== STAN (RHS) — sampling, please wait... ===\n")
stan_mod   <- STANselect(train, test, m0=num_strong, K=K, n=n, nfit=nfit)
stan_preds <- STANgetpredict(stan_mod)
rmse_stan  <- median(RMSE_bayes(test$y, stan_preds[, -(1:nfit)]))
stan_bp    <- STANbetapost(stan_mod, sd_x = apply(train[,-1], 2, sd))
stan_conf  <- STANconfusionRates(stan_bp, ts_raw)

cat("RMSE:", round(rmse_stan, 4), "\n")
cat("TPR:", round(stan_conf$beta_TPR, 4), "  TNR:", round(stan_conf$beta_TNR, 4), "\n")
cat("FPR:", round(stan_conf$beta_FPR, 4), "  FNR:", round(stan_conf$beta_FNR, 4), "\n")

print(stan_bp[, c("mean", "low", "high")])

# ---- Coverage rates ----
cat("\n=== COVERAGE RATES (95%) ===\n")
cov_rates <- coverageRates(ts_raw, train, aic_mod, stan_mod, glm_mod, K=K)
print(round(cov_rates, 4))

cat("\n=== DONE ===\n")
