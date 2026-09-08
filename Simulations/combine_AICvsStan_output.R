# Combine AICvsStan array-job output (.rds files) with the config file
# and write a single tidy .csv, with one row per simulation trial.
#
# Usage:
#   Rscript Simulations/combine_AICvsStan_output.R <data_dir> [config_file] [out_file]
#
#   <data_dir>    directory holding the per-array-job .rds output files
#                 (e.g. Data/AICvsStan_data/)
#   [config_file] path to the config file used to launch the array jobs
#                 (default: Simulations/AICvsStanConfig.txt)
#   [out_file]    path to write the combined .csv
#                 (default: Simulations/AICvsStanResults.csv)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript combine_AICvsStan_output.R <data_dir> [config_file] [out_file]")
}

data_dir <- sub("/+$", "", args[1])
config_file <- if (length(args) >= 2) args[2] else "Simulations/AICvsStanConfig.txt"
out_file <- if (length(args) >= 3) args[3] else "Simulations/AICvsStanResults.csv"

config <- read.table(config_file, header = TRUE, stringsAsFactors = FALSE)

# unlist a mapply-style matrix of list-mode values (e.g. STANconfusion) into
# a plain numeric matrix, keeping its dim/dimnames
as_numeric_matrix <- function(x) {
  m <- matrix(as.numeric(unlist(x)), nrow = nrow(x), ncol = ncol(x))
  dimnames(m) <- dimnames(x)
  m
}

results_list <- vector("list", nrow(config))

for (i in seq_len(nrow(config))) {

  nameID <- config$nameID[i]

  files <- c(
    RMSEAIC        = file.path(data_dir, paste0("RMSEAIC_", nameID, ".rds")),
    RMSEGLM        = file.path(data_dir, paste0("RMSEGLM_", nameID, ".rds")),
    RMSESTAN       = file.path(data_dir, paste0("RMSESTAN_", nameID, ".rds")),
    RMSEmodelAvg   = file.path(data_dir, paste0("RMSEmodelAvg_", nameID, ".rds")),
    AICconfusion   = file.path(data_dir, paste0("AICconfusion_", nameID, ".rds")),
    GLMconfusion   = file.path(data_dir, paste0("GLMconfusion_", nameID, ".rds")),
    modelAvgConfusion = file.path(data_dir, paste0("modelAvgConfusion_", nameID, ".rds")),
    STANconfusion  = file.path(data_dir, paste0("STANconfusion_", nameID, ".rds")),
    coverage       = file.path(data_dir, paste0("coverage_", nameID, ".rds"))
  )

  missing <- files[!file.exists(files)]
  if (length(missing) > 0) {
    warning(sprintf("nameID '%s': missing %d/%d output file(s), skipping.",
                     nameID, length(missing), length(files)))
    next
  }

  RMSEAIC <- readRDS(files["RMSEAIC"])
  RMSEGLM <- readRDS(files["RMSEGLM"])
  RMSESTAN <- readRDS(files["RMSESTAN"])
  RMSEmodelAvg <- readRDS(files["RMSEmodelAvg"])

  AICconfusion <- readRDS(files["AICconfusion"])
  GLMconfusion <- readRDS(files["GLMconfusion"])
  modelAvgConfusion <- readRDS(files["modelAvgConfusion"])
  STANconfusion <- as_numeric_matrix(readRDS(files["STANconfusion"]))

  coverage <- readRDS(files["coverage"])

  ntrials <- length(RMSEAIC)
  lengths_ok <- c(
    length(RMSEGLM), length(RMSESTAN), length(RMSEmodelAvg),
    ncol(AICconfusion), ncol(GLMconfusion), ncol(modelAvgConfusion),
    ncol(STANconfusion), nrow(coverage)
  ) == ntrials
  if (!all(lengths_ok)) {
    warning(sprintf("nameID '%s': output lengths do not match across files, skipping.",
                     nameID))
    next
  }

  trial_df <- data.frame(
    ArrayTaskID = config$ArrayTaskID[i],
    nameID = nameID,
    trial = seq_len(ntrials),

    RMSE_AIC = as.numeric(RMSEAIC),
    RMSE_GLM = as.numeric(RMSEGLM),
    RMSE_STAN = as.numeric(RMSESTAN),
    RMSE_modelAvg = as.numeric(RMSEmodelAvg),

    AIC_TPR = AICconfusion[1, ],
    AIC_TNR = AICconfusion[2, ],

    GLM_TPR = GLMconfusion[1, ],
    GLM_TNR = GLMconfusion[2, ],

    modelAvg_TPR = modelAvgConfusion["TPR", ],
    modelAvg_TNR = modelAvgConfusion["TNR", ],

    STAN_beta_TPR = STANconfusion["beta_TPR", ],
    STAN_beta_TNR = STANconfusion["beta_TNR", ],
    STAN_beta_FPR = STANconfusion["beta_FPR", ],
    STAN_beta_FNR = STANconfusion["beta_FNR", ],

    stringsAsFactors = FALSE
  )

  trial_df <- cbind(trial_df, coverage)

  # attach the config parameters for this array task (repeated for every trial)
  config_row <- config[i, setdiff(names(config), c("ArrayTaskID", "nameID")), drop = FALSE]
  rownames(config_row) <- NULL
  trial_df <- cbind(trial_df, config_row[rep(1, ntrials), , drop = FALSE])

  results_list[[i]] <- trial_df
}

results <- do.call(rbind, results_list)

cat(sprintf("Combined %d of %d array tasks (%d rows) into %s\n",
            length(Filter(Negate(is.null), results_list)), nrow(config),
            nrow(results), out_file))

write.csv(results, file = out_file, row.names = FALSE)
