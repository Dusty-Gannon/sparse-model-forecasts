# Directory contents

-   `AIC_modelSelection.R`: Basic code to make a set of time series and do stepwise AIC based model selection on them *** maybe not needed any more

-   `AICvsStanBatch.R`: Run a bunch of time series through AIC vs Horseshoe prior methods and recover the RMSE values, set up to run on an HPC

-   `AICvsStanConfig.R`: File to make the config file needed in AICvsStanBatch.R, enter the different parameters of interest eg time series length, number of predictors, etc here

-   `AICvsStanConfig.txt`: The config file output from AICvsStanConfig.R

-   `AICvsSTANcorrResults.csv`: The result table for the case of constant correlations, based on runs listed in AICvsStanConfig.txt and defined in AICvsStanConfig.R

-   `AICvsSTANdecorrelation.csv`: The result table for the case of covariate shift / decorrelation, based on runs listed in DecorrConfig.txt and definted in DecorrConfig.R

-   `AICvsStanRMSE.R`: File of R functions we used to run the AIC vs Horsehoe analyses found in AICvsStanBatch.R *** maybe should be put in the R folder?

-   `AR-p_beta_p_model_sims.R`: Amy is unsure what this is

-   `combine_seasonal_ARp_sim_output.R`: Appears to be a file to summarize results from the HPC

-   `DecorrConfig.R`: File to make the config file needed in AICvsStanBatch.R for the decorrelation cases, enter the different parameters of interest eg time series length, number of predictors, etc here

-   `DecorrConfig.txt`: The config file output from DecorrConfig.R

-   `fit_AR-p_seasonal_model_sims.R`: Amy is unsure what this is

-   `run_AICvsSTAN.sh`: Runs AICvsStanBatch.R on the HPC, requires the config file AICvsStanConfig.txt

-   `run_ARp_sim_compiler.sh`: Amy is unsure why there are so many compiler versions and if they all need to be kept

-   `run_ARp_sim_compiler2.sh`: Amy is unsure why there are so many compiler versions and if they all need to be kept

-   `run_ARp_sim_compiler3.sh`: Amy is unsure why there are so many compiler versions and if they all need to be kept

-   `run_ARp_sim_compiler4.sh`: Amy is unsure why there are so many compiler versions and if they all need to be kept

-   `run_Decorr.sh`: Runs AICvsStanBatch.R on the HPC, requires the config file DecorrConfig.txt

-   `run_seasonal_ARp_model.sh`: Appears to run the aquatic models




# Archived - still needs to be listed

-   `fit_AR-p_seasonal_model_sims2.R`: 

-   `lnorm_Ricker_model_sims_freq_x_int_dist_s4t.R`: 

-   `lnorm_Ricker_model_sims_freq_x_int_dist.R`: 

-   `lnorm_Ricker_model_sims_freq_x_nsp_thin_s4t.R`: 

-   `lnorm_Ricker_model_sims_freq_x_nsp_thin.R`:

-   `lnorm_Ricker_model_sims_thin.R`:

-   `lnorm_Ricker_model_sims_vrtests_1param.R`:

-   `lnorm_Ricker_model_sims_vrtests.R`:

-   `lottery_model_sims.R`:

-   `Ricker_model_sims.R`:

-   `run_ARp_beta_model_wrapper.sh`:

-   `run_ARp_beta_model.sh`:

-   `run_lottery_model_1.sh`:

-   `run_lottery_model_2.sh`:

-   `run_seasonal_ARp_model2.sh`:

-   `run_sims_freq_x_int_dist_s4t.sh`:

-   `run_sims_freq_x_int_dist.sh`:

-   `run_sims_freq_x_nsp_thin_s4t.sh`:

-   `split_sim_files.R`:

-   `teton_lottery_model_sims.R`:
