# Analysis code

* `AICvsSTANfigures.R`: This R script generates figures comparing the performance of stepwise AIC model selection and Bayesian variable selection using the horseshoe prior (via STAN) in time series regression settings. It includes:
    - Figure 1: A simulated example illustrating prediction intervals and coefficient recovery for both AIC and STAN-based models.
    
    - Figure 2: Violin plots comparing true positive rate (TPR), true negative rate (TNR), and RMSE across varying levels of predictor correlation.

    - Figure 3: Violin plots comparing prediction RMSE under covariate shift and under assumptions of stationarity.

    - The script relies on precomputed simulation results that can be replicated using scripts in the `Simulations/` directory. Outputs are saved as PDF figures in the `Figures/` directory.

-   `AICvsSTANfigures.R`: This R script generates figures comparing the performance of AIC (Akaike Information Criterion) and a Bayesian sparse regression in model selection and prediction accuracy.

    -   **Author**: Amy Patterson

    -   **Upstream scripts**

        -   `/Simulations/AICvsStanRMSE.R` (used to create panel 1)

    -   **Data files**

        -   `/Simulations/AICvsSTANcorrResults.csv` (on MedBow)

        -   `/Simulations/AICvsSTANdecorrelation.csv` (on MedBow)

-   `data_ex.R`: Data exploration and detrending of empirical tree growth data.

    -   **Author**: Kaitlyn McKnight
    
    - **Data files**
    
        - ca726-rwl-noaa.txt (unknown origin)
        
        - ca719-rwl-noaa.txt (unknown origin)

    -   **To do**:

        - [] Let's fix the filepaths in this script so they use relative paths and the `here` package.
        
        - [] Decide if/how to make data available for reproducibility.

-   `datatrends.R`: This script contains investigation of trends, breakpoints, and correlations in the PRISM climate data.

    -   **Author**: Kaitlyn McKnight
    
    - **Data files**
    
        - `/SparseTS_prismdata/ca719_BM_Seq_1800_2012.rds` 
        
        - `/SparseTS_prismdata/BMG_tmin.csv`
        
        - `/SparseTS_prismdata/BMG_tmax.csv`
        
        - `/SparseTS_prismdata/BMG_td_mean.csv`
        
        - `/SparseTS_prismdata/BMG_ppt.csv`
        
        - `/SparseTS_prismdata/BMG_tmean.csv`
        
        - `/SparseTS_prismdata/BMG_vpdmin.csv`
        
        - `/SparseTS_prismdata/BMG_vpdmax.csv`

    -   **To do**:

        - [] Let's fix the filepaths in this script so they use the `here` package.

-   `plot_ARp_simulations.R`: This script plots summary data from AR-p beta-p simulation runs. These simulations iterate through different time series lengths and different standard deviations for the random innovations where there are sparse covariates, lags in one covariate, and sparse AR terms comparing the fits of regularized and not regularized models.

    - **Author**: Alice Carter
    
    - **Data files**:
    
        - `/Data/aquatic_sim_data/ARp_err_sims_02_12_condensed.csv` (MedBow?)
        
        - `/Data/aquatic_sim_data/ARp_err_sims_01_03_condensed.csv` (MedBow?)
        
        - `/Data/aquatic_sim_data/ARp_err_sims_10_31_condensed.csv` (MedBow?)
        
        - `/Data/aquatic_sim_data/ARp_err_sims_04_22_condensed.csv` (MedBow?)
        
    - **To Do**
    
        - [ ] Which version (if either) is current? `plot_ARp_simulations.R` or `plot_ARp_simulations2.R`. Archive one or both?
        
- `plot_individual_ARp_simulation.R`: This code plots parameter estimates comparing models fit with different priors from a single simulation.

    - **Author**: Alice Carter?
    
    - **Data files**
    
        - `/Data/aquatic_sim_data/test/simdat_run12.rds` (MedBow?)
        
    - **To Do**
    
        - [] Is this one current? Archive?
        
- `PrismDat_Agg.R`: This script contains different aggregations of prism data to split each variable into summer (June-Sept), winter(Oct-May), and wateryear (Oct-Sept), then take mean daily mins, maxes, or means.

    - **Author**: Kaitlyn McKnight
    
    - **Date files**: 
    
        - `/SparseTS_prismdata/ca719_BM_Seq_1800_2012.rds` 
        
        - `/SparseTS_prismdata/BMG_tmin.csv`
        
        - `/SparseTS_prismdata/BMG_tmax.csv`
        
        - `/SparseTS_prismdata/BMG_td_mean.csv`
        
        - `/SparseTS_prismdata/BMG_ppt.csv`
        
        - `/SparseTS_prismdata/BMG_tmean.csv`
        
        - `/SparseTS_prismdata/BMG_vpdmin.csv`
        
        - `/SparseTS_prismdata/BMG_vpdmax.csv`
        
- `tree_dat_forecasts_env_covariates.R`: This script compares forecasts from sparse model fits to those using stepwise AIC. The tree growth data are split into two groups (all years, and just those years after the year identified in the breakpoint analysis in `/Analyses/datatrends.R`). Each covariate created in `/Analyses/PrismDat_Agg.R` is used, along with 5 lags from each covariate. 

    -**Author**: Dusty Gannon
    
    - **Data files**
    
        - `/SparseTS_prismdata/ca719_BM_Seq_1800_2012.rds`
        
        - `/SparseTS_prismdata/prism_wateryear.rds`
        
        - `/SparseTS_prismdata/prism_winter.rds`
        
        - `/SparseTS_prismdata/prism_summer.rds`
        
    - **Upstream scripts**
    
        - `/Stan/sparse_reg_FHS.stan`
        
- `tree_dat_forecasts_fourier.R`: This script compares forecasts using Bayesian smoothing of Fourier Basis terms and sparse selection of autoregressive terms to forecasts using the [Hyndman-Khandakar algorithm](https://www.jstatsoft.org/article/view/v027i03). For these forecasts, it is assumed that no environmental features are available.

    - **Author**: Dusty Gannon
    
    - **Data files** 
    
        - `/SparseTS_prismdata/ca719_BM_Seq_1800_2012.rds`
        
        - `/SparseTS_prismdata/prism_wateryear.rds`
        
        - `/SparseTS_prismdata/prism_winter.rds`
        
        - `/SparseTS_prismdata/prism_summer.rds`
        
    - **Upstream scripts**
    
        - `/Stan/AR-p_err3_FHS_DG.stan`
