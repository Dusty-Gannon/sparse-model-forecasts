# Analysis code

This directory contains R scripts and slurm scripts used to fit models to simulated and empirical data.

* `fit_lnorm_ricker_mods_parallel.R`: Wrapper script meant to be run from the command line on Beartooth. This script does the bulk of the analyses we run for the simulated competitive communities. The command line arguments are ordered and include:
  
    - `<input_data.rds>`: An R data file that is a list with at least the named matrix `N`, which has species in rows and their abundances through time in columns. This is used to generate the response vector (abundances of the first species through time, arbitrarily) and the covariate matrix `t(N[-1, ])`. The data should be stored in `Data/terrestrial_sim_data/lnorm_ricker/`.
  
    - `<start>`: Numeric time step at which to start using the simulated data. This allows for a 'burnin' period for the community to find a dynamic equilibrium.
    
    - `<stop>`: Numeric time step at which to stop using data. This allows the user to control the length of the time series. All simulated datasets go from `t=1,...,500`.
    
    - `<index1>`
    
    - `<index2>`: Together, `<index1>` and `<index2>` allow the user to partition all simulated datasets in `<input_data.rds>` into smaller chunks to reduce memory usage. For example, providing `1` and `10` as arguments will fit the models to the first 10 datasets in `<input_data.rds>` in parallel. An example of splitting all 40,200 simulated communities for the disturbance experiments into smaller chunks and fitting the model using a SLURM array job can be found in the SLURM script `fun_fit_lnorm_ricker_dist_freq_x_nsp_x_int.sh`.
    
    - `<outfile>`: Path and filename for the outfile, a list wih confusion matrices and posterior draws from the shrunk model coefficients. The path should start below `/Data/terrestrial_sim_data/lnorm_ricker/`. For example, providing `results/out.rds` would look for a directory called `results` in `/Data/terrestrial_sim_data/lnorm_ricker/` in which to put `out.rds`.
    
    - OPTIONAL `dist`: For the disturbance simulations, we need to add an indicator covariate for when the focal species was disturbed to attribute reductions in growth due to disturbance. For those simulations, provide `dist` as the seventh and final argument.
  
* `fit_lnorm_ricker_mods_tandem.R`: Wrapper script meant to be run from the command line on Beartooth. This script will fit ricker models to datasets stored in a list. The command line arguments are ordered and include:
  
    - `<input_data.rds>`: An R data file that is a list with at least the named matrix `N`, which has species in rows and their abundances through time in columns. This is used to generate the response vector (abundances of the first species through time, arbitrarily) and the covariate matrix `t(N[-1, ])`. The data should be stored in `Data/terrestrial_sim_data/lnorm_ricker/`.
  
    - `<start>`: Numeric time step at which to start using the simulated data. This allows for a 'burnin' period for the community to find a dynamic equilibrium.
  
    - `<stop>`: Numeric time step at which to stop using data. This allows the user to control the length of the time series. All simulated datasets go from `t=1,...,500`.
  
    - `<outfile>`: Path and filename for the outfile, a list wih confusion matrices and posterior draws from the shrunk model coefficients. The path should start below `/Data/terrestrial_sim_data/lnorm_ricker/`. For example, providing `results/out.rds` would look for a directory called `results` in `/Data/terrestrial_sim_data/lnorm_ricker/` in which to put `out.rds`.

* `fit_lnorm_ricker_vrtests_tandem.R`: R script to fit population growth models to data assuming lognormally-distributed demographic errors. This script will fit models to a list of datasets in tandem.

* `plot_summarize_Ricker_freq_x_nsp_thintests.R`: This script takes `.rds` files from the SLURM array job implemented by `run_fit_lnorm_ricker_thintests_freq_x_nsp.sh` and creates a heatmap of the confusion metrics across the 10 thinning frequency and 10 proportions of the community combinations.

* `plot_summarize_Ricker_thintests.R`: This script takes results from fitting models to simulations in which we only varied the interval between thinning treatments and plots the confusion metrics over those thinning intervals. The fitting is implemented using `run_fit_lnorm_ricker_thintests.sh`, which calls the `fit_lnorm_ricker_mods_tandem.R` script.

* `run_fit_lnorm_ricker_dist_freq_x_nsp_x_int.sh`: This SLURM script takes the 40,200 datasets generated for the disturbance experiments (200 reps x 201 treatments), splits them into 4020 jobs and fits the lnorm ricker model to the 10 datasets in each job in parallel. The results are stored as `.rds` files with suffixes indicating the indexes of the datasets used based on the original input file with all datasets.

* `run_fit_lnorm_ricker_thintests_freq_x_nsp.sh`: This SLURM script takes the 20,100 datasets generated for the disturbance experiments (200 reps x 101 treatments), splits them into 2010 jobs and fits the lnorm ricker model to the 10 datasets in each job in parallel. The results are stored as `.rds` files with suffixes indicating the indexes of the datasets used based on the original input file with all datasets.

* `run_fit_lnorm_ricker_thintests.sh`: Running initial thinning tests in which only the interval between thinning treatments was varied.

* `run_fit_lnorm_ricker_vrtests.sh`: SLURM script to fit the Ricker models to communities for which the ratio of the variance of the focal species demographic stochasticity to variation in heterospecific abundances was manipulated.

* `run_plot_freq_x_nsp_thintests.sh`: SLURM script to create the heatmaps from the thinning interval x proportion of community thinned experiments.


