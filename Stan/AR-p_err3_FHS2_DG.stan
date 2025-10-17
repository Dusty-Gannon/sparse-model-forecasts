/////////////////////////////////////////////////////////////
// This Stan program fits an AR(p) Gaussian time-series model
// with the potential for covariates. Both of the parameter
// vectors, phi and beta, get sparsity-inducing horseshoe
// priors.
/////////////////////////////////////////////////////////////

#include ar_var_approx.stan

data{

  int<lower = 1> N;           // number of observations
  int<lower = 0> P;           // total number of covariates
  int<lower = 0> P_0;         // number of un-regularized covariates
  int<lower = 0> p;           // guess for the max order of the autoregressive process
  vector[N] y;                // vector of responses
  matrix[N, P] X;             // model matrix for shrinking effects

  // prior inputs for regularization
  real<lower = 0> tau0_phi;        // global shrinkage for phi
  real<lower = 0> slab_scl_phi;    // scale for non-zero coefficients in phi
  real<lower = 0> slab_df_phi;     // degrees of freedom for non-zero coefficients
  real<lower = 0> tau0_beta;       // global shrinkage for beta
  real<lower = 0> slab_scl_beta;   // scale for non-zero coefficients in phi
  real<lower = 0> slab_df_beta;    // degrees of freedom for non-zero coefficients

  // forecasting objects
  int<lower = 0> N_new;       // number of observations to forecast
  matrix[N_new, P] X_new;     // model matrix for forecasting

}

transformed data{

  row_vector[p] ones = rep_row_vector(1, p);     // row vector of ones

  // transformations for horseshoe priors
  real slab_scl2_phi = square(slab_scl_phi);
  real half_slab_df_phi = 0.5 * slab_df_phi;
  real slab_scl2_beta = square(slab_scl_beta);
  real half_slab_df_beta = 0.5 * slab_df_beta;

}

parameters{

  vector[P_0] alpha_std;      // unshrunk coefficients
  vector[P - P_0] beta_std;   // standardized coefficients before shrinkage
  vector[p] phi_std;          // autoregression parameters
  real<lower = 0> sigma;      // sd of the innovations

  // parameters for shrinkage priors
  vector<lower = 0>[P - P_0] local_scale_beta;    // non-regularized local scale for beta
  vector<lower = 0>[p] local_scale_phi;           // non-regularized local scale for phi
  real<lower = 0> c2_std_phi;                     // unscaled version of c2
  real<lower = 0> tau_std_phi;                    // unscaled version of tau
  real<lower = 0> c2_std_beta;                    // unscaled version of c2
  real<lower = 0> tau_std_beta;                   // unscaled version of tau

}
transformed parameters{

  vector[N] mu;               // declare vector of means
  vector[N] err;              // vector of residuals
  vector[P] beta;             // all scaled regression coefficients
  vector[p] phi;              // scaled autoregression coefficients

  // scale c2: c2 ~ inv_gamma(half_slab_df, half_slab_df * slab_scl2)
  real c2_phi = slab_scl2_phi * c2_std_phi;
  real c2_beta = slab_scl2_beta * c2_std_beta;

  // tau ~ cauchy(0, tau0)
  real tau_phi = tau0_phi * tau_std_phi;

  // This calculation follows equation 2.8 in Piironen and Vehtari 2017
  vector[p] local_scale_tilde_p =
    sqrt(c2_phi * square(local_scale_phi) ./ (c2_phi + square(tau_phi) * square(local_scale_phi)));

  // scale phis
  phi = tau_phi * local_scale_tilde_p .* phi_std;

  // same process for betas, but error var depends on phi
  real tau_beta = tau0_beta * tau_std_beta * sqrt(ar_var(phi, pow(sigma, 2), 20));

  vector[P - P_0] local_scale_tilde_b =
    sqrt(c2_beta * square(local_scale_beta) ./ (c2_beta + square(tau_beta) * square(local_scale_beta)));

  beta[1:P_0] = alpha_std * 2.5;
  beta[(P_0 + 1):P] = tau_beta * local_scale_tilde_b .* beta_std;

  // define means and error
  mu = X * beta;
  err = y - mu;

}

model{

  // define random innovations
  vector[N] epsilon;

  // define reverse order to make AR coefficient order match tradition
  vector[p] err_rev;

  // priors
  beta_std ~ std_normal();
  alpha_std ~ std_normal();
  phi_std ~ std_normal();

  tau_std_phi ~ cauchy(0, 1);
  tau_std_beta ~ cauchy(0, 1);
  local_scale_beta ~ cauchy(0, 1);
  local_scale_phi ~ cauchy(0, 1);
  c2_std_phi ~ inv_gamma(half_slab_df_phi, half_slab_df_phi);
  c2_std_beta ~ inv_gamma(half_slab_df_beta, half_slab_df_beta);

  sigma ~ normal(0, 2);

  // complete the AR process
  epsilon[1:p] = rep_vector(0, p);
  for(t in (p + 1):N){
    for(i in 1:p){
      err_rev[i] = err[t - i];
    }
    epsilon[t] = err[t] - err_rev' * phi;
  }

  // likelihood
  epsilon[(p + 1):N] ~ normal(0, sigma);
}


generated quantities{

  // post. pred. sampling plus forecasts
  real y_rep[N + N_new];
  vector[N + N_new] err_rep;                 // extend the error term for the forecasting
  vector[p] err_rep_rev = rep_vector(0, p);  // initialize reverse order for consistency

  err_rep[1:N] = err;

  // fitted values
  for(t in 1:p){
    y_rep[t] = mu[t];
  }
  for(t in (p + 1):N){
    for(i in 1:p){
      err_rep_rev[i] = err_rep[t - i];
    }
    y_rep[t] = mu[t] + err_rep_rev' * phi;
  }

  // forecast
  if(N_new > 0){
    for(t in (N + 1):(N + N_new)){
      for(i in 1:p){
        err_rep_rev[i] = err_rep[t - i];
      }
      err_rep[t] = err_rep_rev' * phi + normal_rng(0, sigma);
      y_rep[t] = X_new[t - N, ] * beta + err_rep[t];
    }
  }

}

