////////////////////////////////////////////////////////////////////
// This Stan program defines a Bayesian sparse regression model with
// regularized Horseshoe priors for the regression coefficients
////////////////////////////////////////////////////////////////////

data {

  int<lower = 0> N;                 // length of the sliced time series
  int<lower = 0> P0;                // number of non-shrinking effects
  int<lower = 0> P;                 // total number of explanatory variables (including intercept)
  vector[N] y;                      // vector of responses
  matrix[N, P] X;                   // standardized model matrix for shrinking effects
  real<lower = 0> tau0;             // scale for global shrinkage parameter
  real<lower = 0> slab_scl;         // scale for non-zero coefficients
  real<lower = 0> slab_df;          // degrees of freedom for non-zero coefficients

  // objects for prediction
  int<lower = 0> N_new;             // number of predictions
  matrix[N_new, P] X_new;           // model matrix for predictions

}

transformed data{

  // transformations for horseshoe priors
  real slab_scl2 = square(slab_scl);
  real half_slab_df = 0.5 * slab_df;

}


parameters{

  vector[P0] alpha_std;                // unregularized effects
  vector[P - P0] beta_std;             // standardized coefficients before shrinkage
  real<lower = 0> sigma;               // demographic stochasticity

  // parameters for shrinkage priors
  vector<lower = 0>[P - P0] local_scale;    // non-regularized local scale
  real<lower = 0> c2_std;              // unscaled version of c2
  real<lower = 0> tau_std;             // unscaled version of tau

}


transformed parameters{

  vector[P] beta;                      // declare vector of scaled regression coefficients
  vector[N] eta;                       // declare vector of linear predictors

  // scale c2: c2 ~ inv_gamma(half_slab_df, half_slab_df * slab_scl2)
  real c2 = slab_scl2 * c2_std;

  // tau ~ cauchy(0, tau0 * sigma)
  real tau = tau0 * tau_std * sigma;

  // This calculation follows equation 2.8 in Piironen and Vehtari 2017
  vector[P - P0] local_scale_tilde =
    sqrt(c2 * square(local_scale) ./ (c2 + square(tau) * square(local_scale)));

  // scale betas
  beta[1:P0] = alpha_std * 5;
  beta[(P0 + 1):P] = tau * local_scale_tilde .* beta_std;

  // construct linear predictors
  eta = X * beta;

}


model{

  // priors
  alpha_std ~ std_normal();
  beta_std ~ std_normal();
  sigma ~ cauchy(0, 2.5);


  tau_std ~ cauchy(0, 1);
  local_scale ~ cauchy(0, 1);
  c2_std ~ inv_gamma(half_slab_df, half_slab_df);

  // likelihood
  y ~ normal(eta, sigma);

}

generated quantities{

  vector[N + N_new] y_rep;

  // posterior predictions with same linear predictor as observed
  for(i in 1:N){
    y_rep[i] = normal_rng(eta[i], sigma);
  }

  // optional predictions to unobserved
  if(N_new > 0){
    for(i in (N + 1):(N + N_new)){
      y_rep[i] = normal_rng(X_new[i - N, ] * beta, sigma);
    }
  }

}
