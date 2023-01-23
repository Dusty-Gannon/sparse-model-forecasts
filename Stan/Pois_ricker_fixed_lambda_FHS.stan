////////////////////////////////////////////////////////////////////////
// This Stan program fits a Ricker competition model with
// regularized horseshoe priors on the effects for species
// assuming conditionally Poisson-distributed demographic noise.
// The growth rate, lambda, is unknown and assumed fixed through time,
// but some prior data exist such that an informative prior can be used.
////////////////////////////////////////////////////////////////////////

data {

  int<lower = 0> N;             // length of the time series
  int<lower = 0> P_h;           // number of heterospecific species
  int<lower = 0> P;             // number of shrinking effects
  real<lower = 0> scl_X;        // scaling factor
  int<lower = 0> y[N];          // vector of responses
  real<lower = 0> a_lambda;     // prior shape parameter for lambda
  real<lower = 0> b_lambda;     // prior rate parameter for lambda
  matrix[N, P_h] X_beta0;       // model matrix for generic effect
  matrix[N, P] X_beta;          // model matrix for shrinking effects
  real<lower = 0> beta0_scl;    // prior scale for generic effects
  real<lower = 0> tau0;         // scale for global shrinkage parameter
  real<lower = 0> slab_scl;     // scale for non-zero coefficients
  real<lower = 0> slab_df;      // degrees of freedom for non-zero coefficients

}

transformed data{

  // transformations for horseshoe priors
  real slab_scl2 = square(slab_scl);
  real half_slab_df = 0.5 * slab_df;

  vector[P_h] ones = rep_vector(1, P_h);

}


parameters{

  real alpha_std;                      // standardized intra-specific competition
  real<upper = 0> beta0_std;           // standardized generic effect
  vector[P] beta_std;                  // standardized coefficients before shrinkage
  real<lower = 0> lambda;              // intrinsic growth

  // parameters for shrinkage priors
  vector<lower = 0>[P] local_scale;    // non-regularized local scale
  real<lower = 0> c2_std;              // unscaled version of c2
  real<lower = 0> tau_std;             // unscaled version of tau

}


transformed parameters{

  vector[N] eta;                // declare vector of linear predictors

  // scale c2: c2 ~ inv_gamma(half_slab_df, half_slab_df * slab_scl2)
  real c2 = slab_scl2 * c2_std;

  // tau ~ cauchy(0, tau0)
  real tau = tau0 * tau_std;

  // This calculation follows equation 2.8 in Piironen and Vehtari 2017
  vector[P] local_scale_tilde =
    sqrt(c2 * square(local_scale) ./ (c2 + square(tau) * square(local_scale)));

  // scale betas
  vector[P] beta = tau * local_scale_tilde .* beta_std;

  // scale alpha
  real alpha = alpha_std * 0.5;

  // scale beta0
  real beta0 = beta0_std * beta0_scl;

  // construct linear predictors
  eta[1] = log(y[1]);
  for(t in 2:N){
    // mean   time-varying growth         intra                     generic                       non-generic          offset
    eta[t] = log(lambda) + (y[t - 1]/scl_X) * alpha + (X_beta0[t - 1, ] * ones) * beta0 + X_beta[t - 1, ] * beta + log(y[t - 1]);

  }

}


model{

  // priors
  alpha_std ~ std_normal();
  beta_std ~ std_normal();
  beta0_std ~ std_normal();
  lambda ~ gamma(a_lambda, b_lambda);


  tau_std ~ cauchy(0, 1);
  local_scale ~ cauchy(0, 1);
  c2_std ~ inv_gamma(half_slab_df, half_slab_df);

  // likelihood
  for(t in 2:N){
   y[t] ~ poisson_log(eta[t]);
  }

}



