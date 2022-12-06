//////////////////////////////////////////////////////////////////////////
// This stan program fits a Poisson GLARMA(1,1) model to an integer-valued
// time series and uses Regularized horseshoe priors for the regression
// coefficients
//////////////////////////////////////////////////////////////////////////

data{

  int<lower = 1> N;           // number of observations
  int<lower = 1> P0;          // number of non-shrinking parameters
  int<lower = 0> P;           // number of shrinking parameters
  int<lower = 0> y[N];        // vector of responses
  real<lower = 0> y_star[N];  // vector of threshold transformed responses
  matrix[N, P0] X_alpha;      // model matrix for unshrunk effects
  matrix[N, P] X_beta;        // model matrix for shrinking effects
  real<lower = 0> tau0;       // scale for global shrinkage parameter
  real<lower = 0> slab_scl;   // scale for non-zero coefficients
  real<lower = 0> slab_df;    // degrees of freedom for non-zero coefficients

}


transformed data{

  // transformations for horseshoe priors
  real slab_scl2 = square(slab_scl);
  real half_slab_df = 0.5 * slab_df;

}


parameters{

  vector[P0] alpha_std;                      // unshrunk coefficients (self-limitation and intercept)
  vector[P] beta_std;                  // standardized coefficients before shrinkage
  real<lower = -1, upper = 1> phi;      // autoregression parameter
  real<lower = -1, upper = 1> theta;    // MA parameter

  // parameters for shrinkage priors
  vector<lower = 0>[P] local_scale;    // non-regularized local scale
  real<lower = 0> c2_std;              // unscaled version of c2
  real<lower = 0> tau_std;             // unscaled version of tau

}


transformed parameters{

  vector[N] mu;                // declare vector of means
  vector[N] z;                 // declare vector of ARMA variables

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
  vector[P0] alpha = alpha_std * 5;

  // construct linear predictors
  mu[1] = y_star[1];
  z[1] = 0;
  for(i in 2:N){

    z[i] = phi * (log(y_star[i-1]) - X_alpha[i-1, ] * alpha - X_beta[i-1, ] * beta) + theta * log(y_star[i-1] / mu[i-1]);
    mu[i] = exp(X_alpha[i, ] * alpha + X_beta[i, ] * beta + z[i]);

  }

}


model{

  // priors
  alpha_std ~ std_normal();
  beta_std ~ std_normal();

  tau_std ~ cauchy(0, 1);
  local_scale ~ cauchy(0, 1);
  c2_std ~ inv_gamma(half_slab_df, half_slab_df);

  // likelihood
  for(i in 2:N){
   y[i] ~ poisson(mu[i]);
  }

}


generated quantities{

  int y_rep[N];     // post pred draws
  y_rep[1] = y[1];
  for(i in 2:N){
    y_rep[i] = poisson_rng(mu[i]);
  }

}

