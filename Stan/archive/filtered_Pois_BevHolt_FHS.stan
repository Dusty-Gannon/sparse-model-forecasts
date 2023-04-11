////////////////////////////////////////////////////////////////////
// This Stan program fits a Beverton-Holt model in which the
// total fecundity for year t is filtered by survival and germination
// to result in the abundace in year t + 1.
//
// Regularized horseshoe priors are put on the effects for species
// assuming conditionally Poisson-distributed demographic noise
/////////////////////////////////////////////////////////////////////

data {

  int<lower = 0> N;             // length of the time series
  int<lower = 0> P0;            // number of non-shrinking effects
  int<lower = 0> P;             // number of shrinking effects
  int<lower = 0> fec[N];        // vector of seed sets over time
  int<lower = 0> y[N];          // vector of cover over time
  vector<lower = 0>[N] y_star;  // augmented data
  matrix[N, P0] X_alpha;        // model matrix for unshrunk effects
  matrix[N, P] X_beta;          // model matrix for shrinking effects
  real<lower = 0> tau0;         // scale for global shrinkage parameter
  real<lower = 0> slab_scl;     // scale for non-zero coefficients
  real<lower = 0> slab_df;      // degrees of freedom for non-zero coefficients

}

transformed data{

  // transformations for horseshoe priors
  real slab_scl2 = square(slab_scl);
  real half_slab_df = 0.5 * slab_df;

}


parameters{

  vector[P0] alpha_std;                // unshrunk coefficients (self-limitation and intercept)
  vector[P] beta_std;                  // standardized coefficients before shrinkage
  real<lower = 0, upper = 1> p;        // transition probability from seed to adult

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
  vector[P0] alpha = alpha_std * 5;

  // construct linear predictors
  eta[1] = log(y_star[1]);
  for(t in 2:N){

    eta[t] = X_alpha[t,] * alpha + X_beta[t, ] * beta + log(y_star[t - 1]);

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
  for(t in 2:N){
   fec[t] ~ normal(eta[t], sigma) * y[t - 1];
   y[t] ~ poisson(p * fec[t]);
  }

}


// generated quantities{
//
//   int<lower = 0> y_rep[N];
//   for(t in 1:N){
//     y_rep[t] = poisson_log_rng(eta[t]);
//   }
//
// }

