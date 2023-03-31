///////////////////////////////////////////////////////////////
// This Stan program fits a Ricker competition model with
// regularized horseshoe priors on the effects for species
// assuming log-normal demographic stochasticity
///////////////////////////////////////////////////////////////

data {

  int<lower = 0> N;             // length of the time series
  int<lower = 0> P0;            // number of non-shrinking effects
  int<lower = 0> P;             // number of heterospecific species effects
  vector<lower = 0>[N] y;       // vector of responses
  matrix[N, P0] X_alpha;        // model matrix for non-shrinking effects
  matrix[N, P] X_beta;          // standardized model matrix for shrinking effects
  real<lower = 0> error_scl;    // prior guess for the scale of demographic stochasticity
  real<lower = 0> tau0;         // scale for global shrinkage parameter
  real<lower = 0> slab_scl;     // scale for non-zero coefficients
  real<lower = 0> slab_df;      // degrees of freedom for non-zero coefficients

}

transformed data{

  // transformations for horseshoe priors
  real slab_scl2 = square(slab_scl);
  real half_slab_df = 0.5 * slab_df;

  // convert counts to growth rates
  vector[N - 1] r = log(y[2:N] ./ y[1:(N - 1)]);

}


parameters{

  vector[P0] alpha_std;                // standardized intra-specific competition
  vector[P] beta_std;                  // standardized coefficients before shrinkage
  real<lower = 0> lambda;              // intrinsic growth of the focal species
  real<lower = 0> sigma;               // demographic stochasticity

  // parameters for shrinkage priors
  vector<lower = 0>[P] local_scale;    // non-regularized local scale
  real<lower = 0> c2_std;              // unscaled version of c2
  real<lower = 0> tau_std;             // unscaled version of tau

}


transformed parameters{

  vector[N - 1] eta;                // declare vector of linear predictors

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
  vector[P0] alpha = alpha_std * 1;

  // construct linear predictors
  for(t in 1:(N - 1)){
    // mean   //intrinsic growth  //intra       //inter
    eta[t] = log(lambda) + X_alpha[t, ] * alpha + X_beta[t, ] * beta;

  }

}


model{

  // priors
  alpha_std ~ std_normal();
  beta_std ~ std_normal();
  lambda ~ gamma(2, 2);
  sigma ~ normal(0, error_scl);


  tau_std ~ cauchy(0, 1);
  local_scale ~ cauchy(0, 1);
  c2_std ~ inv_gamma(half_slab_df, half_slab_df);

  // likelihood
   for(t in 1:(N - 1)){
     r[t] ~ normal(eta[t], sigma / sqrt(y[t]));
   }

}



