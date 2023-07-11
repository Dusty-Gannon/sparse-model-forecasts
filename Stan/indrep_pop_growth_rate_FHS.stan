////////////////////////////////////////////////////////////////////
// This Stan program fits a Ricker competition model with
// regularized horseshoe priors on the effects for species
// assuming log-normal demographic stochasticity. The response,
// y[t], is equal to log(N[t+1] / N[t]), and the predictors are
// indexed by t. Such that t = 1, 2, ..., T - 1, where T is the
// length of the original time series. We treat stretches of 0
// density as "missing data" to handle the undefined steps of
// infinite growth
////////////////////////////////////////////////////////////////////

data {

  int<lower = 0> N;                   // length of the sliced time series
  int<lower = 0> P0;                  // number of non-shrinking effects
  int<lower = 0> P;                   // number of heterospecific species effects
  int<lower = 0> K;                   // number of replicate time series
  vector[N] y;                        // vector of responses
  int<lower = 0, upper = 1> z[N];     // indicator vector for whether an observation is non-missing
  vector<lower = 0>[N] dens_foc;      // density of the focal at time t
  matrix[N, P0] X_alpha;              // model matrix for non-shrinking effects
  matrix[N, P] X_beta;                // standardized model matrix for shrinking effects
  int<lower = 0, upper = K> site[N];  // Indicator variable for the site / replicate ts
  real<lower = 0> error_scl;          // prior guess for the scale of demographic stochasticity
  real<lower = 0> tau0;               // scale for global shrinkage parameter
  real<lower = 0> slab_scl;           // scale for non-zero coefficients
  real<lower = 0> slab_df;            // degrees of freedom for non-zero coefficients

}

transformed data{

  // transformations for horseshoe priors
  real slab_scl2 = square(slab_scl);
  real half_slab_df = 0.5 * slab_df;

}


parameters{

  vector[P0] alpha_std;                // standardized intra-specific competition
  vector[P] beta_std;                  // standardized coefficients before shrinkage
  real<lower = 0> lambda;              // intrinsic growth of the focal species
  real<lower = 0> sigma;               // demographic stochasticity
  vector[K] u_std;                     // unscaled site-level effects
  real<lower = 0> sigma_u;             // scale of site-level effects

  // parameters for shrinkage priors
  vector<lower = 0>[P] local_scale;    // non-regularized local scale
  real<lower = 0> c2_std;              // unscaled version of c2
  real<lower = 0> tau_std;             // unscaled version of tau

}


transformed parameters{

  vector[N] eta;                       // declare vector of linear predictors

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
  for(t in 1:N){
    // mean   //intrinsic growth  //intra       //inter
    eta[t] = log(lambda) + X_alpha[t, ] * alpha + X_beta[t, ] * beta +
      u_std[site[t]] * sigma_u;

  }

}


model{

  // priors
  alpha_std ~ std_normal();
  beta_std ~ std_normal();
  u_std ~ std_normal();
  lambda ~ gamma(2, 2);
  sigma ~ normal(0, error_scl);
  sigma_u ~ normal(0, 1);


  tau_std ~ cauchy(0, 1);
  local_scale ~ cauchy(0, 1);
  c2_std ~ inv_gamma(half_slab_df, half_slab_df);

  // likelihood
   for(t in 1:N){
     if(z[t] == 1){
       target += normal_lpdf(y[t] | eta[t], sigma / sqrt(dens_foc[t]));
     }
   }

}



