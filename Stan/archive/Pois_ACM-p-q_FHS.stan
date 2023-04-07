//////////////////////////////////////
// Stan model to fit an autoregressive
// conditional mean model of order p,q
// with Poisson errors
//////////////////////////////////////


data{

  int<lower = 1> N;           // number of observations
  int<lower = 1> P0;          // number of non-shrinking parameters
  int<lower = 0> P;           // number of shrinking parameters
  int<lower = 0> p;           // order of the AR process
  int<lower = 0> q;           // order of the MA process
  int<lower = 0> y[N];        // vector of responses
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

  int<lower=1> m = max(p, q);          // maximum order of the process

}


parameters{

  real<lower = 0> a;                   // intercept term
  vector[P0] alpha_std;                // unshrunk coefficients (self-limitation and intercept)
  vector[P] beta_std;                  // standardized coefficients before shrinkage
  vector<lower = 0>[p] phi;            // autoregression parameter
  vector<lower = 0>[q] theta;          // MA parameter

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
  z[1:m] = rep_vector(a, m);
  mu[1:m] = to_vector(y[1:m]);

  for(t in (m + 1):N){
    z[t] = a + phi' * to_vector(y[(t - p):(t - 1)]) + theta' * mu[(t - p):(t - 1)];
  }

  for(t in (m + 1):N){
    mu[t] = z[t] * exp(X_alpha[t, ] * alpha + X_beta[t, ] * beta);
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
  for(i in (m+1):N){
   y[i] ~ poisson(mu[i]);
  }

}


generated quantities{

  int<lower = 0> y_rep[N];
  for(t in 1:N){
    y_rep[t] = poisson_rng(mu[t]);
  }

}
