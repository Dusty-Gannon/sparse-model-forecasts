/////////////////////////////////////////////////////////////
// This Stan program fits an AR(p) Gaussian time-series model
// with a model selection component on both the order of the
// AR process (i.e., which of the phi_1, phi_2,..., phi_p)
// lags are important) and on the covariates using Finnish
// Horseshoe priors.
/////////////////////////////////////////////////////////////


data{

  int<lower = 1> N;           // number of observations
  int<lower = 1> P0;          // number of non-shrinking parameters
  int<lower = 0> P;           // number of shrinking parameters
  int<lower = 0> p;           // guess for the max order of the autoregressive process
  vector[N] y;                // vector of responses
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

  row_vector[p] ones = rep_row_vector(1, p);     // row vector of ones

}


parameters{

  vector[P0] alpha_std;                    // unshrunk coefficients (self-limitation and intercept)
  vector[P] beta_std;                      // standardized coefficients before shrinkage
  vector[p] phi_std;                           // autoregression parameters
  real<lower = 0> sigma;                   // sd of the innovations

  // parameters for shrinkage priors
  vector<lower = 0>[P] local_scale_beta;    // non-regularized local scale for beta
  vector<lower = 0>[p] local_scale_phi;     // non-regularized local scale for phi
  real<lower = 0> c2_std;                   // unscaled version of c2
  real<lower = 0> tau_std;                  // unscaled version of tau

}


transformed parameters{

  vector[N] mu;               // declare vector of means
  vector[N] err;              // vector of residuals
  vector[N-p] epsilon;
  // scale c2: c2 ~ inv_gamma(half_slab_df, half_slab_df * slab_scl2)
  real c2 = slab_scl2 * c2_std;

  // tau ~ cauchy(0, tau0)
  real tau = tau0 * tau_std;

  // This calculation follows equation 2.8 in Piironen and Vehtari 2017
  vector[P] local_scale_tilde_b =
    sqrt(c2 * square(local_scale_beta) ./ (c2 + square(tau) * square(local_scale_beta)));
  vector[p] local_scale_tilde_p =
    sqrt(c2 * square(local_scale_phi) ./ (c2 + square(tau) * square(local_scale_phi)));

  // scale betas
  vector[P] beta = tau * local_scale_tilde_b .* beta_std;

  // scale phis
  vector[p] phi = tau * local_scale_tilde_p .* phi_std;

  // scale alpha
  vector[P0] alpha = alpha_std * 5;

  // assume no error for first p observations
  mu = X_alpha * alpha + X_beta * beta;
  err = y - mu;

  // complete the AR process
  for(t in (p + 1):N){
    epsilon[t-p] = err[t] - err[(t-p):(t-1)]' *phi;
  }

}


model{

  // priors
  beta_std ~ std_normal();
  alpha_std ~ std_normal();
  phi_std ~ std_normal();

  tau_std ~ cauchy(0, 1);
  local_scale_beta ~ cauchy(0, 1);
  local_scale_phi ~ cauchy(0, 1);
  c2_std ~ inv_gamma(half_slab_df, half_slab_df);

  sigma ~ normal(0, 2);

  // likelihood
  epsilon ~ normal(0, sigma);
}


generated quantities{

  // post. pred. sampling
  real y_rep[N - p] = normal_rng(err[(1 + p):N] + mu[(1 + p):N], sigma);

  // residuals
  vector[N - p] resid = y[(p + 1):N] - mu[(p + 1):N];

}

