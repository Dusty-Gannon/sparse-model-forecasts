

data{

  int<lower = 1> N;           // number of observations
  int<lower = 1> P0;          // number of non-shrinking parameters
  int<lower = 0> P;           // number of shrinking parameters
  int<lower = 0> p;           // order of the autoregressive process
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
  vector<lower = -1, upper = 1>[p] phi;    // autoregression parameter
  real<lower = 0> sigma;                   // sd of the innovations

  // parameters for shrinkage priors
  vector<lower = 0>[P] local_scale;    // non-regularized local scale
  real<lower = 0> c2_std;              // unscaled version of c2
  real<lower = 0> tau_std;             // unscaled version of tau

}


transformed parameters{

  vector[N] mu;               // declare vector of means

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

  mu[1:p] = rep_vector(0, p);

  for(t in (p + 1):N){
    mu[t] = X_alpha[t, ] * alpha + X_beta[t, ] * beta + (y[(t - p):(t - 1)])' * phi;
  }

}


model{

  // priors
  beta_std ~ std_normal();
  alpha_std ~ std_normal();

  tau_std ~ cauchy(0, 1);
  local_scale ~ cauchy(0, 1);
  c2_std ~ inv_gamma(half_slab_df, half_slab_df);

  phi ~ uniform(-1, 1);
  sigma ~ cauchy(0, 1);

  // likelihood
  y[(p + 1):N] ~ normal(mu[(p + 1):N], sigma);


}


generated quantities{

  // post. pred. sampling
  real y_rep[N - p] = normal_rng(mu[(1 + p):N], sigma);

  // residuals
  vector[N - p] resid = y[(p + 1):N] - mu[(p + 1):N];

}

