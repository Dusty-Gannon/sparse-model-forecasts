

data{

  int<lower = 1> N;           // number of observations
  int<lower = 1> P0;          // number of non-shrinking parameters
  int<lower = 0> P;           // number of shrinking parameters
  int<lower = 0> p;           // order of the autoregressive process
  int<lower = 0> q;           // order of the MA process
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

  int m = max(p, q);                             // max order of the ARMA process
  row_vector[p] ones = rep_row_vector(1, p);     // row vector of ones

}


parameters{

  vector[P0] alpha_std;                    // unshrunk coefficients (self-limitation and intercept)
  vector[P] beta_std;                   // standardized coefficients before shrinkage
  vector<lower = -1, upper = 1>[p] phi;    // autoregression parameter
  vector<lower = -1, upper = 1>[q] theta;  // MA parameter
  real<lower = 0> sigma;                   // sd of the innovations

  // parameters for shrinkage priors
  vector<lower = 0>[P] local_scale;    // non-regularized local scale
  real<lower = 0> c2_std;              // unscaled version of c2
  real<lower = 0> tau_std;             // unscaled version of tau

}


transformed parameters{

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

}


model{

  // declare some extra variables for convenience
  vector[q] epsilon;
  real epsilon_t;
  vector[N] mu;

  // priors
  beta_std ~ std_normal();
  alpha_std ~ std_normal();

  tau_std ~ cauchy(0, 1);
  local_scale ~ cauchy(0, 1);
  c2_std ~ inv_gamma(half_slab_df, half_slab_df);

  phi ~ uniform(-1, 1);
  theta ~ uniform(-1, 1);
  sigma ~ cauchy(0, 1);

  mu = X_alpha * alpha + X_beta * beta;

  // assume no error for the initial values
  epsilon = rep_vector(0, q);

  for (t in (m + 1):N) {
    epsilon_t = y[t] - mu[t] - phi' * (y[(t - p):(t - 1)]) - theta' * epsilon;
    epsilon_t ~ normal(0, sigma);
    epsilon[2:q] = epsilon[1:(q - 1)];
    epsilon[1] = epsilon_t;
  }


}



