functions {
  /** ARR2 prior
   *
   * @param phi_z Real
   * @param psi Simplex
   * @param R2 Real
   * @param sigma Real
   * @param sigma_sd Real
   * @param mean_R2 Real
   * @param prec_R2 Real
   * @param cons Vector
   * @param var_y Real
   */
  real arr2_ncp_lpdf(vector phi_z, vector psi, real R2, real sigma,
		 data real sigma_sd, data real mean_R2,
		 data real prec_R2, data vector cons, data real var_y) {
    return normal_lpdf(phi_z | 0, 1) +
      beta_lpdf(R2 | mean_R2 * prec_R2, (1 - mean_R2) * prec_R2) +
      normal_lpdf(sigma | 0, sigma_sd) +
      dirichlet_lpdf(psi | cons);
  }
}

data {

  int<lower = 1> N;             // number of time points
  int<lower = 0> p;             // AR order
  int<lower = 0> P;             // number of regression coefficients
  vector[N] y;                  // observations
  matrix[N, P] X;               // model matrix

  // data for regularization of regression coefficients //
  real<lower = 0> tau0_beta;       // global shrinkage for beta
  real<lower = 0> slab_scl_beta;   // scale for non-zero coefficients in phi
  real<lower = 0> slab_df_beta;    // degrees of freedom for non-zero coefficients

  vector<lower=0>[p] cons;         // concentration vector of the Dirichlet prior

  // data for the R2D2 prior
  real<lower=0> mean_R2;           // mean of the R2 prior
  real<lower=0> prec_R2;           // precision of the R2 prior
  real<lower=0> sigma_sd;          // sd of sigma prior

  // // Data for forecasting
  // int<lower = 0> N_new;
  // matrix[N_new, K] X_new;

}

transformed data {

  // Variance estimate of y
  real<lower = 0> var_y = variance(y);

  // transformations for horseshoe priors
  real slab_scl2_beta = square(slab_scl_beta);
  real half_slab_df_beta = 0.5 * slab_df_beta;

}

parameters {

  simplex[p] psi;                 // decomposition simplex
  real<lower=0, upper=1> R2;      // coefficient of determination
  real<lower=0> sigma;            // observation model sd
  real alpha;                     // intercept
  vector[p] phi_std;

  vector[P] beta_std;             // standardized betas
  vector[P] local_scale_beta;       // unscaled version of local shrinkage priors
  real<lower = 0> c2_std_beta;    // unscaled version of c2
  real<lower = 0> tau_std_beta;   // unscaled version of tau

}

transformed parameters {

  vector[p] phi;
  vector[P] beta;
  vector[N] mu;
  vector[N] eta;
  vector[N] err;

  // scale c2: c2 ~ inv_gamma(half_slab_df, half_slab_df * slab_scl2)
  real c2_beta = slab_scl2_beta * c2_std_beta;

  // tau ~ cauchy(0, tau0)
  real tau_beta = tau0_beta * tau_std_beta;

  // This calculation follows equation 2.8 in Piironen and Vehtari 2017
  vector[P] local_scale_tilde_b =
    sqrt(c2_beta * square(local_scale_beta) ./ (c2_beta + square(tau_beta) * square(local_scale_beta)));

  // scale phi
  phi = sqrt(sigma^2 / var_y * (R2 / (1 - R2)) * psi) .* phi_std;

  // scaling beta
  beta = tau_beta * local_scale_tilde_b .* beta_std;

   // linear predictor
  mu = alpha + X * beta;
  eta[1:p] = y[1:p] - mu[1:p];
  err[1:p] = rep_vector(0, p);
  vector[p] phi_rev = reverse(phi);       // reverse order phi
  for(t in (p + 1):N){
    eta[t] = err[(t - p):(t - 1)]' * phi_rev;
    err[t] = y[t] - mu[t] - eta[t];
  }
}

model {

  // priors //
  beta_std ~ std_normal();
  alpha ~ std_normal();
  tau_std_beta ~ cauchy(0, 1);
  local_scale_beta ~ cauchy(0, 1);
  c2_std_beta ~ inv_gamma(half_slab_df_beta, half_slab_df_beta);
  target += arr2_ncp_lpdf(phi_std | psi, R2, sigma, sigma_sd, mean_R2, prec_R2, cons, var_y);

  // likelihood
  err[(p + 1):N] ~ normal(0, sigma);

}
