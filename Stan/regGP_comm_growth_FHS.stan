
#include GP_functions.stan

data{

  int<lower = 0> S;             // number of species in each community
  int<lower = 0> N;             // Length of time series minus 1 (number of steps)
  int<lower = 0> K;             // number of sites

  matrix[S * K, N] Ntm1;        // abundunce of each species at time t - 1, stacked over sites
  matrix[S * K, N] y;           // log growth from t - 1 to t

  matrix[S * K, N] X_beta;      // matrix of dummy variables to include disturbance effects
  real<lower = 0> tau0;         // scale for global shrinkage parameter
  real<lower = 0> slab_scl;     // scale for non-zero coefficients
  real<lower = 0> slab_df;      // degrees of freedom for non-zero coefficients

}

transformed data{

  matrix[S * K, S * K] Wtm1[N];                       // diagonal weights matrix
  matrix[K, K] I_K = diag_matrix(rep_vector(1, K));   // Identity matrix
  vector[S] ones = rep_vector(1, S);                  // S-vector of ones

  // transformations for horseshoe priors
  real slab_scl2 = square(slab_scl);
  real half_slab_df = 0.5 * slab_df;

  for(t in 1:N){
    Wtm1[t] = diag_matrix(inv(sqrt(Ntm1[,t])));
  }

}

parameters{

  real<lower = 0> sigma_lambda;         // scale of site-level random effects
  vector<lower = 0>[S] sigma_s;         // scale of errors
  real a0_std;                          // standardized log-generic competitive effect

  vector[S] alpha_std;                  // log intra-specific effect
  vector[S * (S - 1)] a_tilde_std;      // log-deviation from generic effect
  vector<lower = 0>[S] r;               // log growth for each species
  real beta;                            // effect of disturbance on growth
  vector[K] u;                          // random effects for the site

  // parameters for shrinkage priors
  vector<lower = 0>[S * (S - 1)] local_scale;    // non-regularized local scale
  real<lower = 0> c2_std;                        // unscaled version of c2
  real<lower = 0> tau_std;                       // unscaled version of tau

}

transformed parameters{

  vector[S * K] Lambda;                       // declare stacked lambda vector
  matrix[S * K, N] eta;                       // declare matrix of linear predictors
  vector[S * K] sigma_sk;                     // declare repeated vector of sigmas

  vector[S] alpha = alpha_std * 0.1;         // scale log intra-specific effects
  real a0 = a0_std * 0.1;                    // scale log generic effect

  // scale c2: c2 ~ inv_gamma(half_slab_df, half_slab_df * slab_scl2)
  real c2 = slab_scl2 * c2_std;

  // tau ~ cauchy(0, tau0)
  real tau = tau0 * tau_std;

  // This calculation follows equation 2.8 in Piironen and Vehtari 2017
  vector[S * (S - 1)] local_scale_tilde =
    sqrt(c2 * square(local_scale) ./ (c2 + square(tau) * square(local_scale)));

  // scale the deviations from the generic effect
  vector[S * (S - 1)] a_tilde = tau * local_scale_tilde .* a_tilde_std;

  matrix[S, S] Alpha = diag_matrix(alpha);         // diagonal matrix of intra-specific effects
  matrix[S, S] A0 = fill_off_diag(                 // repeating matrix of log-generic effects
    rep_vector(0, S),
    rep_vector(a0, S * (S - 1))
  );

  // matrix of log-non-generic effects with 0s along the diagonal
  matrix[S, S] A_tilde = fill_off_diag(rep_vector(0, S), a_tilde);

  // construct stacked lambda and sigma vectors
  for(k in 1:K){
    Lambda[(S * (k - 1) + 1):(S * k)] = exp(r + ones * u[k] * sigma_lambda);
    sigma_sk[(S * (k - 1) + 1):(S * k)] = sigma_s;
  }

  // construct linear predictor
  for(t in 1:N){
    eta[, t] = Lambda + X_beta[, t] * beta +
    exp(
      kronecker_prod(I_K, Alpha) +
      kronecker_prod(I_K, A0) +
      kronecker_prod(I_K, A_tilde)
      ) * Ntm1[, t];
  }

}

model{

  // priors //
  u ~ std_normal();
  alpha_std ~ std_normal();
  a_tilde_std ~ std_normal();
  r ~ normal(0, 2);
  sigma_lambda ~ normal(0, 2);
  sigma_s ~ normal(0, 1);
  a0_std ~ std_normal();

  // priors for horseshoe hierarchy
  tau_std ~ cauchy(0, 1);
  local_scale ~ cauchy(0, 1);
  c2_std ~ inv_gamma(half_slab_df, half_slab_df);

  // likelihood //
  for(t in 1:N){

    y[, t] ~ multi_normal(eta[, t], Wtm1[t] * diag_matrix(sigma_sk));

  }

}

generated quantities{

  matrix[S, S] A_mat = Alpha + exp(A0 + A_tilde);

}

