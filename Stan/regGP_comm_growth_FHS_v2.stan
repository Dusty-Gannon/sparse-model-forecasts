
#include GP_functions.stan

data{

  int<lower = 0> S;                          // number of species in each community
  int<lower = 0> N;                          // Length of time series minus 1 (number of steps)
  int<lower = 0> P;                          // additional weakly regularized effects

  matrix[S, N] Ntm1;                         // abundunce of each species at time t - 1, stacked
  vector[S * N] y;                           // log growth from t - 1 to t, stacked for each species

  matrix[S * N, P] X_beta;                   // matrix of non-shrinking effects
  int<lower = 0, upper = S> sp_id[N * S];    // species index
  real<lower = 0> tau0;                      // scale for global shrinkage parameter
  real<lower = 0> slab_scl;                  // scale for non-zero coefficients
  real<lower = 0> slab_df;                   // degrees of freedom for non-zero coefficients

}

transformed data{

  vector[S * N] Wtm1;                       // demographic stochasticity weights
  matrix[S, S * N] Ntm1_rep;                // repeating columns of Ntm1
  // matrix[K, K] I_K = diag_matrix(rep_vector(1, K));   // Identity matrix
  // vector[S] ones = rep_vector(1, S);                  // S-vector of ones

  // transformations for horseshoe priors
  real slab_scl2 = square(slab_scl);
  real half_slab_df = 0.5 * slab_df;

  for(t in 1:N){
    Wtm1[(S * (t - 1) + 1):(S * t)] = inv(sqrt(Ntm1[,t]));
    for(s in 1:S){
      Ntm1_rep[, S * (t - 1) + s] = Ntm1[, t];
    }
  }

}

parameters{

  vector<lower = 0>[S] sigma_s;         // scale of errors
  real a0_std;                          // standardized log-generic competitive effect

  vector[S] alpha_std;                  // log intra-specific effect
  vector[S * (S - 1)] a_tilde_std;      // log-deviation from generic effect
  vector<lower = 0>[S] r;               // log growth for each species
  vector[P] beta_std;                   // additional covariate effects

  // parameters for shrinkage priors
  vector<lower = 0>[S * (S - 1)] local_scale;    // non-regularized local scale
  real<lower = 0> c2_std;                        // unscaled version of c2
  real<lower = 0> tau_std;                       // unscaled version of tau

}

transformed parameters{

  vector[S * N] eta;                          // declare vector of linear predictors
  matrix[S, S] A;                             // compiled interaction matrix

  vector[S] alpha = alpha_std;         // scale log intra-specific effects
  vector[P] beta = beta_std * 2;       // scale betas
  real a0 = a0_std;                    // scale log generic effect

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

  // compile interaction matrix
  A = exp(Alpha + A0 + A_tilde);

  // construct linear predictor
  for(i in 1:(S * N)){
    eta[i] = r[sp_id[i]] + X_beta[i, ] * beta -
      A[sp_id[i], ] * Ntm1_rep[, i];
  }

}

model{

  // priors //
  alpha_std ~ std_normal();
  a_tilde_std ~ std_normal();
  beta_std ~ std_normal();
  r ~ normal(0, 2);
  sigma_s ~ normal(0, 1);
  a0_std ~ std_normal();

  // priors for horseshoe hierarchy
  tau_std ~ cauchy(0, 1);
  local_scale ~ cauchy(0, 1);
  c2_std ~ inv_gamma(half_slab_df, half_slab_df);

  // likelihood //
  for(i in 1:N){

    y[i] ~ normal(eta[i], Wtm1[i] * sigma_s[sp_id[i]]);

  }

}


