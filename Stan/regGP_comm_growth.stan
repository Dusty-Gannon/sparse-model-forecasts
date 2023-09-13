
#include GP_functions.stan

data{

  int<lower = 0> S;             // number of species in each community
  int<lower = 0> N;             // Length of time series minus 1 (number of steps)
  int<lower = 0> K;             // number of sites

  matrix[S * K, N] Ntm1;        // abundunce of each species at time t - 1
  matrix[S * K, N] y;           // log growth from t - 1 to t

}

transformed data{

  matrix[S * K, S * K] Wtm1[N];                       // diagonal weights matrix
  matrix[K, K] I_K = diag_matrix(rep_vector(1, K));   // Identity matrix
  vector[S] ones = rep_vector(1, S);                  // S-vector of ones

  for(t in 1:N){
    Wtm1[t] = diag_matrix(inv(sqrt(Ntm1[,t])));
  }

}

parameters{

  real<lower = 0> sigma_lambda;         // scale of site-level random effects
  real<lower = 0> sigma_e;              // scale of errors
  real a0;                              // log-generic competitive effect

  vector<lower = 0>[S] alpha_std;       // intra-specific effect
  vector[S * (S - 1)] a_tilde_std;      // log-deviation from generic effect
  vector[S] r;                          // log growth
  vector[K] u;                          // random effects for the site

}

transformed parameters{

  vector[S * K] Lambda;                       // declare stacked lambda vector
  matrix[S * K, N] eta;                       // declare matrix of linear predictors

  vector[S] alpha = alpha_std * 0.25;         // scale intra-specific effects

  // scale the deviations from the generic effect
  vector[S * (S - 1)] a_tilde = a_tilde_std * 0.1;

  matrix[S, S] Alpha = diag_matrix(alpha);    // diagonal matrix of intra-specific effects
  matrix[S, S] A0 = fill_off_diag(            // repeating matrix of log-generic effects
    rep_vector(0, S),
    rep_vector(a0, S * (S - 1))
  );

  // matrix of log-non-generic effects with 0s along the diagonal
  matrix[S, S] A_tilde = fill_off_diag(rep_vector(0, S), a_tilde);

  // construct stacked lambda vector
  for(k in 1:K){
    Lambda[(S * (k - 1) + 1):(S * k)] = exp(r + ones * u[k] * sigma_lambda);
  }

  // construct linear predictor
  for(t in 1:N){
    eta[, t] = Lambda + (kronecker_prod(I_K, Alpha) +
      exp(kronecker_prod(I_K, A0) + kronecker_prod(I_K, A_tilde))) * Ntm1[, t];
  }

}

model{

  // priors //
  u ~ std_normal();
  alpha_std ~ std_normal();
  a_tilde_std ~ std_normal();
  r ~ normal(0, 2);
  sigma_lambda ~ normal(0, 2);
  sigma_e ~ normal(0, 2);
  a0 ~ normal(0, 5);

  // likelihood //
  for(t in 1:N){

    y[, t] ~ multi_normal(eta[, t], sigma_e * Wtm1[t]);

  }

}

generated quantities{

  matrix[S, S] A_mat = Alpha + exp(A0 + A_tilde);

}

