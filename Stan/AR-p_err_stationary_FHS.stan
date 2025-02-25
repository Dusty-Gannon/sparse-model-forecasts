

data {

  int<lower = 0> N;            // Length of time series
  int<lower = 0> P;            // number of covariates
  int<lower = 0> p;            // AR order
  int<lower = 0> q;            // MA order
  vector[N] y;                 // time series

  // prior inputs for regularization
  real<lower = 0> tau0_phi;        // global shrinkage for phi
  real<lower = 0> slab_scl_phi;    // scale for non-zero coefficients in phi
  real<lower = 0> slab_df_phi;     // degrees of freedom for non-zero coefficients
  // real<lower = 0> tau0_beta;       // global shrinkage for beta
  // real<lower = 0> slab_scl_beta;   // scale for non-zero coefficients in phi
  // real<lower = 0> slab_df_beta;    // degrees of freedom for non-zero coefficients

  // Data for forecasting
  int<lower = 0> N_new;
}


transformed data {

  int r = max(p, q);
  int q_plus = q + 1;

  // transformations for horseshoe priors
  real slab_scl2_phi = square(slab_scl_phi);
  real half_slab_df_phi = 0.5 * slab_df_phi;
  // real slab_scl2_beta = square(slab_scl_beta);
  // real half_slab_df_beta = 0.5 * slab_df_beta;

}


parameters {

  vector[P] beta;
  vector<lower = -1, upper = 1>[p] r_phi;      // partial correlations based on phi
  vector<lower = -1, upper = 1>[q] r_theta;    // partial correlations based on theta
  real<lower = 0> sigma;                       // standard deviation of random innovations

  vector<lower = 0>[p + q] local_scale_phi;       // non-regularized local scale for phi
  real<lower = 0> c2_std_phi;                     // unscaled version of c2
  real<lower = 0> tau_std_phi;                    // unscaled version of tau

  // real<lower = 0> c2_std_beta;                    // unscaled version of c2
  // real<lower = 0> tau_std_beta;                   // unscaled version of tau

}


transformed parameters{

  vector[p] phi_std;
  vector[p] phi;
  vector[q] theta_std;
  vector[q] theta;

  // tracking the linear predictor
  vector[N] mu;
  vector[N] eta;
  vector[N] err;


  matrix[p, p] P = diag_matrix(r_phi);        // matrix to track recursions
  matrix[q, q] Q = diag_matrix(r_theta);      // matrix to track recursions

  // scale c2: c2 ~ inv_gamma(half_slab_df, half_slab_df * slab_scl2)
  real c2_phi = slab_scl2_phi * c2_std_phi;
  // real c2_beta = slab_scl2_beta * c2_std_beta;

  // tau ~ cauchy(0, tau0)
  real tau_phi = tau0_phi * tau_std_phi;
  // real tau_beta = tau0_beta * tau_std_beta;

  // This calculation follows equation 2.8 in Piironen and Vehtari 2017
  vector[p + q] local_scale_tilde_p =
    sqrt(c2_phi * square(local_scale_phi) ./ (c2_phi + square(tau_phi) * square(local_scale_phi)));
  // vector[P - P_0] local_scale_tilde_b =
  //   sqrt(c2_beta * square(local_scale_beta) ./ (c2_beta + square(tau_beta) * square(local_scale_beta)));

  // Recursions to compute phi and theta
  for(k in 2:p){
    for(i in 1:(k - 1)){
      P[i, k] = P[i, (k - 1)] - r_phi[k] * P[(k - i), (k - 1)];
    }
  }
  phi_std = P[, p];

  for(k in 2:q){
    for(i in 1:(k - 1)){
      Q[i, k] = Q[i, (k - 1)] - r_theta[k] * Q[(k - i), (k - 1)];
    }
  }
  theta_std = Q[, q];

  // scaling phi and theta
  phi = tau_phi * local_scale_tilde_p[1:p] .* phi_std;

  theta = tau_phi * local_scale_tilde_p[(p + 1):(p + q)] .* theta_std;

  // linear predictor
  mu = X * beta;
  eta[1:r] = rep_vector(0, r);
  err[1:r] = (y[1:r] - mu[1:r]) - eta[1:r];
  vector[p] phi_rev = reverse(phi);       // reverse order phi
  vector[q] theta_rev = reverse(theta);   // reverse order theta
  for(t in (r + 1):N){
    eta[t] = (eta[(t-p):(t-1)]' * phi_rev + err[(t - q):(t - 1)]' * theta_rev);
    err[t] = y[t] - mu[t] - eta[t];
  }
}


model {

  real J_phi = 0;          // Jacobian adjustment
  real J_theta = 0;        // Jacobian adjustment

  // priors //
  sigma ~ normal(0, 2.5);
  beta ~ normal(0, 2.5);

  // beta_std ~ std_normal();
  // alpha_std ~ std_normal();
  phi_std ~ std_normal();

  tau_std_phi ~ cauchy(0, 1);
  // tau_std_beta ~ cauchy(0, 1);
  // local_scale_beta ~ cauchy(0, 1);
  local_scale_phi ~ cauchy(0, 1);
  c2_std_phi ~ inv_gamma(half_slab_df_phi, half_slab_df_phi);
  // c2_std_beta ~ inv_gamma(half_slab_df_beta, half_slab_df_beta);

  // Jacobian adjustments based on Jones (1987) JRSS Series C
  for(k in 2:p){
    J_phi += 0.5 * k * log(1 - r_phi[k]) + 0.5 * (k - 1) * log(1 + r_phi[k]);
  }
  for(k in 2:q){
    J_theta += 0.5 * k * log(1 - r_theta[k]) + 0.5 * (k - 1) * log(1 + r_theta[k]);
  }

  target += J_phi + J_theta;

  // likelihood
  err[(r + 1):N] ~ normal(0, sigma);

}

generated quantities{

  // post. pred. sampling plus forecasts
  vector[N + N_new] y_rep;
  vector[N + N_new] err_rep;
  vector[N + N_new] eta_rep;
  vector[N + N_new] mu_rep;

  // fitted values
  y_rep[1:N] = mu + eta;
  eta_rep[1:N] = eta;
  err_rep[1:N] = err;
  err_rep[(N + 1):(N + N_new)] = rep_vector(0, N_new);

  // forecast
  if(N_new > 0){
    for(t in (N + 1):(N + N_new)){
      eta_rep[t] = eta_rep[(t-p):(t-1)]' * phi_rev;
      y_rep[t] = mu_rep[t] + eta_rep[t] + err_rep[(t-q):(t-1)]' * theta_rev +
        normal_rng(0, sigma);
    }
  }

}


