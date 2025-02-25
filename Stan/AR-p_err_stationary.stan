

data {

  int<lower = 0> N;            // Length of time series
  int<lower = 0> p;            // AR order
  int<lower = 0> q;            // MA order
  vector[N] y;                 // time series

}


transformed data {

  int r = max(p, q);
  int q_plus = q + 1;

}


parameters {

  real mu;
  vector<lower = 0, upper = 1>[p] r_phi_std;      // partial correlations based on phi
  vector<lower = 0, upper = 1>[q] r_theta_std;    // partial correlations based on theta
  // vector<lower = -1, upper = 1>[p] phi;
  // vector<lower = -1, upper = 1>[q] theta;
  real<lower = 0> sigma;                          // standard deviation of random innovations

}


transformed parameters{

  vector[p] phi;
  vector[q] theta;
  vector[p] r_phi = 2 * r_phi_std - 1;
  vector[q] r_theta = 2 * r_theta_std - 1;

  matrix[p, p] P = diag_matrix(r_phi);        // matrix to track recursions
  matrix[q, q] Q = diag_matrix(r_theta);      // matrix to track recursions

  // Recursions to compute phi and theta
  for(k in 2:p){
    for(i in 1:(k - 1)){
      P[i, k] = P[i, (k - 1)] - r_phi[k] * P[(k - i), (k - 1)];
    }
  }
  phi = P[, p];

  for(k in 2:q){
    for(i in 1:(k - 1)){
      Q[i, k] = Q[i, (k - 1)] - r_theta[k] * Q[(k - i), (k - 1)];
    }
  }
  theta = Q[, q];

}


model {

  vector[N] err;      // error term
  vector[N] eta;           // mean-centered model

  // priors //
  sigma ~ normal(0, 2.5);
  mu ~ normal(0, 2.5);

  for(i in 1:p){
    r_phi_std[i] ~ beta(0.5 * (i + 1), 0.5 * i + 1);
  }

  for(i in 1:q){
    r_theta_std[i] ~ beta(0.5 * (i + 1), 0.5 * i + 1);
  }

  // likelihood //
  eta = y - rep_vector(mu, N);
  err[1:r] = rep_vector(0, r);
  vector[p] phi_rev = reverse(phi);
  vector[q] theta_rev = reverse(theta);
  for(t in (r + 1):N){
    err[t] = eta[t] - (eta[(t-p):(t-1)]' * phi_rev + err[(t - q):(t - 1)]' * theta_rev);
  }

  // likelihood
  err[(r + 1):N] ~ normal(0, sigma);

}


