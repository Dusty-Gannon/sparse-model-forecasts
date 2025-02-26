

data {

  int<lower = 0> N;            // Length of time series
  int<lower = 0> p;            // AR order
  // int<lower = 0> K;
  vector[N] y;                 // time series
  // matrix[N, K] X;

}


parameters {

  // vector[K] beta;
  vector<lower = 0, upper = 1>[p] r_phi_std;      // partial correlations based on phi
  real<lower = 0> sigma;                          // standard deviation of random innovations

}


transformed parameters{

  vector[N] eta;
  vector[N] err;
  vector[p] phi;

  vector[p] r_phi = 2 * r_phi_std - 1;

  matrix[p, p] P = diag_matrix(r_phi);        // matrix to track recursions

  // Recursions to compute phi and theta
  for(k in 2:p){
    for(i in 1:(k - 1)){
      P[i, k] = P[i, (k - 1)] - r_phi[k] * P[(k - i), (k - 1)];
    }
  }
  phi = P[, p];


  // mu = X * beta;
  eta[1:p] = y[1:p];
  err[1:p] = rep_vector(0, p);
  vector[p] phi_rev = reverse(phi);
  // vector[q] theta_rev = reverse(theta);
  for(t in (p + 1):N){
    eta[t] = err[(t - p):(t - 1)]' * phi_rev;
    err[t] = y[t] - eta[t];
  }

}


model {

  // priors //
  sigma ~ normal(0, 2.5);
  // beta ~ normal(0, 2.5);

  for(i in 1:p){
    r_phi_std[i] ~ beta(0.5 * (i + 1), 0.5 * i + 1);
  }

  // likelihood
  err[(p + 1):N] ~ normal(0, sigma);

}


