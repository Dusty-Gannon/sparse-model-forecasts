

data {

  int<lower = 0> N;            // Length of time series
  int<lower = 0> p;            // AR order
  vector[N] y;                 // data

}

parameters{

  vector<lower = 0, upper = 1>[p] r_std;    // partial correlations
  real<lower = 0> sigma;                    // scale of random innovations

}

transformed parameters{

  vector[p] phi;                                // AR parameters
  vector[p] r = 2 * r_std - rep_vector(1, p);   // scale and shift to be between -1 and 1
  matrix[p, p] P = diag_matrix(r);              // matrix to store recursions

  for(k in 2:p){
    for(i in 1:(k - 1)){
      P[i, k] = P[i, (k - 1)] - r[k] * P[(k - i), (k - 1)];
    }
  }

  phi = P[, p];

}

model{

  // additional objects
  vector[p] y_rev;

  // priors
  for(i in 1:p){
    r_std[i] ~ beta(0.5 * (i + 1), 0.5 * i + 1);
  }

  for(t in (p + 1):N){
    for(i in 1:p){
      y_rev[i] = y[t - i];
    }
    y[t] ~ normal(y_rev' * phi, sigma);
  }

}
