

data {

  int<lower = 0> N;            // Length of time series
  int<lower = 0> p;            // AR order
  vector[N] y;                 // data

}

parameters{

  vector<lower = -1, upper = 1>[p] r;       // partial correlations
  real<lower = 0> sigma;                    // scale of random innovations

}

transformed parameters{

  vector[p] phi;                                // AR parameters
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
  real J;           // Jacobian adjustment

  // priors
  phi ~ normal(0, 0.1);
  J = 0;
  for(k in 2:p){
    J += 0.5 * k * log(1 - r[k]) + 0.5 * (k - 1) * log(1 + r[k]);
  }
  target += J;

  for(t in (p + 1):N){
    for(i in 1:p){
      y_rev[i] = y[t - i];
    }
    y[t] ~ normal(y_rev' * phi, sigma);
  }

}
