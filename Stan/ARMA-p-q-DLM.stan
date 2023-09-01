////////////////////////////////////////////////////////////
// This Stan program fits stationary ARMA models using a DLM
// parameterization and the Kalman Filter
////////////////////////////////////////////////////////////

#include DLM_functions.stan

data {

  int<lower = 0> N;            // Length of time series
  int<lower = 0> p;            // AR order
  int<lower = 0> q;            // MA order
  real m_0;                    // initial state mean
  real<lower = 0> C_0;         // initial state variance
  vector[N] y;                 // data

}

transformed data{

  // calculate r
  int r = max(p, q + 1);

  // compute F
  row_vector[r] F = make_F(r);

  // make V with no observation error
  real<lower = 0> V = 0;

  // turn m_0 into a matrix
  vector[r] m0 = rep_vector(m_0, r);

  // turn C_0 into covariance matrix
  matrix[r, r] C0 = diag_matrix(rep_vector(C_0, r));

}


parameters {

  vector<lower = 0, upper = 1>[p] r_phi_std;      // partial correlations based on phi
  vector<lower = 0, upper = 1>[q] r_theta_std;    // partial correlations based on theta
  real<lower = 0> sigma;                          // standard deviation of random innovations

}

transformed parameters{

  vector[p] phi;
  vector[q] theta;
  vector[p] r_phi = 2 * r_phi_std - 1;
  vector[q] r_theta = 2 * r_theta_std - 1;

  matrix[p, p] P = diag_matrix(r_phi);        // matrix to track recursions
  matrix[q, q] Q = diag_matrix(r_theta);      // matrix to track recursions

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

  // model //
  // transition matrix
  matrix[r, r] G = make_G(augment_phi(r, phi));

  // R vector
  vector[r] R = make_R(augment_theta(r, theta));

  // Covariance matrix
  matrix[r, r] W = R * R' * square(sigma);

  // // state vectors
  // vector[r] x[N];
  //
  // // state covariance
  // matrix[r, r] C[N];

  // Kalman filter objects //
  vector[r] x_t = m0;        // filtered state prediction
  matrix[r, r] C_t = C0;     // filtered covariance
  vector[r] x_pred;          // step ahead state predictions
  matrix[r, r] C_pred;       // step ahead state covariance
  vector[N] y_pred;          // step ahead observation prediction
  vector[N] S_pred;          // step ahead observation covariance
  vector[r] K_t;             // Kalman gain

  for(t in 1:N){

    // one-step ahead state prediction and variance
    x_pred = G * x_t;
    C_pred = G * C_t * G' + W;

    // one step ahead observation prediction and variance
    y_pred[t] = F * x_pred;
    S_pred[t] = F * C_pred * F' + V;

    // Update with Kalman filter
    K_t = (C_pred * F') / S_pred[t];
    x_t = x_pred + K_t * (y[t] - y_pred[t]);
    C_t = C_pred - K_t * F * C_pred;

    // // store for sampling statement in generated quantities block
    // x[t] = x_t;
    // C[t] = C_t;

  }

  // priors //
  sigma ~ normal(0, 2.5);

  for(i in 1:p){
    r_phi_std[i] ~ beta(0.5 * (i + 1), 0.5 * i + 1);
  }

  for(i in 1:q){
    r_theta_std[i] ~ beta(0.5 * (i + 1), 0.5 * i + 1);
  }

  // likelihood //
  for(t in 1:N){
    y[t] ~ normal(y_pred[t], sqrt(S_pred[t]));
  }

}

