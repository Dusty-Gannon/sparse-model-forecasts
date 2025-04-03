
data {

  int<lower = 0> p;            // AR order

}

parameters {

  vector<lower = 0, upper = 1>[p] r_phi_std;      // partial correlations based on phi

}

transformed parameters {

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

}

model {

  for(i in 1:p){
    r_phi_std[i] ~ beta(0.5 * (i + 1), 0.5 * i + 1);
  }

}

