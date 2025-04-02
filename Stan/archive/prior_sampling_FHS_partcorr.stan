
data {

  int<lower = 0> p;
  real<lower = 0> tau0_phi;        // global shrinkage for phi
  real<lower = 0> slab_scl_phi;    // scale for non-zero coefficients in phi
  real<lower = 0> slab_df_phi;     // degrees of freedom for non-zero coefficients

}

transformed data {

  real slab_scl2_phi = square(slab_scl_phi);
  real half_slab_df_phi = 0.5 * slab_df_phi;

}

parameters {

  vector<lower = -1, upper = 1>[p] r_phi_std;      // partial correlations based on phi

  vector<lower = 0>[p] local_scale_phi;       // non-regularized local scale for phi
  real<lower = 0> c2_std_phi;                     // unscaled version of c2
  real<lower = 0> tau_std_phi;                    // unscaled version of tau

}

transformed parameters {

  vector[p] r_phi;
  vector[p] phi;

  // scale c2: c2 ~ inv_gamma(half_slab_df, half_slab_df * slab_scl2)
  real c2_phi = slab_scl2_phi * c2_std_phi;
  // real c2_beta = slab_scl2_beta * c2_std_beta;

  // tau ~ cauchy(0, tau0)
  real tau_phi = tau0_phi * tau_std_phi;

  // This calculation follows equation 2.8 in Piironen and Vehtari 2017
  vector[p] local_scale_tilde_p =
    sqrt(c2_phi * square(local_scale_phi) ./ (c2_phi + square(tau_phi) * square(local_scale_phi)));

  // scaling the partial correlations
  r_phi = tau_phi * local_scale_tilde_p .* r_phi_std;

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

  r_phi_std ~ std_normal();
  tau_std_phi ~ cauchy(0, 1);
  local_scale_phi ~ cauchy(0, 1);
  c2_std_phi ~ inv_gamma(half_slab_df_phi, half_slab_df_phi);

}


