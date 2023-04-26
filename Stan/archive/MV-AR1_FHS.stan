

data{

  int<lower = 1> N;           // number of observations
  int<lower = 1> S;           // number of species in the community
  int<lower = 1> P;           // number of non-shrinking parameters
  matrix[S, N] y;             // matrix of responses
  matrix[S, P] X_alpha[N];    // model matrix for unshrunk effects
  real<lower = 0> tau0;       // scale for global shrinkage parameter
  real<lower = 0> slab_scl;   // scale for non-zero coefficients
  real<lower = 0> slab_df;    // degrees of freedom for non-zero coefficients

}


transformed data{

  // transformations for horseshoe priors
  real slab_scl2 = square(slab_scl);
  real half_slab_df = 0.5 * slab_df;

  // vector of ones to combine alphas
  vector[S] ones_alpha = rep_vector(1, S);

}


parameters{

  matrix[P, S] alpha_std;                  // unshrunk coefficients (self-limitation and intercept)
  matrix[S, S] beta_std;                   // standardized coefficients before shrinkage
  // vector<lower = 0, upper = 1>[S] phi;     // AR params (loosely interpreted as growth rates)
  vector<lower = 0>[S] sigma;              // sd of the innovations for each species

  // parameters for shrinkage priors
  matrix<lower = 0>[S, S] local_scale;    // non-regularized local scale
  real<lower = 0> c2_std;                 // unscaled version of c2
  real<lower = 0> tau_std;                // unscaled version of tau

}


transformed parameters{

  // scale c2: c2 ~ inv_gamma(half_slab_df, half_slab_df * slab_scl2)
  real c2 = slab_scl2 * c2_std;

  // tau ~ cauchy(0, tau0)
  real tau = tau0 * tau_std;

  // This calculation follows equation 2.8 in Piironen and Vehtari 2017
  matrix[S, S] local_scale_tilde =
    sqrt(c2 * square(local_scale) ./ (c2 + square(tau) * square(local_scale)));

  // scale betas
  matrix[S, S] beta = tau * local_scale_tilde .* beta_std;

  // scale alpha
  matrix[P, S] alpha = alpha_std * 5;

}


model{

  // declare some extra variables for convenience
  matrix[S, N] epsilon;

  // priors
  to_vector(beta_std) ~ std_normal();
  to_vector(alpha_std) ~ std_normal();

  tau_std ~ cauchy(0, 1);
  to_vector(local_scale) ~ cauchy(0, 1);
  c2_std ~ inv_gamma(half_slab_df, half_slab_df);

  sigma ~ cauchy(0, 1);

  // likelihood
  epsilon[, 1] = rep_vector(0, S);

  for(t in 2:N){

    epsilon[, t] = y[, t] - diagonal(X_alpha[t] * alpha) - beta * y[, (t - 1)];
    epsilon[, t] ~ multi_normal(rep_vector(0, S), diag_matrix(sigma));

  }

}


generated quantities{

  // post. pred draws
  vector[S] y_rep[N];

  y_rep[1] = y[, 1];
  for(t in 2:N){
    y_rep[t] = multi_normal_rng(diagonal(X_alpha[t] * alpha) + beta * y[, (t - 1)], diag_matrix(sigma));
  }

}



