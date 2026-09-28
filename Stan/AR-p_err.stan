/////////////////////////////////////////////////////////////
// This Stan program fits an AR(p) Gaussian time-series model
// with the potential for covariates. The AR(p) coefficients,
// phi, are reparameterized as partial correlations, r_phi, and
// mapped back to phi via the Levinson-Durbin recursion, which
// guarantees a stationary process. Both r_phi and the
// regularized subset of the regression coefficients, beta, get
// sparsity-inducing horseshoe priors; the global shrinkage of the beta
// horseshoe is scaled by the implied innovation variance of the
// AR(p) process.
/////////////////////////////////////////////////////////////

#include ar_var_approx.stan

data {

  int<lower = 1> N;           // number of observations
  int<lower = 0> P;           // total number of covariates
  int<lower = 0, upper = P> P_0;  // number of un-regularized covariates (first P_0 columns of X)
  int<lower = 1> p;           // guess for the max order of the autoregressive process
  vector[N] y;                // vector of responses
  matrix[N, P] X;              // model matrix (un-regularized columns first)

  // prior inputs for regularization of phi (via partial correlations)
  real<lower = 0> tau0_phi;
  real<lower = 0> slab_scl_phi;
  real<lower = 0> slab_df_phi;

  // prior inputs for regularization of beta
  real<lower = 0> tau0_beta;
  real<lower = 0> slab_scl_beta;
  real<lower = 0> slab_df_beta;

  // forecasting objects
  int<lower = 0> N_new;       // number of observations to forecast
  matrix[N_new, P] X_new;     // model matrix for forecasting

}

transformed data {

  real slab_scl2_phi    = square(slab_scl_phi);
  real half_slab_df_phi = 0.5 * slab_df_phi;
  real slab_scl2_beta    = square(slab_scl_beta);
  real half_slab_df_beta = 0.5 * slab_df_beta;

}

parameters {

  vector[P_0] alpha_std;                  // unshrunk coefficients (raw scale)
  vector[P - P_0] beta_std;               // standardized coefficients before shrinkage
  real<lower = 0> sigma;                  // innovation SD

  vector[p] r_phi_std;                    // unconstrained auxiliary (non-centered) for partial correlations
  vector<lower = 0>[p] local_scale_phi;   // local shrinkage scales for phi
  real<lower = 0> c2_std_phi;             // unscaled slab variance for phi
  real<lower = 0> tau_std_phi;            // unscaled global shrinkage scale for phi

  vector<lower = 0>[P - P_0] local_scale_beta;  // local shrinkage scales for beta
  real<lower = 0> c2_std_beta;                  // unscaled slab variance for beta
  real<lower = 0> tau_std_beta;                 // unscaled global shrinkage scale for beta

}

transformed parameters {

  // ---- phi: regularized horseshoe on partial correlations ----

  // tau ~ cauchy(0, tau0_phi)
  real tau_phi = tau0_phi * tau_std_phi;

  // c2 ~ inv_gamma(half_slab_df, half_slab_df * slab_scl2)
  real c2_phi  = slab_scl2_phi * c2_std_phi;

  // Piironen & Vehtari (2017) eq. 2.8
  vector[p] local_scale_tilde_phi =
    sqrt(c2_phi * square(local_scale_phi) ./
         (c2_phi + square(tau_phi) * square(local_scale_phi)));

  // tanh maps the horseshoe scale to (-1, 1) for partial correlations
  vector[p] r_phi = tanh(tau_phi * local_scale_tilde_phi .* r_phi_std);

  // Levinson-Durbin recursion: partial correlations -> AR coefficients
  vector[p] phi;
  {
    matrix[p, p] Pmat = diag_matrix(r_phi);
    for(k in 2:p){
      for(i in 1:(k - 1)){
        Pmat[i, k] = Pmat[i, (k - 1)] - r_phi[k] * Pmat[(k - i), (k - 1)];
      }
    }
    phi = Pmat[, p];
  }

  // ---- beta: regularized horseshoe, scaled by the implied AR(p) innovation variance ----

  real c2_beta  = slab_scl2_beta * c2_std_beta;

  // same process as for phi, but the global scale accounts for the
  // marginal variance implied by the (stationary) AR(p) process
  real tau_beta = tau0_beta * tau_std_beta * sqrt(ar_var(phi, square(sigma), 20));

  vector[P - P_0] local_scale_tilde_beta =
    sqrt(c2_beta * square(local_scale_beta) ./
         (c2_beta + square(tau_beta) * square(local_scale_beta)));

  vector[P] beta;
  beta[1:P_0] = alpha_std * 2.5;
  beta[(P_0 + 1):P] = tau_beta * local_scale_tilde_beta .* beta_std;

  vector[N] mu  = X * beta;
  vector[N] err = y - mu;

}

model {

  // define random innovations
  vector[N] epsilon;

  // define reverse order to make AR coefficient order match tradition
  vector[p] err_rev;

  // priors: phi
  r_phi_std ~ std_normal();
  local_scale_phi ~ cauchy(0, 1);
  tau_std_phi ~ cauchy(0, 1);
  c2_std_phi ~ inv_gamma(half_slab_df_phi, half_slab_df_phi);

  // priors: beta
  alpha_std ~ std_normal();
  beta_std  ~ std_normal();
  local_scale_beta ~ cauchy(0, 1);
  tau_std_beta ~ cauchy(0, 1);
  c2_std_beta ~ inv_gamma(half_slab_df_beta, half_slab_df_beta);

  sigma ~ normal(0, 2);

  // complete the AR process
  epsilon[1:p] = rep_vector(0, p);
  for(t in (p + 1):N){
    for(i in 1:p){
      err_rev[i] = err[t - i];
    }
    epsilon[t] = err[t] - err_rev' * phi;
  }

  // likelihood
  epsilon[(p + 1):N] ~ normal(0, sigma);

}

generated quantities {

  // post. pred. sampling plus forecasts
  real y_rep[N + N_new];
  vector[N + N_new] err_rep;                 // extend the error term for forecasting
  vector[p] err_rep_rev = rep_vector(0, p);  // initialize reverse order for consistency

  err_rep[1:N] = err;

  // fitted values
  for(t in 1:p){
    y_rep[t] = mu[t];
  }
  for(t in (p + 1):N){
    for(i in 1:p){
      err_rep_rev[i] = err_rep[t - i];
    }
    y_rep[t] = mu[t] + err_rep_rev' * phi;
  }

  // forecast
  if(N_new > 0){
    for(t in (N + 1):(N + N_new)){
      for(i in 1:p){
        err_rep_rev[i] = err_rep[t - i];
      }
      err_rep[t] = err_rep_rev' * phi + normal_rng(0, sigma);
      y_rep[t] = X_new[t - N, ] * beta + err_rep[t];
    }
  }

}
