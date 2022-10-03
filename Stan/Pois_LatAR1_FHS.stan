//
// This Stan program defines an auto-regressive Poisson model with a latent AR(1)
//   process. The regression coefficients are shrunk using regularized horsheshoe priors
//

functions{

  // create AR(1) correlation matrix
  matrix chol_ar1_corrmat(int N, real rho){
    // declare variables
    matrix[N, N] Q;
    matrix[N, N] L;

    // calculations
    for(i in 1:N){
      for(j in 1:N){
        Q[i, j] = pow(rho, fabs(i - j));
      }
    }
    L = cholesky_decompose(Q);
    return(L);
  }

}

data{

  int<lower = 1> N;           // number of observations
  int<lower = 1> P;           // number of non-shrinking parameters
  int<lower = 0> S;           // number of competing species
  int<lower = 0> y[N];        // vector of responses
  matrix[N, P] X_alpha;       // model matrix for unshrunk effects
  matrix[N, S] X_beta;        // model matrix for competing species
  real<lower = 0> tau0;       // scale for global shrinkage parameter
  real<lower = 0> slab_scl;   // scale for non-zero coefficients
  real<lower = 0> slab_df;    // degrees of freedom for non-zero coefficients

}

transformed data{

  // transformations for horseshoe priors
  real slab_scl2 = square(slab_scl);
  real half_slab_df = 0.5 * slab_df;

}

parameters{

  vector[P] alpha_std;                 // unshrunk coefficients (self-limitation and intercept)
  vector[S] beta_std;                  // standardized coefficients for competing species effects
  real<lower = 0, upper = 1> phi;      // autocorrelation parameter
  real<lower = 0> sigma_e;             // sd of random innovations
  vector[N] gamma_std;                 // standardized latent AR(1) variables

  // parameters for shrinkage priors
  vector<lower = 0>[S] local_scale;    // non-regularized local scale
  real<lower = 0> c2_std;              // unscaled version of c2
  real<lower = 0> tau_std;             // unscaled version of tau

}

transformed parameters{

  vector[N] eta;               // declare vector of means
  vector[N] gamma;             // declare vector of latent AR(1) variables

  // scale c2: c2 ~ inv_gamma(half_slab_df, half_slab_df * slab_scl2)
  real c2 = slab_scl2 * c2_std;

  // tau ~ cauchy(0, tau0)
  real tau = tau0 * tau_std;

  // This calculation follows equation 2.8 in Piironen and Vehtari 2017
  vector[S] local_scale_tilde =
    sqrt(c2 * square(local_scale) ./ (c2 + square(tau) * square(local_scale)));

  // scale betas
  vector[S] beta = tau * local_scale_tilde .* beta_std;

  // scale alpha
  vector[P] alpha = alpha_std;

  // scale gamma
  gamma = (sigma_e)/sqrt(1 - square(phi)) * chol_ar1_corrmat(N, phi) * gamma_std;

  // construct linear predictors
  eta = X_alpha * alpha + X_beta * beta + gamma;
}

model{

  // priors
  alpha_std ~ std_normal();
  beta_std ~ std_normal();
  gamma_std ~ std_normal();
  sigma_e ~ normal(0, 3);

  tau_std ~ cauchy(0, 1);
  local_scale ~ cauchy(0, 1);
  c2_std ~ inv_gamma(half_slab_df, half_slab_df);

  // likelihood
  y ~ poisson_log(eta);

}

