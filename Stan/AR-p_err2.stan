/////////////////////////////////////////////////////////////
// This Stan program fits an AR(p) Gaussian time-series model
// with the potential for covariates
/////////////////////////////////////////////////////////////
data{
  int<lower = 1> N;           // number of observations
  int<lower = 0> P;           // number of covariates
  int<lower = 0> p;           // guess for the max order of the autoregressive process
  vector[N] y;                // vector of responses
  matrix[N, P] X;             // model matrix for shrinking effects
}
transformed data{
  row_vector[p] ones = rep_row_vector(1, p);     // row vector of ones
}
parameters{
  vector[P] beta_std;         // standardized coefficients before shrinkage
  vector[p] phi_std;          // autoregression parameters
  real<lower = 0> sigma;      // sd of the innovations
}
transformed parameters{
  vector[N] mu;               // declare vector of means

  vector[N] err;              // vector of residuals

  vector[N - p] epsilon;      // random innovations
  // scale betas
  vector[P] beta = beta_std * 2.5;
  // scale phis
  vector[p] phi = phi_std * 0.25;

  mu = X * beta;
  err = y - mu;

  for(t in (p + 1):N){
    epsilon[t-p] = err[t] - err[(t-p):(t-1)]' * phi;
  }
}

model{
  // priors
  beta_std ~ std_normal();
  phi_std ~ std_normal();

  sigma ~ cauchy(0, 1);

  epsilon ~ normal(0, sigma);
}


generated quantities{

  // post. pred. sampling
  real y_rep[N - p] = normal_rng(err[(1 + p):N] + mu[(1 + p):N], sigma);

  // residuals
  vector[N - p] resid = y[(p + 1):N] - mu[(p + 1):N];

}

