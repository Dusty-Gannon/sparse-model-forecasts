/////////////////////////////////////////////////////////////
// This Stan program fits an AR(p) Gaussian time-series model
// with the potential for covariates
/////////////////////////////////////////////////////////////
data{
  int<lower = 1> N;           // number of observations
  int<lower = 0> P;           // total number of covariates
  int<lower = 0> P_0;         // number of un-regularized covariates
  int<lower = 0> p;           // guess for the max order of the autoregressive process
  vector[N] y;                // vector of responses
  matrix[N, P] X;             // model matrix for shrinking effects

  // forecasting objects
  int<lower = 0> N_new;       // number of observations to forecast
  matrix[N_new, P] X_new;     // model matrix for forecasting

}

transformed data{

  row_vector[p] ones = rep_row_vector(1, p);     // row vector of ones

}

parameters{

  vector[P_0] alpha;          // unshrunk coefficients
  vector[P - P_0] beta_std;   // standardized coefficients before shrinkage
  vector[p] phi_std;          // autoregression parameters
  real<lower = 0> sigma;      // sd of the innovations

}
transformed parameters{

  vector[N] mu;               // declare vector of means
  vector[N] err;              // vector of residuals
  vector[P] beta;             // all scaled regression coefficients

  // scale phis
  vector[p] phi = phi_std * 0.25;

  // scale betas
  beta[1:P_0] = alpha;
  beta[(P_0 + 1):P] = beta_std * 0.5;

  mu = X * beta;
  err = y - mu;

}

model{

  // define random innovations
  vector[N] epsilon;

  // define reverse order to make AR coefficient order match tradition
  vector[p] err_rev;

  // priors
  beta_std ~ std_normal();
  phi_std ~ std_normal();

  sigma ~ cauchy(0, 1);

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


generated quantities{

  // post. pred. sampling plus forecasts
  real y_rep[N + N_new];
  vector[N + N_new] err_rep;                 // extend the error term for the forecasting
  vector[p] err_rep_rev = rep_vector(0, p);  // initialize reverse order for consistency

  err_rep[1:N] = err;

  // fitted values
  for(t in 1:N){
    y_rep[t] = mu[t] + err[t] + normal_rng(0, sigma);
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

