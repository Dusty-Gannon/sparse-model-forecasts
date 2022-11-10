


data{

  int<lower = 1> N;           // number of observations
  int<lower = 0> P;           // number of regression coefficients
  int<lower = 0> y[N];        // vector of responses
  real<lower = 0> y_star[N];  // vector of threshold transformed responses
  matrix[N, P] X;             // model matrix

}


parameters{

  vector[P] beta_std;                  // standardized coefficients before scaling
  real<lower = 0, upper = 1> phi;      // autoregression parameter
  real<lower = 0, upper = 1> theta;    // MA parameter

}


transformed parameters{

  vector[N] mu;                // declare vector of means
  vector[N] z;                 // declare vector of ARMA variables

  // scale betas
  vector[P] beta = beta_std * 5;

  // construct linear predictors
  mu[1] = y_star[1];
  z[1] = 0;
  for(i in 2:N){

    z[i] = phi * (log(y_star[i-1]) - X[i-1, ] * beta) + theta * log(y_star[i-1] / mu[i-1]);
    mu[i] = exp(X[i, ] * beta + z[i]);

  }

}


model{

  // priors
  beta_std ~ std_normal();

  // likelihood
  for(i in 1:N){
   y[i] ~ poisson(mu[i]);
  }

}



