/////////////////////////////////////
// This Stan model fits a Poisson
// GLARMA(1,1) model with traditional
// horseshoe priors on the regression
// coefficients
/////////////////////////////////////

data{

  int<lower = 2> N;            // length of the time series
  int<lower = 0> P0;           // number of non-shrinking coefficients
  int<lower = 0> P;            // number of shrinking coefficients
  int<lower = 0> y[N];         // conditional Poisson responses
  vector[N] y_star;            // augmented responses
  matrix[N, P0] X_alpha;       // model matrix for non-shrinking coefficients
  matrix[N, P] X_beta;         // model matrix for shrinking coefficients
  real<lower = 0> tau0;        // scale of the non-zero effects global shrinkage

}


parameters{

  // sparse regression components
  vector[P0] alpha_std;                 // standardized un-shrunk coefficients
  vector[P] beta_std;                   // standardized shrunk coefficients
  real<lower = 0> tau;                  // global shrinkage parameter
  vector<lower = 0>[P] lambda;          // local shrinkage parameters

  // time series components
  real<lower = 0, upper = 1> phi;       // AR coefficient
  real<lower = 0, upper = 1> theta;     // MA coefficient
}


transformed parameters{

  vector[P] beta;             // scaled coefficients
  vector[N] mu;               // declare vector of means
  vector[N] z;                // declare vector of ARMA variables

  //scale the alphas
  vector[P0] alpha = alpha_std * 5;

  // scale the betas
  beta = beta_std .* lambda * tau;

  // construct linear predictors
  mu[1] = y_star[1];
  z[1] = 0;
  for(i in 2:N){

    z[i] = phi * (log(y_star[i-1]) - X_alpha[i-1, ] * alpha - X_beta[i-1, ] * beta) + theta * log(y_star[i-1] / mu[i-1]);
    mu[i] = exp(X_alpha[i, ] * alpha + X_beta[i, ] * beta + z[i]);

  }

}


model{

  // priors
  alpha_std ~ std_normal();
  beta_std ~ std_normal();

  tau ~ cauchy(0, tau0);
  lambda ~ cauchy(0, 1);

  // likelihood
  for(i in 1:N){
   y[i] ~ poisson(mu[i]);
  }

}





