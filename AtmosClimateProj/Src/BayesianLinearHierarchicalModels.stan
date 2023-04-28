//Standard Hierarchical Linear Model
//Notes:
data {
  int<lower=0> N;   // number of observations items through time, across regions
  int<lower=0> NP;   // number of time periods
  int<lower=0> P[N];   // pointer vector for number of periods
  int<lower=0> NR;   // number of regions
  int<lower=0> R[N];   // pointer vector for number of regions
  int<lower=0> K;   // number of predictors (only 1 for this analysis)
  matrix[N, K] x;   // predictor matrix - i.e., SST, etc.
  vector[N] y;      // outcome vector - first pass this is upwelling
}
parameters {
  real alphaP[NP];  // intercepts 
  real betaP[NP];  // coefficients for predictors
//  real alphaR[NR];  // intercepts 
  real betaR[NR];  // coefficients for predictors
  real<lower=0> sigma;  // error scale
}
model {
    // priors
  for(p in 1:NP) {
  alphaP[p] ~ normal(0, 10);
  betaP[p] ~ normal(0, 10);
  }
  for(r in 1:NR) {
 // alphaR[r] ~ normal(0, 10);
  betaR[r] ~ normal(0, 10);
  }
  sigma ~ normal(0, 10);
  
  //LIKELIHOOD
  for(n in 1:N){
    y[n] ~ normal(x[n] * betaP[P[n]]+x[n] * betaR[R[n]] + alphaP[P[n]], sigma);  // likelihood
  }
}