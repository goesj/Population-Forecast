// AR1 Model for TFR
data {
  int<lower=0> T;
  int<lower = 0> R;

  matrix[R,T] y; // TFR values as real matrix
  int<lower = 0> H; // Forecast Horizon
}
parameters {
  vector[R] mu_raw;
  real m; // global hyperparaemter for group specific intercept
  real<lower = 0> s; // global hyperparameter for group specific variance

  // AR(1) parameter
  vector<lower=-1, upper = 1>[R] a_r;
  vector<lower=0>[R] sigma;
}
transformed parameters{
  vector[R] mu = m + s*mu_raw; // non-central parameterization
}
model {
  // Priors for mean (hierarchical)
  target += std_normal_lpdf(mu_raw); // non-central parameterization
  target += normal_lpdf(m|0,5);

  // priors of variances
  target += normal_lpdf(s|0,5);
  target += normal_lpdf(sigma|0,5);

  // prior of AR parameter
  target += normal_lpdf(a_r|0, 2);


  for(r in 1:R){ // for all regions
    target += normal_lpdf(y[r,1]|mu[r], sigma[r]);
    target += normal_lpdf(y[r,2:T]|mu[r]+a_r[r]*(y[r,1:(T-1)]-mu[r]), sigma[r]);
  }
}
generated quantities{
  vector[R*T] loglike; // vector of log likelihood for WAIC
  matrix[R, T] y_rep; // Y_Rep Matrix for PPC
  matrix[R, H] y_for; // forecated y
  int pos = 1;

   // In-Sample Measures
   for(r in 1:R){
     loglike[pos] = normal_lpdf(y[r, 1] | mu[r], sigma[r]);
     y_rep[r, 1] = mu[r] + std_normal_rng()*sigma[r];
    pos += 1; //pos = pos + 1
      for(t in 2:T){
        loglike[pos] = normal_lpdf(y[r, t]|mu[r]+a_r[r]*(y[r,t-1]-mu[r]), sigma[r]);
        y_rep[r, t] = mu[r]+a_r[r]*(y[r,t-1]-mu[r])+std_normal_rng()*sigma[r];
        pos += 1; //pos = pos + 1
      }
    }


  // Forecasted Values
  for(r in 1:R){
    y_for[r,1] = mu[r]+a_r[r]*(y[r,T]-mu[r]) + std_normal_rng()*sigma[r];
    for(h in 2:H){
      y_for[r,h] = mu[r]+a_r[r]*(y_for[r,h-1]-mu[r]) + std_normal_rng()*sigma[r];
    }
  }
} 
