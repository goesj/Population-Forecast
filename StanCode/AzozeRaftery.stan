// Some functions
data {
  int<lower=0> T;
  int<lower = 0> R; 
  
  matrix[R,T] y; // Net Migration matrix
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
  //target += normal_lpdf(mu|m,s); //alternative
  target += normal_lpdf(m|0,5); 
  
  // priors of variances
  target += normal_lpdf(s|0,2);
  target += normal_lpdf(sigma|0,2); 
  
  // prior of AR parameter
  target += normal_lpdf(a_r|0, 2);
  
  for(r in 1:R){ // for all regions
    target += normal_lpdf(y[r,1]|mu[r]*(1-a_r[r]), sigma[r]);
    for(t in 2:T){
      target += normal_lpdf(y[r,t]|a_r[r]*y[r,t-1]+mu[r]*(1-a_r[r]), sigma[r]);
    }
  }
}

generated quantities{
  vector[R*T] loglike; // vector of log likelihood for WAIC
  matrix[R, T] y_rep; // Y_Rep Matrix for PPC
  matrix[R, H] y_for; // forecated y 
  matrix[R, H] mu_for; // forecasted mean value
  int pos = 1;
   
   // In-Sample Measures
   for(r in 1:R){
     loglike[pos] = normal_lpdf(y[r, 1] | mu[r]*(1-a_r[r]), sigma[r]);
     
     y_rep[r, 1] = mu[r]*(1-a_r[r]) + normal_rng(0, sigma[r]); 
     pos += 1; //pos = pos + 1
      for(t in 2:T){
        loglike[pos] = normal_lpdf(y[r, t]|a_r[r]*y[r,t-1]+mu[r]*(1-a_r[r]), sigma[r]);
                                                  
        y_rep[r, t] = mu[r]*(1-a_r[r])+a_r[r]*y[r,t-1]+normal_rng(0, sigma[r]); 
        pos += 1; //pos = pos + 1
      }
    }
    
    
  // Forecasted Values
  for(r in 1:R){
    y_for[r,1] = mu[r]*(1-a_r[r])+a_r[r]*y[r,T] + normal_rng(0, sigma[r]);
    mu_for[r, 1] = mu[r]*(1-a_r[r])+a_r[r]*y[r,T];
    for(h in 2:H){
      y_for[r,h] = mu[r]*(1-a_r[r])+a_r[r]*y_for[r,h-1] + normal_rng(0, sigma[r]);
      mu_for[r, h] = mu[r]*(1-a_r[r])+a_r[r]*y_for[r,h-1]; 
    }
}

}

