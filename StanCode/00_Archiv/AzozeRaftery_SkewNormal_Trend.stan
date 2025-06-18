// Some functions
functions {
  real normal_skewed_lpdf(real x, real mu, real sigma, real gamma_param) {
    real z = x - mu;
    real normalization = log(2 / (gamma_param + inv(gamma_param)));
    real log_sd = -0.5 * log(2 * pi()) - log(sigma);
    real scale_term;
    
    if (z >= 0) {
      scale_term = -0.5 * square(z / sigma) * inv_square(gamma_param);
    } else {
      scale_term = -0.5 * square(z / sigma) * square(gamma_param);
    }
    
    return normalization + log_sd + scale_term;
  }
  
  real normal_skewed_pdf(real x, real mu, real sigma, real gamma_param) {
    return exp(normal_skewed_lpdf(x | mu, sigma, gamma_param));
  }
  
  real normal_skewed_cdf(real x, real location, real scale, real gamma) {
    real z = x - location;
    real result;
    if (z < 0) {
      result = 2 / (square(gamma) + 1) * normal_cdf(z*gamma, 0 , scale); 
    } else {
      result = 1 / (square(gamma) + 1) +
               2 / (1 + inv_square(gamma)) *
               (normal_cdf(z/gamma, 0, scale) - 0.5); 
    }
    return result;
  }
  
  real normal_skewed_quantile(real p, real location, real scale, real gamma) {
    real probzero = normal_skewed_cdf(0, 0, scale, gamma); // Location is handled at end
    real result;
    
    if (p < probzero) {
      result = inv(gamma) *
               inv_Phi(((square(gamma) + 1) * p) / 2) * scale;
    } else {
      result = gamma * 
               inv_Phi((1 + inv_square(gamma)) / 2 * (p - probzero) + 0.5)*scale; 
              
    }

    return result + location;
  }
  
  real normal_skewed_rng(real location, real scale, real gamma) {
    real u = uniform_rng(0, 1);
    real z = normal_skewed_quantile(u, location, scale, gamma); 
    return(z);
  }
  
  real mean_normal_skewed(real location, real scale, real gamma){
      // see Fernandez, Steel 1998 (Formula 5)
     real absMoment = sqrt(2/pi()); // E(|X|) central t
     
     real scaleFactor = (square(gamma)+(-1/square(gamma)))/(gamma+1/gamma); 
     
     real Mean_scale_location = ((absMoment*scaleFactor)*scale)+location; 
     
     return(Mean_scale_location); // return mean of skewed distribution
  }

}

// The input data is a vector 'y' of length 'N'.
data {
  int<lower=0> T;
  int<lower = 0> R; 
  
  matrix[R,T] y; // Net Migration matrix
  vector[T] X_t; // linear trend component
  int<lower = 0> H; // Forecast Horizon
}

parameters {
  vector[R] mu_raw;
  real m; // global hyperparaemter for group specific intercept
  real<lower = 0> s; // global hyperparameter for group specific variance
  
  // AR(1) parameter
  vector<lower=-1, upper = 1>[R] a_r; 
  vector<lower=0>[R] sigma;
  
  // Skewness Parameter globals
  real <lower = 0> gamma; 
  
  vector[R] beta;  
}
transformed parameters{
  vector[R] mu = m + s*mu_raw; // non-central parameterization
  matrix[R, T] Eps; 
  for(r in 1:R){
    Eps[r, ] = y[r, ] - (X_t*beta[r])'; // subtract linear trend 
  }
}
model {
  
  // Priors for mean (hierarchical)
  target += std_normal_lpdf(mu_raw); // non-central parameterization
  //target += normal_lpdf(mu|m,s); //alternative
  target += normal_lpdf(m|0,1); 
  
  // priors of variances
  target += normal_lpdf(s|0,1);
  target += normal_lpdf(sigma|0,2); 
  
  // prior of AR parameter
  target += normal_lpdf(a_r|0, 2);
  
  // Skewnes Parameter (global)
  target += normal_lpdf(gamma| 1, 1); 
   
  target += normal_lpdf(beta|0, 5); 
  
  for(r in 1:R){ // for all regions
    target += normal_skewed_lpdf(Eps[r,1]|mu[r]*(1-a_r[r]), sigma[r], gamma);
    for(t in 2:T){
      target += normal_skewed_lpdf(Eps[r,t]|a_r[r]*Eps[r,t-1]+mu[r]*(1-a_r[r]), 
                                          sigma[r], gamma);
    }
  }
}

generated quantities{
  vector[R*T] loglike; // vector of log likelihood for WAIC
  matrix[R, T] y_rep; // Y_Rep Matrix for PPC
  matrix[R, H] y_for; // forecated y 
  matrix[R, H] Eps_for; 
  vector[H] X_For;
  
  int pos = 1;
   
   // In-Sample Measures
   for(r in 1:R){
     loglike[pos] = normal_skewed_lpdf(y[r, 1] | X_t[1]*beta[r]+mu[r]*(1-a_r[r]), sigma[r], gamma);
     
     y_rep[r, 1] = X_t[1]*beta[r]+mu[r]*(1-a_r[r]) + 
                   normal_skewed_rng(0, sigma[r], gamma); 
     pos += 1; //pos = pos + 1
      for(t in 2:T){
     
        loglike[pos] = normal_skewed_lpdf(y[r, t]| X_t[t]*beta[r]+a_r[r]*Eps[r,t-1]+
                                                   mu[r]*(1-a_r[r]), sigma[r], gamma);
                                                  
        y_rep[r, t] = X_t[t]*beta[r]+mu[r]*(1-a_r[r])+a_r[r]*Eps[r,t-1]+
                      normal_skewed_rng(0, sigma[r], gamma); 
        pos += 1; //pos = pos + 1
      }
    }
    
    
  // Forecasted Values
  for(r in 1:R){
    X_For[1] = T+1; 
    Eps_for[r,1] = mu[r]*(1-a_r[r])+a_r[r]*y[r,T] + 
                 normal_skewed_rng(0, sigma[r], gamma);
  
    y_for[r, 1] = Eps_for[r, 1] + X_For[1]*beta[r]; 
    for(h in 2:H){
      X_For[h] = T+h; 
      Eps_for[r,h] = mu[r]*(1-a_r[r])+a_r[r]*y_for[r,h-1] + 
                   normal_skewed_rng(0, sigma[r], gamma);
     
      y_for[r, h] = Eps_for[r, h] + X_For[h]*beta[r]; 
    }
}

}

