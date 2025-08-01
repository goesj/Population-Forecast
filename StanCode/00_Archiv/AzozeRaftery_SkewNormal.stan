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
  int<lower = 0> H; // Forecast Horizon
}

parameters {
  vector[R] mu_raw;
  real m; // global hyperparaemter for group specific intercept
  real<lower = 0> s; // global hyperparameter for group specific variance
  
  
  real sig_mean; 
  real<lower = 0> sig_sd;
  // AR(1) parameter
  vector<lower=-1, upper = 1>[R] a_r; 
  vector<lower=0>[R] sigma;
  
  // Skewness Parameter globals
  vector <lower = 0>[R] gamma; 
}
transformed parameters{
  vector[R] mu = m + s*mu_raw; // non-central parameterization
}
model {
  
  // Priors for mean (hierarchical)
  target += normal_lpdf(mu_raw| 0, 1);
  //target += normal_lpdf(mu|m,s); //alternative
  target += normal_lpdf(m|0, 0.1); //weakly informative priors
  
  // priors of variances
  target += normal_lpdf(s|0, 0.5);
  target += normal_lpdf(sigma|sig_mean, sig_sd); 
  
  target += normal_lpdf(sig_mean| 0 , 0.1);
  target += normal_lpdf(sig_sd| 0, 0.1); 
  
  // prior of AR parameter
  target += normal_lpdf(a_r|0, 2);
  
  // Skewnes Parameter (global)
  target += normal_lpdf(gamma| 1, 1); 
   
  
  for(r in 1:R){ // for all regions
    target += normal_skewed_lpdf(y[r,1]|mu[r]*(1-a_r[r]), sigma[r], gamma[r]);
    for(t in 2:T){
      target += normal_skewed_lpdf(y[r,t]|a_r[r]*y[r,t-1]+mu[r]*(1-a_r[r]), 
                                          sigma[r], gamma[r]);
    }
  }
}

generated quantities{
  vector[R*T] loglike; // vector of log likelihood for WAIC
  matrix[R, T] y_rep; // Y_Rep Matrix for PPC
  matrix[R, H] y_for; // forecated y 
  matrix[R, H] mu_for; // for exact log score of future OOS values
  matrix[R, H] error_for; // forecasting of skewed error
  
  int pos = 1;
   
   // In-Sample Measures
   for(r in 1:R){
     loglike[pos] = normal_skewed_lpdf(y[r, 1] | mu[r]*(1-a_r[r]), sigma[r], gamma[r]);
     
     y_rep[r, 1] = mu[r]*(1-a_r[r]) + 
                   normal_skewed_rng(0, sigma[r], gamma[r]); 
     pos += 1; //pos = pos + 1
      for(t in 2:T){
        loglike[pos] = normal_skewed_lpdf(y[r, t]|a_r[r]*y[r,t-1]+mu[r]*(1-a_r[r]), 
                                                  sigma[r], gamma[r]);
                                                  
        y_rep[r, t] = mu[r]*(1-a_r[r])+a_r[r]*y[r,t-1]+
                      normal_skewed_rng(0, sigma[r], gamma[r]); 
        pos += 1; //pos = pos + 1
      }
    }
    
    
  // Forecasted Values
  for(r in 1:R){
    // y_for[r,1] = mu[r]*(1-a_r[r])+a_r[r]*y[r,T] + 
    //              normal_skewed_rng(0, sigma[r], gamma);
                 
                 
    mu_for[r, 1] = mu[r]*(1-a_r[r])+a_r[r]*y[r,T] ; 
    error_for[r, 1] = normal_skewed_rng(0, sigma[r], gamma[r]);
    
    y_for[r, 1] = mu_for[r, 1] + error_for[r, 1]; 
  
    // mu_for[r, 1] = mu[r]*(1-a_r[r])+a_r[r]*y[r,T] + 
    //                mean_normal_skewed(0, sigma[r], gamma); 
                 
    for(h in 2:H){
      // y_for[r,h] = mu[r]*(1-a_r[r])+a_r[r]*y_for[r,h-1] + 
      //              normal_skewed_rng(0, sigma[r], gamma);
                   
      mu_for[r, h] = mu[r]*(1-a_r[r])+a_r[r]*y_for[r,h-1] ; 
      error_for[r, h] = normal_skewed_rng(0, sigma[r], gamma[r]);
      
      y_for[r, h] = mu_for[r, h] + error_for[r, h]; 
        // 
      // mu_for[r, h] = mu[r]*(1-a_r[r])+a_r[r]*y_for[r,h-1] +      // for log score
      //             mean_normal_skewed(0, sigma[r], gamma);               
    }
}

}


     