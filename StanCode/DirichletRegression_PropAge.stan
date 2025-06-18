data {
  int<lower=0> T; // Number of Data Points
  int<lower=0> A; // Number of Age Groups/Categories
  int<lower=0> R; // Number of regions
  
  array[R, T] vector[A] y; // vector of data points per region & time
  vector[T] X_t; // vector of linear trend (for linear trend)
  
  int<lower = 0> H; // Forecasting Horizon
}
parameters {
  vector[A] InterceptVec;
  vector[A] betaVec;
  matrix[R-1, A] deltaMat_tilde; // Matrix of  age-region specific effects
}
transformed parameters {
  array[R, T] vector[A] eta; // array of mean vectors
  matrix[R, A] deltaMat; 
  
  deltaMat[, 1] = rep_vector(0, R); // Corner constraint on first age group 
  deltaMat[1:(R-1), ] = deltaMat_tilde; 
  
  for (t in 1:T) {
    for (r in 1:R) {
      for (a in 1:A) {
        eta[r, t][a] = InterceptVec[a] + betaVec[a] * X_t[t] + deltaMat[r,a];
      }
    }
  }
}
model {
  target += normal_lpdf(InterceptVec | 0, 5); 
  target += normal_lpdf(betaVec | 0, 5);
  
  for(a in 1:A){
    target += normal_lpdf(deltaMat_tilde[,a]| 0, 5); 
  }

  
  for (t in 1:T) {
    for (r in 1:R) {
      target += dirichlet_lpdf(y[r, t] | exp(eta[r, t])); // Likelihood with vector input
    }
  }
}

generated quantities{
  vector[R*T] loglike; // vector of log likelihood for WAIC
  array[R, T] vector[A] alpha = exp(eta); // Parameters of Dirichlet for PPC 
  array[R, T] vector[A] Y_Rep; 

  array[R, H] vector[A] Y_For; 
  array[R, H] vector[A] eta_For; 
  vector[H] X_For;
  
  int pos = 1;
  
  for(t in 1:T) for(r in 1:R){
    loglike[pos] = dirichlet_lpdf(y[r, t] | exp(eta[r, t]));
    
    Y_Rep[r,t] = dirichlet_rng(alpha[r, t]); // create new values of Y for PPC
    
    pos += 1; //pos = pos + 1 
  }
  
  for(h in 1:H){
    X_For[h] = T+h; 
    for(r in 1:R){
      for(a in 1:A){
        eta_For[r, h][a] = InterceptVec[a] + betaVec[a] * X_For[h] + 
                                             deltaMat[r,a]; // calculate linear prediction
      }
      Y_For[r, h] = dirichlet_rng(exp(eta_For[r, h]));  // Generate some forecasts
    }
  }
}
