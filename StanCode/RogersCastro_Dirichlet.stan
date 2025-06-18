data {
  int<lower=0> T; // Number of Data Points
  int<lower=0> A; // Number of Age Groups/Categories
  int<lower=0> R; // Number of regions
  
  array[R, T] vector[A] y; // vector of data points per region & time
}
parameters {
  //vector[A] InterceptVec;
  matrix[R, A] gammaMat; 
  
  //matrix[R-1, A-1] gammaMat_tilde; // Matrix of  age-region specific effects
  //vector[R-1] RegionVec_tilde; 

}
transformed parameters {
  array[R, T] vector[A] eta; // array of mean vectors
  //matrix[R, A] gammaMat; 
  //vector[R] RegionVec; 
  
  // RegionVec[R] = 0; // Corner Constraint
  // RegionVec[1:(R-1)] = RegionVec_tilde; 
  
  //gammaMat[, 1] = rep_vector(0, R); // Corner constraint on first age group 
  // gammaMat[R, ] = rep_row_vector(0, A);
  // gammaMat[1:(R-1), ] = gammaMat_tilde;
  // 
  for (t in 1:T) {
    for (r in 1:R) {
      for (a in 1:A) {
        eta[r, t][a] =  gammaMat[r,a]; //InterceptVec[a] + RegionVec[r]; //+ gammaMat[r,a]; 
      }
    }
  }
}
model {
  // target += normal_lpdf(InterceptVec | 0, 5); 
  // target += normal_lpdf(RegionVec_tilde | 0, 5);

  //  for(a in 1:(A-1)){
  //   target += normal_lpdf(gammaMat_tilde[,a]| 0, 5); 
  // }
   for(a in 1:A){
    target += normal_lpdf(gammaMat[,a]| 0, 5);
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
  
  int pos = 1;
  
  for(t in 1:T) for(r in 1:R){
    loglike[pos] = dirichlet_lpdf(y[r, t] | exp(eta[r, t]));
    
    Y_Rep[r,t] = dirichlet_rng(alpha[r, t]); // create new values of Y for PPC
    
    pos += 1; //pos = pos + 1 
  }
}
