data{
  int <lower=1> T; //Time Index
  int <lower=1> A; // Age Index
  int <lower=1> R; //Region Index

  int <lower=0> y[T*A*R]; // count outcomes
  vector<lower=0>[T*A*R] E;// exposure
  int <lower=1> M; // How much wider is Age interval to Time

  int<lower = 1> TFor;// number of forecast years

} transformed data {
  int N = T*A*R; //Total Number of Observations

  vector[N] log_E = log(E); // log of exposure
  int<lower = 1> Pred = TFor*A*R; // size of prediction vector

  int T_2 = T/2; // half way point of time index
  int L_K = (T-4); // length kappa_time_tilde
  int L_K_2 = L_K / 2; // half length of kappa_time_tilde

}
parameters{

  // Age Parameters
  vector[A] alpha_age; // Age Specific Isntercept
  unit_vector[A] beta_age1; //Unit of 1

  // Time Parameter
  vector<lower=0>[T-4] kappa_time_tilde; //
  positive_ordered[3] kappa_tails; // tails of kappa_tilde
  real<lower = 0> drift; // drift Parameter with positivity constraint

  //AgeRegion Interaction Terms
  matrix[R-1, A] deltaMat_tilde; // Matrix of  age-region specific effects

  //Variances
  real<lower=0> sigma_time;    // Variance of Time
  real<lower=0> sigma_eps;     //Variance of epsilon
  vector<lower=0>[A] sigma_age; // vector of variances

  //Error Terms
  vector[N] eps; // Vector of overdispersion effects
}
transformed parameters{

  vector[T] kappa_time;
  matrix[R, A] deltaMat;

  kappa_time[1] = 0;
  kappa_time[2] = kappa_tails[1]*sigma_time + drift;

  kappa_time[3:(T_2)] = kappa_time_tilde[1:L_K_2] * sigma_time + drift;
  kappa_time[T_2+1] = kappa_tails[2]*sigma_time+drift;

  kappa_time[(T_2+2):(T-1)] = kappa_time_tilde[(L_K_2+1):L_K] * sigma_time + drift;
  kappa_time[T] = kappa_tails[3]* sigma_time + drift;

  // Age Region Interaction w/ corner constraint
  deltaMat[R, ] = rep_row_vector(0, A); // Corner Constraint on last region
  deltaMat[1:(R-1), ] = deltaMat_tilde;

}
model {
   // Stan only supports variable definitions at top of the block
  vector[N] mu; // Vector of random effects
  int pos = 1;
  for(t in 1:T) for(r in 1:R) for (a in 1:A){
    mu[pos]=alpha_age[a]+beta_age1[a]*kappa_time[t]+deltaMat[r, a];
    pos += 1; //pos = pos + 1
   }

  //Varianzen
  target += student_t_lpdf(sigma_time|5,0,1);
  target += student_t_lpdf(sigma_eps|5,0,1);
  target += student_t_lpdf(sigma_age|5,0,1);

 //Other Effects
  target += normal_lpdf(drift|0, 1); // prior Intercept (medium vague)
  target += normal_lpdf(eps|0, 1); //Epsilon Effects
  
  //TIME Effect
  target += normal_lpdf(kappa_tails[1]|0,1); // first element of RW
  target += normal_lpdf(kappa_time_tilde[1]|kappa_tails[1],1); // first element of RW

  target += normal_lpdf(kappa_time_tilde[2:(L_K_2)]|kappa_time_tilde[1:(L_K_2-1)],1); // first half of kappa_time_tilde

  target += normal_lpdf(kappa_tails[2]|kappa_time_tilde[L_K_2],1); // middle element of order constraint

  target += normal_lpdf(kappa_time_tilde[L_K_2+1]|kappa_tails[2],1); // first half
  target += normal_lpdf(kappa_time_tilde[(L_K_2+2):(L_K)]|kappa_time_tilde[(L_K_2+1):(L_K-1)],1); // second half of kappa_time_tilde
  target += normal_lpdf(kappa_tails[3]|kappa_time_tilde[L_K],1); // last element of RW


  //Age Effect
  target += normal_lpdf(alpha_age|0, sigma_age); // Prior on alpha_x

  target += normal_lpdf(beta_age1|0, 5); // unit vector by construction
  // Age Region Interaction effect
  for(a in 1:A){
     target += normal_lpdf(deltaMat_tilde[, a]|0, 10);
  }

  target += poisson_log_lpmf(y| log_E + mu + eps*sigma_eps);


} generated quantities {

  vector[N] log_like_y;
  vector[N] mu_rep;
  vector[N] y_rep;
  int pos = 1;

  // Quantites for Forecast
  vector[TFor] kappa_time_pred;
  vector[Pred] mufor; // predicted death rates
  int pos_f = 1;
  // get mu for insample fit
  for(t in 1:T) for(r in 1:R) for (a in 1:A){
    mu_rep[pos]=exp(alpha_age[a]+beta_age1[a]*kappa_time[t]+
                        deltaMat[r, a]+
                        normal_rng(0,1)*sigma_eps);

    if(mu_rep[pos] >100){
      mu_rep[pos] = 0.05;
    }
    y_rep[pos] = poisson_rng(mu_rep[pos]*E[pos]);

    log_like_y[pos] = poisson_log_lpmf(y[pos]|alpha_age[a]+beta_age1[a]*kappa_time[t]+
                                              deltaMat[r, a]+
                                              eps[pos]*sigma_eps+
                                              log_E[pos]);
    pos += 1; //pos = pos + 1
   }


  kappa_time_pred[1] = drift+kappa_time[T]+sigma_time * normal_rng(0,1);

  // Check if Forecast Period is greater than 1
  if(TFor > 1){
  //RW (1) - Drift Cohort Index
   for (h in 2:TFor) kappa_time_pred[h] = drift+kappa_time_pred[h - 1] + sigma_time * normal_rng(0,1);
   }

  //See that it has the same structure as loop in model block (else predictions are not comparable)
  for (h in 1:TFor) for (r in 1:R) for (a in 1:A)  {
    mufor[pos_f] = alpha_age[a] +
                   beta_age1[a] * kappa_time_pred[h]+
                   deltaMat[r, a]+
                   normal_rng(0,1) * sigma_eps;
    pos_f += 1;
  }
} 

