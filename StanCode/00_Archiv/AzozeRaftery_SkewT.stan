// Some functions
functions {
  real student_t_skewed_lpdf(real x, real mu, real sigma, real gamma_param, real df) {
    real z = x - mu;
    real log_c = log(2 / (gamma_param + inv(gamma_param)));
    real log_t = lgamma((df + 1) / 2) - lgamma(df / 2) - 
                 log(sigma * sqrt(df * pi()));
    
    real quad_term;
    if(z >= 0){
      quad_term = log1p( square(z/sigma) * inv_square(gamma_param) / df);
    } else {
      quad_term = log1p( square(z/sigma) * square(gamma_param) / df);
    }
    real log_prob = log_c + log_t - ((df + 1) / 2) * quad_term;

    return log_prob;
  }

  real student_t_skewed_pdf(real x, real mu, real sigma, real gamma_param, real df) {
    return exp(student_t_skewed_lpdf(x| mu, sigma, gamma_param, df));
  }
  
  real student_t_skewed_cdf(real x, real location, real scale, real gamma, real df) {
    real z = x - location;
    real result;
    if (z < 0) {
      result = 2 / (square(gamma) + 1) * student_t_cdf(z*gamma| df, 0 ,scale); 
    } else {
      result = 1 / (square(gamma) + 1) +
               2 / (1 + inv_square(gamma)) *
               (student_t_cdf(z/gamma | df , 0,  scale) - 0.5); 
    }
    return result;
  }
  
  real student_t_qf(real p, real df) { // taken from https://discourse.mc-stan.org/t/student-t-quantile-function/31668
    real eps = 1e-12;
    real d_epsilon = 2.220446e-16;
    real d_max = 1.797693e+308;
    real d_min = 2.225074e-308;
    int d_mant_dig = 53;
    real q;
    
    if (is_nan(p) || is_nan(df)) {
      return p + df;
    }
    
    if (df <= 0) {
      reject("Invalid value for df");
    }
    
    if (df < 1) {
      // Find the upper and lower bounds
      real accu = 1e-13;
      real Eps = 1e-11;
      if (p > 1 - d_epsilon) {
        return positive_infinity();
      }
      real pp = min({1 - d_epsilon, p * (1 + Eps)});
      real ux = 1;
      while (student_t_cdf(ux | df, 0., 1.) < pp) {
        ux = ux * 2;
      }
      pp = p * (1 - Eps);
      real lx = -1;
      while (student_t_cdf(lx | df, 0., 1.) > pp) {
        lx = lx * 2;
      }
      
      // Find the quantile using interval halving
      real nx = 0.5 * (lx + ux);
      int iter = 0;
      while ((ux - lx) / abs(nx) > accu && iter < 1000) {
        iter += 1;
        if (student_t_cdf(nx | df, 0., 1.) > p) {
          ux = nx;
        } else {
          lx = nx;
        }
        nx = 0.5 * (lx + ux);
      }
      return 0.5 * (lx + ux);
    }
    
    if (df > 1e20) {
      return inv_Phi(p);
    }
    
    int neg = p < 0.5 ? 1 : 0;
    int is_neg_lower = neg;
    real P = neg == 1 ? 2 * p : 2 * (0.5 - p + 0.5);
    
    P = min({max({P, 0}), 1});
    
    if (abs(df - 2) < eps) {
      if (P > d_min) {
        if (3 * P < d_epsilon) {
          q = 1 / sqrt(P);
        } else if (P > 0.9) {
          q = (1 - P) * sqrt(2 / (P * (2 - P)));
        } else {
          q = sqrt(2 / (P * (2 - P)) - 2);
        }
      } else {
        q = positive_infinity();
      }
    } else if (df < 1. + eps) {
      if (P == 1.) {
        q = 0;
      } else if (P > 0) {
        q = 1 / tan(pi() * p / 2);
      } else {
        q = negative_infinity();
      }
    } else {
      real x = 0;
      real y;
      real log_P2 = 0;
      real a = 1 / (df - 0.5);
      real b = 48 / (a * a);
      real c = ((20700 * a / b - 98) * a - 16) * a + 96.36;
      real d = ((94.5 / (b + c) - 3) / b + 1) * sqrt(a * pi() / 2) * df;
      
      y = pow((d * P), (2.0 / df));
      int P_ok = y >= d_epsilon ? 1 : 0;
      
      if (P_ok != 1) {
        log_P2 = is_neg_lower == 1 ? log(p) : log1m_exp(p);
        x = (log(d) + log2() + log_P2) / df;
        y = exp(2 * x);
      }
      
      if ((df < 2.1 && P > 0.5) || y > 0.05 + a) {
        if (P_ok == 1) {
          x = inv_Phi(0.5 * P);
        } else {
          x = inv_Phi(log_P2);
        }
        
        y = square(x);
        
        if (df < 5) {
          c += 0.3 * (df - 4.5) * (x + 0.6);
        }
        
        c = (((0.05 * d * x - 5) * x - 7) * x - 2) * x + b + c;
        y = (((((0.4 * y + 6.3) * y + 36) * y + 94.5) / c - y - 3) / b + 1)
            * x;
        y = expm1(a * square(y));
        q = sqrt(df * y);
      } else if (P_ok != 1 && x < -log2() * d_mant_dig) {
        q = sqrt(df) * exp(-x);
      } else {
        y = ((1
              / (((df + 6) / (df * y) - 0.089 * d - 0.822) * (df + 2) * 3)
              + 0.5 / (df + 4))
             * y - 1)
            * (df + 1) / (df + 2) + 1 / y;
        
        q = sqrt(df * y);
      }
      
      if (P_ok == 1) {
      int it = 0;
      while (it < 10) {
        y = exp(student_t_lpdf(q | df, 0., 1.));
        if (y <= 0 || is_inf(y)) {
          break;
        }
        
        real t = (exp(student_t_lccdf(q | df, 0., 1.)) - P / 2) / y;
        if (abs(t) <= 1e-14 * abs(q)) {
          break;
        }
        
        q = q + t * (1. + t * q * (df + 1) / (2 * (q * q + df)));
        
        it += 1;
      }
      }
    }
    
    if (neg) {
      q = -q;
    }
    
    return q;
  }
  
  real student_t_skewed_quantile(real p, real location, real scale, real gamma, real df) {
    real probzero = student_t_skewed_cdf(0, 0, scale, gamma, df); // Location is handled at end
    real result;
    
    if (p < probzero) {
      result = inv(gamma) *
               student_t_qf(((square(gamma) + 1) * p) / 2, df) * scale;
    } else {
      result = gamma * 
               student_t_qf((1 + inv_square(gamma)) / 2 * (p - probzero) + 0.5, 
               df)*scale; 
              
    }

    return result + location;
  }
  
  real student_t_skewed_rng(real location, real scale, real gamma, real df) {
    real u = uniform_rng(0, 1);
    real z = student_t_skewed_quantile(u, location, scale, gamma, df); 
    return(z);
  }
  
  real mean_student_t_skewed(real location, real scale, real gamma, real df){
      // see Fernandez, Steel 1998 (Formula 5)
     real absMoment = sqrt(df/pi())*tgamma((df-1)/2)/tgamma(df/2); // E(|X|) central t
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
  
  // AR(1) parameter
  vector<lower=-1, upper = 1>[R] a_r; 
  vector<lower=0>[R] sigma;
  
  // Skewness Parameter globals
  real <lower = 0> gamma; 
  real <lower = 2> nu; // degrees of freedom parameter
}
transformed parameters{
  vector[R] mu = m + s*mu_raw; // non-central parameterization
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
  target += normal_lpdf(nu| 5, 5); 
   
  
  for(r in 1:R){ // for all regions
    target += student_t_skewed_lpdf(y[r,1]|mu[r]*(1-a_r[r]), sigma[r], gamma, nu);
    for(t in 2:T){
      target += student_t_skewed_lpdf(y[r,t]|a_r[r]*y[r,t-1]+mu[r]*(1-a_r[r]), 
                                          sigma[r], gamma, nu);
    }
  }
}

generated quantities{
  vector[R*T] loglike; // vector of log likelihood for WAIC
  matrix[R, T] y_rep; // Y_Rep Matrix for PPC
  matrix[R, H] y_for; // forecated y 
  matrix[R, H] mu_for; // for exact log score of future OOS values
  
  int pos = 1;
   
   // In-Sample Measures
   for(r in 1:R){
     loglike[pos] = student_t_skewed_lpdf(y[r, 1] | mu[r]*(1-a_r[r]), sigma[r], gamma, nu);
     
     y_rep[r, 1] = mu[r]*(1-a_r[r]) + 
                   student_t_skewed_rng(0, sigma[r], gamma, nu); 
     pos += 1; //pos = pos + 1
      for(t in 2:T){
        loglike[pos] = student_t_skewed_lpdf(y[r, t]|a_r[r]*y[r,t-1]+mu[r]*(1-a_r[r]), 
                                                  sigma[r], gamma, nu);
                                                  
        y_rep[r, t] = mu[r]*(1-a_r[r])+a_r[r]*y[r,t-1]+
                      student_t_skewed_rng(0, sigma[r], gamma, nu); 
        pos += 1; //pos = pos + 1
      }
    }
    
    
  // Forecasted Values
  for(r in 1:R){
    y_for[r,1] = mu[r]*(1-a_r[r])+a_r[r]*y[r,T] + 
                 student_t_skewed_rng(0, sigma[r], gamma, nu);
                 
    mu_for[r, 1] = mu[r]*(1-a_r[r])+a_r[r]*y[r,T] + 
                   mean_student_t_skewed(0, sigma[r], gamma, nu); 
                 
    for(h in 2:H){
      y_for[r,h] = mu[r]*(1-a_r[r])+a_r[r]*y_for[r,h-1] + 
                   student_t_skewed_rng(0, sigma[r], gamma, nu);
                   
      mu_for[r, h] = mu[r]*(1-a_r[r])+a_r[r]*y_for[r,h-1] +      // for log score
                  mean_student_t_skewed(0, sigma[r], gamma, nu); 
    }
}

}

