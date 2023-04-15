functions {
  real icar_normal_lpdf(vector phi, int N, int[] node1, int[] node2) {
    // Soft sum-to-zero constraint
    return -0.5 * dot_self(phi[node1] - phi[node2]) + normal_lpdf(sum(phi) | 0, 0.001*N);
 }
  real ZIP_partial_lpmf(int[] y_slice, int start, int end, vector pi, vector lambda) {
    int slice_size = size(y_slice) - 1;
    real a; real b; real c;
    real ret = 0;
    
    for (i in 1:slice_size) {
      a = bernoulli_logit_lpmf(1 | pi[start+i]);
      b = bernoulli_logit_lpmf(0 | pi[start+i]);
      c = poisson_log_lpmf(y_slice[i] | lambda[start+i]);
        
      if (y_slice[i] == 0) {
        ret += log_sum_exp(a, b+c);
      } else {
        ret += b + c;
      }
    }
    return ret;
  }
}
data {
  // Number of municipalities
  int<lower=0> N; 
  // Number of years
  int<lower=0> T; 
  // Number of adjacent edges
  int<lower=0> N_edges;
  // node1[i] adjacent to node2[i]
  int<lower=1, upper=N> node1[N_edges];  
  // and node1[i] < node2[i]
  int<lower=1, upper=N> node2[N_edges];  
  // count outcomes
  int y[T,N];           
  // Population exposure
  int E[T,N];
  // Scaling factor-- scales variance of spatial effects
  real<lower=0> scaling_factor;
  // Grainsize
  int grainsize;
  
}
transformed data {
  // Logged population
  vector[N] log_E[T];
  
  for (t in 1:T) {
    log_E[t] = to_vector(log(E[t,1:N]));
  }
}
parameters {
  // Intercepts
  vector[T] beta0; 
  vector[T] alpha0;
  
  // Spatial (phi) and aspaital (theta) error terms
  // Bernoulli spatial error 
  vector[N] phi_pi[T];        
  // Bernoulli aspatial error
  vector[N] theta_pi[T];    
  // Scaling of total bernoulli error 
  real<lower=0.00001> sigma_pi[T];
  // Proportion of spatial error
  real<lower=0, upper=1> rho_pi[T];
  
  // Poisson error
  vector[N] u_lambda[T];    
  
}
transformed parameters{
  vector[N] pi[T];   // Bernoulli GLM term 
  vector[N] lambda[T];  // Poisson GLM term
  
  /*(sqrt(1 - rho_lambda) * theta_lambda[t] + sqrt(rho_lambda/scaling_factor) * phi_lambda[t]);*/
  for (t in 1:T) {
    lambda[t] = log_E[t] + beta0[t] + u_lambda[t];
    pi[t] = alpha0[t] + sigma_pi[t] * (sqrt(1 - rho_pi[t]) * theta_pi[t] + sqrt(rho_pi[t]/scaling_factor) * phi_pi[t]);
  }
}
model {
  real a;
  real b;
  real c;
  
  // AR Priors for intercepts
  // Year 1
  beta0 ~ normal(-10, 10);//normal(0.0, 1.0);
  alpha0 ~ normal(-5, 10);//normal(0.0, 1.0);
  
  // Priors on rho and sigma
  rho_pi ~ uniform(0,1); 
  sigma_pi ~ std_normal();
  
  
  for (t in 1:T){
    // Priors on spatial errors 
    phi_pi[t] ~ icar_normal(N, node1, node2);
    
    // Priors on aspaial errors
    theta_pi[t] ~ std_normal();
    u_lambda[t] ~ std_normal(); 
  }
  
  for (t in 1:T) {
    // OLD: Regular Poisson
    // y[t] ~ poisson_log(log_E[t] + lambda[t]);
    
    // ZIP model:
    target += reduce_sum(ZIP_partial_lpmf, y[t], grainsize, pi[t], lambda[t]);
    // for (n in 1:N){
    //   // ZIP
    //   a = bernoulli_logit_lpmf(1 | pi[t,n]);
    //   b = bernoulli_logit_lpmf(0 | pi[t,n]);
    //   c = poisson_log_lpmf(y[t,n] | lambda[t,n]);
    //     
    //   if (y[t,n] == 0) {
    //     target += log_sum_exp(a, b+c);
    //   } else {
    //     target += b + c;
    //   }
    // }
  }
}
generated quantities {
  vector[N] E_y[T]; // Expected count 
  
  for (t in 1:T) {
    E_y[t] = (1-inv_logit(pi[t])) .* exp(lambda[t]);
  }
  
}
