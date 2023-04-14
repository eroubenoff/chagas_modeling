functions {
  real icar_normal_lpdf(vector phi, int N, int[] node1, int[] node2) {
    return -0.5 * dot_self(phi[node1] - phi[node2]);
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
  real beta0; 
  real alpha0;
  
  // Spatial (phi) and aspaital (theta) error terms
  // vector[N] phi_lambda[T];        
  vector[N] phi_pi[T];        
  vector[N] theta_lambda[T];    
  vector[N] theta_pi[T];    
  real<lower=0> sigma_lambda;
  real<lower=0> sigma_pi;
  
  // real<lower=0, upper=1> rho_lambda;
  real<lower=0, upper=1> rho_pi;
}
transformed parameters{
  vector[N] pi[T];   // Bernoulli GLM term 
  vector[N] lambda[T];  // Poisson GLM term
  
  vector[N] epsilon[T]; // Poisson error
  vector[N] eta[T]; // Bernoulli error
  
  for (t in 1:T) {
    // Poisson error
    epsilon[t] = sigma_lambda * theta_lambda[t]; /*(sqrt(1 - rho_lambda) * theta_lambda[t] + sqrt(rho_lambda/scaling_factor) * phi_lambda[t]);*/
    // Bernoulli error
    eta[t] = sigma_pi * (sqrt(1 - rho_pi) * theta_pi[t] + sqrt(rho_pi/scaling_factor) * phi_pi[t]);
  }
  
  for (t in 1:T) {
    lambda[t] = log_E[t] + beta0 + epsilon[t];
    pi[t] = alpha0 + eta[t];
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
  //rho_lambda ~ uniform(0,1); 
  rho_pi ~ uniform(0,1); 
  
  sigma_lambda ~ std_normal(); 
  sigma_pi ~ std_normal(); 
  
  for (t in 1:T){
    // Priors on spatial errors 
    //phi_lambda[t] ~ icar_normal(N, node1, node2);
    phi_pi[t] ~ icar_normal(N, node1, node2);
    // Priors on aspaial errors
    theta_lambda[t] ~ std_normal();
    theta_pi[t] ~ std_normal();
    
    // soft sum-to-zero constraint on phi
    // more efficient than mean(phi) ~ normal(0, 0.001)
    // NOT NEEDED with scaling factor-- makes rho unidentifiable
    //sum(phi_lambda[t]) ~ normal(0, 0.001 * N);
    //sum(phi_pi[t]) ~ normal(0, 0.001 * N);
  }
  
  for (t in 1:T) {
    // OLD: Regular Poisson
    // y[t] ~ poisson_log(log_E[t] + lambda[t]);
    
    // ZIP model:
    for (n in 1:N){
      // ZIP
      a = bernoulli_logit_lpmf(1 | pi[t,n]);
      b = bernoulli_logit_lpmf(0 | pi[t,n]);
      c = poisson_log_lpmf(y[t,n] | lambda[t,n]);
        
      if (y[t,n] == 0) {
        target += log_sum_exp(a, b+c);
      } else {
        target += b + c;
      }
    }
  }
}
generated quantities {
  vector[N] E_y[T]; // Expected count 
  
  for (t in 1:T) {
    E_y[t] = (1-inv_logit(pi[t])) .* exp(lambda[t]);
  }
  
}
