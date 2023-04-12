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
  
  // Zeros and nonzero structure-- refer to main code
  // All are indexed by T and N
  int zeros[T,N]; 
  int max_zeros[T];
  int nonzeros[T,N];
  int max_nonzeros[T];

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
  
  // Temporally AR params
  /* vector[T] real delta_psi_1;*/
  /* vector[T] real delta_pi_1; */
  /* vector[T] real delta_psi_2;*/
  /* vector[T] real delta_pi_2; */
  
  // Spatial (phi) and aspaital (theta) error terms
  vector[N] phi_psi[T];        
  vector[N] phi_pi[T];        
  vector[N] theta_psi[T];    
  vector[N] theta_pi[T];    
  real<lower=0, upper=1> rho_psi[T];
  real<lower=0, upper=1> rho_pi[T];
  real<lower=0> sigma_psi[T];
  real<lower=0> sigma_pi[T];
}
transformed parameters{
  vector[N] pi[T];   // Bernoulli GLM term 
  vector[N] psi[T];  // Poisson GLM term
  
  vector[N] epsilon[T]; // Poisson error
  vector[N] eta[T]; // Bernoulli error
  
  for (t in 1:T) {
    // Poisson error
    epsilon[t] = sigma_psi[t] * (sqrt(1 - rho_psi[t]) *theta_psi[t] + sqrt(rho_psi[t]/scaling_factor) * phi_psi[t]);
    // Bernoulli error
    eta[t] = sigma_pi[t] * (sqrt(1 - rho_pi[t]) *theta_pi[t] + sqrt(rho_pi[t]/scaling_factor) * phi_pi[t]);
  }
  
  // Now construct the GLm eqns
  psi[1] = beta0 + epsilon[1];
  pi[1] = beta0 + eta[1];
  psi[2] = beta0 + /*delta_psi_1 * psi[1]*/ + epsilon[2];
  pi[2] = beta0 + /*delta_pi_1 * pi[1]*/ + eta[2];
  for (t in 3:T) {
    psi[t] = beta0 + /*delta_psi_1 * psi[t-1] + delta_psi_2*psi[t-2] */ + epsilon[t];
    pi[t] = alpha0 + /*delta_pi_1 * pi[t-1] + delta_pi_2*pi[t-2] */+ eta[t];
  }
  
}
model {
  real a;
  real b;
  real c;
  // Priors for intercepts
  beta0 ~ uniform(-100, 100);//normal(0.0, 1.0);
  alpha0~ uniform(-100, 100);//normal(0.0, 1.0);
  
  // Priors for temporal AR
  /*
  delta_psi_1 ~ std_normal();
  delta_pi_1 ~ std_normal();
  delta_psi_2 ~ std_normal();
  delta_pi_2 ~ std_normal();
  */
  
  // Priors on rho and sigma
  rho_psi ~ uniform(0,1);
  rho_pi ~ uniform(0,1);
  sigma_psi ~ std_normal();
  sigma_pi ~ std_normal();
  
  for (t in 1:T){
    // Priors on spatial errors 
    phi_psi[t] ~ icar_normal(N, node1, node2);
    phi_pi[t] ~ icar_normal(N, node1, node2);
    // Priors on aspaial errors
    theta_psi[t] ~ std_normal();
    theta_pi[t] ~ std_normal();
  }
  
  for (t in 1:T) {
    // OLD: Regular Poisson
    // y[t] ~ poisson_log(log_E[t] + psi[t]);
    
    for (n in 1:N){
      /*
      // ZIP
      real a = bernoulli_logit_lpmf(1 | pi[t,n]);
      real b = bernoulli_logit_lpmf(0 | pi[t,n]);
      real c = poisson_log_lpmf(y[t,n] | log_E[t,n] + psi[t,n]);
        
      if (y[t,n] == 0) {
        target += log_sum_exp(a, b+c);
      }
      else {
        target += b + c;
      }
    */
    
      // Hurdle model 
      /*
      if (y[t,n] == 0) {
        target += log_inv_logit(pi[t,n]);
      }
      else {
        real a = log1m_inv_logit(pi[t,n]);
        real b = poisson_log_lpmf(y[t,n] | log_E[t,n] + psi[t,n]);
        real c = log1m_exp(-exp(log_E[t,n] + psi[t,n])); // Prev: log1m_exp. Return to this...
        target +=  a +  b + c;
      }
      */
    }
    
    // Vectorized hurdle model:
    // First, handle zeros, which are indexed as:
    // zeros[t, 1:max_zeros[t]]
    target += sum(log_inv_logit(pi[t, zeros[t, 1:max_zeros[t]]]));
    
    // Then handle nonzeros:
    a = sum(log1m_inv_logit(pi[t,nonzeros[t, 1:max_nonzeros[t]]]));
    b = poisson_log_lpmf(y[t,nonzeros[t, 1:max_nonzeros[t]]] | log_E[t,nonzeros[t, 1:max_nonzeros[t]]] + psi[t,nonzeros[t, 1:max_nonzeros[t]]]);
    c = sum(log1m_exp(-exp(log_E[t,nonzeros[t, 1:max_nonzeros[t]]] + psi[t,nonzeros[t, 1:max_nonzeros[t]]]))); // Prev: log1m_exp. Return to this...
    target +=  a +  b + c;
    
    
    // soft sum-to-zero constraint on phi
    // more efficient than mean(phi) ~ normal(0, 0.001)
    sum(phi_psi[t]) ~ normal(0, 0.001 * N);
    sum(phi_pi[t]) ~ normal(0, 0.001 * N);
  }
}
