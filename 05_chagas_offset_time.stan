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
  int zeros_array[N] = rep_array(0, N);
  int ones_array[N] = rep_array(1, N);
  
  // Logged population
  vector[N] log_E[T];
  vector[N] log_E_zero[T];
  vector[N] log_E_nonzero[T];
  
  for (t in 1:T) {
    log_E[t] = to_vector(log(E[t,1:N]));
  }
  
  // Partition log_E conditional on y == 0
  for (t in 1:T){
    log_E_zero[t, 1:max_zeros[t]] = log_E[t, zeros[t,1:max_zeros[t]]];
    log_E_nonzero[t, 1:max_nonzeros[t]] = log_E[t, nonzeros[t,1:max_nonzeros[t]]];
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
  
  vector[N] pi_zero[T];
  vector[N] pi_nonzero[T];
  vector[N] psi_zero[T];
  vector[N] psi_nonzero[T];
  
  for (t in 1:T) {
    // Poisson error
    epsilon[t] = sigma_psi[t] * (sqrt(1 - rho_psi[t]) *theta_psi[t] + sqrt(rho_psi[t]/scaling_factor) * phi_psi[t]);
    // Bernoulli error
    eta[t] = sigma_pi[t] * (sqrt(1 - rho_pi[t]) *theta_pi[t] + sqrt(rho_pi[t]/scaling_factor) * phi_pi[t]);
  }
  
  // Now construct the GLm eqns
  psi[1] = beta0 + epsilon[1];
  pi[1] = alpha0 + eta[1];
  psi[2] = beta0 + /*delta_psi_1 * psi[1]*/ + epsilon[2];
  pi[2] = alpha0 + /*delta_pi_1 * pi[1]*/ + eta[2];
  for (t in 3:T) {
    psi[t] = beta0 + /*delta_psi_1 * psi[t-1] + delta_psi_2*psi[t-2] */ + epsilon[t];
    pi[t] = alpha0 + /*delta_pi_1 * pi[t-1] + delta_pi_2*pi[t-2] */+ eta[t];
  }
  
  // Partition pi and  psi into separate objects conditioned on y==0
  profile("partition"){
  for (t in 1:T){
    psi_zero[t, 1:max_zeros[t]] = psi[t, zeros[t, 1:max_zeros[t]]];
    pi_zero[t, 1:max_zeros[t]] = pi[t, zeros[t, 1:max_zeros[t]]];
    
    psi_nonzero[t, 1:max_nonzeros[t]] = psi[t, nonzeros[t, 1:max_nonzeros[t]]];
    pi_nonzero[t, 1:max_nonzeros[t]] = pi[t, nonzeros[t, 1:max_nonzeros[t]]];
    
  }
  }
  
  
  
}
model {
  real a;
  real b;
  real c;
  real d;
  real e;
  // Priors for intercepts
  profile("prior alpha beta") {
  beta0 ~ normal(0, 10);//normal(0.0, 1.0);
  alpha0~ normal(0, 10);//normal(0.0, 1.0);
  }
  
  // Priors for temporal AR
  /*
  delta_psi_1 ~ std_normal();
  delta_pi_1 ~ std_normal();
  delta_psi_2 ~ std_normal();
  delta_pi_2 ~ std_normal();
  */
  
  // Priors on rho and sigma
  profile("prior rho") {
  rho_psi ~ normal(0.5, 2); //uniform(0,1);
  rho_pi ~ normal(0.5, 2); //uniform(0,1);
  }
  profile("prior sigma") {
  sigma_psi ~ normal(0, 5); //std_normal();
  sigma_pi ~ normal(0, 5); //std_normal();
  }
  
  profile("priors 2") {
  for (t in 1:T){
    // Priors on spatial errors 
    phi_psi[t] ~ icar_normal(N, node1, node2);
    phi_pi[t] ~ icar_normal(N, node1, node2);
    // Priors on aspaial errors
    theta_psi[t] ~ std_normal();
    theta_pi[t] ~ std_normal();
  }
  }
  
  for (t in 1:T) {
    // OLD: Regular Poisson
    // y[t] ~ poisson_log(log_E[t] + psi[t]);
    
    for (n in 1:N){
      // ZIP
      profile("a") {
      a = bernoulli_logit_lpmf(1 | pi[t,n]);
      }
      profile("b"){
      b = bernoulli_logit_lpmf(0 | pi[t,n]);
      }
      profile("c"){
      c = poisson_log_lpmf(y[t,n] | log_E[t,n] + psi[t,n]);
      }
        
      if (y[t,n] == 0) {
        target += log_sum_exp(a, b+c);
      }
      else {
        target += b + c;
      }
    
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
    
    // Vectorized ZIP model
    /*
    // Zeros
    profile("a"){
      a= bernoulli_logit_lpmf(head(ones_array, max_zeros[t]) | head(pi_zero[t], max_zeros[t]));
    }
    profile("b"){
      b = bernoulli_logit_lpmf(head(zeros_array, max_zeros[t]) | head(pi_zero[t], max_zeros[t])); 
    }
    profile("c"){
      c = poisson_log_lpmf(head(zeros_array, max_zeros[t]) | head(log_E_zero[t], max_zeros[t]) + head(psi_zero[t], max_zeros[t]));
    }
    profile("zero likelihood"){
    target += log_sum_exp(a, b+c);
    }
    
    // Nonzeros
    profile("d"){
      d = bernoulli_logit_lpmf(head(zeros_array, max_nonzeros[t]) | head(pi_nonzero[t], max_nonzeros[t]));
    }
    profile("e"){
      e = poisson_log_lpmf(y[t, nonzeros[t, 1:max_nonzeros[t]]] | head(log_E_nonzero[t], max_nonzeros[t]) + head(psi_nonzero[t], max_nonzeros[t]));
    }
    profile("nonzero likelihood"){
    target += d + e;
    }
    */
    
    // Vectorized hurdle model:
    // First, handle zeros, which are indexed as:
    // zeros[t, 1:max_zeros[t]]
    
    // The statement we want to evaluate is 1 ~ bernoulli(.)
    // However that is equivalent to target += log(.)
    // Since pi is on the logit scale, need to do target  += log_inv_logit(pi)
    // And then need to sum them: target += log_inv_logit(pi)
    // Is there any way to simplify this?
    // target += sum(log_inv_logit(pi[t, zeros[t, 1:max_zeros[t]]]));
    
    /*
    // Trying instead:
    profile("bernoulli") {
      rep_array(0, max_zeros[t]) ~ bernoulli_logit(pi[t, zeros[t, 1:max_zeros[t]]]);
    }
    // Possible that this statement is equivalent, and faster:
    rep_array(0, max_zeros[t]) ~bernoulli_logit_glm((sqrt(1 - rho_pi[t]) *theta_pi[t] + sqrt(rho_pi[t]/scaling_factor) * phi_pi[t]), alpha0, sigma_pi[t])
    
    // Then handle nonzeros:
    profile("a") {
    a = sum(log1m_inv_logit(pi[t,nonzeros[t, 1:max_nonzeros[t]]]));
    }
    profile("b") {
    b = poisson_log_lpmf(y[t,nonzeros[t, 1:max_nonzeros[t]]] | log_E[t,nonzeros[t, 1:max_nonzeros[t]]] + psi[t,nonzeros[t, 1:max_nonzeros[t]]]);
    }
    profile("c") {
    c = sum(log1m_exp(-exp(log_E[t,nonzeros[t, 1:max_nonzeros[t]]] + psi[t,nonzeros[t, 1:max_nonzeros[t]]]))); // Prev: log1m_exp. Return to this...
    }
    profile("poisson_likelihood") {
    target +=  a +  b + c;
    }
    */
    
    
    // soft sum-to-zero constraint on phi
    // more efficient than mean(phi) ~ normal(0, 0.001)
    sum(phi_psi[t]) ~ normal(0, 0.001 * N);
    sum(phi_pi[t]) ~ normal(0, 0.001 * N);
  }
}
