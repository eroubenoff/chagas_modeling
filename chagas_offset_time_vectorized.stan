functions {
  real icar_normal_lpdf(vector phi, int N, int[] node1, int[] node2) {
    // Soft sum-to-zero constraint
    return -0.5 * dot_self(phi[node1] - phi[node2]) + normal_lpdf(sum(phi) | 0, 0.001*N);
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
  // indices of zero counts
  int zero_idx[T,N];
  // Max number of zero counts
  int zero_max[T];
  // indices of nonzero counts
  int nonzero_idx[T,N];
  // max number of nonzero counts
  int nonzero_max[T];
  
  for (t in 1:T) {
    log_E[t] = to_vector(log(E[t,1:N]));
  }
  
  // Loop through the y data and tally zeros and nonzeros appropriately
  for (t in 1:T) {
    zero_max[t] = 0;
    nonzero_max[t] = 0;
    for (n in 1:N){
      if (y[t,n] == 0) {
        zero_max[t] += 1;
        zero_idx[t, zero_max[t]] = n;
      }
      else {
        nonzero_max[t] += 1;
        nonzero_idx[t, nonzero_max[t]] = n;
      }
    }
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
  
  for (t in 1:T) {
    lambda[t] = log_E[t] + beta0[t] + u_lambda[t];
    pi[t] = alpha0[t] + sigma_pi[t] * (sqrt(1 - rho_pi[t]) * theta_pi[t] + sqrt(rho_pi[t]/scaling_factor) * phi_pi[t]);
  }
}
model {
  
  vector[N] pi_inv_logit; 
  vector[N] lambda_exp;
  
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
    
    // Vectorized ZIP
    pi_inv_logit = inv_logit(pi[t]);
    lambda_exp = exp(lambda[t]);
    
    // Zeros
    target += sum(log(
           pi_inv_logit[zero_idx[t, 1:zero_max[t]]] + 
        (1-pi_inv_logit[zero_idx[t, 1:zero_max[t]]]) .* 
        exp(-lambda_exp[zero_idx[t, 1:zero_max[t]]])
      ));
    
    // Nonzeros
    target += bernoulli_lpmf(
                rep_array(0, nonzero_max[t]) | 
                pi_inv_logit[nonzero_idx[t,1:nonzero_max[t]]]
          ) + poisson_lpmf(
                y[t, nonzero_idx[t, 1:nonzero_max[t]]] | 
                lambda_exp[nonzero_idx[t, 1:nonzero_max[t]]]
          );
      
  }
}
generated quantities {
  vector[N] E_y[T]; // Expected count 
  
  for (t in 1:T) {
    E_y[t] = (1-inv_logit(pi[t])) .* exp(lambda[t]);
  }
  
}
