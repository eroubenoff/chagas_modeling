functions {
  real icar_normal_lpdf(vector phi, int N, array[] int node1, array[] int node2) {
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
  array[N_edges] int<lower=1, upper=N> node1;  
  // and node1[i] < node2[i]
  array[N_edges] int<lower=1, upper=N> node2;  
  // count outcomes
  array[T,N] int y;           
  // Population exposure
  array[T,N] int E;
  // Scaling factor-- scales variance of spatial effects
  real<lower=0> scaling_factor;
  // indices of zero counts
  array[T,N] int zero_idx;
  // Max number of zero counts
  array[T] int zero_max;
  // indices of nonzero counts
  array[T,N] int nonzero_idx;
  // max number of nonzero counts
  array[T] int nonzero_max;
}
transformed data {
  // Logged population
  array[T] vector[N] log_E;
  
  for (t in 1:T) {
    log_E[t] = to_vector(log(E[t,1:N]));
  }
}
parameters {
  // Bernoulli part: Knorr-Held model
  // Intercept
  real mu_pi;
  real mu_lambda;
  
  // Structured temporal trend
  vector[T] alpha_pi;
  vector[T] alpha_lambda;
  real<lower=1e-10, upper=10> sigma_alpha_pi;
  real<lower=1e-10, upper=10> sigma_alpha_lambda;
  
  // Structured spatial pattern 
  vector[N] phi_pi;
  // vector[N] phi_lambda;
  
  // Unstructured spatial pattern
  vector[N] theta_pi;
  vector[N] theta_lambda;
  
  // Proportion of spatial/aspatial error
  real<lower=0, upper=1> rho_pi;
  // real<lower=0, upper=1> rho_lambda;
  real<lower=1e-10, upper=10> sigma_convolved_pi;
  real<lower=1e-10, upper=10> sigma_convolved_lambda;
  
  // Knorr-Held Type I spatio-temporal interaction
  array[T] vector[N] delta_pi;
  array[T] vector[N] delta_lambda;
  real<lower=1e-10, upper=10> sigma_delta_lambda;
  real<lower=1e-10, upper=10> sigma_delta_pi;
  
  
}
transformed parameters{
  array[T] vector[N] pi;   // Bernoulli GLM term 
  array[T] vector[N] lambda;  // Poisson GLM term
  
  for (t in 1:T) {
    pi[t] = inv_logit(mu_pi + 
            alpha_pi[t] +
            sigma_convolved_pi * (
              sqrt(rho_pi/scaling_factor) * phi_pi + sqrt(1-rho_pi)*theta_pi
            ) +
            sigma_delta_pi * delta_pi[t]);
    lambda[t] = exp(log_E[t] + mu_lambda + 
              alpha_lambda[t] +
              sigma_convolved_lambda * (
                // sqrt(rho_lambda/scaling_factor) * 
                // phi_lambda 
                theta_lambda
                // sqrt(1-rho_lambda)*theta_lambda
              ) +
              sigma_delta_lambda * delta_lambda[t]);
  }
}
model {  
  // Intercepts
  mu_pi ~ normal(-10, 10);
  mu_lambda ~ normal(-5, 10);
  
  // Structured temporal trend
  alpha_pi[1] ~ normal(0, sigma_alpha_pi);
  alpha_pi[2:T] ~ normal(alpha_pi[1:(T-1)], sigma_alpha_pi);
  sigma_alpha_pi ~ gamma(2, 1);
  
  alpha_lambda[1] ~ normal(0, sigma_alpha_lambda);
  alpha_lambda[2:T] ~ normal(alpha_lambda[1:(T-1)], sigma_alpha_lambda);
  sigma_alpha_lambda ~ gamma(2, 1);
  
  // Structured spatial patten
  phi_pi ~ icar_normal(N, node1, node2);
  // phi_lambda ~ icar_normal(N, node1, node2);
  
  // Unstructured spatial error
  theta_pi ~ std_normal();
  theta_lambda ~ std_normal();
  
  // Prior on Rho
  rho_pi ~ beta(.5, .5);
  // rho_lambda ~ beta(.5, .5);
  
  // Convolved variance
  sigma_convolved_pi ~ gamma(2,1);
  sigma_convolved_lambda ~ gamma(2,1);
  
  for (t in 1:T){
    // Interaction
    delta_pi[t] ~ std_normal();
    delta_lambda[t] ~ std_normal();
    
  }
  sigma_delta_pi ~ gamma(2, 1);
  sigma_delta_lambda ~ gamma(2, 1);
  
  
  
  // Likelihood
  for (t in 1:T) {
    
    // Vectorized ZIP
    // Zeros
    if (zero_max[t] > 0) {
      target += log(
             pi[t, zero_idx[t, 1:zero_max[t]]] +
          (1 - pi[t, zero_idx[t, 1:zero_max[t]]]) .*
          exp(-lambda[t, zero_idx[t, 1:zero_max[t]]])
        );
    }
    
    // Nonzeros
    if (nonzero_max[t] > 0) {
      target += bernoulli_lpmf(
                  rep_array(0, nonzero_max[t]) |
                  pi[t, nonzero_idx[t,1:nonzero_max[t]]]
            );
      target += poisson_lpmf(
                  y[t, nonzero_idx[t, 1:nonzero_max[t]]] |
                  lambda[t, nonzero_idx[t, 1:nonzero_max[t]]]
            );
    }
  }
}
