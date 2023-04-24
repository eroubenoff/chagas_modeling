functions {
  real icar_normal_lpdf(vector phi, int N, array[] int node1, array[] int node2) {
    // Soft sum-to-zero constraint
    return -0.5 * dot_self(phi[node1] - phi[node2]) + normal_lpdf(sum(phi) | 0, 0.001*N);
 }
  real knorr_held_type4_lpdf(vector delta_t, vector delta_tm1, int N, array[] int node1, array[] int node2) {
    // Soft sum-to-zero constraint
    return -0.5 * dot_self(delta_t[node1] - delta_t[node2] - delta_tm1[node1] + delta_tm1[node2]) + 
               normal_lpdf(sum(delta_t) | 0, 0.001*N);
    
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
  
  // Structured temporal trend
  vector[T] alpha_pi;
  real<lower=1e-10, upper=10> sigma_alpha_pi;
  
  // Structured spatial pattern 
  vector[N] phi_pi;
  real<lower=1e-10, upper=10> sigma_phi_pi;
  
  // Unstructured spatial pattern
  vector[N] theta_pi;
  real<lower=1e-10, upper = 10> sigma_theta_pi;
  
  // Knorr-Held Type 4 spatio-temporal interaction
  array[T] vector[N] delta_pi;
  real<lower=1e-10, upper=10> sigma_delta_pi; 
  
  // Poisson part: global mean and s-t random effect
  // Intercepts
  real mu_lambda;
  
  // Unstructured error terms
  array[T] vector[N] u;
  real<lower=1e-10, upper = 10> sigma_u;
  
  
}
transformed parameters{
  array[T] vector[N] pi;   // Bernoulli GLM term 
  array[T] vector[N] lambda;  // Poisson GLM term
  
  profile("transformed_params"){
  for (t in 1:T) {
    pi[t] = mu_pi + 
            sigma_alpha_pi * alpha_pi[t] + 
            sigma_phi_pi * phi_pi + 
            sigma_theta_pi * theta_pi + 
            sigma_delta_pi * delta_pi[t];
    pi[t] = inv_logit(pi[t]);
    lambda[t] = mu_lambda + sigma_u*u[t];
    lambda[t] = exp(log_E[t] + lambda[t]);
  }
  }
}
model {  
  // Inverse terms
  vector[N] pi_inv_logit; 
  vector[N] lambda_exp;
  
  // Intercepts
  profile("mu"){
  mu_pi ~ normal(-10, 10);
  mu_lambda ~ normal(-5, 10);
  }
  
  // Structured temporal trend
  // alpha_pi ~ normal(0, 1);
  alpha_pi[1] ~ normal(0, sigma_alpha_pi);
  alpha_pi[2:T] ~ normal(alpha_pi[1:(T-1)], sigma_alpha_pi);
  // // alpha_pi[2:T] ~ normal(0, sigma_alpha_pi);
  sigma_alpha_pi ~ gamma(2, 1);
  
  // Structured spatial patten
  phi_pi ~ icar_normal(N, node1, node2);
  sigma_phi_pi ~ gamma(2, 1);
  
  // Unstructured spatial error
  theta_pi ~ std_normal();
  sigma_theta_pi ~ gamma(2,1);
  
  for (t in 1:T){
    // Lambda: Unstructured heterogeneity
    profile("u"){
    u[t] ~ std_normal();
    }
    
    // Interaction error Type 4
    // if (t == 1) {
    //   delta_pi[1] ~  icar_normal(N, node1, node2);
    // }
    if (t >= 2) {
      delta_pi[t] ~ knorr_held_type4(delta_pi[t-1], N, node1, node2);
    }
  }
  sigma_delta_pi ~ gamma(2, 1);
  sigma_u ~ gamma(2, 1);
  
  
  // Likelihood
  for (t in 1:T) {
    
    // Vectorized ZIP
    
    // Zeros
    if (zero_max[t] > 0) {
      target += sum(log(
             pi[t, zero_idx[t, 1:zero_max[t]]] +
          (1 - pi[t, zero_idx[t, 1:zero_max[t]]]) .*
          exp(-lambda[t, zero_idx[t, 1:zero_max[t]]])
        ));
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
generated quantities {
  profile("generated") {
  array[T] vector[N] E_y; // Expected count 
  
  for (t in 1:T) {
    E_y[t] = (1-inv_logit( pi[t])) .* exp(log_E[t] + lambda[t]);
  }
  }
  
}
