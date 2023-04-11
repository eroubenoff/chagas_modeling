functions {
  real icar_normal_lpdf(vector phi, int N, int[] node1, int[] node2) {
    return -0.5 * dot_self(phi[node1] - phi[node2]);
 }
}
data {
  int<lower=0> N;
  int<lower=0> T; // Years
  int<lower=0> N_edges;
  int<lower=1, upper=N> node1[N_edges];  // node1[i] adjacent to node2[i]
  int<lower=1, upper=N> node2[N_edges];  // and node1[i] < node2[i]

  int y[T,N];            // count outcomes
  int E[T,N];           // exposure
  real<lower=0> scaling_factor; // scales the variance of the spatial effects

}
transformed data {
  vector[N] log_E[T];
  for (t in 1:T) {
    log_E[t] = to_vector(log(E[t,1:N]));
  }
}
parameters {
  vector[T] beta0;                 // intercept
  
  real<lower=0> sigma[T];        // overall standard deviation
  real<lower=0, upper=1> rho[T]; // proportion unstructured vs. spatially structured variance

  vector[N] phi[T];         // spatial effects
  vector[N] theta[T];       // aspatial effects
}
transformed parameters{
  // variance of each component should be approximately equal to 1
  vector[N] convoluted_re[T];
  vector[N] psi[T];
  for (t in 1:T) {
    profile("convoluted_re") {
    convoluted_re[t] = sqrt(1 - rho[t]) * theta[t] + sqrt(rho[t] / scaling_factor) * phi[t];
    }
    profile("convoluted_psi"){
    psi[t] = beta0[t] + convoluted_re[t]*sigma[t];
    }
    // TODO: add the lagged AR term
  }
}
model {
  profile("beta0_prior"){
  beta0 ~ uniform(-20, 20)//normal(0.0, 1.0);
  }
  profile("sigma_prior"){
  sigma ~ normal(0.0, 1.0);
  }
  profile("rho_prior") {
  rho ~ uniform(0,1);
  }
  for (t in 1:T) {
    profile("theta_prior"){
    theta[t] ~ normal(0.0, 1.0);
    }
    profile("phi_prior"){
    phi[t] ~ icar_normal(N, node1, node2);
    }
  }
  for (t in 1:T) {
    profile("y_likelihood"){
    y[t] ~ poisson_log(log_E[t] + psi[t]);
    }
    // soft sum-to-zero constraint on phi
    // more efficient than mean(phi) ~ normal(0, 0.001)
    profile("phi_sum0"){
    sum(phi[t]) ~ normal(0, 0.001 * N);
    }
  }
}
