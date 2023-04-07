functions {
  real icar_normal_lpdf(vector phi, int N, int[] node1, int[] node2) {
    return -0.5 * dot_self(phi[node1] - phi[node2]);
 }
}
data {
  int<lower=0> N;
  int<lower=0> N_edges;
  int<lower=1, upper=N> node1[N_edges];  // node1[i] adjacent to node2[i]
  int<lower=1, upper=N> node2[N_edges];  // and node1[i] < node2[i]

  int<lower=0> y[N];              // count outcomes
  // vector<lower=0>[N] x;           // coefficient
  vector<lower=0>[N] E;           // exposure
  real<lower=0> scaling_factor; // scales the variance of the spatial effects

}
transformed data {
  vector[N] log_E = log(E);
}
parameters {
  real beta0;             // intercept
  
  real<lower=0> sigma;        // overall standard deviation
  real<lower=0, upper=1> rho; // proportion unstructured vs. spatially structured variance

  vector[N] phi;         // spatial effects
  vector[N] theta;         // aspatial effects
}
transformed parameters{
  // variance of each component should be approximately equal to 1
  vector[N] convoluved_re = sqrt(1 - rho) * theta + sqrt(rho / scaling_factor) * phi;
  vector[N] psi = beta0 + convoluved_re*sigma;
}
model {
  rho ~ beta(0.5, 0.5);
  beta0 ~ normal(0.0, 1.0);
  sigma ~ normal(0.0, 1.0);
  theta ~ normal(0.0, 1.0);
  phi ~ icar_normal_lpdf(N, node1, node2);
  y ~ poisson_log(log_E + psi);
  
  // soft sum-to-zero constraint on phi
  // more efficient than mean(phi) ~ normal(0, 0.001)
  sum(phi) ~ normal(0, 0.001 * N);
  
}
