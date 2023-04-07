functions {
  real icar_normal_lpdf(vector phi, int[] node1, int[] node2) {
    return -0.5 * dot_self(phi[node1] - phi[node2]);
 }
}
data {
  int<lower=0> UF;                       // Number of UFs
  int<lower=0> N_max;                   // Maximum number of municipios per UF
  int<lower=0> N[UF];                    // Number of municipios within UF u 
  int<lower=0> Y;                       // Year
  int<lower=0> chagas[UF, N_max, Y];     // Array of counts
  // int<lower=0> pop[UF, N_max, Y];        // Array of Population
  vector[N_max] pop[UF,Y];
  
  // Vectors containing adjacency information 
  int n_mu_edges[UF];                    // Number of municipality edges in each UF
  int n_mu_edges_max;                   // Maximum number of municipality edges
  int<lower=1> m_node1[UF, n_mu_edges_max];  // node1[i] adjacent to node2[i]
  int<lower=1> m_node2[UF, n_mu_edges_max];  // and node1[i] < node2[i]
  
  int n_u_edges;
  int<lower=1> u_node1[n_u_edges];  // node1[i] adjacent to node2[i]
  int<lower=1> u_node2[n_u_edges];  // and node1[i] < node2[i]
}
transformed data{
  vector[N_max] log_pop[UF,Y];
  for (uf in 1:UF) {
    for (y in 1:Y) {
      log_pop[uf, y] = log(pop[uf,y]);
    }
  }
}
parameters {
  vector[N_max] r[UF,Y]; // Municipality/year specific rates
  vector[N_max] u[UF,Y]; // IAR error term
  vector[UF] mu[Y];      // UF/year specific averages
  real tau_u;        // Variances for IAR dist
  real tau_mu;       // 
  
}
model {
  // Priors on tau
  tau_u ~ uniform(0, 1);
  tau_mu ~ uniform(0, 1);
  
  mu[1] ~ icar_normal_lpdf(u_node1, u_node2);
  for (y in 2:Y) {
    mu[y] ~ icar_normal_lpdf(u_node1, u_node2);
  }
  
  // Priors on r
  for (y in 1:Y) {
    for (uf in 1:UF) {
      chagas[uf,y] ~ poisson_log(log_pop[uf,y] + mu[y] + u[uf,y]);
      u[uf,y] ~ icar_normal_lpdf(m_node1[uf], m_node1[uf]);
    }
  }
}



