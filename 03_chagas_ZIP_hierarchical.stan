data {
  int<lower=0> U;                       // Number of UFs
  int<lower=0> N_max;                   // Maximum number of municipios
  int<lower=0> N[U];                    // Number of municipios within UF u 
  int<lower=0> Y;                       // Year
  int<lower=0> M;                       // Month
  int<lower=0> chagas[U, N_max, Y, M];  // Array of counts
  int<lower=0> pop[U, N_max, Y];        // Array of Population
}
parameters {
  vector<lower=0, upper=1>[U] theta;     // Probability of 0 count
  vector<lower=0>[U] lambda;             // Poisson mean P(chagas count | theta == 0)
}
model {
  theta ~ uniform(0, 1);            // Bounded by 0 and 1
  lambda ~ uniform(0, 100);         // Greater than 0
  
  for (u in 1:U) {
    // Lambda and theta within each UF
    // theta[u] ~ uniform(0, 1);
    // lambda[u] ~ uniform(0, 100);
    
    for (i in 1:N[u]) {
      for (y in 1:Y) {
        for (m in 1:M) {
          if (chagas[u, i, y, m] == 0) 
                  target += log_sum_exp(bernoulli_lpmf(1 | theta[u]),
                                  bernoulli_lpmf(0 | theta[u])
                                + poisson_lpmf(chagas[u, i, y, m] | lambda[u] * pop[u, i, y]));
          else
                  target += bernoulli_lpmf(0 | theta[u])
                             + poisson_lpmf(chagas[u, i, y, m] | lambda[u] * pop[u, i, y]);             
        }
      }
    }
  }
}