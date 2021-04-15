// This is a simple ZIP model for Chagas.
// This model assumes no hierarchical structure of the data,
// and there is only one of each parameter theta and lambda.
// Adapted from: https://mc-stan.org/docs/2_20/stan-users-guide/zero-inflated-section.html

// Data
data {
  int<lower=0> U;                       // Number of UFs
  int<lower=0> N_max;                   // Maximum number of municipios
  int<lower=0> N[U];                    // Number of municipios within UF u 
  int<lower=0> Y;                       // Year
  int<lower=0> M;                       // Month
  int<lower=0> chagas[U, N_max, Y, M];  // Array of counts
  int<lower=0> pop[U, N_max, Y];        // Array of Population
}

// Parameters
parameters {
  real<lower=0, upper=1> theta;     // Probability of 0 count
  real<lower=0> lambda;             // Poisson mean P(chagas count | theta == 0)
}

// Model code
model {
  theta ~ uniform(0, 1);            // Bounded by 0 and 1
  lambda ~ uniform(0, 100);         // Greater than 0
  
  for (u in 1:U) {
    for (i in 1:N[u]) {
      for (y in 1:Y) {
        for (m in 1:M) {
          
          // The "Zero-Inflated" part
          if (chagas[u, i, y, m] == 0) 
                  target += log_sum_exp(bernoulli_lpmf(1 | theta),
                                  bernoulli_lpmf(0 | theta)
                                + poisson_lpmf(chagas[u, i, y, m] | lambda * pop[u, i, y]));
          // The "Poisson" part
          else
                  target += bernoulli_lpmf(0 | theta)
                             + poisson_lpmf(chagas[u, i, y, m] | lambda * pop[u, i, y]);             
        }
      }
    }
  }
}
