data {
  int<lower=1> N; // Number of data
  int<lower=1> M; // Number of covariates // verifies that you have a legit number of observations
  matrix[M, N] X;
  real y[N];
}

parameters {
  vector[M] beta;
  real alpha;
  real<lower=0> sigma;
}


model {
  vector[N] mu = X' * beta + alpha;
  
  for (m in 1:M)
    beta[m] ~ normal(0, 10);
  
  alpha ~ normal(0, 10);
  sigma ~ normal(0, 5);

  for (n in 1:N)
    y[n] ~ normal(mu[n], sigma);

}


generated quantities{
  real y_ppc[N];
  {
    vector[N] mu = X' * beta + alpha;
    for (n in 1:N)
      y_ppc[n] = normal_rng(mu[n], sigma);
  }
}
