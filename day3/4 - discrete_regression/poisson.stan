data {
  int<lower=1> N;
  int<lower=1> M;
  matrix[M,N] X;
  int y[N];
}

parameters {
  vector[M] beta;
  real alpha;
}

model {
  beta ~ normal(0,10);
  alpha ~ normal(0,10);
  
  y ~ poisson_log(X' * beta + alpha);
}

generated quantities {
  int y_ppc[N];
  for (n in 1:N)
    y_ppc[n] = poisson_log_rng(X[1:M,n]' * beta + alpha);
}
