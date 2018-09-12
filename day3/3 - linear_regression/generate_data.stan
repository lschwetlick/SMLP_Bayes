transformed data {
  int<lower=0> M = 3;
  int<lower=0> N = 500;
  vector[M] beta = [5,-3, 2]';
  real alpha = 10;
  real sigma = 1;
}

generated quantities {
  matrix[M, N] X;
  real y[N];

  for (n in 1:N) {
    real x = uniform_rng(-1, 1);
    
    // cubic regression, these are the covariates
    X[1, n] = x;
    X[2, n] = x*x;
    X[3, n] = x*x*x;
    
    y[n] = normal_rng(X[1:M,n]' * beta + alpha, sigma);
  }
}
