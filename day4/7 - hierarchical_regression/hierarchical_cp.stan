data{
  int<lower=0> N;
  vector[N] y;
  real<lower=0> sigma;
}

parameters{
  real mu;
  real<lower=0> tau;
  vector[N] theta;
}

model{
  mu ~ normal(0,5);
  tau ~ normal(0,5);
  theta ~ normal(mu, tau);
  y ~ normal(theta, sigma);
}

generated quantities{
  vector[N] y_ppc;
  for (n in 1:N){
    y_ppc[n] = normal_rng(theta[n], sigma);
  }
}
