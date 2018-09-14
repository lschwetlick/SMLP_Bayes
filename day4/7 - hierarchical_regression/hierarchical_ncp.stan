data{
  int<lower=0> N;
  vector[N] y;
  real<lower=0> sigma;
}

parameters{
  real mu;
  real<lower=0> tau;
  vector[N] theta_tilde;
}

transformed parameters{
  vector[N] theta = mu + tau * theta_tilde;
}

model{
  mu ~ normal(0,5);
  tau ~ normal(0,5);
  theta_tilde ~ normal(0,1);
  y ~ normal(theta, sigma);
}

generated quantities{
  vector[N] y_ppc;
  for (n in 1:N){
    y_ppc[n] = normal_rng(theta[n], sigma);
  }
}
