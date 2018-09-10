// Excersize 1 - generate_poisson.stan

data{
  real l;
}

generated quantities{
  int x = poisson_rng(l);
}
