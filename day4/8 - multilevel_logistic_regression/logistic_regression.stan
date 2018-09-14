data{
  int<lower=1> N;
  int<lower=1> N_location;
  int<lower=1, upper=N_location> location[N];

  int<lower=1> N_method;
  int<lower=1, upper=N_method> method[N];

  vector[N] x;
  int y[N];
}

parameters{
  real beta;
  real alpha;
}

model{
  beta ~ normal(0,5);
  alpha ~ normal(0,5);

  y ~ bernoulli_logit(beta*x+alpha);
}

generated quantities{
  real p_hat_ppc = 0;
  vector[N_location] p_hat_location_ppc = rep_vector(0, N_location);
  vector[N_method] p_hat_method_ppc = rep_vector(0, N_method);

  {
    vector[N_location] location_counts = rep_vector(0, N_location);
    vector[N_method] method_counts = rep_vector(0, N_method);
    
    for(n in 1:N){
      if (bernoulli_logit_rng(beta*x[n]+alpha)){
        p_hat_ppc = p_hat_ppc + 1;
        p_hat_location_ppc[location[n]] =  p_hat_location_ppc[location[n]] + 1;
        p_hat_method_ppc[method[n]] =  p_hat_method_ppc[method[n]] + 1;
        
      }
      location_counts[location[n]]=location_counts[location[n]]+1;
      method_counts[method[n]]=method_counts[method[n]]+1;
    
    }
    p_hat_ppc =p_hat_ppc/N;
    p_hat_location_ppc =p_hat_location_ppc ./ location_counts;
    p_hat_method_ppc =p_hat_method_ppc ./ method_counts;
    
  }
}
