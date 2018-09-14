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
  real mu_alpha;
  real<lower=0> sigma_alpha;
  vector[N_location] alpha_location_tilde;
  
  real mu_beta;
  real<lower=0> sigma_beta;
  vector[N_location] beta_location_tilde;
}

transformed parameters{
  vector[N_location] alpha_location = mu_alpha + sigma_alpha * alpha_location_tilde;
  vector[N_location] beta_location = mu_beta + sigma_beta * beta_location_tilde;
}

model{
  mu_alpha~normal(0,2.5);
  sigma_alpha~normal(0,2.5);
  alpha_location_tilde~normal(0,1);
  
  mu_beta~normal(0,2.5);
  sigma_beta~normal(0,2.5);
  beta_location_tilde~normal(0,1);
  
  y~bernoulli_logit(beta_location[location] .* x+alpha_location[location]);
}

generated quantities{
  real p_hat_ppc = 0;
  vector[N_location] p_hat_location_ppc= rep_vector(0, N_location);
  vector[N_method] p_hat_method_ppc= rep_vector(0, N_method);
  
  {
    vector[N_location] location_counts= rep_vector(0, N_location);
    vector[N_method] method_counts= rep_vector(0, N_method);
    
    vector[N] alpha_indiv = alpha_location[location];
    vector[N] beta_indiv = beta_location[location];
    
    for (n in 1:N){
      real logit_theta = beta_indiv[n]* x[n] + alpha_indiv[n];
      if (bernoulli_logit_rng(logit_theta)){
        p_hat_ppc = p_hat_ppc + 1;
        p_hat_location_ppc[location[n]] = p_hat_location_ppc[location[n]] + 1;
        p_hat_method_ppc[method[n]] = p_hat_method_ppc[method[n]] + 1;
      }
      location_counts[location[n]] = location_counts[location[n]] +1;
      method_counts[method[n]] = method_counts[method[n]] + 1;
    }
    p_hat_ppc = p_hat_ppc / N;
    p_hat_location_ppc = p_hat_location_ppc ./ location_counts;
    p_hat_method_ppc = p_hat_method_ppc ./ method_counts;
  }
}
