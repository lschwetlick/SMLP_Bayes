data { 
  int<lower=1> N_trials;
  int<lower=1,upper=5> w_ans[N_trials];
  int item[N_trials];
  real complexity[N_trials];
}
parameters {
  real<lower=0,upper=1> a;
  real<lower=0,upper=1> t;
  real<lower=0,upper=1> c;
  
  real f_alpha;
  real f_beta;
} 
transformed parameters {
  simplex[5] theta[N_trials];
  
  for(n in 1:N_trials){
    real f = inv_logit(f_alpha  + complexity[n] * f_beta);
    
    theta[n, 1] = 1 - a; //Pr_NR
    theta[n, 2] = a * (1 - t) * (1 - f) * (1 - c) + a * t * (1 - f) * (1 - c); //Pr_Neologism
    theta[n, 3] = a * (1 - t) * (1 - f) * c +  a * t * (1 - f) * c;  //Pr_Formal
    theta[n, 4] = a * (1 - t) * f; //Pr_Mixed
    theta[n, 5] = a * t * f; //Pr_Correct
  }
  
}
model {
  a ~ beta(2, 2);  
  t ~ beta(2, 2);  
  c ~ beta(2, 2);
  
  f_alpha ~ normal(0, 2);  
  f_beta ~ normal(0, 2);  
  
  for(n in 1:N_trials)
    w_ans[n] ~ categorical(theta[n]);
}
generated quantities{
	int pred_w_ans[N_trials];
  for(n in 1:N_trials)
    pred_w_ans[n] = categorical_rng(theta[n]);
}
