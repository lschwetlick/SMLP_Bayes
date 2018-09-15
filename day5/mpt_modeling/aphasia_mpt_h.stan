data { 
  int<lower = 1> N_trials;
  int<lower = 1,upper = 5> w_ans[N_trials];
  real complexity[N_trials];
  int<lower = 1> N_subject;
  int<lower = 1, upper = N_subject> subject[N_trials];
}
parameters {
  real<lower = 0, upper = 1> t;
  real<lower = 0, upper = 1> c;

  real a_alpha;
  real<lower = 0> sigma_a_alpha_subj;
  vector[N_subj] t_alpha_subj_tilde;
  
  real f_alpha;
  real f_beta;


} 
transformed parameters {
  simplex[5] theta[N_trials];
  
  for (n in 1:N_trials){
    real a = inv_logit(t_alpha + sigma_a_subja_item * a_alpha_subj_tilde[item[n]]);
    real f = inv_logit(f_alpha + complexity[n] * f_beta);

    theta[n, 1] = 1 - a; //Pr_NR
    theta[n, 2] = a * (1 - t) * (1 - f) * (1 - c) + a * t * (1 - f) * (1 - c); //Pr_Neologism
    theta[n, 3] = a * (1 - t) * (1 - f) * c +  a * t * (1 - f) * c;  //Pr_Formal
    theta[n, 4] = a * (1 - t) * f; //Pr_Mixed
    theta[n, 5] = a * t * f; //Pr_Correct
  }
}
model {
  t ~ beta(2, 2);  
  c ~ beta(2, 2);  
  f_alpha ~ normal(0, 2);  
  f_beta ~ normal(0, 2);  
  a_alpha ~ normal(0, 2);  
  a_alpha_subject_tilde ~ normal(0, 1);

  for(n in 1:N_trials)
    w_ans[n] ~ categorical(theta[n]);
}
generated quantities{
  int<lower = 1,upper = 5> pred_w_ans[N_trials];
  for(n in 1:N_trials) 
    pred_w_ans[n] = categorical_rng(theta[n]);
}
