data { 
  int<lower=1> N_trials;
  int<lower=1,upper=5> w_ans[N_trials];
}
parameters {
  simplex[5] theta;
} 
model {
  theta ~ beta(2, 2);  
  for(n in 1:N_trials)
    w_ans[n] ~ categorical(theta);
}
generated quantities{
  int pred_w_ans[N_trials];
  for(n in 1:N_trials)
    pred_w_ans[n] = categorical_rng(theta);
}
