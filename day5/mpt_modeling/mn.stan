data { 
  int<lower=1> N_trials;
  int<lower=0,upper=N_trials> ans[5];
}
parameters {
  simplex[5] theta;
} 
model {
  theta[5] ~ beta(15, 10);  
  theta[1:4] ~ beta(2, 2);  


  ans ~ multinomial(theta);
}
generated quantities{
	int pred_ans[5];
  pred_ans = multinomial_rng(theta, 5);
}
