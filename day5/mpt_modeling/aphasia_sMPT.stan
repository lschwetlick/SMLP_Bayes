data { 
  int<lower=1> N_trials;
  int<lower=0,upper=N_trials> ans[5];
}
parameters {
  real<lower=0,upper=1> a;
  real<lower=0,upper=1> t;
  real<lower=0,upper=1> f;
  real<lower=0,upper=1> c;
} 
transformed parameters {
  simplex[5] theta;
  theta[1] = 1 - a; //Pr_NR
  theta[2] = a * (1 - t) * (1 - f) * (1 - c) + a * t * (1 - f) * (1 - c); //Pr_Neologism
  theta[3] = a * (1 - t) * (1 - f) * c +  a * t * (1 - f) * c;  //Pr_Formal
  theta[4] = a * (1 - t) * f; //Pr_Mixed
  theta[5] = a * t * f; //Pr_Correct
}
model {
  a ~ beta(5, 2);  
  t ~ beta(2, 2);  
  f ~ beta(2, 2);  
  c ~ beta(2, 2);  
  ans ~ multinomial(theta);
}
generated quantities{
	int pred_ans[5];
  pred_ans = multinomial_rng(theta, 5);
}
