setwd('/Users/lisa/Documents/SMLP/material/day5/mpt_modeling/')
## ----setup13B, include=FALSE,cache = FALSE-------------------------------
knitr::opts_chunk$set(tidy = TRUE,cache=TRUE, autodep=TRUE,tidy.opts=list(width.cutoff=59))
options(htmltools.dir.version = FALSE)
options(scipen=999, digits=3)
library(ggplot2)
theme_set(theme_classic())

## ---- message = FALSE----------------------------------------------------
set.seed(42)
library(dplyr)
library(extraDistr)
library(rstan)
# Save compiled models:
rstan_options(auto_write = TRUE)
# Parallelize the chains using all the cores:
options(mc.cores = parallel::detectCores())  
# options(mc.cores = parallel::detectCores() - 1) # If you want to have an extra core free
library(bayesplot)

## ------------------------------------------------------------------------
true_theta <- tibble(theta_NR = .2, 
					 theta_Neologism = .1,
					 theta_Formal = .2,
					 theta_Mixed = .08,
					 theta_Correct =  1 - (theta_NR + theta_Neologism + theta_Formal + theta_Mixed))
true_theta
sum(true_theta)


## ------------------------------------------------------------------------
N_trials <- 200
(ans_mn <- rmultinom(1, N_trials, true_theta))

## ----mn, tidy = TRUE, comment="", echo=FALSE-----------------------------
cat(readLines("mn.stan"), sep = "\n")  

## ---- results = "hide", message = "FALSE"--------------------------------
# Create a list:
# c(ans_mn) makes a vector out of the "tibble"
data_mn <-  list(N_trials = N_trials,
			   ans = c(ans_mn)) 

# We fit the model with the default values of number of chains and iterations:
# chains = 4,    iter = 2000
fit_mn <- stan(file = 'mn.stan', data = data_mn)

## ------------------------------------------------------------------------
print(fit_mn, pars = c("theta"), probs = c(.0275,.5,.975))

## ---- message = FALSE----------------------------------------------------
mcmc_hist(as.array(fit_mn), regex_pars = "theta") + 
        geom_vline(data = tibble(Parameter = paste0("theta[",1:5,"]"), value = unlist(true_theta)),
                   aes(xintercept = value)) +
  scale_x_continuous(limits = c(0, 1))


## ----cat, tidy = TRUE, comment="", echo=FALSE----------------------------
cat(readLines("cat.stan"), sep = "\n")  

## ---- results = "hide", message = "FALSE"--------------------------------
ans_cat <- rmultinom(N_trials, 1, true_theta) 
head(ans_cat[,1:20])

# We need to change the answer to make it compatible to what Stan expects (a response between 1-5 that corresponds to our categories).

data_cat <-  list(N_trials = N_trials,
			   w_ans = c(t(ans_cat) %*% 1:5)) # This makes matrix into a vector of responses between 1 and 5

## ---- message = FALSE, results = "hide"----------------------------------
fit_cat <- stan(file = 'cat.stan', data = data_cat)

## ------------------------------------------------------------------------
print(fit_cat, pars = c("theta"), probs = c(.0275,.5,.975))

## ---- message = FALSE----------------------------------------------------
mcmc_hist(as.array(fit_cat), regex_pars = "theta") + 
        geom_vline(data = tibble(Parameter = paste0("theta[",1:5,"]"), value = unlist(true_theta)),
                   aes(xintercept = value)) +
  scale_x_continuous(limits = c(0, 1))


## ---- message=F, warning=F-----------------------------------------------
# Probabilities of different answers
Pr_NR <- function(a, t, f, c){
  1 - a
}  
Pr_Neologism <- function(a, t, f, c){
 a * (1 - t) * (1 - f) * (1 - c) + a * t * (1 - f) * (1 - c)
} 
Pr_Formal <- function(a, t, f, c){
 a * (1 - t) * (1 - f) * c +  a * t * (1 - f) * c
} 
Pr_Mixed <- function(a, t, f, c){
  a * (1 - t) * f
}
Pr_Correct <- function(a, t, f, c){
  a * t * f
}

# Fake true underlying values
a_true <- .75
t_true <- .9
f_true <- .8
c_true <- .1

# Prob of the different answers:
(theta_NR <- Pr_NR(a_true, t_true, f_true, c_true))
(theta_Neologism <- Pr_Neologism(a_true, t_true, f_true, c_true))
(theta_Formal <- Pr_Formal(a_true, t_true, f_true, c_true))
(theta_Mixed <- Pr_Mixed(a_true, t_true, f_true, c_true))
(theta_Correct <- Pr_Correct(a_true, t_true, f_true, c_true))

N_trials <- 100

ans <- rmultinom(1, N_trials, c(theta_NR,
	theta_Neologism,
	theta_Formal,
	theta_Mixed,
	theta_Correct
	))


## ---- tidy = TRUE, comment="", echo=FALSE--------------------------------
cat(readLines("aphasia_sMPT.stan"), sep = "\n")  

## ---- results = "hide", message = "FALSE"--------------------------------
data_sMPT <-  list(N_trials = N_trials,
			   ans = c(ans)) 
fit_sMPT <- stan(file = 'aphasia_sMPT.stan', data = data_sMPT)  

## ------------------------------------------------------------------------
print(fit_sMPT, pars=c("a","t","f","c"), probs = c(.0275,.5,.975))

## ---- message = FALSE----------------------------------------------------
mcmc_hist(as.array(fit_sMPT), pars = c("a","t","f","c")) + 
        geom_vline(data = tibble(Parameter = c("a","t","f","c"), value = c(a_true, t_true, f_true, c_true)),
                   aes(xintercept = value)) +
  scale_x_continuous(limits = c(0, 1))


## ------------------------------------------------------------------------
print(fit_sMPT, pars=c("theta"), probs = c(.0275,.5,.975))

## ------------------------------------------------------------------------
N_trials <- 1000
complexity <- rnorm(N_trials, 0, 2)

f_alpha <- .3
f_beta <- -.05

# This is where the magic happens:
f_true <- plogis(f_alpha + complexity * f_beta)

## ------------------------------------------------------------------------
theta_NR_v <- rep(Pr_NR(a_true, t_true, f_true, c_true), N_trials)  
theta_Neologism_v <- Pr_Neologism(a_true, t_true, f_true, c_true)
theta_Formal_v <- Pr_Formal(a_true, t_true, f_true, c_true)
theta_Mixed_v <- Pr_Mixed(a_true, t_true, f_true, c_true)
theta_Correct_v <- Pr_Correct(a_true, t_true, f_true, c_true)

theta_item <- matrix(
				  c(theta_NR_v,
				  theta_Neologism_v,
				  theta_Formal_v,
				  theta_Mixed_v,
				  theta_Correct_v),ncol=5)

fake_data_cx <- tibble(item = 1:N_trials,
	   complexity =  complexity,
	   w_ans = c(rmnom(N_trials, 1, theta_item) %*% 1:5))

## ----apha, tidy = FALSE, comment="", echo=FALSE--------------------------
cat(readLines("aphasia_sMPT_cat.stan"), sep = "\n")  

## ------------------------------------------------------------------------
N_item <- 20 
N_subj <- 30 
N_trials <- N_item * N_subj 
subject <- rep(1:N_subj, each = N_item)
item <- rep(1:N_item, time = N_subj)
complexity <- rnorm(N_trials, 0, 2)

exp_fake <- tibble(subject = subject,
						       item = item,
                   complexity = complexity)

sigma_a_alpha_subject <- 1.1
a_alpha_subject_tilde <- rnorm(N_subj, 0, sigma_a_alpha_subject)
a_true_h <- plogis(qlogis(a_true) + sigma_a_alpha_subject * a_alpha_subject_tilde[subject])

f_true <- plogis(f_alpha + complexity * f_beta)

theta_NR_v_h <- Pr_NR(a_true_h, t_true, f_true, c_true) 
theta_Neologism_v_h <- Pr_Neologism(a_true_h, t_true, f_true, c_true)
theta_Formal_v_h <- Pr_Formal(a_true_h, t_true, f_true, c_true)
theta_Mixed_v_h <- Pr_Mixed(a_true_h, t_true, f_true, c_true)
theta_Correct_v_h <- Pr_Correct(a_true_h, t_true, f_true, c_true)

theta_h <- matrix(
          c(theta_NR_v_h,
          theta_Neologism_v_h,
          theta_Formal_v_h,
          theta_Mixed_v_h,
          theta_Correct_v_h),ncol=5)

fake_data_h <- mutate(exp_fake,
                 w_ans = c(rmnom(N_trials, 1, theta_h) %*% 1:5))

## ----h, tidy = FALSE, comment="", echo=FALSE-----------------------------
cat(readLines("aphasia_MPT_h.stan"), sep = "\n")   

## ---- results = "hide", message = "FALSE"--------------------------------
fake_list_h <-  list(N_trials = nrow(fake_data_h),
         w_ans = fake_data_h$w_ans,
         N_subject = max(fake_data_h$subject),
         subject = fake_data_h$subject,
         complexity = fake_data_h$complexity) 

fit_h <- stan(file = 'aphasia_MPT_h.stan', data = fake_list_h, control = list(adapt_delta = .9))  

## ------------------------------------------------------------------------
print(fit_h, pars=c("t", "c", "sigma_a_alpha_subject", "a_alpha", "f_alpha", "f_beta"), probs = c(.0275,.5,.975))

## ---- message = FALSE----------------------------------------------------
mcmc_hist(as.array(fit_h), pars = c("sigma_a_alpha_subject","a_alpha", "t","f_alpha","f_beta","c")) + 
        geom_vline(data = tibble(Parameter = c("sigma_a_alpha_subject","a_alpha", "t","f_alpha","f_beta","c"), value = c(sigma_a_alpha_subject, qlogis(a_true),  t_true, f_alpha, f_beta, c_true)),
                   aes(xintercept = value)) +
  scale_x_continuous(limits = c(0, 2))


## ------------------------------------------------------------------------
gen_data <- rstan::extract(fit_h)$pred_w_ans
ppc_bars(fake_list_h$w_ans, gen_data) + ggtitle ("Hierarchical model")

## ---- fig.height=8-------------------------------------------------------
ppc_bars_grouped(fake_list_h$w_ans, gen_data, group = subject) + ggtitle ("Individual subjects in the hierarchical model") 

## ---- message = FALSE, results = "hide"----------------------------------
fit_sh <- stan(file = 'aphasia_sMPT_cat.stan', data = fake_list_h) 

## ------------------------------------------------------------------------
gen_data_sMPT <- rstan::extract(fit_sh)$pred_w_ans
ppc_bars(fake_list_h$w_ans, gen_data_sMPT) + ggtitle ("Non hierarchical model")

## ---- fig.height=8-------------------------------------------------------
ppc_bars_grouped(fake_list_h$w_ans, gen_data_sMPT, group = subject) + ggtitle ("Individual subjects in the non hierarchical model")

