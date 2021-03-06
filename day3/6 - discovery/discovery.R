############################################################
#
# Initial Setup
#
############################################################

library(rstan)
rstan_options(auto_write = TRUE)

library(foreach)
library(doParallel)

c_light <- c("#DCBCBC")
c_light_highlight <- c("#C79999")
c_mid <- c("#B97C7C")
c_mid_highlight <- c("#A25050")
c_dark <- c("#8F2727")
c_dark_highlight <- c("#7C0000")

############################################################
# Calibrating Discovery Claims
############################################################

# First we simulate an ensemble from the Bayesian joint distribution
R <- 500

fit <- stan(file='generate_data.stan',
            iter=R, warmup=0, chains=1, refresh=R,
            seed=4838282, algorithm="Fixed_param")

N <- 50
M <- 3

simu_alpha_t <- extract(fit)$alpha_treatment
simu_Xs <- extract(fit)$X
simu_ys <- extract(fit)$y
simu_treatments <- extract(fit)$treatment

# Compile the posterior fit model
fit_model = stan_model(file='logistic_regression.stan')

# Given a fit to a particular observaiton we can use posterior 
# probabilities relative to a _scale of practical equivalence_
# to decide whether or not to make a discovery claim
sope <- 0.2

input_data <- list("N" = N, "M" = M, 
                   "y" = simu_ys[12,], "X" = simu_Xs[12,,], 
                   "treatment" = simu_treatments[12,],
                   "sope" = sope)

fit <- sampling(fit_model, data=input_data, seed=4938483)
params <- extract(fit)

hist(params$alpha_treatment, breaks=seq(-1, 1, 0.1),
     main="", xlab="alpha_treatment", yaxt='n', ylab="",
     col="white", border=c_dark_highlight)

hist(params$alpha_treatment[params$alpha_treatment > sope], breaks=seq(-1, 1, 0.1),
     col=c_dark, border=c_dark_highlight, add=TRUE)

abline(v=sope, col=c_light, lty=1, lw=2)

# Let's calibrate the decision making process for P = 0.25
tryCatch({
  registerDoParallel(makeCluster(detectCores()))

  discovery_stats <- foreach(r=1:R,
                             .combine='cbind') %dopar% {
    # Truth
    true_discovery <- 0
    if (simu_alpha_t[r] >= sope)
      true_discovery <- 1
    
    # Fit the simulated observation
    input_data <- list("N" = N, "M" = M, 
                       "y" = simu_ys[r,], "X" = simu_Xs[r,,], 
                       "treatment" = simu_treatments[r,],
                       "sope" = sope)
    
    capture.output(library(rstan))
    capture.output(fit <- sampling(fit_model, data=input_data, seed=4938483))
    
    # Inferred Claim
    inferred_discovery <- 0
    s <- summary(fit, probs = c(), pars='discovery')$summary
    if (s[,1] > 0.25)
      inferred_discovery <- 1
    
    c(true_discovery, inferred_discovery)
  }
}, finally={ stopImplicitCluster() })
  
# Compute False Discovery Rate Posterior
theta <- seq(0, 1, 0.001)

false_disc <- sum(discovery_stats[2,discovery_stats[1,] == 0])
N_no_disc <- length(discovery_stats[2,discovery_stats[1,] == 0])

dens <- dbeta(theta, 1 + false_disc, 1 + N_no_disc)
plot(theta, dens, type="l", col=c_dark_highlight, lwd=2,
     xlab="False Discovery Rate", xlim=c(0, 1), yaxt='n')

# Compute True Discovery Rate Posterior
true_disc <- sum(discovery_stats[2,discovery_stats[1,] == 1])
N_disc <- length(discovery_stats[2,discovery_stats[1,] == 1])

dens <- dbeta(theta, 1 + true_disc, 1 + N_disc)
plot(theta, dens, type="l", col=c_dark_highlight, lwd=2,
     xlab="True Discovery Rate", xlim=c(0, 1), yaxt='n')

