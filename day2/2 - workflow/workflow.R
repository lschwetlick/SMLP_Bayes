setwd('')

############################################################
#
# Initial Setup
#
############################################################

library(rstan)
rstan_options(auto_write = TRUE)

library(foreach)
library(doParallel)

util <- new.env()
source('stan_utility.R', local=util)

c_light <- c("#DCBCBC")
c_light_highlight <- c("#C79999")
c_mid <- c("#B97C7C")
c_mid_highlight <- c("#A25050")
c_dark <- c("#8F2727")
c_dark_highlight <- c("#7C0000")

############################################################
#
# Workflow
#
############################################################

############################################################
# PRIOR TO OBSERVATION
############################################################

############################################################
# 1. Conceptual Analysis
############################################################

# We are working with a suite of detectors that each record
# the same source over a given time interval. The source
# strength and detector response are not expected to vary
# significantly in time.  Each detector is identical and
# returns discrete counts.

############################################################
# 2. Define Observations
############################################################

# Mathematically our observation takes the form of integer
# counts, y, for each of the N detectors.  In the Stan
# modeling lanugage this is specified as

writeLines(readLines("fit_data.stan", n=4))

############################################################
# 3. Identify Relevant Summary Statistics
############################################################

# There are N components in each observation, one for each
# detector.  We could analyze each component independently,
# but because we assume that the detectors are all identical
# we can analyze their comprehensive responses with a histogram
# of their counts.  In other words we consider the histogram
# of detector counts _as the summary statistic_!

# In this conceptual example assume that our conceptual
# domain expertise informs us that 25 counts in a detector
# would be an extreme but not impossible observation.

############################################################
# 4. Build a Generative Model
############################################################

# The constant source strength and detector responds suggests
# a Poisson observation model for each of the detectors with
# a single source strength, lambda.
#
# Our domain expertise that 25 counts is extreme suggests that
# we want our prior for lambda to keep most of its probability
# mass below lambda = 15, which corresponds to fluctations in
# the observations around 15 + 3 * sqrt(15) ~ 25.
#
# We achieve this with a half-normal prior with standard
# deviation = 6.44787 such that only 1% of the prior probability
# mass is above lambda = 15.

lambda <- seq(0, 20, 0.001)

plot(lambda, dnorm(lambda, 0, 6.44787), type="l", col=c_dark_highlight, lwd=2,
     xlab="lambda", ylab="Prior Density", yaxt='n')

lambda99 <- seq(0, 15, 0.001)
dens <- dnorm(lambda99, 0, 6.44787)
lambda99 <- c(lambda99, 15, 0)
dens <- c(dens, 0, 0)

polygon(lambda99, dens, col=c_dark, border=NA)

# This generative model is implemented in the Stan programs

writeLines(readLines("generative_ensemble.stan"))
writeLines(readLines("fit_data.stan"))

############################################################
# 5. Analyze the Generative Ensemble
############################################################

R <- 1000 # 1000 draws from the Bayesian joint distribution
N <- 100

############################################################
# 5a. Analyze the Prior Predictive Distribution
############################################################

simu_data <- list("N" = N)

fit <- stan(file='generative_ensemble.stan', data=simu_data,
            iter=R, warmup=0, chains=1, refresh=R,
            seed=4838282, algorithm="Fixed_param")

simu_lambdas <- extract(fit)$lambda
simu_ys <- extract(fit)$y

# Plot aggregated summary histogram for simulated observations
B <- 40
counts <- sapply(1:R, function(r) hist(simu_ys[r,], breaks=(0:(B + 1))-0.5, plot=FALSE)$counts)
probs = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
cred <- sapply(1:(B + 1), function(b) quantile(counts[b,], probs=probs))

idx <- rep(0:B, each=2)
x <- sapply(1:length(idx), function(b) if(b %% 2 == 0) idx[b] + 0.5 else idx[b] - 0.5)
pad_cred <- do.call(cbind, lapply(idx, function(n) cred[1:9,n + 1]))

plot(1, type="n", main="Prior Predictive Distribution",
     xlim=c(-0.5, B + 0.5), xlab="y", ylim=c(0, max(cred[9,])), ylab="")

polygon(c(x, rev(x)), c(pad_cred[1,], rev(pad_cred[9,])),
        col = c_light, border = NA)
polygon(c(x, rev(x)), c(pad_cred[2,], rev(pad_cred[8,])),
        col = c_light_highlight, border = NA)
polygon(c(x, rev(x)), c(pad_cred[3,], rev(pad_cred[7,])),
        col = c_mid, border = NA)
polygon(c(x, rev(x)), c(pad_cred[4,], rev(pad_cred[6,])),
        col = c_mid_highlight, border = NA)
lines(x, pad_cred[5,], col=c_dark, lwd=2)

abline(v=25, col="white", lty=1, lw=2.5)
abline(v=25, col="black", lty=1, lw=2)

# We see a very small prior predictive probability above the
# extreme observation scale from our domain expertise

length(simu_ys[simu_ys > 25]) / length(simu_ys)

############################################################
# 5b. Fit the simulated observations and evaluate
############################################################

tryCatch({
  registerDoParallel(makeCluster(detectCores()))

  simu_list <- t(data.matrix(data.frame(simu_lambdas, simu_ys)))

  # Compile the posterior fit model
  fit_model = stan_model(file='fit_data.stan')

  ensemble_output <- foreach(simu=simu_list,
                             .combine='cbind') %dopar% {
    simu_lambda <- simu[1]
    simu_y <- simu[2:(N + 1)];

    # Fit the simulated observation
    input_data <- list("N" = N, "y" = simu_y)

    capture.output(library(rstan))
    capture.output(fit <- sampling(fit_model, data=input_data, seed=4938483))

    # Compute diagnostics
    util <- new.env()
    source('stan_utility.R', local=util)

    warning_code <- util$check_all_diagnostics(fit, quiet=TRUE)

    # Compute rank of prior draw with respect to thinned posterior draws
    sbc_rank <- sum(simu_lambda < extract(fit)$lambda[seq(1, 4000 - 8, 8)])

    # Compute posterior sensitivities
    s <- summary(fit, probs = c(), pars='lambda')$summary
    post_mean_lambda <- s[,1]
    post_sd_lambda <- s[,3]

    prior_sd_lambda <- 6.44787

    z_score <- (post_mean_lambda - simu_lambda) / post_sd_lambda
    shrinkage <- 1 - (post_sd_lambda / prior_sd_lambda)**2

    c(warning_code, sbc_rank, z_score, shrinkage)
  }
}, finally={ stopImplicitCluster() })

# Check for fit diagnostics
warning_code <- ensemble_output[1,]
if (sum(warning_code) != 0) {
  print ("Some simulated posterior fits in the generative ensemble encountered problems!")
  for (r in 1:R) {
    if (warning_code[r] != 0) {
      print(sprintf('Replication %s of %s', r, R))
      util$parse_warning_code(warning_code[r])
      print(sprintf('Simulated lambda = %s', simu_lambdas[r]))
      print(" ")
    }
  }
} else {
  print ("No posterior fits in the generative ensemble encountered problems!")
}

# Check SBC histogram
sbc_rank <- ensemble_output[2,]
sbc_hist <- hist(sbc_rank, seq(0, 500, 25) - 0.5, col=c_dark, border=c_dark_highlight)
plot(sbc_hist, main="", xlab="Prior Rank", yaxt='n', ylab="")

low <- qbinom(0.005, R, 1 / 20)
mid <- qbinom(0.5, R, 1 / 20)
high <- qbinom(0.995, R, 1 / 20)
bar_x <- c(-10, 510, 500, 510, -10, 0, -10)
bar_y <- c(high, high, mid, low, low, mid, high)

polygon(bar_x, bar_y, col=c("#DDDDDD"), border=NA)
segments(x0=0, x1=500, y0=mid, y1=mid, col=c("#999999"), lwd=2)

plot(sbc_hist, col=c_dark, border=c_dark_highlight, add=T)

# Plot posterior sensitivities
z_score <- ensemble_output[3,]
shrinkage <- ensemble_output[4,]

plot(shrinkage, z_score, col=c("#8F272720"), lwd=2, pch=16, cex=0.8,
     xlim=c(0, 1), xlab="Posterior Shrinkage", ylim=c(-5, 5), ylab="Posterior z-Score")

############################################################
# POSTERIOR TO OBSERVATION
############################################################

############################################################
# 6. Fit the observations and evaluate
############################################################

input_data <- read_rdump('workflow.data.R')
fit <- stan(file='fit_data_ppc.stan', data=input_data, seed=4938483)

# Check diagnostics
util$check_all_diagnostics(fit)

# Plot marginal posterior
params = extract(fit)

hist(params$lambda, main="", xlab="lambda", yaxt='n', ylab="",
     col=c_dark, border=c_dark_highlight)

############################################################
# 7. Analyze the Posterior Predictive Distribution
############################################################

B <- 30

obs_counts <- hist(input_data$y, breaks=(0:(B + 1))-0.5, plot=FALSE)$counts

idx <- rep(0:B, each=2)
x <- sapply(1:length(idx), function(b) if(b %% 2 == 0) idx[b] + 0.5 else idx[b] - 0.5)
pad_obs <- do.call(cbind, lapply(idx, function(n) obs_counts[n + 1]))

counts <- sapply(1:4000, function(n) hist(params$y_ppc[n,], breaks=(0:(B + 1))-0.5, plot=FALSE)$counts)
probs = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
cred <- sapply(1:(B + 1), function(b) quantile(counts[b,], probs=probs))
pad_cred <- do.call(cbind, lapply(idx, function(n) cred[1:9,n + 1]))

plot(1, type="n", main="Posterior Predictive Distribution",
     xlim=c(-0.5, B + 0.5), xlab="y",
     ylim=c(0, max(c(obs_counts, cred[9,]))), ylab="")

polygon(c(x, rev(x)), c(pad_cred[1,], rev(pad_cred[9,])),
        col = c_light, border = NA)
polygon(c(x, rev(x)), c(pad_cred[2,], rev(pad_cred[8,])),
        col = c_light_highlight, border = NA)
polygon(c(x, rev(x)), c(pad_cred[3,], rev(pad_cred[7,])),
        col = c_mid, border = NA)
polygon(c(x, rev(x)), c(pad_cred[4,], rev(pad_cred[6,])),
        col = c_mid_highlight, border = NA)
lines(x, pad_cred[5,], col=c_dark, lwd=2)

lines(x, pad_obs, col="white", lty=1, lw=2.5)
lines(x, pad_obs, col="black", lty=1, lw=2)

