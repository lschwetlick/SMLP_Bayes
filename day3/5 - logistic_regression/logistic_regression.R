setwd('/Users/lisa/Documents/SMLP/material/day3/5 - logistic_regression')

############################################################
# Initial setup
############################################################

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

util <- new.env()
source('stan_utility.R', local=util)

c_light <- c("#DCBCBC")
c_light_highlight <- c("#C79999")
c_mid <- c("#B97C7C")
c_mid_highlight <- c("#A25050")
c_dark <- c("#8F2727")
c_dark_highlight <- c("#7C0000")

############################################################
# Ungrouped logistic regression 1
############################################################

input_data <- read_rdump('logistic_regression.data.R')
fit <- stan(file='logistic_regression1.stan', data=input_data, seed=4938483)

# Check diagnostics
util$check_all_diagnostics(fit)

# Plot marginal posteriors
params = extract(fit)

par(mfrow=c(2, 3))

hist(params$beta[,1], main="", xlab="beta[1]", yaxt='n', ylab="",
     col=c_dark, border=c_dark_highlight)

hist(params$beta[,2], main="", xlab="beta[2]", yaxt='n', ylab="",
     col=c_dark, border=c_dark_highlight)

hist(params$beta[,3], main="", xlab="beta[3]", yaxt='n', ylab="",
     col=c_dark, border=c_dark_highlight)

hist(params$beta[,4], main="", xlab="beta[4]", yaxt='n', ylab="",
     col=c_dark, border=c_dark_highlight)

hist(params$beta[,5], main="", xlab="beta[5]", yaxt='n', ylab="",
     col=c_dark, border=c_dark_highlight)

hist(params$alpha, main="", xlab="alpha", yaxt='n', ylab="",
     col=c_dark, border=c_dark_highlight)

# Now let's look at the PPCs
par(mfrow=c(1, 3))

# Aggegrate
breaks_delta <- 1.0 / length(input_data$y)
breaks <- seq(- 0.5 * breaks_delta, 1 + 0.5 * breaks_delta, breaks_delta)

hist(params$p_hat_ppc, breaks=breaks, main="", xlab="p_hat_ppc", yaxt='n', ylab="",
     col=c_dark, border=c_dark_highlight)
abline(v=sum(input_data$y) / input_data$N, col=c_light, lty=1, lw=2)

# Left-Handed Group
breaks_delta <- 1.0 / length(input_data$y[input_data$h == 1])
breaks <- seq(- 0.5 * breaks_delta, 1 + 0.5 * breaks_delta, breaks_delta)

hist(params$p_hat_left_ppc, breaks=breaks, 
     col=c_dark, border=c_dark_highlight, 
     main="", xlab="p_hat_left_ppc")
abline(v=sum(input_data$y * (input_data$h == 1)) / sum(input_data$h == 1),
       col=c_light, lty=1, lw=2)
    
# Right-handed Group   
breaks_delta <- 1.0 / length(input_data$y[input_data$h == 2])
breaks <- seq(- 0.5 * breaks_delta, 1 + 0.5 * breaks_delta, breaks_delta)

hist(params$p_hat_right_ppc, breaks=breaks, 
     col=c_dark, border=c_dark_highlight, 
     main="", xlab="p_hat_right_ppc")
abline(v=sum(input_data$y * (input_data$h == 2)) / sum(input_data$h == 2),
       col=c_light, lty=1, lw=2)

# Aggregate fit looks okay, but individual groups
# are very different from aggregrate response

############################################################
# Ungrouped logistic regression 2
############################################################

# Now let's run again with the more efficient implentation
# of the logistic regression and confirm that we get an
# equivalent fit

fit <- stan(file='logistic_regression2.stan', data=input_data, seed=4938483)

# Check diagnostics
util$check_all_diagnostics(fit)

# Plot marginal posteriors
params = extract(fit)

par(mfrow=c(2, 3))

hist(params$beta[,1], main="", xlab="beta[1]", yaxt='n', ylab="",
     col=c_dark, border=c_dark_highlight)

hist(params$beta[,2], main="", xlab="beta[2]", yaxt='n', ylab="",
     col=c_dark, border=c_dark_highlight)

hist(params$beta[,3], main="", xlab="beta[3]", yaxt='n', ylab="",
     col=c_dark, border=c_dark_highlight)

hist(params$beta[,4], main="", xlab="beta[4]", yaxt='n', ylab="",
     col=c_dark, border=c_dark_highlight)

hist(params$beta[,5], main="", xlab="beta[5]", yaxt='n', ylab="",
     col=c_dark, border=c_dark_highlight)

hist(params$alpha, main="", xlab="alpha", yaxt='n', ylab="",
     col=c_dark, border=c_dark_highlight)

# Now let's look at the PPCs
par(mfrow=c(1, 3))

# Aggegrate
breaks_delta <- 1.0 / length(input_data$y)
breaks <- seq(- 0.5 * breaks_delta, 1 + 0.5 * breaks_delta, breaks_delta)

hist(params$p_hat_ppc, breaks=breaks, main="", xlab="p_hat_ppc", yaxt='n', ylab="",
     col=c_dark, border=c_dark_highlight)
abline(v=sum(input_data$y) / input_data$N, col=c_light, lty=1, lw=2)

# Left-Handed Group
breaks_delta <- 1.0 / length(input_data$y[input_data$h == 1])
breaks <- seq(- 0.5 * breaks_delta, 1 + 0.5 * breaks_delta, breaks_delta)

hist(params$p_hat_left_ppc, breaks=breaks, 
     col=c_dark, border=c_dark_highlight, 
     main="", xlab="p_hat_left_ppc")
abline(v=sum(input_data$y * (input_data$h == 1)) / sum(input_data$h == 1),
       col=c_light, lty=1, lw=2)

# Right-handed Group   
breaks_delta <- 1.0 / length(input_data$y[input_data$h == 2])
breaks <- seq(- 0.5 * breaks_delta, 1 + 0.5 * breaks_delta, breaks_delta)

hist(params$p_hat_right_ppc, breaks=breaks, 
     col=c_dark, border=c_dark_highlight, 
     main="", xlab="p_hat_right_ppc")
abline(v=sum(input_data$y * (input_data$h == 2)) / sum(input_data$h == 2),
       col=c_light, lty=1, lw=2)

############################################################
# Grouped logistic regression
############################################################

fit <- stan(file='grouped_logistic_regression.stan', data=input_data, seed=4938483)

# Check diagnostics
util$check_all_diagnostics(fit)

# Plot marginal posteriors
params = extract(fit)

par(mfrow=c(2, 4))

hist(params$beta[,1], main="", xlab="beta[1]", yaxt='n', ylab="",
     col=c_dark, border=c_dark_highlight)

hist(params$beta[,2], main="", xlab="beta[2]", yaxt='n', ylab="",
     col=c_dark, border=c_dark_highlight)

hist(params$beta[,3], main="", xlab="beta[3]", yaxt='n', ylab="",
     col=c_dark, border=c_dark_highlight)

hist(params$beta[,4], main="", xlab="beta[4]", yaxt='n', ylab="",
     col=c_dark, border=c_dark_highlight)

hist(params$beta[,5], main="", xlab="beta[5]", yaxt='n', ylab="",
     col=c_dark, border=c_dark_highlight)

hist(params$alpha[,1], main="", xlab="alpha[1]", yaxt='n', ylab="",
     col=c_dark, border=c_dark_highlight)

hist(params$alpha[,2], main="", xlab="alpha[2]", yaxt='n', ylab="",
     col=c_dark, border=c_dark_highlight)

# Now let's look at the PPCs
par(mfrow=c(1, 3))

# Aggegrate
breaks_delta <- 1.0 / length(input_data$y)
breaks <- seq(- 0.5 * breaks_delta, 1 + 0.5 * breaks_delta, breaks_delta)

hist(params$p_hat_ppc, breaks=breaks, main="", xlab="p_hat_ppc", yaxt='n', ylab="",
     col=c_dark, border=c_dark_highlight)
abline(v=sum(input_data$y) / input_data$N, col=c_light, lty=1, lw=2)

# Left-Handed Group
breaks_delta <- 1.0 / length(input_data$y[input_data$h == 1])
breaks <- seq(- 0.5 * breaks_delta, 1 + 0.5 * breaks_delta, breaks_delta)

hist(params$p_hat_left_ppc, breaks=breaks, 
     col=c_dark, border=c_dark_highlight, 
     main="", xlab="p_hat_left_ppc")
abline(v=sum(input_data$y * (input_data$h == 1)) / sum(input_data$h == 1),
       col=c_light, lty=1, lw=2)

# Right-handed Group   
breaks_delta <- 1.0 / length(input_data$y[input_data$h == 2])
breaks <- seq(- 0.5 * breaks_delta, 1 + 0.5 * breaks_delta, breaks_delta)

hist(params$p_hat_right_ppc, breaks=breaks, 
     col=c_dark, border=c_dark_highlight, 
     main="", xlab="p_hat_right_ppc")
abline(v=sum(input_data$y * (input_data$h == 2)) / sum(input_data$h == 2),
       col=c_light, lty=1, lw=2)

