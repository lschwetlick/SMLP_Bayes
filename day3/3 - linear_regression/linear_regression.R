setwd('/Users/lisa/Documents/SMLP/material/day3/3 - linear_regression/')

############################################################
# Initial setup
############################################################

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores()) # parelleization

#get utility fnct/ diagnistics
util <- new.env() 
source('stan_utility.R', local=util)

c_light <- c("#DCBCBC")
c_light_highlight <- c("#C79999")
c_mid <- c("#B97C7C")
c_mid_highlight <- c("#A25050")
c_dark <- c("#8F2727")
c_dark_highlight <- c("#7C0000")

############################################################
# Create data
############################################################

fit <- stan(file='generate_data.stan', iter=1,
            chains=1, seed=194838, algorithm="Fixed_param")

N <- 500
M <- 3
X <- extract(fit)$X[1,,]
y <- extract(fit)$y[1,]

stan_rdump(c("N", "M", "X", "y"), file="linear_regression.data.R")

############################################################
# Fit initial Stan program
############################################################

input_data <- read_rdump("linear_regression.data.R")
fit <- stan(file='linear_regression1.stan', data=input_data, seed=4938483)

# Check diagnostics one by one
util$check_n_eff(fit)
util$check_rhat(fit)
util$check_div(fit)
util$check_treedepth(fit)
util$check_energy(fit)

# Or all at once
util$check_all_diagnostics(fit)

# Plot marginal posteriors
params = extract(fit)

par(mfrow=c(3, 2))

hist(params$sigma, main="", xlab="sigma", yaxt='n', ylab="",
     col=c_dark, border=c_dark_highlight)
abline(v=1, col=c_light, lty=1, lw=2)

hist(params$alpha, main="", xlab="alpha",yaxt='n', ylab="",
     col=c_dark, border=c_dark_highlight)
abline(v=10, col=c_light, lty=1, lw=2)

hist(params$beta[,1], main="", xlab="beta[1]", yaxt='n', ylab="",
     col=c_dark, border=c_dark_highlight)
abline(v=5, col=c_light, lty=1, lw=2)

hist(params$beta[,2], main="", xlab="beta[2]", yaxt='n', ylab="",
     col=c_dark, border=c_dark_highlight)
abline(v=-3, col=c_light, lty=1, lw=2)

hist(params$beta[,3], main="", xlab="beta[3]", yaxt='n', ylab="",
     col=c_dark, border=c_dark_highlight)
abline(v=2, col=c_light, lty=1, lw=2)

# Perform a posterior predictive check by plotting
# the posterior predictive distribution of various
# components against the observed data
par(mfrow=c(2, 2))

hist(params$y_ppc[,5], main="", xlab="y[5]", yaxt='n', ylab="",
     col=c_dark, border=c_dark_highlight)
abline(v=input_data$y[5], col=c_light, lty=1, lw=2)

hist(params$y_ppc[,10], main="", xlab="y[10]", yaxt='n', ylab="",
     col=c_dark, border=c_dark_highlight)
abline(v=input_data$y[10], col=c_light, lty=1, lw=2)

hist(params$y_ppc[,15], main="", xlab="y[15]", yaxt='n', ylab="",
     col=c_dark, border=c_dark_highlight)
abline(v=input_data$y[15], col=c_light, lty=1, lw=2)

hist(params$y_ppc[,20], main="", xlab="y[20]", yaxt='n', ylab="",
     col=c_dark, border=c_dark_highlight)
abline(v=input_data$y[20], col=c_light, lty=1, lw=2)

# Or consider the posterior predictive distribution
# of a particular histogram of the data
# Plot aggregated summary histogram for simulated observations
B <- 10

min_y <- min(params$y_ppc)
max_y <- max(params$y_ppc)
breaks <- seq(min_y, max_y, (max_y - min_y) / B)

idx <- rep(1:B, each=2)
xs <- sapply(1:length(idx), function(b) if(b %% 2 == 0) idx[b] + 0.5 else idx[b] - 0.5)

obs <- hist(input_data$y, breaks=breaks, plot=FALSE)$counts
pad_obs <- do.call(cbind, lapply(idx, function(n) obs[n]))

ppc <- sapply(1:4000, function(n) hist(params$y_ppc[n,], breaks=breaks, plot=FALSE)$counts)
probs = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
cred <- sapply(1:B, function(b) quantile(ppc[b,], probs=probs))
pad_cred <- do.call(cbind, lapply(idx, function(n) cred[1:9, n]))

plot(1, type="n", main="Posterior Predictive Distribution",
     xlim=c(min_y, max_y), xlab="y",
     ylim=c(0, max(c(obs, cred[9,]))), ylab="")

polygon(c(xs, rev(xs)), c(pad_cred[1,], rev(pad_cred[9,])),
        col = c_light, border = NA)
polygon(c(xs, rev(xs)), c(pad_cred[2,], rev(pad_cred[8,])),
        col = c_light_highlight, border = NA)
polygon(c(xs, rev(xs)), c(pad_cred[3,], rev(pad_cred[7,])),
        col = c_mid, border = NA)
polygon(c(xs, rev(xs)), c(pad_cred[4,], rev(pad_cred[6,])),
        col = c_mid_highlight, border = NA)
lines(xs, pad_cred[5,], col=c_dark, lwd=2)

lines(xs, pad_obs, col="white", lty=1, lw=2.5)
lines(xs, pad_obs, col="black", lty=1, lw=2)

############################################################
# Fit with vectorized Stan program
############################################################

fit <- stan(file='linear_regression2.stan', data=input_data, seed=4938483)

# Check diagnostics
util$check_all_diagnostics(fit)

# Plot marginal posteriors
params = extract(fit)

par(mfrow=c(3, 2))

hist(params$sigma, main="", xlab="sigma", yaxt='n', ylab="",
     col=c_dark, border=c_dark_highlight)
abline(v=1, col=c_light, lty=1, lw=2)

hist(params$alpha, main="", xlab="alpha", yaxt='n', ylab="",
     col=c_dark, border=c_dark_highlight)
abline(v=10, col=c_light, lty=1, lw=2)

hist(params$beta[,1], main="", xlab="beta[1]", yaxt='n', ylab="",
     col=c_dark, border=c_dark_highlight)
abline(v=5, col=c_light, lty=1, lw=2)

hist(params$beta[,2], main="", xlab="beta[2]", yaxt='n', ylab="",
     col=c_dark, border=c_dark_highlight)
abline(v=-3, col=c_light, lty=1, lw=2)

hist(params$beta[,3], main="", xlab="beta[3]", yaxt='n', ylab="",
     col=c_dark, border=c_dark_highlight)
abline(v=2, col=c_light, lty=1, lw=2)

# Perform a posterior predictive check by plotting
# the posterior predictive distribution of various
# components against the observed data
par(mfrow=c(2, 2))

hist(params$y_ppc[,5], main="", xlab="y[5]", yaxt='n', ylab="",
     col=c_dark, border=c_dark_highlight)
abline(v=input_data$y[5], col=c_light, lty=1, lw=2)

hist(params$y_ppc[,10], main="", xlab="y[10]", yaxt='n', ylab="",
     col=c_dark, border=c_dark_highlight)
abline(v=input_data$y[10], col=c_light, lty=1, lw=2)

hist(params$y_ppc[,15], main="", xlab="y[15]", yaxt='n', ylab="",
     col=c_dark, border=c_dark_highlight)
abline(v=input_data$y[15], col=c_light, lty=1, lw=2)

hist(params$y_ppc[,20], main="", xlab="y[20]", yaxt='n', ylab="",
     col=c_dark, border=c_dark_highlight)
abline(v=input_data$y[20], col=c_light, lty=1, lw=2)

# Or consider the posterior predictive distribution
# of a particular histogram of the data
# Plot aggregated summary histogram for simulated observations
B <- 10

min_y <- min(params$y_ppc)
max_y <- max(params$y_ppc)
breaks <- seq(min_y, max_y, (max_y - min_y) / B)

idx <- rep(1:B, each=2)
xs <- sapply(1:length(idx), function(b) if(b %% 2 == 0) idx[b] + 0.5 else idx[b] - 0.5)

obs <- hist(input_data$y, breaks=breaks, plot=FALSE)$counts
pad_obs <- do.call(cbind, lapply(idx, function(n) obs[n]))

ppc <- sapply(1:4000, function(n) hist(params$y_ppc[n,], breaks=breaks, plot=FALSE)$counts)
probs = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
cred <- sapply(1:B, function(b) quantile(ppc[b,], probs=probs))
pad_cred <- do.call(cbind, lapply(idx, function(n) cred[1:9, n]))

plot(1, type="n", main="Posterior Predictive Distribution",
     xlim=c(min_y, max_y), xlab="y",
     ylim=c(0, max(c(obs, cred[9,]))), ylab="")

polygon(c(xs, rev(xs)), c(pad_cred[1,], rev(pad_cred[9,])),
        col = c_light, border = NA)
polygon(c(xs, rev(xs)), c(pad_cred[2,], rev(pad_cred[8,])),
        col = c_light_highlight, border = NA)
polygon(c(xs, rev(xs)), c(pad_cred[3,], rev(pad_cred[7,])),
        col = c_mid, border = NA)
polygon(c(xs, rev(xs)), c(pad_cred[4,], rev(pad_cred[6,])),
        col = c_mid_highlight, border = NA)
lines(xs, pad_cred[5,], col=c_dark, lwd=2)

lines(xs, pad_obs, col="white", lty=1, lw=2.5)
lines(xs, pad_obs, col="black", lty=1, lw=2)

