setwd('')

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

c_mid_trans <- c("#B97C7C88")
c_mid_highlight_trans <- c("#A2505088")
c_dark_trans <- c("#8F272788")
c_dark_highlight_trans <- c("#7C000088")

############################################################
# Fit Poisson model
############################################################

input_data <- read_rdump('discrete_regression.data.R')
fit <- stan(file='poisson.stan', data=input_data, seed=4938483)

# Check diagnostics
util$check_all_diagnostics(fit)

# Plot marginal posteriors
params = extract(fit)

par(mfrow=c(2, 2))

hist(params$beta[,1], main="", xlab="beta[1]", yaxt='n', ylab="",
     col=c_dark, border=c_dark_highlight)
hist(params$beta[,2], main="", xlab="beta[2]", yaxt='n', ylab="",
     col=c_dark, border=c_dark_highlight)
hist(params$beta[,3], main="", xlab="beta[3]", yaxt='n', ylab="",
     col=c_dark, border=c_dark_highlight)
hist(params$alpha, main="", xlab="alpha", yaxt='n', ylab="",
     col=c_dark, border=c_dark_highlight)

# Perform a posterior predictive check by plotting
# the posterior predictive distribution of various
# components against the observed data
B <- 50

idx <- rep(0:B, each=2)
x <- sapply(1:length(idx), function(b) if(b %% 2 == 0) idx[b] + 0.5 else idx[b] - 0.5)

obs_counts <- hist(input_data$y, breaks=(0:(B + 1))-0.5, plot=FALSE)$counts
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
