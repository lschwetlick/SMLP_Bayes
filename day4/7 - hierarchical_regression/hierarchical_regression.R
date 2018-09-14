setwd('/Users/lisa/Documents/SMLP/material/day4/7 - hierarchical_regression/')

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

c_light_trans <- c("#DCBCBC80")
c_dark_trans <- c("#8F272780")
c_green_trans <- c("#00FF0080")

############################################################
#
# Strongly Informative Data
#
############################################################

data <- read_rdump('strong.data.R')

############################################################
# Centered parameterization
############################################################

cp_fit <- stan(file='hierarchical_cp.stan', data=data, seed=4938483,
               iter=11000, warmup=1000, refresh=11000)

# Check diagnostics
util$check_all_diagnostics(cp_fit)

# Plot marginal posteriors
cp_params = extract(cp_fit)

hist(cp_params$tau, main="", xlab="tau", yaxt='n', ylab="",
     col=c_dark, border=c_dark_highlight)

N <- data$N + 1
idx <- rep(1:N, each=2)
x <- sapply(1:length(idx), function(m) if(m %% 2 == 0) idx[m] + 0.5 else idx[m] - 0.5)

probs = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
cred <- sapply(1:data$N, function(n) quantile(cp_params$theta[,n], probs=probs))
cred <- cbind(cred, quantile(cp_params$mu, probs=probs))
pad_cred <- do.call(cbind, lapply(idx, function(n) cred[1:9,n]))

plot(1, type="n", main="", 
     xlim=c(0.5, N + 0.5), xlab="", 
     ylim=c(min(pad_cred[1,]), max(pad_cred[9,])), ylab="Marginal Posteriors")
abline(v=N-0.5, col="gray80", lwd=2, lty=3)

polygon(c(x, rev(x)), c(pad_cred[1,], rev(pad_cred[9,])),
       col = c_light, border = NA)
polygon(c(x, rev(x)), c(pad_cred[2,], rev(pad_cred[8,])),
       col = c_light_highlight, border = NA)
polygon(c(x, rev(x)), c(pad_cred[3,], rev(pad_cred[7,])),
       col = c_mid, border = NA)
polygon(c(x, rev(x)), c(pad_cred[4,], rev(pad_cred[6,])),
       col = c_mid_highlight, border = NA)
lines(x, pad_cred[5,], col=c_dark, lwd=2)

# Check PPCs
N <- data$N
idx <- rep(1:N, each=2)
x <- sapply(1:length(idx), function(m) if(m %% 2 == 0) idx[m] + 0.5 else idx[m] - 0.5)

probs = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
cred <- sapply(1:data$N, function(n) quantile(cp_params$y_ppc[,n], probs=probs))
pad_cred <- do.call(cbind, lapply(idx, function(n) cred[1:9,n]))

plot(1, type="n", main="", 
     xlim=c(0.5, N + 0.5), xlab="", 
     ylim=c(min(pad_cred[1,]), max(pad_cred[9,])), ylab="Marginal Posteriors")

polygon(c(x, rev(x)), c(pad_cred[1,], rev(pad_cred[9,])),
        col = c_light, border = NA)
polygon(c(x, rev(x)), c(pad_cred[2,], rev(pad_cred[8,])),
        col = c_light_highlight, border = NA)
polygon(c(x, rev(x)), c(pad_cred[3,], rev(pad_cred[7,])),
        col = c_mid, border = NA)
polygon(c(x, rev(x)), c(pad_cred[4,], rev(pad_cred[6,])),
        col = c_mid_highlight, border = NA)
lines(x, pad_cred[5,], col=c_dark, lwd=2)

pad_obs <- do.call(cbind, lapply(idx, function(n) data$y[n]))

lines(x, pad_obs, lwd=1.5, col="white")
lines(x, pad_obs, lwd=1.25, col="black")

############################################################
# Non-centered parameterization
############################################################

ncp_fit <- stan(file='hierarchical_ncp.stan', data=data, seed=4938483,
                iter=11000, warmup=1000, refresh=11000)

# Check diagnostics
util$check_all_diagnostics(ncp_fit)

# Let's look at some marginal posterior correlations
# to see what the non-centered posterior looks like
ncp_params = extract(ncp_fit)

plot(ncp_params$theta[,1], log(ncp_params$tau),
     col=c_dark_trans, pch=16, cex=0.8, 
     xlab="theta[1]", ylab="log(tau)")

plot(ncp_params$theta_tilde[,1], log(ncp_params$tau),
     col=c_dark_trans, pch=16, cex=0.8, 
     xlab="theta_tilde[1]", ylab="log(tau)")

plot(ncp_params$"mu", log(ncp_params$tau),
     col=c_dark_trans, pch=16, cex=0.8, 
     xlab="mu", ylab="log(tau)")

# Now let's compare performance
get_sampler_params(cp_fit, inc_warmup=FALSE)[[1]][,'stepsize__'][1]
get_sampler_params(ncp_fit, inc_warmup=FALSE)[[1]][,'stepsize__'][1]

cp_steps <- do.call(rbind, get_sampler_params(cp_fit, inc_warmup=FALSE))[,'n_leapfrog__']
hist(cp_steps, breaks=0:700-0.5, main="", 
     col=c_dark, border=c_dark_highlight, 
     xlab="Number of Leapfrog Steps", yaxt='n', ylab="")

ncp_steps <- do.call(rbind, get_sampler_params(ncp_fit, inc_warmup=FALSE))[,'n_leapfrog__']
hist(ncp_steps, breaks=0:700-0.5,
     col=c_light, border=c_light_highlight, add=T)

############################################################
#
# Weakly Informative Data
#
############################################################

data <- read_rdump('weak.data.R')

############################################################
# Centered parameterization
############################################################

cp_fit <- stan(file='hierarchical_cp.stan', data=data, seed=4938483,
               iter=11000, warmup=1000, refresh=11000)

# Check diagnostics
util$check_all_diagnostics(cp_fit)

# Divergences concentrate where the chain is
# getting stuck at small values of sigma_beta
traceplot(cp_fit, pars=c("tau"))

partition <- util$partition_div(cp_fit)
div_params <- partition[[1]]
nondiv_params <- partition[[2]]

plot(nondiv_params$"theta[1]", log(nondiv_params$tau),
     col=c_dark_trans, pch=16, cex=0.8, 
     xlab="theta[1]", ylab="log(tau)")
points(div_params$"theta[1]", log(div_params$tau),
       col=c_green_trans, pch=16, cex=0.8)

# Any better with a higher adapt_delta? A little --
# we explore deeper into the funnel but not completely.
cp_fit <- stan(file='hierarchical_cp.stan', data=data, seed=4938483,
            control=list(adapt_delta=0.99))

util$check_all_diagnostics(cp_fit)

partition <- util$partition_div(cp_fit)
div_params <- partition[[1]]
nondiv_params <- partition[[2]]

plot(nondiv_params$"theta[1]", log(nondiv_params$tau),
     col=c_dark_trans, pch=16, cex=0.8, 
     xlab="theta[1]", ylab="log(tau)")
points(div_params$"theta[1]", log(div_params$tau),
       col=c_green_trans, pch=16, cex=0.8)

cp_params <- extract(cp_fit)

############################################################
# Non-centered parameterization
############################################################

ncp_fit <- stan(file='hierarchical_ncp.stan', data=data, seed=4938483,
            iter=11000, warmup=1000, refresh=11000)

# Check diagnostics
util$check_all_diagnostics(ncp_fit)

# Plot marginal posteriors
ncp_params = extract(ncp_fit)

hist(ncp_params$tau, main="", xlab="tau", yaxt='n', ylab="",
     col=c_dark, border=c_dark_highlight)

N <- data$N + 1
idx <- rep(1:N, each=2)
x <- sapply(1:length(idx), function(m) if(m %% 2 == 0) idx[m] + 0.5 else idx[m] - 0.5)

probs = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
cred <- sapply(1:data$N, function(n) quantile(ncp_params$theta[,n], probs=probs))
cred <- cbind(cred, quantile(ncp_params$mu, probs=probs))
pad_cred <- do.call(cbind, lapply(idx, function(n) cred[1:9,n]))

plot(1, type="n", main="", 
     xlim=c(0.5, N + 0.5), xlab="", 
     ylim=c(min(pad_cred[1,]), max(pad_cred[9,])), ylab="Marginal Posteriors")
abline(v=N-0.5, col="gray80", lwd=2, lty=3)

polygon(c(x, rev(x)), c(pad_cred[1,], rev(pad_cred[9,])),
        col = c_light, border = NA)
polygon(c(x, rev(x)), c(pad_cred[2,], rev(pad_cred[8,])),
        col = c_light_highlight, border = NA)
polygon(c(x, rev(x)), c(pad_cred[3,], rev(pad_cred[7,])),
        col = c_mid, border = NA)
polygon(c(x, rev(x)), c(pad_cred[4,], rev(pad_cred[6,])),
        col = c_mid_highlight, border = NA)
lines(x, pad_cred[5,], col=c_dark, lwd=2)

# Check PPCs
N <- data$N
idx <- rep(1:N, each=2)
x <- sapply(1:length(idx), function(m) if(m %% 2 == 0) idx[m] + 0.5 else idx[m] - 0.5)

probs = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
cred <- sapply(1:data$N, function(n) quantile(ncp_params$y_ppc[,n], probs=probs))
pad_cred <- do.call(cbind, lapply(idx, function(n) cred[1:9,n]))

plot(1, type="n", main="", 
     xlim=c(0.5, N + 0.5), xlab="", 
     ylim=c(min(pad_cred[1,]), max(pad_cred[9,])), ylab="Marginal Posteriors")

polygon(c(x, rev(x)), c(pad_cred[1,], rev(pad_cred[9,])),
        col = c_light, border = NA)
polygon(c(x, rev(x)), c(pad_cred[2,], rev(pad_cred[8,])),
        col = c_light_highlight, border = NA)
polygon(c(x, rev(x)), c(pad_cred[3,], rev(pad_cred[7,])),
        col = c_mid, border = NA)
polygon(c(x, rev(x)), c(pad_cred[4,], rev(pad_cred[6,])),
        col = c_mid_highlight, border = NA)
lines(x, pad_cred[5,], col=c_dark, lwd=2)

pad_obs <- do.call(cbind, lapply(idx, function(n) data$y[n]))

lines(x, pad_obs, lwd=1.5, col="white")
lines(x, pad_obs, lwd=1.25, col="black")

# Let's compare the centered and non-centered fits
plot(ncp_params$theta[,1], log(ncp_params$tau),
     col=c_dark_trans, pch=16, cex=0.8, 
     xlab="theta[1]", ylab="log(tau)")
points(cp_params$theta[,1], log(cp_params$tau),
       col=c_light_trans, pch=16, cex=0.8)

