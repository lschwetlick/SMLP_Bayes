setwd('')

############################################################
#
# Initial setup
#
############################################################

library(rstan)

c_light <- c("#DCBCBC")
c_light_highlight <- c("#C79999")
c_mid <- c("#B97C7C")
c_mid_highlight <- c("#A25050")
c_dark <- c("#8F2727")
c_dark_highlight <- c("#7C0000")

c_light_trans <- c("#DCBCBC80")
c_light_highlight_trans <- c("#C7999980")
c_mid_trans <- c("#B97C7C80")
c_mid_highlight_trans <- c("#A2505080")
c_dark_trans <- c("#8F272780")
c_dark_highlight_trans <- c("#7C000080")

############################################################
#
# Setting a Foundation
#
############################################################

# Consider the two intervals over the integers,
A1 <- c(4, 5, 6, 7)
A2 <- c(7, 8, 9)

# What is the complement of A1 conceptually?
# A1c <- c(1, 2, 3, 8, 9, ...)

# What is the union of A1 and A2?
A_union <- union(A1, A2)
A_union

# What is the intersection of A1 and A2?
A_inter <- intersect(A1, A2)
A_inter

# Now consider the intervals over the real numbers
# given by the lower and upper bounds
B1_min <- -2
B1_max <- 1

B2_min <- 0
B2_max <- 3

# What is the complement of B1 conceptually?
# B1c = (-\infty, -2) \cup (1, \infty)

# What are the lower and upper bounds of the union of B1 and B2?
B_union_min <- ?
B_union_max <- ?

# What are the lower and upper bounds of the intersection of B1 and B2?
B_inter_min <- ?
B_inter_max <- ?

############################################################
#
# Probabilty Mass Functions
#
############################################################

# Consider the Poisson distribution with intensity l = 5.
# The probability mass function is defined in scipy as
l = 5
dpois(2, l)

# The cumulative distribution function is similary defined as
ppois(2, l)

# Let's plot the probability mass function as a function of x
plot_poisson <- function(l) {
  p <- hist(0, breaks=0:21-0.5, plot=FALSE)
  p$counts <- dpois(0:20, l)

  par(mar = c(8, 6, 0, 0.5))
  plot(p, main="", col="white", border=c_dark_highlight,
       xlab="x", xlim=c(-0.5, 20.5), 
       ylab="Probability Mass", ylim=c(0, 0.2), yaxt='n',
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
}

plot_poisson(l)

# What is the probabilty of A1 by explicit summation?
sum(sapply(A1, function(x) dpois(x, l)))

# We can plot that, too
plot_poisson_probs <- function(A, l) {
  bin_edges <- c(A, A[length(A)] + 1) - 0.5
  p_sum <- hist(A[1], breaks=bin_edges, plot=FALSE)
  p_sum$counts <- dpois(A, 5)

  plot(p_sum, col=c_dark, border=c_dark_highlight, add=T)
}

plot_poisson(l)
plot_poisson_probs(A1, l)

# What is the probabilty of A1 by subtracting cumulative distribution functions?
poisson_prob <- function(A, l) {
  return(ppois(A[length(A)], l) - ppois(A[1] - 1, l))
}

poisson_prob(A1, l)

# What is the probabilty of the complement of A1?

# We can compute this with cumulative distibution functions
?

# or by using the sum rule for complements, P[A^c] = 1 - P[A],
?

# We can also plot it
plot_poisson(l)
plot_poisson_probs(0:(A[1] - 1), l)
plot_poisson_probs((A[length(A)] + 1):20, l)

# What is the probabilty of the union of A1 and A2?

# We can compute this using cumulative distribution functions
poisson_prob(A_union, l)

# or by using the general sum rule,
# P[A \cup B] = P[A] + P[B] - P[A \cap B]
poisson_prob(A1, l) + poisson_prob(A2, l) - poisson_prob(A_inter, l)

# In pictures
plot_poisson(l)
plot_poisson_probs(A_union, l)

# What is the probability of the intersection of A1 and A2?

# We can compute this using cumulative distribution functions
poisson_prob(A_inter, l)

# or by using the general sum rule,
# P[A \cap B] = P[A] + P[B] - P[A \cup B]
poisson_prob(A1, l) + poisson_prob(A2, l) - poisson_prob(A_union, l)

# In pictures,
plot_poisson(l)
plot_poisson_probs(A_inter, l)

# What is the mean of this Poisson distribution?
# We can compute it by explicitly by approximating the expectation
# value of the embedding funcion with a finite sum over a bounded
# interval containing most of the probability mass.

# First we have to define the embedding function, which in R
# is an explicit cast from the integers to floating point numbers
iota <- function(x) {
  return(as.double(x))
}
# R will do this casting implicitly, but it's important
# to recognize what's going on in the calculation!

# The mean is then estimated by
sum(sapply(0:100, function(x) iota(x) * dpois(x, l)))

# We can also look at how the mean estimate convergences as
# we consider longer and longer intervals
xs <- 1:100
partial_means <- sapply(xs, function(x) 
                        sum(sapply(0:x, function(y) iota(y) * dpois(y, l))))

plot(xs, partial_means, type="l", col=c_dark_highlight, lwd=2,
     xlab="x", ylab="Partial Mean",
     cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5, yaxt='n')

# In the same way we can compute the variance
sum(sapply(0:100, function(x) (iota(x) - 5)**2 * dpois(x, l)))

# What is the probability mass function for
# the pushforward distribution of y = sqrt(x)?
bin_edges <- c(-0.5, sapply(1:21, function(x) sqrt(x - 0.5)))
p <- hist(0, breaks=bin_edges, plot=FALSE)
p$density <- dpois(0:20, l)

par(mar = c(8, 6, 0, 0.5))
plot(p, main="", col=c_dark, border=c_dark_highlight,
     xlab="y = sqrt(x)", yaxt='n', ylab="Probability Mass",
     cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)

# What is the mean of the pushforward distribution?
emp_mean <- sum(sapply(0:100, function(x) sqrt(iota(x)) * dpois(x, l)))
emp_mean

# What is the variance of the pushforward distribution?
sum(sapply(0:100, function(x) (sqrt(iota(x)) - emp_mean)**2 * dpois(x, l)))

# What is the variance of the pushforward distribution
# for various choices of the intensity?
emp_mean <- sum(sapply(0:100, function(x) sqrt(iota(x)) * dpois(x, 1)))
sum(sapply(0:100, function(x) (sqrt(iota(x)) - emp_mean)**2 * dpois(x, 1)))

emp_mean <- sum(sapply(0:100, function(x) sqrt(iota(x)) * dpois(x, 2.5)))
sum(sapply(0:100, function(x) (sqrt(iota(x)) - emp_mean)**2 * dpois(x, 2.5)))

emp_mean <- sum(sapply(0:100, function(x) sqrt(iota(x)) * dpois(x, 5)))
sum(sapply(0:100, function(x) (sqrt(iota(x)) - emp_mean)**2 * dpois(x, 5)))

emp_mean <- sum(sapply(0:100, function(x) sqrt(iota(x)) * dpois(x, 10)))
sum(sapply(0:100, function(x) (sqrt(iota(x)) - emp_mean)**2 * dpois(x, 10)))

emp_mean <- sum(sapply(0:100, function(x) sqrt(iota(x)) * dpois(x, 25)))
sum(sapply(0:100, function(x) (sqrt(iota(x)) - emp_mean)**2 * dpois(x, 25)))
  
emp_mean <- sum(sapply(0:100, function(x) sqrt(iota(x)) * dpois(x, 50)))
sum(sapply(0:100, function(x) (sqrt(iota(x)) - emp_mean)**2 * dpois(x, 50)))
                    
# Notice how the variance of the pushforward distribution is essentially
# independent of the intensity.  The sqrt function is known as a
# variance stabilizing transformation for the Poisson distribution.
                      
############################################################
#
# Probability Density Functions
#
############################################################

# Consider the Gaussian distribution with intensity mu = 1 and sigma = 1.25
# The probability density function is defined in scipy as
mu <- 1
sigma <- 1.25
dnorm(2, mu, sigma);

# The cumulative distribution function is similary defined as
pnorm(2, mu, sigma);

# Let's plot the probability density function as a function of x
plot_norm <- function(mu, sigma) {
  x <- seq(-8, 8, 0.001)

  plot(x, dnorm(x, mu, sigma), type="l", col=c_dark_highlight, lwd=2,
       xlab="x", ylab="Probability Density",
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5, yaxt='n')
}

plot_norm(mu, sigma)

# What is the probability of B1?
norm_prob <- function(B_min, B_max, mu, sigma) {
    return(pnorm(B_max, mu, sigma) - pnorm(B_min, mu, sigma))
}

norm_prob(B1_min, B1_max, mu, sigma)

# We can also visualize this
plot_norm_probs <- function(mu, sigma, B_min, B_max) {
  x_int <- seq(B_min, B_max, 0.001)
  x <- c(x_int, B_max, B_min)
  y <- c(dnorm(x_int, mu, sigma), 0, 0)

  polygon(x, y, col=c_dark, border=NA)
}

plot_norm(mu, sigma)
plot_norm_probs(mu, sigma, B1_min, B1_max)

# What is the probabilty of the complement of A1?

# We can compute this with cumulative distibution functions
(1 - pnorm(B1_max, mu, sigma)) + pnorm(B1_min, mu, sigma)

# or by using the sum rule for complements, P[A^c] = 1 - P[A],
1 - norm_prob(B1_min, B1_max, mu, sigma)

# In pictures,
plot_norm(mu, sigma)
plot_norm_probs(mu, sigma, -8, B1_min)
plot_norm_probs(mu, sigma, B1_max, 8)

# What is the probabilty of the union of A1 and A2?

# We can compute this using cumulative distribution functions
norm_prob(B_union_min, B_union_max, mu, sigma)

# or by using the general sum rule,
# P[A \cup B] = P[A] + P[B] - P[A \cap B]
norm_prob(B1_min, B1_max, mu, sigma) + 
norm_prob(B2_min, B2_max, mu, sigma) - 
norm_prob(B_inter_min, B_inter_max, mu, sigma)

# In pictures
plot_norm(mu, sigma)
plot_norm_probs(mu, sigma, B_union_min, B_union_max)

# What is the probability of the intersection of A1 and A2?
                      
# We can compute this using cumulative distribution functions
?

# or by using the general sum rule,
# P[A \cap B] = P[A] + P[B] - P[A \cup B]
?

# In pictures
?

# What is the probability density function for
# the pushforward distribution of y = x^2?
y_inv <- function(y) {
  return(sqrt(y))
}

# First we have to compute the Jacobian, dx / dy (y)
abs_J <- function(y) {
  ?
}

# Then we can plot the pushforward density as
ys <- seq(0, 16, 0.001)
pdfs <- sapply(ys, function(y) dnorm(y_inv(y), mu, sigma) * abs_J(y))

plot(ys, pdfs, type="l", col=c_dark_highlight, lwd=2,
     xlab="x", ylab="Probability Density",
     cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5, yaxt='n')

############################################################
#
# Samples!
#
############################################################

# Now let's see how to approximate expectations with Monte Carlo estimators

############################################################
# Poisson Distribution
############################################################

# We can generate exact samples using R's pseudo random number generator
set.seed(8675309)

r_samples <- rpois(1000, l)

# Or we can generate exact samples using Stan's pseudo random number generator
writeLines(readLines("generate_poisson.stan"))

simu_data <- list("l" = l)

fit <- stan(file='generate_poisson.stan', data=simu_data, 
            seed=194838, algorithm="Fixed_param",
            iter=1000, warmup=0, chains=1)

stan_samples <- extract(fit)$x[]

# We can compare these two samples to each other, and the true Poisson
# probability mass function, by constructing their empirical histograms
plot_poisson(l)
hist(r_samples, breaks=0:21-0.5, 
     col=c_dark_trans, border=c_dark_highlight_trans, probability=T, add=T)
hist(stan_samples, breaks=0:21-0.5, 
     col=c_mid_trans, border=c_mid_highlight_trans, probability=T, add=T)


# To facilitate the computation of Monte Carlo estimators let's define a
# Welford accumulator that accurately computes empirical means and variances
# of a sample in a single pass
welford_summary <- function(x) {
  summary = c(0, 0)
  for (n in 1:length(x)) {
    delta <- x[n] - summary[1]
    summary[1] <- summary[1] + delta / (n + 1)
    summary[2] <- summary[2] + delta * (x[n] - summary[1])
  }
  summary[2] <- summary[2] / (length(x) - 1)
  return(summary)
}

# We can then use the Welford accumulator to compute the Monte Carlo
# estimator of a function and an estiamte of its Monte Carlo Standard Error
compute_mc_stats <- function(x) {
  summary <- welford_summary(x)
  return(c(summary[1], sqrt(summary[2] / length(x))))
}

# In order to estimate probabilities we estiamte the expecation value of
# an indicator function
indicator <- function(x, A) {
  return(ifelse(A[1] <= x & x <= A[length(A)], 1, 0))
}

# Finally we can estimate the probability of A1 as
pushforward_samples = sapply(stan_samples, function(x) indicator(x, A1))
compute_mc_stats(pushforward_samples)

# Which is consisent with the exact value,
poisson_prob(A1, l)

# Even better we can plot how this Monte Carlo estimator converges
# to the exact value as the size of the sample increases.  The bands
# in red show the variation expected from the Monte Carlo Central Limit Theorem
iter <- 2:1000
mc_stats <- sapply(iter, function(n) compute_mc_stats(pushforward_samples[0:n]))
     
plot_mc_evo <- function(iter, mc_stats, truth) {
  plot(1, type="n", main="", 
       xlab="Iteration", xlim=c(0, max(iter)),
       ylab="Monte Carlo Estimator",
       ylim=c(min(mc_stats[1,] - 3 * mc_stats[2,]), max(mc_stats[1,] + 3 * mc_stats[2,])))
  
  polygon(c(iter, rev(iter)), 
          c(mc_stats[1,] - 3 * mc_stats[2,], 
            rev(mc_stats[1,] + 3 * mc_stats[2,])),
          col = c_light_highlight, border = NA)
  polygon(c(iter, rev(iter)),
          c(mc_stats[1,] - 2 * mc_stats[2,], 
            rev(mc_stats[1,] + 2 * mc_stats[2,])),
          col = c_mid, border = NA)
  polygon(c(iter, rev(iter)), 
          c(mc_stats[1,] - 1 * mc_stats[2,], 
            rev(mc_stats[1,] + 1 * mc_stats[2,])),
          col = c_mid_highlight, border = NA)
  lines(iter, mc_stats[1,], col=c_dark, lwd=2)
  abline(h=truth, col="grey", lty="dashed", lw=2)
}

plot_mc_evo(iter, mc_stats, poisson_prob(A1, l))

# What is the Monte Carlo estimator of the probability of the complement of A1?
pushforward_samples = sapply(stan_samples, function(x) 1 - indicator(x, A1))
compute_mc_stats(pushforward_samples)

1 - poisson_prob(A1, l)

# What is the Monte Carlo estimator of the probability of the union of A1 and A2?
pushforward_samples = sapply(stan_samples, function(x) indicator(x, A_union))
compute_mc_stats(pushforward_samples)

poisson_prob(A_union, l)

# # What is the Monte Carlo estimator of probability of the intersection of A1 and A2?
pushforward_samples = sapply(stan_samples, function(x) indicator(x, A_inter))
compute_mc_stats(pushforward_samples)

poisson_prob(A_inter, l)

# What is the Monte Carlo estimator of the mean of the Poisson distribution?
pushforward_samples = sapply(stan_samples, function(x) iota(x))
iter <- 2:1000
mc_stats <- sapply(iter, function(n) compute_mc_stats(pushforward_samples[0:n]))
plot_mc_evo(iter, mc_stats, l)

# What is the Monte Carlo estimator of the variance of the Poisson distribution?
pushforward_samples = sapply(stan_samples, function(x) (iota(x) - 5)**2)
compute_mc_stats(pushforward_samples)

# What is the Monte Carlo estimator of the mean of the
# pushforward distribution of y = sqrt(x)?
pushforward_samples = sapply(stan_samples, function(x) sqrt(x))
compute_mc_stats(pushforward_samples)

sum(sapply(0:100, function(x) sqrt(iota(x)) * dpois(x, l)))

# What is the Monte Carlo estimator of the pushforward probability mass function?
bin_edges <- c(-0.5, sapply(1:21, function(x) sqrt(x - 0.5)))

ps <- hist(pushforward_samples, breaks=bin_edges, freq=T)
ps$density <- ps$counts / 1000
plot(ps, col=c_mid_trans, border=c_mid_highlight_trans, 
     main="", xlab="y = sqrt(x)", yaxt='n', ylab="Probability Mass",
     cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
           
p <- hist(0, breaks=bin_edges, plot=FALSE)
p$density <- dpois(0:20, l)
plot(p, col="000000", border=c_dark_highlight, add=T)

############################################################
# Gaussian Distribution
############################################################

# We can generate exact samples using scipy's psuedo random number generator
r_samples <- rnorm(1000, mu, sigma)

# Or we can generate exact samples using Stan's pseudo random number generator
writeLines(readLines("generate_normal.stan"))

simu_data <- list("mu" = mu, "sigma" = sigma)

fit <- stan(file='generate_normal.stan', data=simu_data, 
            seed=194838, algorithm="Fixed_param",
            iter=1000, warmup=0, chains=1)

stan_samples <- extract(fit)$x[]

# We can compare these two samples to each other, and the true Gaussian
# probability density function, by constructing their empirical histograms
plot_norm(mu, sigma)
hist(r_samples, breaks=seq(-8, 8, 0.5), 
     col=c_dark_trans, border=c_dark_highlight_trans, probability=T, add=T)
hist(stan_samples, breaks=seq(-8, 8, 0.5), 
     col=c_mid_trans, border=c_mid_highlight_trans, probability=T, add=T)

# What is the Monte Carlo estimator of the probability of B1?
indicator <- function(x, B_min, B_max) {
  return(ifelse(B_min <= x & x <= B_max, 1, 0))
}

pushforward_samples = ?
compute_mc_stats(pushforward_samples)

# Which compares nicely to the exact value
norm_prob(B1_min, B1_max, mu, sigma)

iter <- 2:1000
mc_stats <- sapply(iter, function(n) compute_mc_stats(pushforward_samples[0:n]))
plot_mc_evo(iter, mc_stats, norm_prob(B1_min, B1_max, mu, sigma))

# What is the Monte Carlo estimator of the probability of the complement of B1?
pushforward_samples = ?
compute_mc_stats(pushforward_samples)

1 - norm_prob(B1_min, B1_max, mu, sigma)

# What is the Monte Carlo estimator of the probability of the union of B1 and B2?
pushforward_samples = ?
compute_mc_stats(pushforward_samples)

norm_prob(B_union_min, B_union_max, mu, sigma)

# What is the Monte Carlo estimator of the probability of the intersection of B1 and B2?
pushforward_samples = ?
compute_mc_stats(pushforward_samples)

norm_prob(B_inter_min, B_inter_max, mu, sigma)

# What is the Monte Carlo estimator of the mean?
identity <- function(x) {
  return(x)
}

pushforward_samples = sapply(stan_samples, function(x) identity(x))
compute_mc_stats(pushforward_samples)

iter <- 2:1000
mc_stats <- sapply(iter, function(n) compute_mc_stats(pushforward_samples[0:n]))
plot_mc_evo(iter, mc_stats, mu)

# What is the Monte Carlo estimator of the variance?
pushforward_samples = ?
compute_mc_stats(pushforward_samples)

# What is the Monte Carlo estimator of the mean
# of the pushforward distribution of y = x^2?
pushforward_samples = ?
compute_mc_stats(pushforward_samples)

# What is the Monte Carlo estimator of the pushforward histogram?
ys <- seq(0, 25, 0.001)
pdfs <- sapply(ys, function(y) dnorm(y_inv(y), mu, sigma) * abs_J(y))

plot(ys, pdfs, type="l", col=c_dark_highlight, lwd=2,
     xlab="y = x^2", ylab="Probability Density",
     cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5, yaxt='n')
     
hist(pushforward_samples, breaks=seq(0, 25, 0.5), probability=T,
     col=c_mid_trans, border=c_mid_highlight_trans, add=T)
