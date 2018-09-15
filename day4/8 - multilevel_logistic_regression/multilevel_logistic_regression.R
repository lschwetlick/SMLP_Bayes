setwd('/Users/lisa/Documents/SMLP/material/day4/8 - multilevel_logistic_regression/')

############################################################
# Initial setup
############################################################

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

util <- new.env()
source('stan_utility.R', local=util)
source('multilevel_utility.R', local=util)

c_light <- c("#DCBCBC")
c_light_highlight <- c("#C79999")
c_mid <- c("#B97C7C")
c_mid_highlight <- c("#A25050")
c_dark <- c("#8F2727")
c_dark_highlight <- c("#7C0000")

############################################################
#    A simple logistic regression
############################################################

data <- read_rdump('multilevel_logistic_regression.data.R')

fit <- stan(file='logistic_regression.stan', data=data, seed=4938483)

# Check diagnostics
util$check_all_diagnostics(fit)

# Now let's look at the PPCs
params = extract(fit)

# Aggregate fit looks okay, but...
util$plot_ppc(params, data)

# Individual groups are very different from aggregrate response
util$plot_group_ppc(params, data, "location", 5)

util$plot_group_ppc(params, data, "method", 4)

############################################################
# Fit a grouped logistic regression
############################################################

fit <- stan(file='grouped_logistic_regression.stan', data=data, seed=4938483)

# Check diagnostics
util$check_all_diagnostics(fit)

# Let's look at some marginal posteriors
grouped_params = extract(fit)

util$plot_group_params(grouped_params, "alpha_location", 5)

util$plot_group_params(grouped_params, "beta_location", 5)

# Now let's look at the PPCs

# Aggregate fit looks okay
util$plot_ppc(grouped_params, data)

# And now the heterogeneity in the located groups looks like
# it's being captured well, at the cost of large uncertainties
util$plot_group_ppc(grouped_params, data, "location", 5)

# Keep in mind that because we haven't yet modeled the herterogeneity
# in the the method groups those PPCs still look bad
util$plot_group_ppc(grouped_params, data, "method", 4)

############################################################
#    One-level, centered parameterization
############################################################

fit <- stan(file='one_level_cp.stan', data=data, seed=4938483)

# Check diagnostics
util$check_all_diagnostics(fit)

# Divergences concentrate where the chain is
# getting stuck at small values of sigma_beta
traceplot(fit, pars=c("sigma_beta"))

partition <- util$partition_div(fit)
div_params <- partition[[1]]
nondiv_params <- partition[[2]]

c_dark_trans <- c("#8F272780")
c_green_trans <- c("#00FF0080")

plot(nondiv_params$"beta_location[1]", log(nondiv_params$sigma_beta),
     col=c_dark_trans, pch=16, cex=0.8, 
     xlab="beta_location[1]", xlim=c(-11, 1), 
     ylab="log(sigma_beta)", ylim=c(-3, 2))
points(div_params$"beta_location[1]", log(div_params$sigma_beta),
       col=c_green_trans, pch=16, cex=0.8)

# Any better with a higher adapt_delta? A little --
# we explore deeper into the funnel but not completely.
fit <- stan(file='one_level_cp.stan', data=data, seed=4938483,
            control=list(adapt_delta=0.99))

util$check_all_diagnostics(fit)

partition <- util$partition_div(fit)
div_params <- partition[[1]]
nondiv_params <- partition[[2]]

plot(nondiv_params$"beta_location[1]", log(nondiv_params$sigma_beta),
     col=c_dark_trans, pch=16, cex=0.8,     
     xlab="beta_location[1]", xlim=c(-11, 1), 
     ylab="log(sigma_beta)", ylim=c(-3, 2))
points(div_params$"beta_location[1]", log(div_params$sigma_beta),
       col=c_green_trans, pch=16, cex=0.8)

############################################################
#    One-level, non-centered parameterization
############################################################

fit <- stan(file='one_level_ncp.stan', data=data, seed=4938483,
            control=list(adapt_delta=0.85))

# Check diagnostics
util$check_all_diagnostics(fit)

# Let's look at some marginal posteriors
hier_params = extract(fit)

util$plot_group_params_comp(grouped_params, hier_params, "alpha_location", 5)

util$plot_group_params_comp(grouped_params, hier_params, "beta_location", 5)

# Now let's look at the PPCs

# Aggregate fit looks good, with slightly smaller
# variances than in the no pooling model
util$plot_ppc_comp(grouped_params, hier_params, data)

# Same story for the location groups
util$plot_group_ppc_comp(grouped_params, hier_params, data, "location", 5)

# Alas, we still haven't modeled the herterogeneity in the method groups
util$plot_group_ppc_comp(grouped_params, hier_params, data, "method", 4)

############################################################
#    Two-level, non-centered parameterization
############################################################

fit <- stan(file='two_level_ncp.stan', data=data, seed=4938483,
            control=list(adapt_delta=0.85))

# Check diagnostics
util$check_all_diagnostics(fit)

# Now let's look at the PPCs
hier_params = extract(fit)

# Aggregate PPCs looks good
util$plot_ppc_comp(grouped_params, hier_params, data)

# As do those grouped by location,
util$plot_group_ppc_comp(grouped_params, hier_params, data, "location", 5)

# And now those grouped by method, too!
util$plot_group_ppc_comp(grouped_params, hier_params, data, "location", 5)
