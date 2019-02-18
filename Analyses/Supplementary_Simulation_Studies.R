# This file contains code to perform some supplementary simulation analyses for
# J. Bracher (2019): An INARMA(1, 1) model with Poisson marginals.

#######################################################
# Simulation studies
#######################################################

setwd("/home/johannes/Documents/underreporting/Theory_Equivalence/Proceedings_INARMA")
source("Analyses/Functions.R")

#######
### 1. Likelihood evaluation for of hidden INAR(1) model and INARMA(1, 1) model

# set true parameters of hidden INAR(1):
nu <- 2
alpha <- 0.3
q <- 0.5
lgt <- 3

# vector for which to evaluate likelihood:
vect <- c(1, 2, 1)

# get many simulations from the model (short chains, otherwise too much computation effort):
sim_hi <- sim_many_hidden_inar(nu, alpha, q, lgt = 3, n_sim = 100000)
# evaluate likelihood empirically:
emp_prob(sim_hi$Y, vect)
# evaluate analytically:
llik_hidden_inar(vect, nu, alpha, q, log = FALSE)
# works

# set true parameters of INARMA(1, 1) model:
tau <- 1
phi <- 0.3
kappa <- 0.2

# get many simulations (short chains of length 3, otherwise too much computation effort):
sim_ia <- sim_many_inarma(tau, phi, kappa, lgt = 3, n_sim = 100000)

# vector (of length 3) for which to evaluate likelihood:
vect <- c(1, 2, 4)
# evaluate likelihood empirically:
emp_prob(sim_ia$X, vect)
# evaluate analytically:
llik_inarma(vect, tau, phi, kappa, log = FALSE)
# works

#######
### 2. Fitting an INARMA(1, 1) model

# set true parameters of INARMA(1, 1) model:
tau2 <- 1
phi2 <- 0.7
kappa2 <- 0.8
lgt2 <- 10000

# simulate some (very) long chains from the model:
sim_ia2 <- sim_many_inarma(tau2, phi2, kappa2, lgt2, n_sim = 10)

# fit the model to one of these series:
fit_ia2 <- fit_inarma(sim_ia2$X[1, ])
fit_ia2$par
# seems to work.


#######
### 3. One-step-ahead forecasts from INARMA(1, 1) model

# test one-step-ahead function:
tau3 <- 1
phi3 <- 0.5
kappa3 <- 0.4

sim_ia3 <- sim_inarma(tau3, phi3, kappa3, lgt = 1000000)
# vector of observations for t - 3, t - 2, t - 1 for which to obtain one-
# step-ahead forecast
pattern <- c(1, 3, 2)
# get all observations following the defined pattern:
inds_following_pattern <- numeric(0)
for(i in 5:length(sim_ia3$X)){
  if(all(sim_ia3$X[i - (length(pattern):1)] == pattern)){
    inds_following_pattern <- c(inds_following_pattern, i)
  }
}
# take a look at their empirical distribution:
table(sim_ia3$X[inds_following_pattern])/length(inds_following_pattern)

# use analytical functions:
pred_distr_ia3 <- osa_inarma0(pattern, ind_forecast = length(pattern) + 1, tau3, phi3, kappa3)

# plot analytical forecast probabilities:
plot(0:(length(pred_distr_ia3) - 1), pred_distr_ia3, ylim = c(0, 0.5), xlim = c(0, 10))
# and add empirical ones:
lines(table(sim_ia3$X[inds_following_pattern])/length(inds_following_pattern))
# works


#######
### 4. Check some analytical expressions

# choose true parameter values:
tau4 <- 5
kappa4 <- 0.6
phi4 <- 0.3
xi4 <- 1 - phi4*(1 - kappa4)

# INGARCH:
sim4 <- sim_ingarch(tau = tau4, phi = phi4, kappa = kappa4, lgt = 100000)
# mean:
mean(sim4)
tau4/(1 - kappa4)
# variance:
var(sim4)
(1 - xi4^2 + kappa4^2*phi4^2)/(1 - xi4^2)*tau4/(1 - kappa4)
# autocorrelation function:
acf(sim4, lag.max = 5)$acf
(1 - xi4^2 + kappa4^2*phi4^2 + kappa4*phi4*(1 - phi4))/(1 - xi4^2 + kappa4^2*phi4^2)*kappa4*phi4*xi4^(0:4)
# works

# INARMA:
sim4.2 <- sim_inarma(tau = tau4, phi = phi4, kappa = kappa4, lgt = 100000)
# mean and variance:
mean(sim4.2$X)
var(sim4.2$X)
tau4/(1 - kappa4)
# autocorrelation function:
acf(sim4.2$X, lag.max = 5)$acf
kappa4*phi4*xi4^(0:4)
# works

# Verify equivalence of INARMA(1, 1) and thinned INAR(1)
params_hidden_inar4 <- reparam_inarma_to_hidden_inar(tau = tau4, phi = phi4, kappa = kappa4)
sim4.3 <- sim_hidden_inar(nu = params_hidden_inar4$nu,
                          alpha = params_hidden_inar4$alpha,
                          q = params_hidden_inar4$q, lgt = 100000)
# means:
mean(sim4.2$X)
mean(sim4.3$Y)
# variances:
var(sim4.2$X)
var(sim4.3$Y)
# acf:
(acf(sim4.2$X))
(acf(sim4.3$Y))
