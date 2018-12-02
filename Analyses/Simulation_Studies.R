#######################################################
# Simulation studies
#######################################################

setwd("/home/johannes/Documents/underreporting/Theory_Equivalence/Proceedings_INARMA")
source("Analyses/Functions.R")

#######
### 1. valuation of likelihood in hidden INAR(1)

# set true parameters of hidden INAR(1):
nu <- 2
alpha <- 0.3
q <- 0.5
lgt <- 3

# vector for which to evaluate likelihood:
vect <- c(1, 2, 1)

# get many simulations:
sim_hi <- sim_many_hidden_inar(nu, alpha, q, lgt = 3, n_sim = 100000)
# evaluate likelihood empirically:
emp_prob(sim_hi$Y, vect)
# evaluate analytically:
lik_hidden_inar(vect, nu, alpha, q, log = FALSE)
# works

### Evaluation of likelihood in INARMA(1, 1) model:
tau <- 1
phi <- 0.3
kappa <- 0.2

# vector for which to evaluate likelihood:
vect <- c(1, 2, 2)

# get many simulations:
sim_ia <- sim_many_inarma(tau, phi, kappa, lgt = 3, n_sim = 100000)
# evaluate likelihood empirically:
emp_prob(sim_ia$X, vect)
# evaluate analytically:
lik_inarma(vect, tau, phi, kappa, log = FALSE)
# works

#######
### 2. fitting an INARMA(1, 1) model
tau2 <- 1
phi2 <- 0.7
kappa2 <- 0.7
lgt2 <- 10000

sim_ia2 <- sim_many_inarma(tau2, phi2, kappa2, lgt2, n_sim = 10)

fit_lik_inarma(sim_ia2$X[1, ])
# seems to work.


#######
### 3. One-step-ahead forecasts from INARMA

# test one-step-ahead function:
tau <- 1
phi <- 0.5
kappa <- 0.4

sim <- sim_inarma(tau, phi, kappa, lgt = 1000000)
pattern <- c(1, 3, 1)
inds_following_pattern <- numeric(0)
for(i in 5:length(sim$X)){
  if(all(sim$X[i - (length(pattern):1)] == pattern)) inds_following_pattern <- c(inds_following_pattern, i)
}
table(sim$X[inds_following_pattern])

pred_distr <- osa_inarma0(pattern, ind_forecast = length(pattern) + 1, tau, phi, kappa)
plot(0:(length(pred_distr) - 1), pred_distr, ylim = c(0, 0.5), xlim = c(0, 10))
lines(table(sim$X[inds_following_pattern])/length(inds_following_pattern))
# works


#######
### 4. Check re-parameterized formulas for moments
tau <- 5
kappa <- 0.6
phi <- 0.3
xi <- 1 - phi + phi*kappa

# INGARH:
sim1 <- sim_ingarch(tau = tau, phi = phi, kappa = kappa, lgt = 100000)
mean(sim1)
tau/(1 - kappa)
var(sim1)
(1 - xi^2 + kappa^2*phi^2)/(1 - xi^2)*tau/(1 - kappa)
acf(sim1, lag.max = 5)$acf
(1 - xi^2 + kappa^2*phi^2 + kappa*phi*(1 - phi))/(1 - xi^2 + kappa^2*phi^2)*kappa*phi*xi^(0:4)
# works

# INARMA:
sim2 <- sim_inarma(tau = tau, phi = phi, kappa = kappa, lgt = 100000)
mean(sim2$X)
var(sim2$X)
tau/(1 - kappa)

acf(sim2$X, lag.max = 5)$acf
kappa*phi*xi^(0:4)
# works

# verify equivalence with thinned INAR(1)
sim3a <- sim_inarma(tau = tau*xi/kappa, phi = 1, kappa = xi)
sim3b <- rbinom(length(sim3a$X), sim3a$X, phi*kappa/xi)
mean(sim3b)
var(sim3b)
tau/(1 - kappa)

acf(sim3b, lag.max = 5)$acf
kappa*phi*xi^(0:4)
