# This file contains functions to reproduce the analyses from
# J. Bracher (2019): An INARMA(1, 1) model with Poisson marginals.


#######################################################
# Functions for likelihood evaluation and fitting
#######################################################

# Transition matrix of an INAR(1) process
#
# Arguments:
# nu: immigration parameter of the INAR(1)
# alpha: autoregressive parameter of the INAR(1)
# upper limit: an upper limit at which to truncate the state space of the INAR(1)
# Return: the transition matrix
get_transition_matrix <- function(nu, alpha, upper_limit){
  n_states <- upper_limit + 1
  immigration_probs <- dpois(0:upper_limit, nu)

  matr <- matrix(nrow = n_states, ncol = n_states)

  for(from in 0:upper_limit){
    for(to in 0:upper_limit){
      matr[from + 1, to + 1] = sum(immigration_probs*dbinom(to - 0:upper_limit, size = from, prob = alpha))
    }
  }
  return(matr)
}

# Evaluate the log-likelihood of an INAR(1) model
#
# Arguments:
# vect: a vector containing the observed time series
# nu, alpha: the parameters of the INAR(1)
# return_fitted: should a vector of fitted values be returned?
# Return: the log-likelihood or (if return_fitted) a list containing the
# log-likelihood and a vector of fitted values
llik_inar <- function(vect, nu, alpha, return_fitted = FALSE){
  upper_limit <- max(qpois(0.999, nu/(1 - alpha)), round(1.5*max(vect)))
  transition_matrix <- get_transition_matrix(nu = nu, alpha = alpha, upper_limit = upper_limit)
  llik_contributions <- fitted_vals <- fitted_vars <- numeric(length(vect))
  llik_contributions[1] <- dpois(vect[1], nu/(1 - alpha), log = TRUE)
  fitted_vals[1] <- fitted_vars[1] <- nu/(1 - alpha)
  for(t in 2:length(vect)){
    llik_contributions[t] <- log(transition_matrix[vect[t - 1] + 1, vect[t] + 1])
    if(return_fitted){
      fitted_vals[t] <- sum(transition_matrix[vect[t - 1] + 1, ]*seq(from = 0, length.out = nrow(transition_matrix)))
      fitted_vars[t] <- sum(transition_matrix[vect[t - 1] + 1, ]*seq(from = 0, length.out = nrow(transition_matrix))^2) - fitted_vals[t]^2
    }
  }
  llik <- sum(llik_contributions)
  if(return_fitted){
    return(list(value = llik, fitted_vals = fitted_vals, fitted_vars = fitted_vars))
  }else{
    return(llik)
  }
}

# Evaluate log-likelihood of hidden INAR(1) model:
#
# Arguments:
# vect: a vector containing the observed time series
# nu, alpha, q: the parameters of the hidden INAR(1)
# log: should log-likelihood or untransformed likelihood be returned?
# return_fp: should forward probabilities be returned?
# Return: the log-likelihood or (if return_fp) a list containing the
# log-likelihood and a matrix of forward probabilities
llik_hidden_inar <- function(vect, nu, alpha, q, log = TRUE, return_fp = FALSE){
  lgt <- length(vect)
  upper_limit <- max(qpois(0.999, nu/(1 - alpha)), round(1.5*max(vect)/q))

  support <- 0:upper_limit
  transition_matrix <- get_transition_matrix(nu = nu, alpha = alpha, upper_limit = upper_limit)

  fp <- matrix(nrow = lgt, ncol = upper_limit + 1)
  normalizing <- numeric(lgt)

  # idea: normalize after ever step and keep normalizing constants
  probs_temp <- dpois(support, nu/(1 - alpha))*dbinom(vect[1], support, q)
  normalizing[1] <- sum(probs_temp)
  fp[1, ] <- probs_temp/normalizing[1]

  for(i in 2:lgt){
    probs_temp <- (fp[i - 1, ] %*% transition_matrix)*dbinom(vect[i], support, q)
    normalizing[i] <- sum(probs_temp)
    fp[i, ] <- probs_temp/normalizing[i]
  }

  return_value <- ifelse(log, sum(log(normalizing)), prod(normalizing))
  if(return_fp){
    return(list(value = return_value, fp = fp))
  }else{
    return(return_value)
  }
}

# Reparameterization from INARMA(1, 1) to hidden INAR(1):
# Arguments:
# tau, phi, kappa: parameters of INARMA(1, 1) model (compare manuscript)
# Result: a named list of the parameters nu, alpha, q of the hidden INAR(1)
reparam_inarma_to_hidden_inar <- function(tau, phi, kappa){
  q_hidden <- phi*kappa/(1 - phi + phi*kappa)
  nu_hidden <- tau*(1 - phi + phi*kappa)/kappa
  alpha_hidden <- 1 - phi + phi*kappa
  return(list(nu = nu_hidden,
              alpha = alpha_hidden,
              q = q_hidden))
}

# Evaluate likelihood of INARMA(1, 1) model
#
# Arguments:
# vect: a vector containing the observed time series
# tau, phi, kappa: parameters of INARMA(1, 1) model (compare manuscript)
# log, return_fp: passed on to lik_hidden_inar (see there)
llik_inarma <- function(vect, tau, phi, kappa, log = TRUE, return_fp = FALSE){
  # re-parametrize to hidden INAR(1)
  pars_hidden_inar <- reparam_inarma_to_hidden_inar(tau = tau, phi = phi, kappa = kappa)
  # evaluate likelihood of corresponding hidden INAR(1) model
  llik_hidden_inar(vect = vect, nu = pars_hidden_inar$nu, alpha = pars_hidden_inar$alpha,
                   q = pars_hidden_inar$q, log = log, return_fp = return_fp)
}

# Fit INARMA(1, 1) model
#
# Arguments:
# vect: a vector containing the observed time series
# ...: additional arguments passed to optim
# Return: a named list containing the parameter estimates, log-likelihood, model dimension
# and object returned by the call to optim
fit_inarma <- function(vect, ...){
  nllik <- function(pars){
    -llik_inarma(vect, tau = exp(pars["log_tau"]),
                 phi = exp(pars["logit_phi"])/(1 + exp(pars["logit_phi"])),
                 kappa = exp(pars["logit_kappa"])/(1 + exp(pars["logit_kappa"])))
  }
  # run optimization:
  opt <- optim(c(log_tau = 2, logit_phi = 0.5, logit_kappa = 0.5), nllik, ...)
  # structure results
  llik <- -opt$value
  par <- c(tau = exp(opt$par["log_tau"]),
           phi = exp(opt$par["logit_phi"])/(1 + exp(opt$par["logit_phi"])),
           kappa = exp(opt$par["logit_kappa"])/(1 + exp(opt$par["logit_kappa"])))
  names(par) <- c("tau", "phi", "kappa")
  return(list(par = par, llik = llik, dim = 3, opt = opt))
}

# Fit hidden INAR(1) model
#
# Arguments:
# vect: a vector containing the observed time series
# ...: additional arguments passed to optim
# Return: a named list containing the parameter estimates, log-likelihood, model dimension
# and object returned by the call to optim
fit_hidden_inar <- function(vect, ...){
  # negative log-likelihood as function of parameter vector
  nllik <- function(pars){
    -llik_hidden_inar(vect, nu = exp(pars["log_nu"]),
                      alpha = exp(pars["logit_alpha"])/(1 + exp(pars["logit_alpha"])),
                      q = exp(pars["logit_q"])/(1 + exp(pars["logit_q"])))
  }
  # run optimization
  opt <- optim(c(log_nu = 2, logit_alpha = 0.5, logit_q = 0.5), nllik, ...)
  # structure results
  llik <- -opt$value
  par <- c(tau = exp(opt$par["log_nu"]),
           alpha = exp(opt$par["logit_alpha"])/(1 + exp(opt$par["logit_alpha"])),
           beta = exp(opt$par["logit_beta"])/(1 + exp(opt$par["logit_beta"])))
  names(par) <- c("nu", "alpha", "beta")
  return(list(par = par, llik = llik, dim = 3, opt = opt))
}

# Fit INAR(1) model:
#
# Arguments:
# vect: a vector containing the observed time series
# ...: additional arguments passed to optim
# Return: a named list containing the parameter estimates, log-likelihood, model dimension
# and object returned by the call to optim
fit_inar <- function(vect, ...){
  # negative log-likelihood as function of parameter vector
  nllik <- function(pars){
    -llik_hidden_inar(vect, nu = exp(pars["log_nu"]),
                      alpha = exp(pars["logit_alpha"])/(1 + exp(pars["logit_alpha"])),
                      q = 1)
  }
  # run optimization:
  opt <- optim(c(log_nu = 2, logit_alpha = 0.5), nllik, ...)
  # structure results:
  par <- c(tau = exp(opt$par["log_nu"]),
           alpha = exp(opt$par["logit_alpha"])/(1 + exp(opt$par["logit_alpha"])))
  llik <- -opt$value
  names(par) <- c("nu", "alpha")
  return(list(par = par, llik = llik, dim = 2, opt = opt))
}

# Evaluate log-likelihood of an INARCH(1) model
#
# Arguments:
# vect: a vector containing the observed time series
# nu, alpha, lambda1: parameters of the INARCH(1) model
# return fitted: should fitted values be returned?
# Return: log-likelihood, or (if return_fitted) a list containing the log-likelihood
# and the fitted values
llik_inarch <- function(vect, nu, alpha, lambda1, log = TRUE, return_fitted = FALSE){
  lgt <- length(vect)
  lambda <- numeric(lgt)
  lambda[1] <- lambda1
  lambda[2:lgt] <- nu + alpha*vect[1:(lgt - 1)]
  llik <- sum(dpois(vect, lambda, log = log))
  if(return_fitted){
    return(list(value = llik, fitted = lambda))
  }else{
    return(llik)
  }
}

# Fitting an INARCH(1) model
#
# Arguments:
# vect: a vector containing the observed time series
# ...: additional arguments passed to optim
# Return: a named list containing the parameter estimates, log-likelihood, model dimension
# and object returned by the call to optim
fit_inarch <- function(vect, ...){
  # negative log-likelihood as function of parameter vector:
  nllik <- function(pars){
    -llik_inarch(vect, nu = exp(pars["log_nu"]),
                 alpha = exp(pars["logit_alpha"])/(1 + exp(pars["logit_alpha"])),
                 lambda1 = exp(pars["log_lambda1"]))
  }
  # run optmization:
  opt <- optim(c(log_nu = 2, logit_alpha = 0.5, log_lambda1 = 0.5), nllik,...)
  # structure results
  llik <- -opt$value
  par <- c(nu = exp(opt$par["log_nu"]),
           alpha = exp(opt$par["logit_alpha"])/(1 + exp(opt$par["logit_alpha"])),
           kappa = exp(opt$par["log_lambda1"]))
  names(par) <- c("nu", "alpha", "lambda1")
  return(list(par = par, llik = llik, dim = 3, opt = opt))
}

# Evaluate log-likelihood for an INGARCH(1, 1) model
#
# Arguments:
# vect: a vector containing the observed time series
# tau, phi, kappa, S1: parameters of the INGARCH(1) model
# return fitted: should fitted values be returned?
# Return: log-likelihood, or (if return_fitted) a list containing the log-likelihood
# and the fitted values
llik_ingarch <- function(vect, tau, phi, kappa, S1, log = TRUE, return_fitted = FALSE){
  lgt <- length(vect)
  # re-parameterize:
  nu <- phi*tau
  alpha <- phi*kappa
  beta <- 1 - phi
  lambda <- numeric(lgt)
  lambda[1] <- (1 - beta)*S1 + nu/(1 - beta)
  for(i in 2:lgt){
    lambda[i] <- nu + alpha*vect[i - 1] + beta*lambda[i - 1]
  }
  llik <- sum(dpois(vect, lambda, log = log))
  if(return_fitted){
    return(list(value = llik, fitted = lambda))
  }else{
    return(llik)
  }
}

# Fitting an INGARCH(1, 1) model
#
# Arguments:
# vect: a vector containing the observed time series
# ...: additional arguments passed to optim
# Return: a named list containing the parameter estimates, log-likelihood, model dimension
# and object returned by the call to optim
fit_ingarch <- function(vect, ...){
  # negative log-likelihood as function of parameter vector:
  nllik <- function(pars){
    -llik_ingarch(vect, tau = exp(pars["log_tau"]),
                  phi = exp(pars["logit_phi"])/(1 + exp(pars["logit_phi"])),
                  kappa = exp(pars["logit_kappa"])/(1 + exp(pars["logit_kappa"])),
                  S1 = exp(pars["log_S1"]))
  }
  opt <- optim(c(log_tau = 2, logit_phi = 0.5, logit_kappa = 0.5, log_S1 = 0.5), nllik,...)
  # structure results
  llik <- -opt$value
  par <- c(tau = exp(opt$par["log_tau"]),
           phi = exp(opt$par["logit_phi"])/(1 + exp(opt$par["logit_phi"])),
           kappa = exp(opt$par["logit_kappa"])/(1 + exp(opt$par["logit_kappa"])),
           kappa = exp(opt$par["log_S1"]))
  names(par) <- c("tau", "phi", "kappa", "S1")
  return(list(par = par, llik = llik, dim = 4, opt = opt))
}

#######################################################
# Simulation functions
#######################################################

# Simulate from a hidden INAR(1) model
#
# Arguments:
# nu, alpha, q: parameters of the hidden INAR(1)
# lgt: the length of the simulated time series
# Return: a named list with the simulated hidden and observed processes
sim_hidden_inar <- function(nu, alpha, q, lgt = 10){
  X <- numeric(lgt)
  X[1] <- rpois(1, nu/(1 - alpha))
  for(i in 2:lgt){
    X[i] <- rpois(1, nu) + rbinom(1, X[i - 1], alpha)
  }
  Y <- rbinom(lgt, X, q)
  return(list(X = X, Y = Y))
}

# Simulate from an INARMA(1, 1) model
#
# Arguments:
# tau, phi, kappa: parameters of the INARMA(1, 1)
# lgt: the length of the simulated time series
# Return: a named list with the simulated processes X, S and I
sim_inarma <- function(tau, phi, kappa, lgt = 10){
  lgt_total <- lgt + 15
  I <- rpois(lgt_total, tau)
  X <- rep(NA, lgt_total); X[1] <- 0
  S <- X
  for(t in 2:lgt_total){
    # version that works:
    S[t] <- S[t - 1] + I[t - 1] - rbinom(1, X[t - 1], 1 - kappa)
    X[t] <- I[t] + rbinom(1, S[t], phi)
  }
  return(list(X = tail(X, lgt), S = tail(S, lgt), I = tail(I, lgt)))
}

# Simulate from an INGARCH(1, 1) model
#
# Arguments:
# tau, phi, kappa: parameters of the INGARCH(1, 1)
# lgt: length of the simulated time series
# Return: The simulated process
sim_ingarch <- function(tau, phi, kappa, lgt){
  lgt_total <- lgt + 100
  X <- S <- numeric(lgt_total)
  for(t in 2:lgt_total){
    S[t] <- kappa*X[t - 1] + (1 - phi)*S[t - 1]
    X[t] <- rpois(1, phi*S[t] + tau)
  }
  return(tail(X, lgt))
}


#######################################################
# Functions for one-step-ahead forecasts
#######################################################

# Get one one-step-ahead-forecast from an INAR(1) model with fixed parameters
#
# Arguments:
# vect: a vector containing the observed time series
# ind_forecast: for which time period shall a forecast be issued?
# nu, alpha: parameters of the INAR(1)
# Return: a vector of forecast probabilities (for 0, ..., round(1.5*max(vect))
osa_inar0 <- function(vect, ind_forecast, nu, alpha){
  upper_limit <- max(qpois(0.999, nu/(1 - alpha)), round(1.5*max(vect)))
  transition_matrix <- get_transition_matrix(nu = nu, alpha, upper_limit = upper_limit)
  return(transition_matrix[vect[ind_forecast - 1] + 1, ])
}

# Get one one-step-ahead-forecast from a hidden INAR(1) model with fixed parameters
#
# Arguments:
# vect: a vector containing the observed time series
# ind_forecast: for which time period shall a forecast be issued?
# nu, alpha, q: parameters of the hidden INAR(1)
# Return: a vector of forecast probabilities (for 0, ..., round(1.5*max(vect))
osa_hidden_inar0 <- function(vect, ind_forecast, nu, alpha, q){

  upper_limit <- max(qpois(0.999, nu/(1 - alpha)), round(1.5*max(vect)/q))
  support <- 0:upper_limit

  # get transistion matrix of latent INAR(1):
  transition_matrix <- get_transition_matrix(nu = nu, alpha = alpha, upper_limit = upper_limit)

  # run through forward pass procedure
  fp <- matrix(nrow = ind_forecast, ncol = upper_limit + 1)
  normalizing <- numeric(ind_forecast - 1)

  # initialize:
  probs_temp <- dpois(support, nu/(1 - alpha))*dbinom(vect[1], support, q)
  normalizing[1] <- sum(probs_temp)
  fp[1, ] <- probs_temp/normalizing[1]

  # run through forward pass scheme:
  if(ind_forecast > 2){
    for(i in 2:(ind_forecast - 1)){
      probs_temp <- (fp[i - 1, ] %*% transition_matrix)*dbinom(vect[i], support, q)
      normalizing[i] <- sum(probs_temp)
      fp[i, ] <- probs_temp/normalizing[i]
    }
  }

  # get forecast for latent INAR(1) process
  forecast_X <- fp[ind_forecast - 1, ] %*% transition_matrix

  # get forecast for observable (i.e. INARMA(1, 1)) process:
  forecast_Y <- numeric(upper_limit + 1)
  for(i in support){
    forecast_Y <- forecast_Y + forecast_X[i + 1]*dbinom(support, size = i, q)
  }

  return(forecast_Y)
}

# Get one-step-ahead forecasts from INAR(1) model for a range of weeks
# The model will be re-fit for each forecast.
#
# Arguments:
# vect: a vector containing the observed time series
# from, to: begin and end of the period to forecast
# Return: a matrix of forecast probabilities (for each week: 0, ..., round(1.5*max(vect))
osa_inar <- function(vect, from, to){
  inds_to_forecast <- from:to
  pred_distr <- list()
  logS_vector <- numeric(length(inds_to_forecast)); names(logS_vector) <- inds_to_forecast

  pb <- txtProgressBar(min = 1, max = length(inds_to_forecast))
  for(i in seq_along(inds_to_forecast)){
    ind <- inds_to_forecast[i]
    # re-fit model:
    fit <- fit_inar(vect[1:(ind - 1)])
    # run one-step-ahead procedure:
    osa_temp <- osa_inar0(vect = vect, ind = ind, nu = fit$par["nu"], alpha = fit$par["alpha"])
    pred_distr[[i]] <- osa_temp
    logS_vector[i] <- log(osa_temp[vect[ind] + 1])
    setTxtProgressBar(pb, i)
  }
  return(list(logS = sum(logS_vector), logS_vector = logS_vector, pred_distr = pred_distr))
}

# Get one one-step-ahead-forecast from an INARMA(1, 1) model with fixed parameters
#
# Arguments:
# vect: a vector containing the observed time series
# ind_forecast: for which time period shall a forecast be issued?
# tau, phi, kappa: parameters of the hidden INAR(1)
# Return: a vector of forecast probabilities (for 0, ..., round(1.5*max(vect))
osa_inarma0 <- function(vect, ind_forecast, tau, phi, kappa){
  # re-parameterize
  par_hidden_inar <- reparam_inarma_to_hidden_inar(tau = tau, phi = phi, kappa = kappa)
  # plug into forecasting function for hidden INAR(1):
  osa_hidden_inar0(vect = vect, ind_forecast = ind_forecast, nu = par_hidden_inar$nu,
                   alpha = par_hidden_inar$alpha, q = par_hidden_inar$q)
}

# Get one-step-ahead forecasts from INARMA(1, 1) model for a range of weeks
# The model will be re-fit for each forecast.
#
# Arguments:
# vect: a vector containing the observed time series
# from, to: begin and end of the period to forecast
# Return: a matrix of forecast probabilities (for each week: 0, ..., round(1.5*max(vect))
osa_inarma <- function(vect, from, to){
  inds_to_forecast <- from:to
  pred_distr <- list()
  logS_vector <- numeric(length(inds_to_forecast)); names(logS_vector) <- inds_to_forecast

  pb <- txtProgressBar(min = 1, max = length(inds_to_forecast))
  for(i in seq_along(inds_to_forecast)){
    ind <- inds_to_forecast[i]
    # re-fit model:
    fit <- fit_inarma(vect[1:(ind - 1)])
    # run one-step-ahead procedure:
    osa_temp <- osa_inarma0(vect = vect, ind_forecast = ind, tau = fit$par["tau"],
                            phi = fit$par["phi"], kappa = fit$par["kappa"])
    pred_distr[[i]] <- osa_temp
    logS_vector[i] <- log(osa_temp[vect[ind] + 1])
    setTxtProgressBar(pb, i)
  }
  return(list(logS = sum(logS_vector), logS_vector = logS_vector, pred_distr = pred_distr))
}

# Get one-step-ahead forecasts from INARCH(1) model for a range of weeks
# The model will be re-fit for each forecast.
#
# Arguments:
# vect: a vector containing the observed time series
# from, to: begin and end of the period to forecast
# Return: a matrix of forecast probabilities (for each week: 0, ..., round(1.5*max(vect))
osa_inarch <- function(vect, from, to){
  inds_to_forecast <- from:to
  logS_vector <- numeric(length(inds_to_forecast)); names(logS_vector) <- inds_to_forecast
  # set up progress bar:
  pb <- txtProgressBar(min = 1, max = length(inds_to_forecast))
  for(i in seq_along(inds_to_forecast)){
    # index the i-th timepoint in from:to takes in vect
    ind <- inds_to_forecast[i]
    # fit model up to ind - 1
    fit <- fit_inarch(vect[1:(ind - 1)])
    # extract parameters:
    nu <- fit$par["nu"]
    alpha <- fit$par["alpha"]
    # compute predicted value and evaluate log score
    logS_vector[i] <- dpois(vect[ind], nu + alpha*vect[ind - 1], log = TRUE)
    setTxtProgressBar(pb, i)
  }
  return(list(logS = sum(logS_vector), logS_vector = logS_vector))
}

# Get one-step-ahead forecasts from INGARCH(1, 1) model for a range of weeks
# The model will be re-fit for each forecast.
#
# Arguments:
# vect: a vector containing the observed time series
# from, to: begin and end of the period to forecast
# Return: a matrix of forecast probabilities (for each week: 0, ..., round(1.5*max(vect))
osa_ingarch <- function(vect, from, to){
  inds_to_forecast <- from:to
  logS_vector <- numeric(length(inds_to_forecast)); names(logS_vector) <- inds_to_forecast
  # set up progress bar:
  pb <- txtProgressBar(min = 1, max = length(inds_to_forecast))
  for(i in seq_along(inds_to_forecast)){
    # index the i-th timepoint in from:to takes in vect
    ind <- inds_to_forecast[i]
    # fit model up to ind - 1
    fit <- fit_ingarch(vect[1:(ind - 1)])
    # re-parameterize:
    nu <- fit$par["phi"]*fit$par["tau"]
    alpha <- fit$par["phi"]*fit$par["kappa"]
    beta <- 1 - fit$par["phi"]
    # compute lambda_ind (truncating lags at 10)
    lambda_ind <- nu/(1 - beta) + sum(alpha*beta^(0:9)*vect[(ind - 1):(ind - 10)])
    # evaluate log score
    logS_vector[i] <- dpois(vect[ind], lambda_ind, log = TRUE)
    setTxtProgressBar(pb, i)
  }
  return(list(logS = sum(logS_vector), logS_vector = logS_vector))
}

#######################################################
# Additional helper functions to for simulation studies
#######################################################

# wrapper to simulate many time series from an INARMA(1, 1) model:
sim_many_inarma <- function(tau, phi, kappa, lgt, n_sim){
  X_samples <- Y_samples <- matrix(ncol = lgt, nrow = n_sim)
  for(i in 1:n_sim){
    sim_temp <- sim_inarma(tau, phi, kappa, lgt)
    X_samples[i, ] <- sim_temp$X
  }
  return(list(X = X_samples))
}

# wrapper to simulate many time series from a hidden INAR(1) model:
sim_many_hidden_inar <- function(nu, alpha, q, lgt, n_sim){
  X_samples <- Y_samples <- matrix(ncol = lgt, nrow = n_sim)
  for(i in 1:n_sim){
    sim_temp <- sim_hidden_inar(nu, alpha, q, lgt)
    X_samples[i, ] <- sim_temp$X
    Y_samples[i, ] <- sim_temp$Y
  }
  return(list(X = X_samples, Y = Y_samples))
}

# evaluate empirical proportion of a given sequence (used in empirical likelihood
# evaluation of short sequences):
emp_prob <- function(Y_samples, vect){
  mean(colSums(t(Y_samples) == vect) == ncol(Y_samples))
}

# get acf in convenient format:
decomp_acf <- function(vect){
  acf <- acf(vect)
  return(c(gamma1 = acf$acf[2], decay = acf$acf[3]/acf$acf[2]))
}