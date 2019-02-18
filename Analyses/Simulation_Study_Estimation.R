# This file contains code to perform the simulation studies underlying Section 6
# of J. Bracher (2019): An INARMA(1, 1) model with Poisson marginals.

setwd("/home/johannes/Documents/underreporting/Proceedings_INARMA")
source("Analyses/Functions.R")

# define true values for the three scenarios:
vals_tau <- c(1, 1, 1)
vals_phi <- c(0.5, 0.8, 0.9)
vals_kappa <- c(0.5, 0.6, 0.8)
vals_lgt <- c(250, 500, 1000) # lengths of simulated time series
n_sim <- 1000 # number of simulation runs

# run simulation:
for(lgt in vals_lgt){ # loop over lengths of time series
  for(s in 1:3){ # loop over scenarios
    print(paste("Started scenario", s))

    # grab true parameter values:
    tau <- vals_tau[s]
    phi <- vals_phi[s]
    kappa <- vals_kappa[s]

    set.seed(111) # set seed
    # simulate all time series at once (wrote custom function for this)
    sim <- sim_many_inarma(tau = tau, phi = phi, kappa = kappa, lgt = lgt, n_sim = n_sim)

    # initializ matrix to store results:
    pars <- ses <- matrix(NA, nrow = n_sim, ncol = 3,
                          dimnames = list(NULL, c("tau", "phi", "kappa")))

    # fir models:
    for(i in 1:n_sim){
      # acf(sim$X[i, ])
      fit_temp <- fit_inarma(sim$X[i, ], hessian = TRUE)
      pars[i, ] <- fit_temp$par # store parameters
      # and store estimated standard errors:
      ses[i, ] <- sqrt(diag(solve(fit_temp$opt$hessian)))
      if(i%%10 == 0 | i < 5) print(i)
    }
    # store (commented out in order not to overwrite reults):
    # print(paste("Storing scenario", s))
    # write.csv(pars, file = paste0("Analyses/Results/means_", tau, "_", phi,
    #                               "_", kappa, "_", lgt, ".csv"))
    # write.csv(ses, file = paste0("Analyses/Results/ses_", tau, "_", phi,
    #                              "_", kappa, "_", lgt, ".csv"))
  }
}


# summarize results:
# initialize lists
results_sim <- list()
removed_est_se <- list()

for(s in 1:3){ # loop over scenarios
  tau <- vals_tau[s]
  phi <- vals_phi[s]
  kappa <- vals_kappa[s]

  scenario_temp <- matrix(NA, nrow = 3, ncol = 13,
                          dimnames = list(c(250, 500, 1000),
                                          c("tau", "phi", "kappa", "T",
                                            "mean_tau", "sd_tau", "est_se_tau",
                                            "mean_phi", "sd_phi", "est_se_phi",
                                            "mean_kappa", "sd_kappa", "est_se_kappa")))
  scenario_temp[1, c("tau", "phi", "kappa")] <- c(tau, phi, kappa)
  scenario_temp[, "T"] <- c(250, 500, 1000)
  removed_est_se_temp <- matrix(NA, ncol = 3, nrow = 3,
                            dimnames = list(c(250, 500, 1000), c("tau", "phi", "kappa")))

  # loop over lengths of time series:
  for(lgt in c(250, 500, 1000)){

    # get parameter estimates:
    pars <- read.csv(file = paste0("Analyses/Results/means_", tau, "_", phi, "_", kappa, "_", lgt, ".csv"))[, c("tau", "phi", "kappa")]
    # fill in means and standard deviations:
    scenario_temp[as.character(lgt), c("mean_tau", "mean_phi", "mean_kappa")] <- colMeans(pars)
    scenario_temp[as.character(lgt), c("sd_tau", "sd_phi", "sd_kappa")] <- apply(pars, 2, sd)

    # get estimated standard errors:
    ses <- read.csv(file = paste0("Analyses/Results/ses_", tau, "_", phi, "_", kappa, "_", lgt, ".csv"))[, c("tau", "phi", "kappa")]
    ses_delta <- NA*ses
    # for tau:
    ses_delta[, "tau"] <- ses[, "tau"]*pars[, "tau"]
    # adapt for phi and kappa
    ses_delta[, "phi"] <- ses[, "phi"]*exp(pars[, "phi"])/(1 + exp(pars[, "phi"]))^2
    ses_delta[, "kappa"] <- ses[, "kappa"]*exp(pars[, "kappa"])/(1 + exp(pars[, "kappa"]))^2

    ses_cleaned <- ses_delta; ses_cleaned[ses_cleaned > 1] <- NA # set numerically instable ses to NA
    scenario_temp[as.character(lgt), c("est_se_tau", "est_se_phi", "est_se_kappa")] <-
      apply(ses_cleaned, 2, mean, na.rm = TRUE)

    # how many times was se estimation not possible or obviously instable?
    removed_est_se_temp[as.character(lgt), ] <- colSums(is.na(ses_cleaned))

  }
  # store results in list:
  results_sim[[s]] <- scenario_temp
  removed_est_se[[s]] <- removed_est_se_temp
}

results_sim
removed_est_se

k <- 1
tau <- vals_tau[k]; phi <- vals_phi[k]; kappa <- vals_kappa[k]; lgt <- 250
pars <- read.csv(file = paste0("Analyses/Results/means_", tau, "_", phi, "_", kappa, "_", lgt, ".csv"))[, c("tau", "phi", "kappa")]

ses <- read.csv(file = paste0("Analyses/Results/ses_", tau, "_", phi, "_", kappa, "_", lgt, ".csv"))[, c("tau", "phi", "kappa")]
ses_delta <- NA*ses
# for tau:
ses_delta[, "tau"] <- ses[, "tau"]*pars[, "tau"]
# adapt for phi and kappa
ses_delta[, "phi"] <- ses[, "phi"]*exp(pars[, "phi"])/(1 + exp(pars[, "phi"]))^2
ses_delta[, "kappa"] <- ses[, "kappa"]*exp(pars[, "kappa"])/(1 + exp(pars[, "kappa"]))^2

ses_cleaned <- ses_delta; ses_cleaned[ses_cleaned > 1] <- NA # set numerically instable ses to NA
scenario_temp[as.character(lgt), c("est_se_tau", "est_se_phi", "est_se_kappa")] <-
  apply(ses_cleaned, 2, mean, na.rm = TRUE)


pars_not_conv <- pars[is.na(rowSums(ses_cleaned)), ]
pars_conv <- pars[!is.na(rowSums(ses_cleaned)), ]

colMeans(pars_conv)
colMeans(pars_not_conv)

ses_delta[is.na(rowSums(ses_cleaned)), ]
