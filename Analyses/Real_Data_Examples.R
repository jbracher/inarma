#######################################################
# Real data study
#######################################################

setwd("/home/johannes/Documents/underreporting/Theory_Equivalence/Proceedings_INARMA")
source("Analyses/Functions.R")

################################################
# Westgren gold particle data:
# get data:
data_gold <- read.csv("Data/gold_westgren.csv")

# Plot:
par(mfrow = 1:2)
plot(data_gold$counts, type = "l", xlab = "", ylab = "particle counts", main  ="Westgren's gold particle data", ylim = c(0, 10))
acf(data_gold$counts, lag.max = 15, ci.col = "black", main = "")
text(10, 0.84, paste("mean =", round(mean(data_gold$counts), 2)))
text(10, 0.61, paste("   var =", round(var(data_gold$counts), 2)))

lgt <- length(data_gold)

# fit four models:
fit_gold_inar <- fit_inar(data_gold$counts)
fit_gold_inarma <- fit_inarma(data_gold$counts)
fit_gold_inarch <- fit_inarch(data_gold$counts)
# # for comparison: fit INARCH(1) using hhh4:
# fit_gold_hhh4 <- hhh4(new("sts", observed = data_gold),
#                       list(ar = list(f = ~1), subset = 2:lgt))
fit_gold_ingarch <- fit_ingarch(data_gold$counts)

# compare parameter estimates:
fit_gold_inar$par
fit_gold_inarma$par
fit_gold_inarch$par
# summary(fit_gold_hhh4, idx2Exp = TRUE)
fit_gold_ingarch$par

# compare model fits:
fit_gold_inar$llik
fit_gold_inarma$llik
fit_gold_inarch$llik
fit_gold_ingarch$llik

# AICs:
-2*fit_gold_inar$llik + 2*2
-2*fit_gold_inarma$llik + 2*3
-2*fit_gold_inarch$llik + 2*3
-2*fit_gold_ingarch$llik + 2*4

# run one-step ahead forecasts:
osa_gold_inar <- osa_inar(data_gold$counts, from = floor(length(data_gold$counts)/2) + 1, to = length(data_gold$counts))
osa_gold_inarma <- osa_inarma(data_gold$counts, from = floor(length(data_gold$counts)/2) + 1, to = length(data_gold$counts))
osa_gold_inarch <- osa_inarch(data_gold$counts, from = floor(length(data_gold$counts)/2) + 1, to = length(data_gold$counts))
osa_gold_ingarch <- osa_ingarch(data_gold$counts, from = floor(length(data_gold$counts)/2) + 1, to = length(data_gold$counts))

mean(osa_gold_inar$logS_vector)
mean(osa_gold_inarma$logS_vector)
mean(osa_gold_inarch$logS_vector)
mean(osa_gold_ingarch$logS_vector)

logS_gold <- cbind(time = (floor(length(data_gold$counts)/2) + 1):length(data_gold$counts),
                   inar = osa_gold_inar$logS_vector,
                   inarma = osa_gold_inarma$logS_vector,
                   inarch = osa_gold_inarch$logS_vector,
                   ingarch = osa_gold_ingarch$logS_vector)

# store results:
# write.csv(logS_gold, file = "Analyses/Results/logS_gold.csv", row.names = FALSE)
# save(osa_gold_inar, osa_gold_inarma, osa_gold_inarch, osa_gold_ingarch,
#      file = "Analyses/Results/osa_gold.rda")
# load("Analyses/Results/osa_gold.rda")


# Analyses for INAR(2) and CINAR(2) models are omitted here as they are based on
# code by C.H. Weiss which is part of the (protected) Supplementary Material of
# C.H. Weiss: An Introduction to Discrete-Valued Time Series, Wiley, 2018.
# (See here: https://www.wiley.com/en-us/An+Introduction+to+Discrete+Valued+Time+Series-p-9781119096962)

################################################
# Bavarian mumps data:
data_mumps <- read.csv("Data/mumps_germany.csv")

# (Note that the data set also contains data from the 15 other states of Germany.
# Just replace Bavaria by e.g. Hamburg).

# Plot:
par(mfrow = 1:2)
plot(data_mumps$time, data_mumps$Bavaria, type = "l", xlab = "", ylab = "case counts",
     main = "Mumps in Bavaria (weekly)")
acf(data_mumps$Bavaria, lag.max = 15, ci.col = "black", main = "")
text(10, 0.84, paste("mean =", round(mean(data_mumps$Bavaria), 2)))
text(10, 0.61, paste("   var =", round(var(data_mumps$Bavaria), 2)))


# fit four models:
fit_mumps_inar <- fit_inar(data_mumps$Bavaria)
fit_mumps_inarma <- fit_inarma(data_mumps$Bavaria)
fit_mumps_inarch <- fit_inarch(data_mumps$Bavaria)
fit_mumps_ingarch <- fit_ingarch(data_mumps$Bavaria)

# AICs:
-2*fit_mumps_inar$llik + 2*2
-2*fit_mumps_inarma$llik + 2*3
-2*fit_mumps_inarch$llik + 2*3
-2*fit_mumps_ingarch$llik + 2*4

# run one-step ahead forecasts:
osa_mumps_inar <- osa_inar(data_mumps$Bavaria, from = floor(length(data_mumps$Bavaria)/2) + 1, to = length(data_mumps$Bavaria))
osa_mumps_inarma <- osa_inarma(data_mumps$Bavaria, from = floor(length(data_mumps$Bavaria)/2) + 1, to = length(data_mumps$Bavaria))
osa_mumps_inarch <- osa_inarch(data_mumps$Bavaria, from = floor(length(data_mumps$Bavaria)/2) + 1, to = length(data_mumps$Bavaria))
osa_mumps_ingarch <- osa_ingarch(data_mumps$Bavaria, from = floor(length(data_mumps$Bavaria)/2) + 1, to = length(data_mumps$Bavaria))

# compare results:
mean(osa_mumps_inar$logS_vector)
mean(osa_mumps_inarma$logS_vector)
mean(osa_mumps_inarch$logS_vector)
mean(osa_mumps_ingarch$logS_vector)


# store results:
# write.csv(logS_mumps, file = "Analyses/Results/logS_mumps.csv", row.names = FALSE)
# save(osa_mumps_inar, osa_mumps_inarma, osa_mumps_inarch, osa_mumps_ingarch,
#      file = "Analyses/Results/osa_mumps.rda")
# load("Analyses/Results/osa_mumps.rda")