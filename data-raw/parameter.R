library(truncdist)
library(inferenceFitnessLandscape)
library(readr)
library(usethis)

#### Custom truncated distributions ####
rgeom_trunc <- function(n, a, b, prob) {rtrunc(n, spec = "geom", a, b, prob)}
rexp_trunc <- function(n, a, b, rate) {rtrunc(n, spec = "exp", a, b, rate)}

#### Generates parameters environment reference ####
data("reference_empirical_fl")
fitness_wt_ref  <- reference_empirical_fl[1, dim(reference_empirical_fl)[2]]

nb_simul    <- 10^6

pn          <- 1/5
rlambda     <- 1/0.2
rmaxfitness <- 1/2
alpha_inf   <- 0.1
alpha_sup   <- 10
Q_inf       <- 0.5
Q_sup       <- 6

n_prior     <- rgeom_trunc(nb_simul, 0, Inf, pn)

reference_fgmrmut_parameter <- cbind(n = n_prior,
                                    lambda = rexp(nb_simul, rlambda),
                                    maxfitness = rexp_trunc(nb_simul, fitness_wt_ref, Inf, rmaxfitness),
                                    alpha = runif(nb_simul, alpha_inf, alpha_sup),
                                    Q = runif(nb_simul, Q_inf, Q_sup),
                                    m = n_prior)

write.table(x = reference_fgmrmut_parameter,
            file = "data-raw/reference_fgmrmut_parameter.csv",
            col.names = TRUE)

use_data(reference_fgmrmut_parameter, overwrite = TRUE)

#### Generates parameters environment new ####
data("new_empirical_fl")
fitness_wt_new  <- new_empirical_fl[1, dim(new_empirical_fl)[2]]

nb_simul    <- 10^6

rlambda     <- 1/0.2
rmaxfitness <- 1/2
alpha_inf   <- 0.1
alpha_sup   <- 10
Q_inf       <- 0.5
Q_sup       <- 6
theta_inf   <- 0
theta_sup   <- pi

new_fgmrmut2env_parameter <- cbind(lambda = rexp(nb_simul, rlambda),
                                    maxfitness = rexp_trunc(nb_simul, fitness_wt_new, Inf, rmaxfitness),
                                    alpha = runif(nb_simul, alpha_inf, alpha_sup),
                                    Q = runif(nb_simul, Q_inf, Q_sup),
                                    theta = runif(nb_simul, theta_inf, theta_sup))

write.table(x = new_fgmrmut2env_parameter,
            file = "data-raw/new_fgmrmut2env_parameter.csv",
            col.names = TRUE)

use_data(new_fgmrmut2env_parameter, overwrite = TRUE)
