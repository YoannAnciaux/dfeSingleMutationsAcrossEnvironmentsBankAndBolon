library(MASS)
library(inferenceFitnessLandscape)
library(usethis)

#### Parameters ####
nb_mut         <- 115
n              <- 2
lambda         <- 0.1
maxfitness     <- 1
alpha          <- 1/2
Q              <- 2
theta          <- pi / 2

fitness_wt_ref <- 0
fitness_wt_new <- 0

# #### Save reference parameters (temporary) ####
# reference_inferred_parameters <- list(fitness_wt_ref = fitness_wt_ref,
#                                       n_ref = n,
#                                       lambda_ref = lambda,
#                                       maxfitness_ref = maxfitness,
#                                       alpha_ref = alpha,
#                                       Q_ref = Q,
#                                       m_ref = n)
# write.table(x = t(as.matrix(reference_inferred_parameters)),
#             file = "data-raw/reference_inferred_parameters.csv",
#             col.names = T)
# use_data(reference_inferred_parameters, overwrite = T)

#### Creates wt and mutations effects ####
lambda_I_n   <- lambda * diag(n)
mut_effects  <- mvrnorm(n = nb_mut, mu = numeric(n), Sigma = lambda_I_n)

pheno_wt <- ftop_fgm_iso(fitness_wt_ref, n, maxfitness, alpha, Q)
mutant  <- mut_effects + matrix(pheno_wt, nb_mut, n, byrow = T)

#### environment reference ####
rmut_fit_ref           <- ptof_fgm_iso(mutant, maxfitness, alpha, Q)
reference_empirical_fl <- cbind(rbind(numeric(nb_mut),
                                      diag(nrow = nb_mut,
                                           ncol = nb_mut)),
                                c(fitness_wt_ref, rmut_fit_ref))

write.table(x = reference_empirical_fl,
            file = "data-raw/fake_reference_empirical_fl.csv",
            row.names = FALSE, col.names = TRUE)

use_data(reference_empirical_fl, overwrite = TRUE)

#### environment new ####
rmut_fit_new     <- ptof_fgm_iso(mutant, maxfitness, alpha, Q,
                                 pheno_opt = c((2 * (maxfitness - fitness_wt_new))^(1/2)* cos(theta),
                                               (2 * (maxfitness - fitness_wt_new))^(1/2)* sin(theta)))
new_empirical_fl <- cbind(rbind(numeric(nb_mut),
                                diag(nrow = nb_mut,
                                     ncol = nb_mut)),
                          c(fitness_wt_new, rmut_fit_new))

write.table(x = new_empirical_fl,
            file= "data-raw/fake_new_empirical_fl.csv",
            row.names = FALSE, col.names = TRUE)

use_data(new_empirical_fl, overwrite = TRUE)
