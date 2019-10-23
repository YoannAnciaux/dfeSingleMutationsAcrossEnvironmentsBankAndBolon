#!/usr/bin/env Rscript

library(abc)

#### Parameters ####
nval <- 100
tols <- c(10^(-3), 10^(-4))
method <- c("rejection", "loclinear", "neuralnet", "ridge")

#### One environment CV ####
data("reference_fgmrmut_parameter")
data("reference_fgmrmut_simulation")

one_environment_cv <- sapply(method,
                             FUN = function(m){tryCatch(cv4abc(param   = reference_fgmrmut_parameter,
                                                               sumstat = reference_fgmrmut_simulation,
                                                               nval    = nval,
                                                               tols    = tols,
                                                               method  = m),
                                                        error = function(e) e)},
                             simplify = F,
                             USE.NAMES = T)
use_data(one_environment_cv, overwrite = T)

summary_one_environment_cv <- sapply(one_environment_cv,
                                     FUN = function(i){s <- summary(i)
                                     if(!is.numeric(s)) s <- NA
                                     round(s, digits = 2)},
                                     simplify = F)
use_data(summary_one_environment_cv, overwrite = T)


#### Two environments CV ####
data("new_fgmrmut2env_parameter")
data("new_fgmrmut2env_simulation")

two_environments_cv <- sapply(method,
                             FUN = function(m){tryCatch(cv4abc(param   = new_fgmrmut2env_parameter,
                                                               sumstat = new_fgmrmut2env_simulation,
                                                               nval    = nval,
                                                               tols    = tols,
                                                               method  = m),
                                                        error = function(e) e)},
                             simplify = F,
                             USE.NAMES = T)
use_data(two_environments_cv, overwrite = T)

summary_two_environments_cv <- sapply(two_environments_cv,
                                     FUN = function(i){s <- summary(i)
                                     if(!is.numeric(s)) s <- NA
                                     round(s, digits = 2)},
                                     simplify = F)
use_data(summary_two_environments_cv, overwrite = T)
