# library(inferenceFitnessLandscape)
# library(here)
# library(abc)
# library(readr)
# library(dplyr)
#
# #### Parameters for simulations ####
# fitness_wt_new <- 0
# pheno_wt_new <- ftop_fgm_iso(fitness_wt_new, 2, 1)
# rexp_trunc <- function(n, a, b, rate) {truncdist::rtrunc(n, spec = "exp", a, b, rate)}
# nb_simul <- 10^4
# n_prior <- rep(2, nb_simul)
# rlambda <- 1/0.2
# rmaxfitness <- 1/2
# alpha_inf <- 1/2; alpha_sup <- 1/2
# Q_inf <- 2; Q_sup <- 2
# theta_inf <- 0; theta_sup <- pi
# matrix_param <- cbind(lambda = rexp(nb_simul, rlambda),
#                       maxfitness = rexp_trunc(nb_simul, fitness_wt_new, Inf, rmaxfitness),
#                       alpha = runif(nb_simul, alpha_inf, alpha_sup),
#                       Q = runif(nb_simul, Q_inf, Q_sup),
#                       theta = runif(nb_simul, theta_inf, theta_sup))
# write.table(x = matrix_param,
#             file = here("inst", "simulation", "115_mutations_fgmrmut2env_10e4_parameters.csv"),
#             col.names = TRUE, sep = "\t")
#
# #### Parameters of reference environment ####
# fun_args <- list(fitness_wt_ref = 0, n_ref = 2, lambda_ref = 0.1, maxfitness_ref = 1,
#                  alpha_ref = 1/2, Q_ref = 2, m_ref = 2)
# save(fun_args, file = here("inst", "simulation", "115_mutations_fgmrmut2env_10e4_reference_environment.Rda"))
#
# #### Fitness landscape ####
# lambda_I_n <- 0.1 * diag(2)
# A = MASS::mvrnorm(n = 115, mu = numeric(2), Sigma = lambda_I_n)
# rmut_fit_ref <- ptof_fgm_iso(A + matrix(pheno_wt_new, 115, 2, byrow = T),
#                              1, 1/2, 2)
# empirical_fl_ref <- cbind(rbind(numeric(115), diag(nrow = 115, ncol = 115)), c(0, rmut_fit_ref))
# empirical_fl_ref <- as_tibble(empirical_fl_ref)
# write_csv(x = empirical_fl_ref,
#           path = "data-raw/fake_empirical_fl_environment_reference.csv")
# test <- read_csv(file = "data-raw/fake_empirical_fl_environment_reference.csv")
# rmut_fit_new <- ptof_fgm_iso(A + matrix(pheno_wt_new, 115, 2, byrow = T),
#                              1, 1/2, 2, pheno_opt = c(2^(1/2)* cos(pi/2), 2^(1/2)* sin(pi/2)))
# empirical_fl_new <- cbind(rbind(numeric(115), diag(nrow = 115, ncol = 115)), c(0, rmut_fit_new))
# empirical_fl_new <- tibble(empirical_fl_new)
# write_csv(x = empirical_fl_ref,
#           path = "data-raw/fake_empirical_fl_environment_new.csv")
#
# #### Simulations ####
# sim <- simulate_fl(parameter = matrix_param, simulation_model = "fgmrmut2env",
#                    empirical_fl = empirical_fl, ncore = 7, fun_args = fun_args,
#                    file_output = here("inst", "simulation", "115_mutations_fgmrmut2env_10e4_simulations.csv"),
#                    multi_file = F, output_args = list(sep = "\t"))
#
# #### Estimation ####
# input <- read_clean_output(file = here("inst", "simulation", "115_mutations_fgmrmut2env_10e4_simulations.csv"),
#                            keep_param = c("lambda", "maxfitness", "theta"),
#                            header = TRUE, sep = "\t")
# target <- cbind(matrix(read.table(file = here("inst", "raw_data", "fake_115_mutations_environment_reference.csv"), skip = 1)[2:116, 116], nrow = 1),
#                 matrix(read.table(file = here("inst", "raw_data", "fake_115_mutations_environment_new.csv"), skip = 1)[2:116, 116], nrow = 1))
#
# abc_out <- list()
# abc_out$rejection <- abc(target = target, param = input$parameter,
#                          sumstat = input$simulation, tol = 0.01,
#                          method ="rejection")
# abc_out$loclinear <- abc(target = target, param = input$parameter,
#                          sumstat = input$simulation, tol = 0.01,
#                          method = "loclinear")
# abc_out$neuralnet <- abc(target = target, param = input$parameter,
#                          sumstat = input$simulation, tol = 0.01,
#                          method = "neuralnet")
# abc_out$ridge <- abc(target = target, param = input$parameter,
#                      sumstat = input$simulation, tol = 0.01,
#                      method = "ridge")
#
# summary(abc_out$rejection)
# summary(abc_out$loclinear)
# summary(abc_out$neuralnet)
# summary(abc_out$ridge)
#
# #### CV ####
# if (is.null(names(files))) {
#   names(files) = sapply(files, FUN = file_path_sans_ext)
# }
#
# cv_out <- list()
# cv_out$rejection <- cv4abc(param = input$parameter,
#                            sumstat = input$simulation, tols = 0.01,
#                            nval = 10, abc.out = abc_out$rejection)
#
# index <- as.vector(sapply(names(input),
#                           FUN = function(n) {rep(n, nrow(input[[n]]$parameter))}))
# sumstat <- do.call(rbind, sapply(names(input),
#                                  FUN = function(n) {input[[n]]$simulation},
#                                  simplify = F))
#
# postpr_out <- postpr(target = target, index = index, sumstat = sumstat,
#                      tol = tol_abc, method = method, ...)
# cv4postpr_out <- cv4postpr(index = index, sumstat = sumstat, nval = nval,
#                            tols = tols_cv, postpr.out = postpr_out)
# f = file()
# sink(file = f)
# parameter_prior <- map(input,
#                        ~ as_tibble(.x$parameter))
#
# parameter_estim <- map(abc_out,
#                        ~ as_tibble(.x$unadj.values)  %>%
#                          setNames(.x$names$parameter.names))
# parameter_error <- list()
# parameter_error$mRSSE <- map(cv4abc_out,
#                              ~ as.data.frame.matrix(summary(.x)) %>%
#                                as_tibble(rownames = 'tol'))
# parameter_error$true_vs_estim <- map(cv4abc_out,
#                                      ~ {
#                                        t <- as_tibble(.x$true)
#                                        map(.x$estim,
#                                            ~ as_tibble(.x) %>%
#                                              setNames(names(t)) %>%
#                                              gather('parameter','estim') %>%
#                                              bind_cols(t %>%
#                                                          gather('parameter','true')
#                                                        %>% select(true))) %>%
#                                          bind_rows(.id = 'tol') %>%
#                                          mutate(tol = as.factor(str_remove(tol, 'tol')))
#                                      })
# model_estim <-
#   tibble(model = postpr_out$values) %>%
#   count(model)
#
# model_misclass <-
#   summary(cv4postpr_out) %>%
#   map(~ map(.x, ~ as.data.frame.matrix(.x) %>%
#               as_tibble(rownames = 'model')) %>%
#         bind_rows(.id = 'tol')) %>%
#   bind_rows(.id = 'misclass') %>%
#   mutate(misclass = as.factor(misclass),
#          tol = as.factor(str_remove(tol, 'tol')),
#          model = as.factor(model))
# sink()
# close(f)
# list(raw = list(input = input,
#                 abc_out = abc_out, cv4abc_out = cv4abc_out,
#                 postpr_out = postpr_out, cv4postpr_out = cv4postpr_out),
#      tidy_summary = list(parameter_prior = parameter_prior,
#                          parameter_estim = parameter_estim,
#                          parameter_error = parameter_error,
#                          model_estim = model_estim,
#                          model_misclass = model_misclass))
# }
#
