library(inferenceFitnessLandscape)
library(here)
library(abc)

#### Parameters for simulations ####
rgeom_trunc <- function(n, a, b, prob) {truncdist::rtrunc(n, spec = "geom", a, b, prob)}
rexp_trunc <- function(n, a, b, rate) {truncdist::rtrunc(n, spec = "exp", a, b, rate)}
fitness_wt <- read.table(file = here("inst", "raw_data", "fake_115_mutations_environment_reference.csv"), skip = 1)[1, 116]
nb_simul <- 10^6
pn <- 1/5
rlambda <- 1/0.2
rmaxfitness <- 1/2
alpha_inf <- 0.1; alpha_sup <- 10
Q_inf <- 0.5; Q_sup <- 6
n_prior <- rgeom_trunc(nb_simul, 0, Inf, pn)
matrix_param <- cbind(n = n_prior,
                      lambda = rexp(nb_simul, rlambda),
                      maxfitness = rexp_trunc(nb_simul, fitness_wt, Inf, rmaxfitness),
                      alpha = runif(nb_simul, alpha_inf, alpha_sup),
                      Q = runif(nb_simul, Q_inf, Q_sup),
                      m = n_prior)
write.table(x = matrix_param,
            file = here("inst", "simulation", "115_mutations_fgmrmut_parameters.csv"),
            col.names = TRUE, sep = "\t")

#### Simulations ####
sim <- simulate_fl(parameter = matrix_param, simulation_model = "fgmrmut",
                   empirical_fl = here("inst", "raw_data", "fake_115_mutations_environment_reference.csv"), ncore = 7, fun_args = fun_args,
                   file_output = here("inst", "simulation", "115_mutations_fgmrmut_simulations.csv"),
                   multi_file = F, output_args = list(sep = "\t"), skip = 1)

#### Estimation ####
input <- read_clean_output(file = here("inst", "simulation", "115_mutations_fgmrmut_simulations.csv"),
                           keep_param = c("n", "lambda", "maxfitness", "alpha", "Q"),
                           header = TRUE, sep = "\t")



target <- matrix(read.table(file = here("inst", "raw_data", "fake_115_mutations_environment_reference.csv"), skip = 1)[2:116, 116], nrow = 1)

abc_out <- list()
abc_out$rejection <- abc(target = target, param = input$parameter,
                         sumstat = input$simulation, tol = 0.01,
                         method ="rejection")
abc_out$loclinear <- abc(target = target, param = input$parameter,
                         sumstat = input$simulation, tol = 0.01,
                         method = "loclinear")
abc_out$neuralnet <- abc(target = target, param = input$parameter,
                         sumstat = input$simulation, tol = 0.01,
                         method = "neuralnet")
abc_out$ridge <- abc(target = target, param = input$parameter,
                     sumstat = input$simulation, tol = 0.01,
                     method = "ridge")

summary(abc_out$rejection)
summary(abc_out$loclinear)
summary(abc_out$neuralnet)
summary(abc_out$ridge)


#### CV ####

cv_out <- list()
cv_out$rejection <- cv4abc(param = input$parameter,
       sumstat = input$simulation, tols = 0.01,
       nval = 10, abc.out = abc_out$rejection)
