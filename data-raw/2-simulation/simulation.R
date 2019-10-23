library(inferenceFitnessLandscape)
library(usethis)

#### One environment (reference environment) ####
data("reference_fgmrmut_parameter")
data("full_raw_empirical_fl")
reference_empirical_fl <- full_raw_empirical_fl$std_rep1
remove(full_raw_empirical_fl)
temp_sim <- simulate_fl(parameter        = reference_fgmrmut_parameter,
                        simulation_model = "fgmrmut",
                        empirical_fl     = reference_empirical_fl,
                        ncore            = 7)

reference_fgmrmut_simulation <- temp_sim$simulation[, -1] #removes wild type which is not useful for the inference

write.table(x = reference_fgmrmut_simulation,
            file = "data-raw/2-simulation/reference_fgmrmut_simulation.csv",
            col.names = TRUE)

use_data(reference_fgmrmut_simulation, overwrite = TRUE)


#### One environment (reference environment) ####
data("reference_fgmrmut_parameters")
data("reference_empirical_fl")
temp_sim <- simulate_fl(parameter        = reference_fgmrmut_parameters,
                        simulation_model = "fgmrmut",
                        empirical_fl     = reference_empirical_fl,
                        ncore            = 7)

reference_fgmrmut_simulation <- temp_sim$simulation[, -1] #removes wild type which is not useful for the inference

write.table(x = reference_fgmrmut_simulation,
            file = "data-raw/reference_fgmrmut_simulation.csv",
            col.names = TRUE)

use_data(reference_fgmrmut_simulation, overwrite = TRUE)

#### Two environments (new environment) ####
data("new_fgmrmut2env_parameters")
data("new_empirical_fl")
data("reference_inferred_parameters")

temp_sim <- simulate_fl(parameter        = new_fgmrmut2env_parameters,
                        simulation_model = "fgmrmut2env",
                        empirical_fl     = empirical_fl_new,
                        ncore            = 7,
                        fun_args         = reference_inferred_parameters)

new_fgmrmut2env_simulation <- temp_sim$simulation[, - c(1, ncol(temp_sim$simulation) / 2 + 1)] #removes wild types which are not useful for the inference

write.table(x = new_fgmrmut2env_simulation,
            file = "data-raw/new_fgmrmut2env_simulation.csv",
            col.names = TRUE)

use_data(new_fgmrmut2env_simulation, overwrite = TRUE)
