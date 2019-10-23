library(tidyverse)
library(devtools)

#### Useful functions ####
raw2fl <- function(data, data_name) {
  # counts the number of mutants and keep the names of the mutants
  data_wo_wt <- data %>%
    filter(Position_AA != "NA_WILDTYPE")
  nb_mut <- data_wo_wt %>%
    nrow(.)
  names_mutants <- data_wo_wt$Position_AA
  # creates the genotype table with associated fitness from the raw data
  empirical_fl <- cbind(rbind(numeric(nb_mut),
                              diag(nrow = nb_mut,
                                   ncol = nb_mut)),
                        data$s)
  colnames(empirical_fl) <- c(names_mutants, "fitness")
  # writes each fl to a csv file
  write.table(x = empirical_fl,
              file = paste0("data-raw/1-empirical_fl/", data_name, "_empirical_fl.csv"),
              row.names = FALSE, col.names = TRUE)
  return(empirical_fl)
}

#### Full raw empirical fl data ####
# Reads the raw data and extract the variable of interest
tbl <- read.csv("data-raw/1-empirical_fl/raw_data_DFE.csv")
tbl <- tbl %>%
  unite(Position_AA, c(Position, AA)) %>%
  transmute(Position_AA,
            Treatment,
            s) %>%
  group_by(Treatment)
list_tbl <-  tbl %>%
  group_split()
names_tbl <- group_keys(tbl)$Treatment

# Creates a list of empirical_fl per treatment
full_raw_empirical_fl <- sapply(names_tbl,
                                FUN = function(data_name) {
                                  raw2fl(list_tbl[[data_name]], data_name)
                                },
                                simplify = F)
names(full_raw_empirical_fl) <- names_tbl
# Saves the data as lazydata (easy loading)
use_data(full_raw_empirical_fl, overwrite = TRUE)


#### Description data ####
# library(tidyverse)
# ggplot() + geom_histogram(aes(x = empirical_fl[, ncol(empirical_fl)], y=..density..))
