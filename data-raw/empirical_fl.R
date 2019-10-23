library(tidyverse)

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
              file = paste0("data-raw/", data_name, "_empirical_fl.csv"),
              row.names = FALSE, col.names = TRUE)
  return(empirical_fl)
}

#### Full raw empirical fl data ####
tbl <- read.csv("data-raw/joint_df.csv")
tbl <- tbl %>%
  unite(Position_AA, c(Position, AA)) %>%
  transmute(Position_AA,
            Treatment,
            s) %>%
  group_by(Treatment)
list_tbl <-  tbl %>%
  group_split()
names_tbl <- group_keys(tbl)$Treatment
names(list_tbl) <- names_tbl

test <- raw2fl(list_tbl[[names_tbl[1]]], names_tbl[1])

full_raw_empirical_fl <- sapply(names_tbl,
                                FUN = function(data_name) {
                                  raw2fl(list_tbl[[data_name]], data_name)
                                },
                                simplify = F,
                                USE.NAMES = T)




use_data(full_raw_empirical_fl, overwrite = TRUE)

#### Separated empirical fl data ####
data(full_raw_empirical_fl)
d = names_tbl[1]



do.call(use_data, list(get(name_fl), overwrite = TRUE))
use_data(get(name_fl), overwrite = TRUE)

#### Description data ####
# library(tidyverse)
# ggplot() + geom_histogram(aes(x = empirical_fl[, ncol(empirical_fl)], y=..density..))
