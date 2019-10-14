#' Read and clean a simulation file produced by \code{\link[inferenceFitnessLandscape]{simulate_fl}}
#'
#' Import simulations, removes simulations with NA and keep only given list of parameters.
#'
#' @param file Name of the simulation file to import.
#' @param keep_param Vector of names corresponding to the names of parameters to keep for the inference.
#' i.e. the parameters varying in the simulations.
#' @inheritParams remove_wt_fitness
#' @param ... Options which will be passed to \code{\link[utils]{read.table}}.
#' @return A list. The first element is a matrix of the simulation raw results
#' (i.e. the fitness of the mutants). Each row correspond to one simulation. The
#' second element is a matrix of the parameters (only the ones named in \code{keep_param}).
#' Each row contains the parameters of the corresponding row of "simulation".
#' @export
read_clean_output <- function(file, keep_param = NULL, two_env = F, ...) {
  result <- as.matrix(read.table(file = file, ...))
  filtered_result <- filter_na(simulation = result[, grepl("^S", colnames(result)), drop = F], parameter = result[, !grepl("^S", colnames(result)), drop = F])
  if (length(filtered_result$idx_na_simulation) > 0) {
    warning(paste(length(filtered_result$idx_na_simulation), "row(s) with NAs removed from simulations and parameters", sep = " "))
  }
  if (is.null(keep_param)) {keep_param = 1:dim(filtered_result$filtered_parameter)[2]}
  list(simulation = remove_wt_fitness(filtered_result$filtered_simulation),
       parameter = filtered_result$filtered_parameter[, keep_param, drop = FALSE])
}

#' Removes lines with NA.
#'
#' Removes the lines with NA in \code{simulation} and remove the corresponding
#' rows in \code(parameter).
#'
#' @param simulation A matrix of real numbers.
#' @param parameter A matrix of real numbers with the same number of rows as
#' \code{simulation}.
#' @return A list containing as first element \code{filtered_simulation} (i.e.
#' without rows with NA), as second element \code{filtered_parameter} (i.e.
#' without the rows with the same index as the ones removed from \code{simulation})
#' and as third element the indexes of the rows of \code{simulation} with NA.
#' @examples
#' simulation <- rnorm(10 * 5)
#' simulation[sample(1:(10 * 5), 5)] <- NA
#' simulation <- matrix(simulation, 10, 5)
#' parameter <- matrix(rnorm(10 * 9), 10, 9)
#' filter_na(simulation = simulation, parameter = parameter)
filter_na <- function(simulation, parameter) {

  #checks : simulation, paramaters
  if (missing(simulation)){
    stop("'simulation' must be supplied.", call.=FALSE)
  }
  stopifnot(is.numeric(simulation), is.matrix(simulation))
  if (missing(parameter)){
    stop("'parameter' must be supplied.", call.=FALSE)
  }
  stopifnot(is.numeric(parameter), is.matrix(parameter), dim(parameter)[1] == dim(simulation)[1])

  idx_na <- which((apply(simulation, MARGIN = 1, FUN=function(s){anyNA(s)})))
  if (length(idx_na) > 0){
    filtered_simulation <- simulation[-idx_na, , drop = FALSE]
    filtered_parameter <- parameter[-idx_na, , drop = FALSE]
  } else {
    filtered_simulation <- simulation
    filtered_parameter <- parameter
  }
  list(filtered_simulation = filtered_simulation,
       filtered_parameter = filtered_parameter,
       idx_na_simulation = idx_na)
}

#' Removes lines with NA.
#'
#' Removes the first value of each simulations if it is the wild type. These values
#' are considered as wild type if they are all identical.
#'
#' @param sumstat A matrix of real numbers. Matrix of summary statistics with each
#' row containing the raw fitness of a simulation.
#' @param two_env Logical. TRUE if teh sumstat correspond to two environments and FALSE if it
#' correspond to one. Default = FALSE.
#' @return The same matrix as \code{sumstat} minus the first column (and the middle column if two envs) if
#' all the values are identical.
remove_wt_fitness <- function(sumstat, two_env = F) {
  to_remove = c()
  if (abs(max(sumstat[, 1]) - min(sumstat[, 1])) < .Machine$double.eps^0.5) {
    to_remove = c(to_remove, 1)
  } else {
    warning("The rows of the first column of sumstat must be all identical to be considered as the wild type.")
  }
  if (two_env) {
    col_wt_new_env = ncol(sumstat) / 2 + 1
    if (abs(max(sumstat[, col_wt_new_env]) - min(sumstat[, col_wt_new_env])) < .Machine$double.eps^0.5) {
      to_remove = c(to_remove, col_wt_new_env)
    } else {
      warning("The rows of the first column of the new env in sumstat must be all identical to be considered as the wild type.")
    }
  }
  sumstat[, -to_remove, drop = FALSE]
}
