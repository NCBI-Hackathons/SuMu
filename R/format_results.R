# Functions to format fit results, depending on object type

#' Format fit result for StanReg object
#' @param fit StanReg fit object
#' @param biodata matrix (with attributes) resulting from prep_biomarker_data
#' @param formula formula relating clinical data to outcomes
#' @param data clinical dataframe
#' @param stan_data dataframe passed to rstanarm fit function
#' @param stan_formula formula passed to rstanarm fit function
#' @param id string name of id column in data & biodata
#' @return structure supporting posterior summary of results
#' @export
format_results.stanreg <- function(fit) {

}

#' Format fit result for Stanfit object
#' @param fit Stanfit object
#' @param biodata matrix (with attributes) resulting from prep_biomarker_data
#' @param formula formula relating clinical data to outcomes
#' @param data clinical dataframe
#' @param stan_data dataframe passed to stan fit function
#' @param stan_formula formula passed to stan fit function
#' @param id string name of id column in data & biodata
#' @return structure supporting posterior summary of results
#' @export
format_results.stanfit <- function(fit) {

}
