# Functions for preparing data for survival analysis
# (in general, supporting functions for prep_data.R)

#' Function to prepare clinical model.matrix for survival analysis
#' @param data dataframe containing clinical data (plus id variable)
#' @param formula formula relating clinical covariates to outcome
#' @param id string indicating ID column of data.frame
#' @import survival
#' @return list containins stan-data (N, M, x, y) & attributes
prep_clin_surv <- function(data, formula, id = NULL) {

}
