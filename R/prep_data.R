# Functions for preparing data for analysis

#' Function to prepare clinical model.matrix for survival analysis
#' @param data dataframe containing clinical data (plus id variable)
#' @param formula formula relating clinical covariates to outcome
#' @param id string indicating ID column of data.frame
#' @import survival
#' @return list containins stan-data (N, M, x, y) & attributes
prep_clin_surv <- function(data, formula, id = NULL) {

}

#' Function to prepare genetic data matrix for analysis
#' @param biomarker_data data.frame or matrix containing biomarker data
#' @param id (optional) string identifying ID column
#' @return model.matrix with attributes
prep_biomarker_data <- function(biomarker_data, id) {

}
