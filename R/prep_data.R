# Functions for preparing data for analysis

#' Function to prepare clinical model.matrix for survival analysis
#' @param data
#' @param formula
#' @param id string indicating ID column of data.frame
#' @import survival
#' @return list containins stan-data (N, M, x, y) & attributes
prep_clin_surv <- function(data, formula, id = NULL) {

}

#' Function to prepare genetic data matrix for analysis
#' @param biomarker_data
#' @param id
#' @return model.matrix with attributes
prep_biomarker_data <- function(biomarker_data, id) {

}
