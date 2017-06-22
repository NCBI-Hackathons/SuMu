# Functions (currently) for working with formula objects



#' Update clinical formula
#' @param biomarker_matrix biomarker data in wide format, as data.frame
#' @param formula given formula for association between clinical & outcome data
#' @param biomarker_placeholder string within clinical_subformula representing biomarkers
#' @param id name of id column
#' @return revised formula object
update_formula <- function(biomarker_matrix,
                           formula,
                           biomarker_placeholder = '__BIOM',
                           id) {
  biomarker_subformula <- get_biomarker_subformula(biomarker_matrix = biomarker_matrix,
                                                   id = id)
  revised_formula <- update_clinical_formula(biomarker_subformula = biomarker_subformula,
                                             clinical_formula = formula,
                                             biomarker_placeholder = biomarker_placeholder
                                             )
  revised_formula
}

#' Update clinical formula
#' @param biomarker_subformula subformula prepared for biomarker association
#' @param clinical_formula given formula for association between clinical & outcome data
#' @param biomarker_placeholder string within clinical_subformula representing biomarkers
#' @return revised formula object
update_clinical_formula <- function(biomarker_subformula,
                                    clinical_formula,
                                    biomarker_placeholder = '__BIOM') {
  str_formula <- stringr::str_c(
    as.character(clinical_formula)[2],
    as.character(clinical_formula)[3],
    sep = as.character(clinical_formula)[1])
  str_formula <- gsub(str_formula,
                      pattern = stringr::str_c('`',biomarker_placeholder, '`'),
                      replacement = biomarker_subformula)
  as.formula(str_formula)
}

#' Get biomarker subformula
#' @param biomarker_matrix biomarker data in wide format, as data.frame
#' @param id name of id column in biomarker_matrix
#' @return string representing biomarker-specific portion of formula
get_biomarker_subformula <- function(biomarker_matrix,
                                id) {

  biomarker_names <- colnames(biomarker_matrix)
  biomarker_names <- biomarker_names[biomarker_names != id]

  biomarker_subformula <- stringr::str_c('`',
                                    stringr::str_c(biomarker_names,
                                                   collapse = '` + `'),
                                    '`')

}


