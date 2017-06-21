# Functions for fitting rstanarm models


#' Helper function to fit an rstanarm::stan_xxx model containing genetic features & clinical data
#'
#' @param data the dataframe containing clinical data
#' @param formula the formula describing association of clinical data to outcome
#' @param biomarker_data a dataframe containing genetic features / potential biomarkers.
#'           Should be provided in denomalized or long format, with one record per subject*marker
#'           An ID column is required with this data format.
#' @param biomarker_matrix a matrix containing genetic features / potential biomarkers.
#'           Should have dimensions NxG where N: number of subjects and G: number features
#'           Can have an ID column (if `id` param provided) or have rownames set to ID values.
#' @param biomarker_formula (optional) formula describing hierarchical structure of association
#'           for biomarker features. Only possible if data are provided in 'long' format.
#'           NOT YET IMPLEMENTED.
#' @param id (optional) name of id column.
#'          If provided, both the clinical & biomarker data/matrix should contain this column.
#'          If not, it is assumed that biomarker data/matrix has rownames, or is sorted in matched order
#'          as the clinical data.
#' @param prior prior on clinical covariates
#' @param biomarker_prior prior on biomarker features
#' @param stanfit_func name of rstanarm::stan_xx function
#' @param ... additional parameters passed to rstanarm fit object (iter, chains, etc)
#'
#' @return stanreg object
#' @import rstanarm
fit_rstanarm <- function(
  formula,
  data,
  biomarker_data = NULL,
  biomarker_matrix = NULL,
  biomarker_formula = NULL,
  stanfit_func,
  id = NULL,
  prior = NULL,
  biomarker_prior = NULL,
  family = NULL,
  ...
) {
  ## prepare genetic data matrix
  biodata <- prep_biomarker_data(biomarker_data = biomarker_data,
                                 biomarker_matrix = biomarker_matrix,
                                 biomarker_formula = biomarker_formula,
                                 data = data,
                                 id = id)

  ## join with clinical data for input to rstanarm

  ## prepare priors to input to rstanarm

  ## update user-supplied formula for rstanarm

  ## prepare call to rstanarm

  ## execute call

  ## format resulting object
}

#' Function to fit an rstanarm::stan_glm model containing genetic features & clinical data
#'
#' @param data the dataframe containing clinical data
#' @param formula the formula describing association of clinical data to outcome
#' @param biomarker_data the dataframe or matrix containing genetic
#'                    features / potential biomarkers. Should have an ID column
#'                    or have rownames set to ID values.
#' @param id (optional) name of id column. Should be common/shared between data & biomarker data
#' @param prior prior on clinical covariates
#' @param biomarker_prior prior on
#' @param ... additional parameters passed to rstanarm fit object (iter, chains, etc)
#'
#' @return stanreg object
#' @import rstanarm
#' @examples
#'
#' # TODO simulate genetic features
#' clin_df <- data(pbcSurv, package = 'rstanarm')
#' gen_df <- data.frame(...)
#'
#' fit <- fit_glm(data = clin_df,
#'                 formula = pfs_90d ~ age + gender + ...,
#'                 biomarker_data = gen_df,
#'                 family = binomial(),
#'                 ...
#'                )
#'
#' @export
fit_glm <- function(
  formula,
  data,
  biomarker_data,
  id = NULL,
  prior = NULL,
  biomarker_prior = NULL,
  family = NULL,
  ...
) {
  ## prepare function inputs to pass to fit_rstanarm()

  ## call fit_rstanarm
}


