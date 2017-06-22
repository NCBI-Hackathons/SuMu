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
#'           for biomarker features. Required if data are provided in 'long' format.
#'           COMPLEX FORMULAS NOT YET IMPLEMENTED.
#' @param id (optional) name of id column.
#'          If provided, both the clinical & biomarker data/matrix should contain this column.
#'          If not, it is assumed that biomarker data/matrix has rownames, or is sorted in matched order
#'          as the clinical data.
#' @param prior prior on clinical covariates
#' @param biomarker_prior prior on biomarker features
#' @param stanfit_func name of rstanarm::stan_xx function
#' @param .fun (optional) function to use when summarizing values, if more than one exists per ID.
#'         defaults to NULL (do not summarize).
#' @param ... additional parameters passed to rstanarm fit object (iter, chains, etc)
#'
#' @return stanreg object
#' @import rstanarm
fit_rstanarm <- function(
  data,
  formula,
  biomarker_data = NULL,
  biomarker_matrix = NULL,
  biomarker_formula = NULL,
  stanfit_func,
  id = NULL,
  prior = rstanarm::hs_plus(),
  biomarker_prior = prior,
  family = NULL,
  .fun = NULL,
  biomarker_placeholder = '__BIOM',
  sparse = TRUE,
  ...
) {

  #if (!assertthat::are_equal(prior, biomarker_prior))
  #  stop('custom biomarker priors not yet implemented.')

  ## prepare genetic data matrix
  prepped <- prep_data(biomarker_data = biomarker_data,
                       biomarker_matrix = biomarker_matrix,
                       biomarker_formula = biomarker_formula,
                       data = data,
                       formula = formula,
                       .fun = .fun,
                       id = id)
  clindata <- prepped$data
  id_colname <- prepped$id
  biomarker_data <- prepped$biomarker_data
  biomarker_matrix <- prepped$biomarker_matrix
  biomarker_formula <- prepped$biomarker_formula

  ## join biomarker with clinical data for input to rstanarm
  analysis_df <- biomarker_matrix %>%
    dplyr::inner_join(clindata, by = id_colname)

  ## update user-supplied formula for rstanarm
  revised_formula <- update_formula(formula = formula,
                                    biomarker_matrix = biomarker_matrix,
                                    id = id,
                                    biomarker_placeholder = biomarker_placeholder
                                    )

  ## execute call to rstanarm
  stanfit <- stanfit_func(
    data = analysis_df,
    formula = revised_formula,
    sparse = sparse,
    prior = prior,
    ...
  )

  ## format resulting object

  format_results(stanfit)
}

#' Function to fit an rstanarm::stan_glm model containing genetic features & clinical data
#'
#' @param data the dataframe containing clinical data
#' @param formula the formula describing association of clinical data to outcome
#' @param biomarker_data the dataframe or matrix containing genetic
#'                    features / potential biomarkers. Should have an ID column
#'                    or have rownames set to ID values.
#' @param biomarker_formula
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
  data,
  formula,
  biomarker_data = NULL,
  biomarker_matrix = NULL,
  biomarker_formula = NULL,
  id = NULL,
  prior = rstanarm::hs_plus(),
  biomarker_prior = prior,
  family = NULL,
  .fun = NULL,
  biomarker_placeholder = '__BIOM',
  sparse = TRUE,
  ...
) {
  updatemc <- mc <- match.call(expand.dots = T)
  updatemc[[1]] <- 'fit_rstanarm' # change function call
  updatemc$fit_func <- rstanarm::stan_glm
  eval(updatemc, parent.frame())
}


