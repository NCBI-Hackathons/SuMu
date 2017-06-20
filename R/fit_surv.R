# Functions for fitting Survival models using rstan::stan



#' Function to fit a survival model containing genetic features & clinical data
#'
#' @param data the dataframe containing clinical data
#' @param formula the formula describing association of clinical data to outcome
#' @param biomarker_data the dataframe or matrix containing genetic
#'                    features / potential biomarkers. Should have an ID column
#'                    or have rownames set to ID values.
#' @param id (optional) name of id column. Should be common/shared between data & biomarker data
#' @param prior prior on clinical covariates
#' @param biomarker_prior prior on
#' @param ... additional parameters passed to stan fit object (iter, chains, etc)
#'
#' @return stanfit object
#' @import rstan
#' @examples
#'
#' # TODO simulate genetic features
#' clin_df <- data(pbcSurv, package = 'rstanarm')
#' gen_df <- data.frame(...)
#'
#' fit <- fit_surv(data = clin_df,
#'                 formula = Surv(event, time) ~ age + gender + ...,
#'                 biomarker_data = gen_df,
#'                 ...
#'                )
#'
#' @export
fit_surv <- function(
  formula,
  data,
  biomarker_data,
  id = NULL,
  prior = NULL,
  biomarker_prior = NULL,
  ...
) {
  ## prepare clinical model.matrix

  ## prepare genetic data matrix

  ## load survival model stan code into memory

  ## fit survival model using rstan::stan
  rstan::stan(...)

  ## format resulting object
}

