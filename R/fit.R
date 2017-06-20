# Functions for fitting models & preparing data


#' Function to fit a model
#'
#'
#' @param dataframe the dataframe containing clinical data
#' @param iter number of iterations
#' @returns fitting stan model object
#' @import rstanarm magrittr
#' @examples
#'
#' data(mtcars)
#' fit(dataframe = mtcars)
#'
#'
#'
#' @export
fit <- function (...) {
  rstanarm::stan_glm(...)
}
