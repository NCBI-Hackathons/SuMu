# functions for preparing clinical data for analysis
# (in general, supporting functions for prep_data.R)

#' Prep clinical data for analysis
#' @param data
#' @param formula
#' @param id
#' @import purrr assertthat
#' @returns data.frame filtered to non-missing observations for terms in the formula
prep_clinical_data <- function(
  data,
  formula,
  id
) {
  ## check input args
  if (is.character(formula))
    formula <- as.formula(formula)
  assertthat::is.string(id)
  assertthat::is.scalar(id)
  assertthat::assert_that(inherits(data, 'data.frame'))

  formula_fields <- as.character(attr(terms(formula), 'variables'))
  formula_fields <- purrr::keep(.x = formula_fields,
                                function(x) x %in% names(data))
  ## identify terms in the formula
  fields <- c(
    formula_fields,
    id
  )

  clin_df <- data %>%
    tidyr::drop_na_(vars = fields)

  clin_df
}


#' Normalize id inputs
#'
#' Depending on input format, users have the option
#' of providing a biomarker_matrix + clinical_data (where rows match)
#' or two data frames with an id variable.
#'
#' Regardless of input, this function creates an id variable
#' in each dataset & converts the input type to a data.frame.
#'
#' @param data clinical data
#' @param id (optional) id string
#' @param biomarker_data (optional) if biomarker data provided as long data frame
#' @param biomarker_matrix (optional) if biomarker data provided as a wide matrix
#' @return list of reconciled data, id, biomarker_data, & biomarker_matrix
normalize_data_inputs <- function(data,
                                  id,
                                  biomarker_data,
                                  biomarker_matrix) {
  ## check for valid inputs
  if (is.null(biomarker_data) && is.null(biomarker_matrix))
    stop('Either biomarker_data or biomarker_matrix must be provided.')
  if (!is.null(biomarker_data) && !is.null(biomarker_matrix))
    stop('Cannot supply both biomarker_data & biomarker_matrix. Please pick one.')
  if (!is.null(biomarker_data) && is.null(id))
    stop('id cannot be NULL if biomarker_data provided.')

  ## convert biomarker_matrix to a data frame, if necessary
  if (!inherits(biomarker_matrix, 'data.frame') &&
      inherits(biomarker_matrix, 'matrix'))
    biomarker_matrix <- as.data.frame(biomarker_matrix)

  ## populate fake id-colname if one doesn't exist already
  if (is.null(id)) {
    id_colname <- '_id'
    data[[id_colname]] <- rownames(data)
    if (!is.null(biomarker_matrix))
      biomarker_matrix[[id_colname]] <- rownames(biomarker_matrix)
  } else {
    id_colname <- id
    assertthat::assert_that(id_colname %in% names(data))
    if (!is.null(biomarker_matrix))
      assertthat::assert_that(id_colname %in% names(biomarker_matrix))
    if (!is.null(biomarker_data))
      assertthat::assert_that(id_colname %in% names(biomarker_data))
  }

  list(biomarker_data = biomarker_data,
       biomarker_matrix = biomarker_matrix,
       id = id_colname,
       data = data
       )
}

