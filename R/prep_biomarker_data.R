# Functions for preparing biomarker data for analysis


#' Function to prepare biomarker data matrix for analysis
#' @param data the dataframe containing clinical data
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
#' @return model.matrix with attributes
prep_biomarker_data <- function(data,
                                biomarker_data = NULL,
                                biomarker_matrix = NULL,
                                biomarker_formula = NULL,
                                id = NULL
) {
  ## check for valid inputs
  if (is.null(biomarker_data) && is.null(biomarker_matrix))
    stop('Either biomarker_data or biomarker_matrix must be provided.')
  if (!is.null(biomarker_data) && !is.null(biomarker_matrix))
    stop('Cannot supply both biomarker_data & biomarker_matrix. Please pick one.')
  if (!is.null(biomarker_data) && is.null(id))
    stop('id cannot be NULL if biomarker_data provided.')
  if (!is.null(biomarker_formula))
    stop('biomarker_formula option not yet implemented.')

  ## populate biomarker_data from biomarker_matrix
  if (!is.null(biomarker_matrix)) {
    biomarker_data <- convert_biom_matrix_to_data(biomarker_matrix, id = id)
    id_colname <- attr(biomarker_data, 'id_colname')
    ## check if ID in clinical data
    if (!id_colname %in% names(data)) {
      ## assume rownames / data sort order of biomarker & clinical data match
      data[id_colname] <- rownames(data)
    }
    biomarker_formula <- attr(biomarker_data, 'biomarker_formula')
  }
  if (!is.null(biomarker_data)) {
    biomarker_matrix <- convert_biom_data_to_matrix(biomarker_data, id = id)
    id_colname <- attr(biomarker_matrix, 'id_colname')
    ## biomarker_formula already given for biomarker_data
  }

  ## filter to non-missing records in clinical data
  ids_in_biomarker_data <- biomarker_data %>%
    dplyr::distinct_(id_colname) %>%
    dplyr::select_(id_colname)
  ids_in_clinical_data <- data %>%
    dplyr::distinct_(id_colname) %>%
    dplyr::select_(id_colname)
  ids_in_both <- ids_in_biomarker_data %>%
    dplyr::semi_join(ids_in_clinical_data, by = id_colname)

  ## filter biomarker_data & biomarker_matrix
  biomarker_data_filtered <- biomarker_data %>%
    dplyr::semi_join(ids_in_both)
  biomarker_matrix_filtered <- biomarker_matrix %>%
    dplyr::semi_join(ids_in_both)
  data_filtered <- data %>%
    dplyr::semi_join(ids_in_both)

  ## return biomarker_matrix
  structure(biomarker_matrix)
}

#' Helper function to convert biomarker_matrix (wide-format) to data (long-format)
#' @param biomarker_matrix a matrix containing genetic features / potential biomarkers.
#'           Should have dimensions NxG where N: number of subjects and G: number features
#'           Can have an ID column (if `id` param provided) or have rownames set to ID values.
#' @param id (optional) name of id column.
#'          If provided, both the clinical & biomarker data/matrix should contain this column.
#'          If not, it is assumed that biomarker data/matrix has rownames, or is sorted in matched order
#'          as the clinical data.
#' @import assertthat tidyr
#' @return long-format data.frame with attributes id_colname & biomarker_formula
convert_biom_matrix_to_data <- function(biomarker_matrix, id) {
  ## populate temporary id column if one not given
  if (is.null(id)) {
    id_colname <- '_id' # name that is unlikely to conflict with existing names
    if (id_colname %in% names(biomarker_matrix))
      id_colname <- '_temp_id'
    biomarker_matrix[id_colname] <- rownames(biomarker_matrix)
  } else {
    id_colname <- id
    assertthat::has_name(biomarker_matrix, id_colname)
  }

  ## reshape data to long-format
  biomarker_data <- biomarker_matrix %>%
    tidyr::gather('biomarker_name', 'biomarker_value', -matches(id_colname))
  biomarker_formula <- formula(biomarker_value ~ biomarker_name)
  structure(biomarker_data, biomarker_formula = biomarker_formula, id_colname = id_colname)
}

#' Helper function to convert biomarker_data (long-format) to biomarker_matrix (wide-format)
#' @param biomarker_data a dataframe containing genetic features / potential biomarkers.
#'           Should be provided in denomalized or long format, with one record per subject*marker
#'           An ID column is required with this data format.
#' @param biomarker_formula (required) formula describing hierarchical structure of association
#'           for biomarker features.
#'           COMPLEX FORMULAS NOT YET IMPLEMENTED.
#' @param id (required) name of id column.
#' @import assertthat tidyr dplyr
#' @return wide-format data.frame with attributes
convert_biom_data_to_matrix <- function(biomarker_data,
                                        id,
                                        biomarker_formula
) {
  ## confirm id input
  assertthat::is.string(id)
  id_colname <- id
  assertthat::has_name(biomarker_data, id_colname)
  ## confirm biomarker_formula
  if (!inherits(biomarker_formula, 'formula'))
    stop('biomarker_formula must be of type `formula`')
  if (!inherits(biomarker_data, 'data.frame'))
    stop('biomarker_data must be of type `data.frame`')

  ## extract value (LHS) & description (RHS) from biomarker_formula
  terms <- terms(biomarker_formula)
  assertthat::assert_that(length(terms) == 3)
  assertthat::assert_that(terms[[1]] == '~')
  value_colname <- terms[[2]]
  biomarker_colname <- terms[[3]]

  ## construct matrix
  biomarker_matrix <- biomarker_data %>%
    dplyr::distinct(one_of(value_colname, id_colname, biomarker_colname)) %>%
    tidyr::spread_(key = biomarker_colname, value = value_colname, fill = 0)

  structure(biomarker_matrix, id_colname = id_colname)
}
