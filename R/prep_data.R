# Functions for preparing data for analysis


#' Function to prepare biomarker data matrix for analysis
#'
#' This function reshapes the biomarker data, normalizes
#' various input formats accepted, and filters to records in the clinical
#' dataset. It does not merge biomarker & clinical data.
#'
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
#' @param .fun (optional) function to use when summarizing values, if more than one exists per ID.
#'         defaults to NULL (do not summarize).
#' @return list containing prepared
#'       biomarker_data, biomarker_matrix, biomarker_formula
#'       data (clinical data), and id
#' @export
prep_data <- function(data,
                      formula,
                      biomarker_data = NULL,
                      biomarker_matrix = NULL,
                      biomarker_formula = NULL,
                      id = NULL,
                      .fun = NULL
) {
  ## check for valid inputs
  if (is.null(biomarker_data) && is.null(biomarker_matrix))
    stop('Either biomarker_data or biomarker_matrix must be provided.')
  if (!is.null(biomarker_data) && !is.null(biomarker_matrix))
    stop('Cannot supply both biomarker_data & biomarker_matrix. Please pick one.')
  if (!is.null(biomarker_data) && is.null(id))
    stop('id cannot be NULL if biomarker_data provided.')
  if (!is.null(biomarker_data) && is.null(biomarker_formula))
    stop('biomarker_formula is required with biomarker_data - should be in the form of `biomarker_value ~ biomarker_name`.')

  ## normalize data inputs
  normalized <- normalize_data_inputs(biomarker_data = biomarker_data,
                                      biomarker_matrix = biomarker_matrix,
                                      id = id,
                                      data = data)
  biomarker_data <- normalized$biomarker_data
  biomarker_matrix <- normalized$biomarker_matrix
  data <- normalized$data
  id_colname <- normalized$id

  if (!is.null(biomarker_data))
    assertthat::assert_that(id_colname %in% names(biomarker_data))
  if (!is.null(biomarker_matrix))
    assertthat::assert_that(id_colname %in% names(biomarker_matrix))
  assertthat::assert_that(id_colname %in% names(data))

  ## prepare clinical data, given clinical formula
  clindata <- prep_clinical_data(data = data,
                                 id = id_colname,
                                 formula = formula
                                 )

  ## prepare biomarker data, given filtered clinical data
  biomarker_data <- prep_biomarker_data(
    data = clindata,
    id = id_colname,
    biomarker_matrix = biomarker_matrix,
    biomarker_data = biomarker_data,
    biomarker_formula = biomarker_formula,
    .fun = .fun
  )

  biomarker_data
}

