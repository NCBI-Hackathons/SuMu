# Functions for importing external data


#' Retrieve clinical matrix for a TCGA cohort.
#'
#' @param cohort four character cohort abbreviation (SKCM, LUAD, BRCA, etc.)
#'
#' @return clinical matrix data frame
#'
#' @import readr
#'
#' @examples
#' get_tcga_clinical(cohort = "SKCM")
#'
#' @export
get_tcga_clinical <- function(cohort) {
  clin_matrix_url = paste0("https://tcga.xenahubs.net/download/TCGA.", cohort, ".sampleMap/", cohort, "_clinicalMatrix")
  clin = readr::read_tsv(clin_matrix_url)
  return(clin)
}

#' Retrieve somatic mutations for a TCGA cohort.
#'
#' @param cohort four character cohort abbreviation (SKCM, LUAD, BRCA, etc.)
#'
#' @return somatic mutation SNPs and small INDELs data frame
#'
#' @import readr
#'
#' @examples
#' get_tcga_somatic_mutations(cohort = "SKCM")
#'
#' @export
get_tcga_somatic_mutations <- function(cohort) {
  mut_matrix_url = paste0("https://tcga.xenahubs.net/download/TCGA.", cohort, ".sampleMap/mutation_broad")
  mut = readr::read_tsv(mut_matrix_url)
  return(mut)
}
