# Functions for importing external (only TCGA so far) data.


#' Retrieve clinical matrix for a TCGA cohort
#'
#' @param cohort four character cohort abbreviation (SKCM, LUAD, BRCA, etc.)
#'
#' @return Clinical matrix data frame.
#'
#' @import readr
#'
#' @examples
#' clinical = get_tcga_clinical(cohort = "SKCM")
#'
#' @export
get_tcga_clinical <- function(cohort) {
  clin_matrix_url = paste0("https://tcga.xenahubs.net/download/TCGA.", cohort, ".sampleMap/", cohort, "_clinicalMatrix")
  clin = readr::read_tsv(clin_matrix_url)
  # fix columns starting with underscores (causes issues with some tools)
  colnames(clin) = gsub("^_", "", colnames(clin))
  return(clin)
}

#' Retrieve somatic mutations for a TCGA cohort
#'
#' @param cohort four character cohort abbreviation ("SKCM", "LUAD", "BRCA", etc.)
#' @param call_type type of call ("broad", "curated_broad", "curated_wustl", "wust", "bcm", etc.)
#'
#' @return Data frame of somatic mutation SNPs and small INDELs.
#'
#' @import readr
#'
#' @examples
#' mutations = get_tcga_somatic_mutations(cohort = "SKCM")
#'
#' @export
get_tcga_somatic_mutations <- function(cohort, call_type = "broad") {
  mut_matrix_url = paste0("https://tcga.xenahubs.net/download/TCGA.", cohort, ".sampleMap/mutation_", call_type)
  mut = readr::read_tsv(mut_matrix_url)
  return(mut)
}

#' Retrieve gene expression data for a TCGA cohort
#'
#' TCGA skin cutaneous melanoma (SKCM) gene expression by RNAseq, mean-normalized (per gene) across all TCGA cohorts.
#' Values in this dataset are generated at UCSC by combining "gene expression RNAseq" values of all TCGA cohorts, values are then mean-centered per gene, then extracting the converted data only belongs to the this cohort.
#' The values are log2 gene expression values mean-normalized (per gene) across all TCGA cohorts.
#'
#' @param cohort four character cohort abbreviation (SKCM, LUAD, BRCA, etc.)
#'
#' @return Data frame of gene expression values.
#'
#' @import readr
#'
#' @examples
#' exp = get_tcga_gene_expression(cohort = "SKCM")
#'
#' @export
get_tcga_gene_expression <- function(cohort) {
  exp_url = paste0("https://tcga.xenahubs.net/download/TCGA.", cohort, ".sampleMap/HiSeqV2_PANCAN")
  exp = readr::read_tsv(exp_url)
  return(exp)
}

#' Retrieve methylation data (Illumina Infinium HumanMethylation450 platform) for a TCGA cohort
#'
#' DNA methylation values, described as beta values, are recorded for each array probe in each sample via BeadStudio software.
#' DNA methylation beta values are continuous variables between 0 and 1, representing the ratio of the intensity of the methylated bead type to the combined locus intensity.
#'
#' @param cohort four character cohort abbreviation (SKCM, LUAD, BRCA, etc.)
#'
#' @return Data frame of beta values.
#'
#' @import readr
#'
#' @examples
#' meth = get_tcga_methylation_450k(cohort = "SKCM")
#'
#' @export
get_tcga_methylation_450k <- function(cohort) {
  meth_url = paste0("https://tcga.xenahubs.net/download/TCGA.", cohort, ".sampleMap/HumanMethylation450")
  meth = readr::read_tsv(meth_url)
  return(meth)
}

#' Retrieve gene-level copy number data for a TCGA cohort
#'
#' Copy number profile was measured experimentally using whole genome microarray at a TCGA genome characterization center.
#' Subsequently, GISTIC2 method was applied using the TCGA FIREHOSE pipeline to produce gene-level copy number estimates.
#' GISTIC2 further thresholded the estimated values to -2,-1,0,1,2, representing homozygous deletion, single copy deletion, diploid normal copy, low-level copy number amplification, or high-level copy number amplification.
#'
#' @param cohort four character cohort abbreviation (SKCM, LUAD, BRCA, etc.)
#'
#' @return Data frame of copy number values (gene-level).
#'
#' @import readr
#'
#' @examples
#' cnv_gene = get_tcga_copy_number_gene(cohort = "SKCM")
#'
#' @export
get_tcga_copy_number_gene <- function(cohort) {
  cnv_url = paste0("https://tcga.xenahubs.net/download/TCGA.", cohort, ".sampleMap/Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes")
  cnv = readr::read_tsv(cnv_url)
  return(cnv)
}
