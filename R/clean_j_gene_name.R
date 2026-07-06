#' Clean J-gene name for CDR3 safety-check logic
#'
#' Standardizes a J-gene name by converting it to uppercase and removing any
#' allele suffix after an asterisk. This helper is used by the CDR3/J-gene
#' safety-check plot when matching special-case J genes, such as TRAJ16,
#' TRAJ61, or TRBJ genes with allele-level names.
#'
#' @param j_gene_name Character vector of J-gene names, such as `"TRAJ42"`,
#'   `"TRAJ42*01"`, `"TRBJ2-5"`, or `"TRBJ2-5*01"`.
#'
#' @return A character vector of cleaned J-gene names with uppercase formatting
#'   and allele suffixes removed. For example, `"TRAJ42*01"` becomes `"TRAJ42"`.
#'
#' @examples
#' clean_j_gene_name("TRAJ42*01")
#' clean_j_gene_name("trbj2-5*01")
#' clean_j_gene_name(c("TRAJ16*01", "TRAJ61", "TRBJ2-5*01"))
#'
#' @importFrom stringr str_remove
#'
#' @export
clean_j_gene_name <- function(j_gene_name) {
  str_remove(toupper(j_gene_name), "\\*.*$")
}
