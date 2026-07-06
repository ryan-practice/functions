#' Infer J-gene type for CDR3 safety-check plotting
#'
#' Classifies a J gene name as alpha-chain TRAJ, beta-chain TRBJ, or a generic
#' J gene label. This helper is used by the CDR3/J-gene safety-check plot to
#' choose accurate row labels, such as "TRAJ gene" versus "TRBJ gene" and
#' "Final alpha chain" versus "Final beta chain".
#'
#' @param j_gene_name Character vector of J gene names, such as `"TRAJ42"`,
#'   `"TRAJ42*01"`, `"TRBJ2-5"`, or `"TRBJ2-5*01"`.
#'
#' @return A character vector with one value per input gene name. Values are
#'   `"TRAJ"` for alpha-chain J genes, `"TRBJ"` for beta-chain J genes, and
#'   `"J"` when the gene type cannot be inferred from the prefix.
#'
#' @examples
#' infer_j_type("TRAJ42")
#' infer_j_type("TRBJ2-5")
#' infer_j_type(c("TRAJ16", "TRBJ1-2", "UNKNOWN"))
#'
#' @importFrom dplyr case_when
#' @importFrom stringr str_starts
#'
#' @export
infer_j_type <- function(j_gene_name) {
  j_gene_name_clean <- toupper(j_gene_name)

  case_when(
    str_starts(j_gene_name_clean, "TRAJ") ~ "TRAJ",
    str_starts(j_gene_name_clean, "TRBJ") ~ "TRBJ",
    TRUE ~ "J"
  )
}





