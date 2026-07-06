#' Plot a CDR3/J-gene reconstruction safety check
#'
#' Creates a visual safety-check plot for a single reconstructed TCR chain. The
#' plot aligns the sequencing-derived CDR3 amino acid sequence, the final
#' reconstructed chain, and the corresponding J-gene amino acid sequence. The
#' CDR3/J-gene overlap is highlighted in green, and the appended J-derived GxG
#' motif is highlighted in yellow.
#'
#' This function is intended to help confirm that amino acids already present in
#' the sequencing-derived CDR3 are not duplicated when the J-gene suffix is
#' appended. It supports both alpha-chain TRAJ inputs and beta-chain TRBJ inputs;
#' plot labels are inferred automatically from `j_gene_name`.
#'
#' @param cdr3 Character string containing one sequencing-derived CDR3 amino acid
#'   sequence.
#' @param j_gene_name Character string containing the J-gene name, such as
#'   `"TRAJ42"`, `"TRAJ42*01"`, `"TRBJ2-5"`, or `"TRBJ2-5*01"`.
#' @param j_gene Character string containing the amino acid sequence of the
#'   corresponding J gene.
#' @param final_chain Character string containing the final reconstructed alpha
#'   or beta chain amino acid sequence.
#' @param plot_title Optional character string to use as a custom plot title. In
#'   the current plotting code, the `title = plot_title` line in `labs()` is
#'   commented out; uncomment that line to display the title.
#'
#' @return A `ggplot` object showing the aligned CDR3, final reconstructed chain,
#'   and J-gene sequence. Green tiles indicate the CDR3/J overlap. Yellow tiles
#'   indicate the selected appended J-derived GxG motif.
#'
#' @details
#' The function first identifies the longest suffix of `cdr3` that appears in
#' `j_gene`. This overlap is used to align the J-gene row beneath the matching
#' region of the CDR3 row. The portion of the J gene after this overlap is
#' treated as the appended J-gene suffix.
#'
#' The GxG motif highlight is chosen from motifs that occur in the appended
#' J-gene suffix, rather than from motifs that are already present in the CDR3
#' overlap. For J genes with multiple GxG motifs, the function relies on
#' `select_appended_gxg()` to choose the motif most consistent with the appended
#' J-derived framework region.
#'
#' Special TRAJ cases with known non-standard GxG counts are labeled in the plot:
#' TRAJ16 and TRAJ61 are expected to have zero GxG motifs, while TRAJ10, TRAJ25,
#' TRAJ28, TRAJ35, and TRAJ45 are expected to have two GxG motifs.
#'
#' @examples
#' plot_cdr3_j_safety_check(
#'   cdr3 = "CAMREGTGANSKLT",
#'   j_gene_name = "TRAJ42",
#'   j_gene = "XYTGANSKLTFGKGITLSVRP",
#'   final_chain = "CAMREGTGANSKLTFGKGITLSVRP"
#' )
#'
#' plot_cdr3_j_safety_check(
#'   cdr3 = "ASSQSPGGTQY",
#'   j_gene_name = "TRBJ2-5",
#'   j_gene = "QETQYFGPGTRLLVL",
#'   final_chain = "ASSQSPGGTQYFGPGTRLLVL"
#' )
#'
#' @importFrom dplyr bind_rows mutate case_when
#' @importFrom tibble tibble
#' @importFrom stringr str_locate fixed
#' @importFrom ggplot2 ggplot aes geom_tile geom_text geom_label
#'   scale_fill_manual labs theme_minimal theme element_blank element_text
#'
#' @export
plot_cdr3_j_safety_check <- function(
    cdr3,
    j_gene_name,
    j_gene,
    final_chain,
    plot_title = NULL
) {
  ...
}
