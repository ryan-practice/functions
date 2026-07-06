#' Select the appended GxG motif from a J-gene sequence
#'
#' Chooses which `G.G` motif should be highlighted as the appended J-derived
#' motif in the CDR3/J-gene safety-check plot. The function first finds all GxG
#' motifs in the supplied J-gene amino acid sequence, then keeps only motifs
#' that begin at or after the first appended J-gene amino acid. If multiple
#' appended motifs are available, it prioritizes motifs immediately preceded by
#' phenylalanine or tryptophan, corresponding to common `FGxG` or `WGxG`
#' framework-region motifs.
#'
#' @param j_gene Character string containing the amino acid sequence of a J gene.
#' @param append_start_j Integer position in `j_gene` where the appended J-gene
#'   suffix begins. This is usually calculated as the J-gene overlap start
#'   position plus the length of the CDR3/J overlap.
#'
#' @return A tibble with one row describing the selected appended GxG motif, with
#'   the following columns:
#'   \describe{
#'     \item{motif_start}{Integer start position of the selected motif in
#'       `j_gene`, or `NA_integer_` if no appended motif is available.}
#'     \item{motif_end}{Integer end position of the selected motif in `j_gene`,
#'       or `NA_integer_` if no appended motif is available.}
#'     \item{motif}{The selected three-amino-acid motif, such as `"GKG"` or
#'       `"GPG"`, or `NA_character_` if no appended motif is available.}
#'     \item{n_gxg_total}{Total number of GxG motifs found anywhere in
#'       `j_gene`.}
#'     \item{gxg_status}{Human-readable status message for labeling the
#'       safety-check plot.}
#'   }
#'
#' @details
#' This function is designed for visual quality control of reconstructed TCR
#' chains. In the safety-check plot, the CDR3/J overlap is highlighted
#' separately from the GxG motif. The selected GxG motif should be the motif
#' contributed by the appended J-gene suffix, not an upstream motif that was
#' already present in the sequencing-derived CDR3.
#'
#' When multiple GxG motifs are present in a J gene, the function selects among
#' motifs in the appended suffix only. It then prefers motifs immediately
#' preceded by `"F"` or `"W"`, because these are more likely to represent the
#' expected `FGxG` or `WGxG` region used in the safety-check visualization.
#'
#' @examples
#' select_appended_gxg(
#'   j_gene = "XYTGANSKLTFGKGITLSVRP",
#'   append_start_j = 11
#' )
#'
#' select_appended_gxg(
#'   j_gene = "QETQYFGPGTRLLVL",
#'   append_start_j = 6
#' )
#'
#' select_appended_gxg(
#'   j_gene = "STDTQYFGPGTRLTVL",
#'   append_start_j = 6
#' )
#'
#' @importFrom dplyr mutate filter case_when n
#' @importFrom tibble tibble
#'
#' @export
select_appended_gxg <- function(j_gene, append_start_j) {
  motifs <- find_gxg_motifs(j_gene) %>%
    mutate(
      n_gxg_total = n(),
      is_appended = motif_start >= append_start_j,
      is_fwgxg_or_wgxg = preceding_aa %in% c("F", "W")
    )

  if (nrow(motifs) == 0) {
    return(tibble(
      motif_start = NA_integer_,
      motif_end = NA_integer_,
      motif = NA_character_,
      n_gxg_total = 0L,
      gxg_status = "0 GxG motifs found in J gene"
    ))
  }

  appended_candidates <- motifs %>%
    filter(is_appended)

  if (nrow(appended_candidates) == 0) {
    return(tibble(
      motif_start = NA_integer_,
      motif_end = NA_integer_,
      motif = NA_character_,
      n_gxg_total = nrow(motifs),
      gxg_status = paste0(
        nrow(motifs),
        " GxG motif(s) found, but none are fully in appended J suffix"
      )
    ))
  }

  selected <- appended_candidates[
    order(
      -as.integer(appended_candidates$is_fwgxg_or_wgxg),
      appended_candidates$motif_start
    )[1],
  ]

  tibble(
    motif_start = selected$motif_start,
    motif_end = selected$motif_end,
    motif = selected$motif,
    n_gxg_total = nrow(motifs),
    gxg_status = case_when(
      nrow(motifs) == 1 ~ paste0("1 GxG motif: highlighted ", selected$motif),
      nrow(motifs) > 1 ~ paste0(
        nrow(motifs),
        " GxG motifs; highlighted appended motif ",
        selected$motif
      ),
      TRUE ~ "GxG motif status unknown"
    )
  )
}
