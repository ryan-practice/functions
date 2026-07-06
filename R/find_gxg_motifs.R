#' Find GxG motifs in a J-gene amino acid sequence
#'
#' Locates all amino acid motifs matching the pattern `G.G` in a J-gene
#' sequence, where the middle position can be any amino acid. This helper is
#' used by the CDR3/J-gene safety-check plot to identify candidate GxG motifs
#' and determine which one belongs to the appended J-gene suffix.
#'
#' @param j_gene Character string containing the amino acid sequence of a J gene.
#'
#' @return A tibble with one row per detected GxG motif and the following
#'   columns:
#'   \describe{
#'     \item{motif_start}{Integer position where the motif starts in `j_gene`.}
#'     \item{motif_end}{Integer position where the motif ends in `j_gene`.}
#'     \item{motif}{The matched three-amino-acid motif, such as `"GKG"` or `"GPG"`.}
#'     \item{preceding_aa}{The amino acid immediately before the motif, or
#'       `NA_character_` if the motif starts at position 1.}
#'   }
#'   If no GxG motifs are found, returns an empty tibble with the same columns.
#'
#' @details
#' This function uses a lookahead regular expression so overlapping GxG motifs
#' can be detected. For example, `"GGGG"` contains two overlapping GxG motifs:
#' positions 1-3 and 2-4.
#'
#' In the CDR3 safety-check workflow, this function only identifies candidate
#' GxG motifs. A downstream helper should decide which motif is actually in the
#' appended J-gene suffix.
#'
#' @examples
#' find_gxg_motifs("XYTGANSKLTFGKGITLSVRP")
#'
#' find_gxg_motifs("QETQYFGPGTRLLVL")
#'
#' find_gxg_motifs("GGGG")
#'
#' @importFrom tibble tibble
#' @importFrom dplyr if_else
#'
#' @export
find_gxg_motifs <- function(j_gene) {
  starts <- gregexpr("(?=G.G)", j_gene, perl = TRUE)[[1]]

  if (length(starts) == 1 && starts[1] == -1) {
    return(tibble(
      motif_start = integer(),
      motif_end = integer(),
      motif = character(),
      preceding_aa = character()
    ))
  }

  tibble(
    motif_start = as.integer(starts),
    motif_end = as.integer(starts) + 2L,
    motif = substr(j_gene, motif_start, motif_end),
    preceding_aa = if_else(
      motif_start > 1,
      substr(j_gene, motif_start - 1, motif_start - 1),
      NA_character_
    )
  )
}
