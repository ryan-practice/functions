#' Find the longest CDR3/J-gene amino acid overlap
#'
#' Identifies the longest suffix of a sequencing-derived CDR3 amino acid
#' sequence that is also present in the supplied J-gene amino acid sequence.
#' This helper is used by the CDR3/J-gene safety-check plot to determine which
#' amino acids are already represented in the CDR3 and therefore should not be
#' duplicated when appending the J-gene suffix.
#'
#' @param cdr3 Character string containing one CDR3 amino acid sequence from
#'   sequencing.
#' @param j_gene Character string containing the amino acid sequence of the
#'   corresponding J gene.
#'
#' @return A character string containing the longest CDR3 suffix found in the
#'   J-gene sequence. Returns `NA_character_` if no suffix of `cdr3` is found
#'   in `j_gene`.
#'
#' @details
#' The function works by generating all suffixes of the CDR3 sequence, testing
#' which suffixes occur in the J-gene sequence, and returning the longest match.
#' In the safety-check plot, this overlap is highlighted to confirm that the
#' reconstructed chain keeps the CDR3 sequence intact while only appending the
#' non-overlapping J-gene sequence.
#'
#' @examples
#' find_longest_cdr3_j_overlap(
#'   cdr3 = "CAMREGTGANSKLT",
#'   j_gene = "XYTGANSKLTFGKGITLSVRP"
#' )
#'
#' find_longest_cdr3_j_overlap(
#'   cdr3 = "ASSQSPGGTQY",
#'   j_gene = "QETQYFGPGTRLLVL"
#' )
#'
#' @importFrom purrr map_chr
#' @importFrom stringr str_detect fixed
#'
#' @export
find_longest_cdr3_j_overlap <- function(cdr3, j_gene) {
  cdr3_len <- nchar(cdr3)

  candidates <- map_chr(seq_len(cdr3_len), function(i) {
    substr(cdr3, i, cdr3_len)
  })

  candidates <- candidates[str_detect(j_gene, fixed(candidates))]

  if (length(candidates) == 0) {
    return(NA_character_)
  }

  candidates[which.max(nchar(candidates))]
}
