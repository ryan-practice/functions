j_suffix_from_gxg <- function(aa_of_j, j_gene = NA, chain = NA) {
  # Expects format like:
  #   "....QY-FGPGTRLLVL"
  # and returns:
  #   "....QYF-GPGTRLLVL"

  if (is.na(aa_of_j) || !grepl("-", aa_of_j, fixed = TRUE)) {
    stop(
      paste0(
        "aa_of_j must contain exactly one dash separator. ",
        "Problem gene: ", j_gene, " chain: ", chain
      )
    )
  }

  parts <- strsplit(aa_of_j, "-", fixed = TRUE)[[1]]

  if (length(parts) != 2) {
    stop(
      paste0(
        "aa_of_j must contain exactly one dash separator. ",
        "Problem gene: ", j_gene, " chain: ", chain,
        ". Value was: ", aa_of_j
      )
    )
  }

  aa_1_old <- parts[1]
  aa_2_old <- parts[2]

  if (nchar(aa_2_old) < 4) {
    stop(
      paste0(
        "aa_of_j suffix is too short to move dash and check GXG. ",
        "Problem gene: ", j_gene, " chain: ", chain,
        ". Suffix was: ", aa_2_old
      )
    )
  }

  residue_moved_past <- substr(aa_2_old, 1, 1)

  aa_1_new <- paste0(aa_1_old, residue_moved_past)
  aa_2_new <- substr(aa_2_old, 2, nchar(aa_2_old))

  # X in GXG means any amino acid, so use regex G.G.
  if (!grepl("^G.G", aa_2_new)) {
    stop(
      paste0(
        "After moving the dash, J suffix does not start with GXG / G.G. ",
        "Problem gene: ", j_gene, " chain: ", chain,
        ". Moved residue: ", residue_moved_past,
        ". New suffix: ", aa_2_new,
        ". Original aa_of_j: ", aa_of_j
      )
    )
  }

  if (!residue_moved_past %in% c("F", "W")) {
    warning(
      paste0(
        "Residue moved past was not F or W. ",
        "Problem gene: ", j_gene, " chain: ", chain,
        ". Residue was: ", residue_moved_past,
        ". New aa_of_j: ", paste0(aa_1_new, "-", aa_2_new)
      )
    )
  }

  paste0(aa_1_new, "-", aa_2_new)
}
