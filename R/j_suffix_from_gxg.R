j_suffix_from_gxg <- function(aa_of_j, j_gene = NA, chain = NA) {
  # Handles both:
  #   10x / already-correct:  "....QYF-GPGTRLLVL"
  #   BD-style shifted split: "....QY-FGPGTRLLVL"
  #
  # Returns a corrected dashed string where the suffix starts with GXG:
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

  aa_1 <- parts[1]
  aa_2 <- parts[2]

  if (nchar(aa_1) < 1 || nchar(aa_2) < 3) {
    stop(
      paste0(
        "aa_of_j has insufficient sequence around dash. ",
        "Problem gene: ", j_gene, " chain: ", chain,
        ". Value was: ", aa_of_j
      )
    )
  }

  # Case 1: already correct, e.g. QYF-GPGTRLLVL
  if (grepl("^G.G", aa_2)) {
    residue_before_dash <- substr(aa_1, nchar(aa_1), nchar(aa_1))

    if (!residue_before_dash %in% c("F", "W")) {
      warning(
        paste0(
          "J suffix already starts with GXG, but residue before dash is not F or W. ",
          "Problem gene: ", j_gene, " chain: ", chain,
          ". Residue before dash was: ", residue_before_dash,
          ". aa_of_j: ", aa_of_j
        )
      )
    }

    return(aa_of_j)
  }

  # Case 2: BD-style shifted split, e.g. QY-FGPGTRLLVL
  # Move first residue of suffix to the left side of dash.
  if (nchar(aa_2) < 4) {
    stop(
      paste0(
        "aa_of_j suffix is too short to move dash and check GXG. ",
        "Problem gene: ", j_gene, " chain: ", chain,
        ". Suffix was: ", aa_2
      )
    )
  }

  residue_moved_past <- substr(aa_2, 1, 1)
  aa_1_new <- paste0(aa_1, residue_moved_past)
  aa_2_new <- substr(aa_2, 2, nchar(aa_2))

  if (!grepl("^G.G", aa_2_new)) {
    stop(
      paste0(
        "J suffix does not start with GXG either before or after moving dash. ",
        "Problem gene: ", j_gene, " chain: ", chain,
        ". Original aa_of_j: ", aa_of_j,
        ". Attempted moved residue: ", residue_moved_past,
        ". Attempted new suffix: ", aa_2_new
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
