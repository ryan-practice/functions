j_suffix_from_gxg <- function(j_aa, gene_name = NA) {
  gxg_hits <- stringr::str_locate_all(j_aa, "G.G")[[1]]

  if (nrow(gxg_hits) == 0) {
    stop(paste0("No G.G motif found for J gene ", gene_name))
  }

  # Your existing logic generally uses the last G.G motif.
  gxg_start <- gxg_hits[nrow(gxg_hits), "start"]

  substr(j_aa, gxg_start, nchar(j_aa))
}
