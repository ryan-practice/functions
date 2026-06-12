#' Tile an amino acid sequence into overlapping peptide tiles
#'
#' Creates peptide tiles from a single amino acid sequence using two tiling tracks:
#' one starting at position 1 and one starting at `offset + 1`. If the sequence
#' ends with partial-length tiles, the function marks those partial tiles as
#' `gets_replaced = "yes"` and adds a full-length terminal replacement tile.
#'
#' For sequences shorter than `tile_size`, the function returns the original
#' short sequence and a repeated/padded version expanded to `tile_size`.
#'
#' @param aa_seq Character string. Amino acid sequence to tile.
#' @param identifier Character string. Identifier to assign to the output
#'   `cluster_id` column.
#' @param tile_size Integer. Desired peptide tile length.
#' @param offset Integer. Offset for the second tiling track. For example, with
#'   `tile_size = 90` and `offset = 45`, the second track starts at amino acid 46.
#' @param keep_partial Logical. Currently unused. Reserved for future behavior
#'   controlling whether partial tiles are retained.
#' @param return_df Logical. Currently unused. Reserved for future output options.
#' @param metadata Optional metadata to append to every output tile. May be
#'   `NULL`, a named list, a named vector, or a one-row data frame. Metadata
#'   column names must not duplicate existing output column names.
#'
#' @return A data frame with one row per tile. Core columns include:
#' \describe{
#'   \item{cluster_id}{Identifier supplied by `identifier`.}
#'   \item{start}{Start coordinate of the tile in the input sequence.}
#'   \item{end}{End coordinate of the tile in the input sequence.}
#'   \item{tile_seq}{Amino acid sequence of the tile.}
#'   \item{tile_length}{Length of `tile_seq`.}
#'   \item{padded}{Whether the tile was added as a padded/replacement tile.}
#'   \item{n_pad_aa}{Number of amino acids added or effectively recovered.}
#'   \item{gets_replaced}{Whether this partial tile is intended to be replaced
#'     by a full-length terminal tile.}
#' }
#'
#' If `metadata` is supplied, its fields are appended as additional columns.
#'
#' @examples
#' tile_aa(
#'   aa_seq = "MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQV",
#'   identifier = "seq_001",
#'   tile_size = 15,
#'   offset = 7
#' )
#'
#' tile_aa(
#'   aa_seq = "MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQV",
#'   identifier = "seq_001",
#'   tile_size = 15,
#'   offset = 7,
#'   metadata = list(
#'     source_file = "HERV_orfs_aa_sequences_v3.fasta",
#'     protein_type = "ORF",
#'     domain_hit = "Gag"
#'   )
#' )
#'
#' @export
tile_aa <- function(
    aa_seq,
    identifier,
    tile_size,
    offset,
    metadata = NULL
){
  full_tiling <- NULL
  aa_seq <- toupper(aa_seq)
  add_metadata <- function(df, metadata) {
    if (is.null(metadata)) {
      return(df)
    }

    if (is.vector(metadata) && !is.list(metadata)) {
      metadata <- as.list(metadata)
    }

    if (is.list(metadata) && !is.data.frame(metadata)) {
      metadata <- as.data.frame(metadata, stringsAsFactors = FALSE)
    }

    if (!is.data.frame(metadata)) {
      stop("metadata must be NULL, a named list, a named vector, or a one-row data.frame")
    }

    if (nrow(metadata) != 1) {
      stop("metadata must contain exactly one row/value per metadata field")
    }

    if (any(names(metadata) %in% names(df))) {
      stop(
        paste0(
          "metadata contains column names already used by tile_aa(): ",
          paste(intersect(names(metadata), names(df)), collapse = ", ")
        )
      )
    }

    cbind(
      df,
      metadata[rep(1, nrow(df)), , drop = FALSE]
    )
  }

  if(nchar(aa_seq) < tile_size){
    full_tiling <- data.frame("cluster_id" = identifier,
                              "start" = 1,
                              "end" = nchar(aa_seq),
                              "tile_seq" = aa_seq,
                              "tile_length" = nchar(aa_seq),
                              "padded" = FALSE,
                              "n_pad_aa" = 0,
                              "gets_replaced" = "yes")
    new_aa_seq <- aa_seq
    while(nchar(new_aa_seq) != tile_size){
      remaining_aa_to_add <- tile_size - nchar(new_aa_seq)
      if(nchar(aa_seq) <= remaining_aa_to_add){
        new_aa_seq <- paste0(new_aa_seq, aa_seq)
      }else{
        new_aa_seq <- paste0(new_aa_seq, substr(new_aa_seq, 1, remaining_aa_to_add))
      }
    }
    add_pad <- data.frame("cluster_id" = identifier,
                          "start" = 1,
                          "end" = nchar(aa_seq),
                          "tile_seq" = new_aa_seq,
                          "tile_length" = nchar(new_aa_seq),
                          "padded" = TRUE,
                          "n_pad_aa" = tile_size - nchar(aa_seq),
                          "gets_replaced" = "no")

    full_tiling <- rbind(full_tiling, add_pad)
    rm(add_pad, new_aa_seq)
  }

  if(nchar(aa_seq) >= tile_size){
    track_one_starts <- seq.int(from = 1, to = nchar(aa_seq), by = tile_size)
    track_one_ends <- pmin(track_one_starts + tile_size -1 , nchar(aa_seq))
    track_one_coords <- data.frame(cluster_id = identifier, start = track_one_starts, end = track_one_ends)
    track_one_coords$tile_seq <- substr(rep(aa_seq, times = nrow(track_one_coords)), track_one_coords[["start"]], track_one_coords[["end"]])
    if(!paste(track_one_coords$tile_seq, collapse = "") == aa_seq){
      stop("track one coordinates and substringing don't recapitulate original input sequence")
    }
    track_one_coords$tile_length <- nchar(track_one_coords$tile_seq)
    track_one_coords$padded <- FALSE
    track_one_coords$n_pad_aa <- 0
    track_one_coords$gets_replaced <- "no"

    track_two_starts <- seq.int(from = (offset+1), to = nchar(aa_seq), by = tile_size)
    track_two_ends <- pmin(track_two_starts + tile_size -1 , nchar(aa_seq))
    track_two_coords <- data.frame(cluster_id = identifier, start = track_two_starts, end = track_two_ends)
    track_two_coords$tile_seq <- substr(rep(aa_seq, times = nrow(track_two_coords)), track_two_coords[["start"]], track_two_coords[["end"]])
    # if(!paste(track_two_coords$tile_seq, collapse = "") == aa_seq){
    #   stop("track two coordinates and substringing don't recapitulate original input sequence")
    # }
    track_two_coords$tile_length <- nchar(track_two_coords$tile_seq)
    track_two_coords$padded <- FALSE
    track_two_coords$n_pad_aa <- 0
    track_two_coords$gets_replaced <- "no"

    full_tiling <- rbind(track_one_coords, track_two_coords)
    full_tiling <- full_tiling[order(full_tiling$start),]

    full_tiling[which(full_tiling$tile_length != tile_size), "gets_replaced"] <- "yes"

    if(any(unique(full_tiling$tile_length) != tile_size)){
      if(length(unique(full_tiling$tile_length)) != 3){
        warning(
          paste0(
            "Found ",
            length(unique(full_tiling$tile_length)),
            " unique tile lengths: ",
            paste(sort(unique(full_tiling$tile_length)), collapse = ", "),
            ". Current code is designed to handle 2 or 3 unique tile lengths."
          )
        )
      }
      if(length(unique(full_tiling$tile_length)) == 2){
        tile_to_replace <- full_tiling[which(full_tiling$tile_length != tile_size & full_tiling$tile_length == min(unique(full_tiling$tile_length))),]
      }
      if(length(unique(full_tiling$tile_length)) == 3){
        tile_to_replace <- full_tiling[which(full_tiling$tile_length != tile_size & full_tiling$tile_length != min(unique(full_tiling$tile_length))),]
      }
      add_pad <- data.frame("cluster_id" = identifier,
                            "start" = (nchar(aa_seq) - tile_size + 1),
                            "end" = nchar(aa_seq),
                            "tile_seq" = substr(aa_seq, (nchar(aa_seq) - tile_size + 1), nchar(aa_seq)),
                            "tile_length" = nchar(substr(aa_seq, (nchar(aa_seq) - tile_size + 1), nchar(aa_seq))),
                            "padded" = TRUE,
                            "n_pad_aa" = tile_size - tile_to_replace$tile_length,
                            "gets_replaced" = "no")
      full_tiling <- rbind(full_tiling, add_pad)
      rm(tile_to_replace, add_pad)
    }

  }

  if (is.null(full_tiling)) {
    stop(
      paste0(
        "tile_aa() failed to create full_tiling for cluster_id = ",
        identifier,
        ". Sequence length = ",
        nchar(aa_seq),
        ", tile_size = ",
        tile_size,
        ", offset = ",
        offset,
        "."
      )
    )
  }

  full_tiling <- add_metadata(full_tiling, metadata)
  return(full_tiling)
}
