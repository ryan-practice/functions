tile_aa <- function(
    aa_seq,
    identifier,
    tile_size,
    offset,
    keep_partial = TRUE,
    return_df = TRUE
){
  full_tiling <- NULL
  aa_seq <- toupper(aa_seq)

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
        warning(paste0(length(unique(full_tiling$tile_length)), " is equal to ", length(unique(full_tiling$tile_length)), ". Currently code should handle length 2 or 3 properly"))
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

  return(full_tiling)
}
