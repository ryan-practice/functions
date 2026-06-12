reduce_similar_tiles <- function(
    tiles,
    tile_col = "tile_seq",
    protected_col = "protected",
    k = 6,
    sketch_size = 50,
    candidate_j_est_threshold = 0.55,
    final_jaccard_threshold = 0.70,
    max_bucket = 200L,
    reduction_method = c("single_linkage", "greedy"),
    na_protected = TRUE,
    return_details = FALSE,
    verbose = TRUE
) {

  # ============================================================
  # Package checks
  # ============================================================

  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package 'data.table' is required.")
  }

  if (!requireNamespace("digest", quietly = TRUE)) {
    stop("Package 'digest' is required.")
  }

  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("Package 'igraph' is required.")
  }

  reduction_method <- match.arg(reduction_method)

  # ============================================================
  # Sanity checks
  # ============================================================

  if (!inherits(tiles, "data.frame")) {
    stop("'tiles' must be a data.frame, data.table, or tibble-like object.")
  }

  if (!tile_col %in% colnames(tiles)) {
    stop("Could not find tile sequence column: ", tile_col)
  }

  if (!protected_col %in% colnames(tiles)) {
    stop("Could not find protected column: ", protected_col)
  }

  n_tiles <- nrow(tiles)

  if (n_tiles == 0L) {
    if (return_details) {
      return(list(
        reduced_tiles = tiles,
        removed_tiles = tiles,
        retained_indices = integer(0),
        removed_indices = integer(0),
        similar_pairs = data.frame(),
        candidates = data.frame(),
        settings = list(
          k = k,
          sketch_size = sketch_size,
          candidate_j_est_threshold = candidate_j_est_threshold,
          final_jaccard_threshold = final_jaccard_threshold,
          max_bucket = max_bucket,
          reduction_method = reduction_method,
          na_protected = na_protected
        )
      ))
    } else {
      return(tiles)
    }
  }

  # ============================================================
  # Identify protected tiles
  # ============================================================
  # Any value other than "no" is considered protected.
  # By default, NA is treated as protected to avoid accidentally
  # discarding rows with missing protection labels.

  protected_clean <- trimws(tolower(as.character(tiles[[protected_col]])))

  if (na_protected) {
    is_protected <- is.na(protected_clean) | protected_clean != "no"
  } else {
    is_protected <- !is.na(protected_clean) & protected_clean != "no"
  }

  # ============================================================
  # Amino-acid k-mer generation
  # ============================================================

  get_aa_kmers <- function(seq, k) {
    seq <- toupper(as.character(seq))
    n <- nchar(seq)

    if (is.na(n) || n < k) {
      return(character(0))
    }

    starts <- seq_len(n - k + 1)
    substring(seq, first = starts, last = starts + k - 1)
  }

  kmers_binary <- lapply(tiles[[tile_col]], get_aa_kmers, k = k)
  kmers_binary <- lapply(kmers_binary, unique)

  # ============================================================
  # Hashing and sketch generation
  # ============================================================

  hash32 <- function(x) {
    hx <- digest::digest(x, algo = "xxhash32", serialize = FALSE)

    # Use first 7 hex characters, 28 bits, to avoid signed integer overflow.
    as.integer(strtoi(substr(hx, 1, 7), base = 16L))
  }

  make_sketch <- function(kmers, sketch_size) {
    if (length(kmers) == 0L) {
      return(integer(0))
    }

    h <- vapply(kmers, hash32, integer(1))
    h <- unique(h)
    h <- h[!is.na(h)]

    if (length(h) == 0L) {
      return(integer(0))
    }

    h <- sort(h)
    h[seq_len(min(sketch_size, length(h)))]
  }

  sketches <- lapply(kmers_binary, make_sketch, sketch_size = sketch_size)
  sketch_lengths <- lengths(sketches)

  # ============================================================
  # Invert sketches: hash -> sequence IDs
  # ============================================================

  sk_dt <- data.table::data.table(
    seq_id = rep.int(seq_along(sketches), sketch_lengths),
    h = unlist(sketches, use.names = FALSE)
  )

  if (nrow(sk_dt) == 0L) {

    pairs <- data.table::data.table(
      i = integer(),
      j = integer(),
      n_shared_sketch = integer(),
      j_est = numeric()
    )

  } else {

    bucket_sizes <- sk_dt[, .N, by = h]

    keep_h <- bucket_sizes[N >= 2L & N <= max_bucket, h]
    sk_dt2 <- sk_dt[h %in% keep_h]

    if (nrow(sk_dt2) == 0L) {

      pairs <- data.table::data.table(
        i = integer(),
        j = integer(),
        n_shared_sketch = integer(),
        j_est = numeric()
      )

    } else {

      pair_hits <- sk_dt2[, {
        v <- sort(unique(seq_id))

        if (length(v) < 2L) {
          data.table::data.table(i = integer(), j = integer())
        } else {
          cmb <- utils::combn(v, 2L)
          data.table::data.table(i = cmb[1L, ], j = cmb[2L, ])
        }
      }, by = h]

      if (nrow(pair_hits) == 0L) {

        pairs <- data.table::data.table(
          i = integer(),
          j = integer(),
          n_shared_sketch = integer(),
          j_est = numeric()
        )

      } else {

        pairs <- pair_hits[, .(n_shared_sketch = .N), by = .(i, j)]

        pairs[, j_est := n_shared_sketch /
                (sketch_lengths[i] + sketch_lengths[j] - n_shared_sketch)]
      }
    }
  }

  # ============================================================
  # Keep likely-similar candidates
  # ============================================================

  candidates <- pairs[j_est >= candidate_j_est_threshold]

  if (nrow(candidates) > 0L) {

    jaccard_scores <- vapply(seq_len(nrow(candidates)), function(idx) {
      a <- kmers_binary[[candidates$i[idx]]]
      b <- kmers_binary[[candidates$j[idx]]]

      union_length <- length(union(a, b))

      if (union_length == 0L) {
        return(0)
      }

      length(intersect(a, b)) / union_length
    }, numeric(1))

    candidates[, jac_kmer := jaccard_scores]

  } else {

    candidates[, jac_kmer := numeric(0)]
  }

  # ============================================================
  # Final similar pairs
  # ============================================================

  keep_pairs <- candidates[
    jac_kmer >= final_jaccard_threshold,
    .(i, j, jac_kmer)
  ]

  # ============================================================
  # Reduction method 1:
  # Single-linkage connected components
  # ============================================================

  if (reduction_method == "single_linkage") {

    if (nrow(keep_pairs) == 0L) {

      rep_indices <- seq_len(n_tiles)

    } else {

      edge_df <- data.frame(
        from = as.character(keep_pairs$i),
        to = as.character(keep_pairs$j)
      )

      vertex_df <- data.frame(
        name = as.character(seq_len(n_tiles))
      )

      g <- igraph::graph_from_data_frame(
        d = edge_df,
        directed = FALSE,
        vertices = vertex_df
      )

      comps <- igraph::components(g)

      membership <- comps$membership[
        as.character(seq_len(n_tiles))
      ]

      component_list <- split(seq_len(n_tiles), membership)

      rep_indices <- unlist(lapply(component_list, function(x) {

        protected_in_component <- x[is_protected[x]]

        if (length(protected_in_component) > 0L) {
          # Keep all protected tiles in this similarity component.
          protected_in_component
        } else {
          # If no protected tiles exist in the component, keep the first tile.
          x[1L]
        }

      }), use.names = FALSE)

      rep_indices <- sort(unique(rep_indices))
    }
  }

  # ============================================================
  # Reduction method 2:
  # Greedy direct-neighbor filtering
  # ============================================================

  if (reduction_method == "greedy") {

    adjacency <- vector("list", n_tiles)

    if (nrow(keep_pairs) > 0L) {

      for (edge_idx in seq_len(nrow(keep_pairs))) {
        a <- keep_pairs$i[edge_idx]
        b <- keep_pairs$j[edge_idx]

        adjacency[[a]] <- c(adjacency[[a]], b)
        adjacency[[b]] <- c(adjacency[[b]], a)
      }

      adjacency <- lapply(adjacency, unique)
    }

    keep_tile <- rep(FALSE, n_tiles)

    # Always retain all protected tiles.
    keep_tile[is_protected] <- TRUE

    # Then process unprotected tiles in existing row order.
    # An unprotected tile is retained only if it is not directly similar
    # to anything already retained.
    for (idx in seq_len(n_tiles)) {

      if (is_protected[idx]) {
        next
      }

      neighbors <- adjacency[[idx]]

      if (length(neighbors) == 0L) {
        keep_tile[idx] <- TRUE
      } else if (!any(keep_tile[neighbors])) {
        keep_tile[idx] <- TRUE
      }
    }

    rep_indices <- which(keep_tile)
  }

  # ============================================================
  # Final output
  # ============================================================

  removed_indices <- setdiff(seq_len(n_tiles), rep_indices)

  reduced_tiles <- tiles[rep_indices, , drop = FALSE]
  removed_tiles <- tiles[removed_indices, , drop = FALSE]

  if (verbose) {
    cat("\nSimilarity filtering complete.\n")
    cat("Reduction method:", reduction_method, "\n")
    cat("Input tiles:", n_tiles, "\n")
    cat("Similar pairs found:", nrow(keep_pairs), "\n")
    cat("Protected tiles in input:", sum(is_protected), "\n")
    cat("Tiles retained:", nrow(reduced_tiles), "\n")
    cat("Tiles removed:", nrow(removed_tiles), "\n")
    cat("Protected tiles retained:", sum(is_protected[rep_indices]), "\n")
    cat("Protected tiles removed:", sum(is_protected[removed_indices]), "\n")
  }

  if (return_details) {

    return(list(
      reduced_tiles = reduced_tiles,
      removed_tiles = removed_tiles,
      retained_indices = rep_indices,
      removed_indices = removed_indices,
      similar_pairs = as.data.frame(keep_pairs),
      candidates = as.data.frame(candidates),
      is_protected = is_protected,
      settings = list(
        tile_col = tile_col,
        protected_col = protected_col,
        k = k,
        sketch_size = sketch_size,
        candidate_j_est_threshold = candidate_j_est_threshold,
        final_jaccard_threshold = final_jaccard_threshold,
        max_bucket = max_bucket,
        reduction_method = reduction_method,
        na_protected = na_protected
      )
    ))

  } else {

    return(reduced_tiles)
  }
}
