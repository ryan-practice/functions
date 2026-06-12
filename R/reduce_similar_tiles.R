#' Reduce peptide tiles by k-mer similarity
#'
#' Reduces a peptide tile table by identifying highly similar peptide sequences
#' using amino-acid k-mer Jaccard similarity. The function first uses a compact
#' hash-sketching step to identify likely similar candidate pairs, then computes
#' exact k-mer Jaccard similarity only for those candidates.
#'
#' Protected tiles are always retained. Protection status is determined from the
#' column specified by `protected_col`. Any value other than `"no"` is treated as
#' protected. By default, missing protection values are also treated as protected
#' to avoid accidentally discarding tiles with missing labels.
#'
#' @param tiles A data frame, data table, or tibble containing peptide tile
#'   sequences and a protection-status column.
#' @param tile_col Character string. Name of the column containing peptide
#'   sequences. Default is `"tile_seq"`.
#' @param protected_col Character string. Name of the column indicating whether
#'   a tile is protected. Default is `"protected"`.
#' @param k Integer. Amino-acid k-mer size used for similarity comparison.
#'   Default is `6`.
#' @param sketch_size Integer. Number of smallest k-mer hashes retained per tile
#'   during sketching. Larger values may improve recall but can increase runtime.
#'   Default is `50`.
#' @param candidate_j_est_threshold Numeric. Minimum approximate sketch-based
#'   Jaccard similarity required for a pair to be evaluated using exact k-mer
#'   Jaccard similarity. Default is `0.55`.
#' @param final_jaccard_threshold Numeric. Minimum exact k-mer Jaccard similarity
#'   required for two tiles to be considered similar or redundant. Default is
#'   `0.70`.
#' @param max_bucket Integer. Maximum number of tiles allowed to share a sketch
#'   hash before that hash is ignored during candidate-pair generation. This
#'   prevents very common hashes from generating excessive candidate pairs.
#'   Default is `200L`.
#' @param reduction_method Character string. Similarity reduction strategy.
#'   Options are `"single_linkage"` and `"greedy"`. Default is
#'   `"single_linkage"`.
#' @param na_protected Logical. If `TRUE`, missing values in `protected_col` are
#'   treated as protected. If `FALSE`, missing values are treated as unprotected.
#'   Default is `TRUE`.
#' @param return_details Logical. If `FALSE`, return only the reduced tile table.
#'   If `TRUE`, return a list with the reduced table, removed table, retained
#'   indices, removed indices, similar pairs, candidate pairs, protection flags,
#'   and settings. Default is `FALSE`.
#' @param verbose Logical. If `TRUE`, print a short filtering summary. Default is
#'   `TRUE`.
#'
#' @details
#' Each peptide sequence is converted into a set of unique overlapping amino-acid
#' k-mers. These k-mers are hashed, and the smallest `sketch_size` hash values
#' are retained as a compact sketch for each tile. Pairs of tiles that share
#' sketch hashes are used as candidate pairs. Exact k-mer Jaccard similarity is
#' then calculated only for those candidate pairs.
#'
#' In `"single_linkage"` mode, similar tiles are represented as an undirected
#' graph. Connected components of similar tiles are reduced together. All
#' protected tiles in each component are retained. If a component contains no
#' protected tiles, the first tile in that component is retained.
#'
#' In `"greedy"` mode, all protected tiles are retained first. Unprotected tiles
#' are then processed in their existing row order. An unprotected tile is retained
#' only if it is not directly similar to any tile already retained. This avoids
#' the transitive chaining behavior of single-linkage clustering.
#'
#' @return
#' If `return_details = FALSE`, returns the retained tile table.
#'
#' If `return_details = TRUE`, returns a list containing the retained tile table,
#' removed tile table, retained row indices, removed row indices, final similar
#' tile pairs, candidate tile pairs, logical protection flags, and parameter
#' settings used for the run.
#'
#' @examples
#' \dontrun{
#' complete_tiles_sim_reduced <- reduce_similar_tiles(complete_tiles)
#'
#' complete_tiles_sim_reduced_greedy <- reduce_similar_tiles(
#'   complete_tiles,
#'   reduction_method = "greedy"
#' )
#'
#' sim_result <- reduce_similar_tiles(
#'   complete_tiles,
#'   reduction_method = "single_linkage",
#'   return_details = TRUE
#' )
#'
#' complete_tiles_sim_reduced <- sim_result$reduced_tiles
#' complete_tiles_sim_removed <- sim_result$removed_tiles
#' keep_pairs <- sim_result$similar_pairs
#' }
#'
#' @export
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

  if (!requireNamespace("digest", quietly = TRUE)) {
    stop("Package 'digest' is required.")
  }

  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("Package 'igraph' is required.")
  }

  reduction_method <- match.arg(reduction_method)

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
      return(tiles)
    }
  }

  protected_clean <- trimws(tolower(as.character(tiles[[protected_col]])))

  if (na_protected) {
    is_protected <- is.na(protected_clean) | protected_clean != "no"
  } else {
    is_protected <- !is.na(protected_clean) & protected_clean != "no"
  }

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

  hash32 <- function(x) {
    hx <- digest::digest(x, algo = "xxhash32", serialize = FALSE)
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

  seq_id_vec <- rep.int(seq_along(sketches), sketch_lengths)
  h_vec <- unlist(sketches, use.names = FALSE)

  if (length(h_vec) == 0L) {

    pairs <- data.frame(
      i = integer(),
      j = integer(),
      n_shared_sketch = integer(),
      j_est = numeric()
    )

  } else {

    unique_h <- unique(h_vec)
    h_group <- match(h_vec, unique_h)
    bucket_n <- tabulate(h_group, nbins = length(unique_h))

    keep_groups <- which(bucket_n >= 2L & bucket_n <= max_bucket)

    if (length(keep_groups) == 0L) {

      pairs <- data.frame(
        i = integer(),
        j = integer(),
        n_shared_sketch = integer(),
        j_est = numeric()
      )

    } else {

      keep_rows <- h_group %in% keep_groups

      bucket_list <- split(
        seq_id_vec[keep_rows],
        h_group[keep_rows],
        drop = TRUE
      )

      pair_keys_list <- lapply(bucket_list, function(v) {
        v <- sort(unique(v))

        if (length(v) < 2L) {
          return(character(0))
        }

        cmb <- utils::combn(v, 2L)
        paste(cmb[1L, ], cmb[2L, ], sep = "\t")
      })

      pair_keys <- unlist(pair_keys_list, use.names = FALSE)

      if (length(pair_keys) == 0L) {

        pairs <- data.frame(
          i = integer(),
          j = integer(),
          n_shared_sketch = integer(),
          j_est = numeric()
        )

      } else {

        pair_tab <- table(pair_keys)
        pair_split <- strsplit(names(pair_tab), "\t", fixed = TRUE)
        pair_mat <- matrix(
          as.integer(unlist(pair_split, use.names = FALSE)),
          ncol = 2,
          byrow = TRUE
        )

        pairs <- data.frame(
          i = pair_mat[, 1],
          j = pair_mat[, 2],
          n_shared_sketch = as.integer(pair_tab),
          stringsAsFactors = FALSE
        )

        pairs$j_est <- pairs$n_shared_sketch /
          (
            sketch_lengths[pairs$i] +
              sketch_lengths[pairs$j] -
              pairs$n_shared_sketch
          )
      }
    }
  }

  candidates <- pairs[
    pairs$j_est >= candidate_j_est_threshold,
    ,
    drop = FALSE
  ]

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

    candidates$jac_kmer <- jaccard_scores

  } else {

    candidates$jac_kmer <- numeric(0)
  }

  keep_pairs <- candidates[
    candidates$jac_kmer >= final_jaccard_threshold,
    c("i", "j", "jac_kmer"),
    drop = FALSE
  ]

  if (reduction_method == "single_linkage") {

    if (nrow(keep_pairs) == 0L) {

      rep_indices <- seq_len(n_tiles)

    } else {

      edge_df <- data.frame(
        from = as.character(keep_pairs$i),
        to = as.character(keep_pairs$j),
        stringsAsFactors = FALSE
      )

      vertex_df <- data.frame(
        name = as.character(seq_len(n_tiles)),
        stringsAsFactors = FALSE
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
          protected_in_component
        } else {
          x[1L]
        }

      }), use.names = FALSE)

      rep_indices <- sort(unique(rep_indices))
    }
  }

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

    keep_tile[is_protected] <- TRUE

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
      similar_pairs = keep_pairs,
      candidates = candidates,
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
