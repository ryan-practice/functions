#' Convert AIRR Rearrangement Records to Paired TCR Records
#'
#' Converts an AIRR-style rearrangement table, where each row represents one
#' receptor chain, into a wide table containing one candidate paired T-cell
#' receptor per row.
#'
#' Within each identifier, every TRA row is paired with every TRB row using a
#' Cartesian product. For example, an identifier with two TRA rows and three
#' TRB rows produces six output rows.
#'
#' No duplicate chain records are removed. Identical input rows are retained
#' and participate independently in the Cartesian pairing. This allows
#' duplicate detection or filtering to be performed downstream.
#'
#' The output preserves the order in which identifiers first appear in the
#' input table. Within each identifier, TRA and TRB records retain their
#' original relative input order.
#'
#' Columns originating from TRA rows are prefixed with `"TRA_"`, while columns
#' originating from TRB rows are prefixed with `"TRB_"`. The identifier column
#' remains unprefixed.
#'
#' @param airr_data A data frame or tibble containing AIRR-style rearrangement
#'   annotations. Each row should describe one receptor chain.
#'
#' @param id_col A single character string naming the column that links chains
#'   belonging to the same cell, sample, or receptor grouping. The default is
#'   `"cell_id"`.
#'
#' @param locus_col A single character string naming the locus column. Values
#'   are expected to include `"TRA"` and `"TRB"`. Matching is
#'   case-insensitive because values are converted to uppercase internally.
#'   The default is `"locus"`.
#'
#' @param include_unpaired Logical. When `FALSE`, identifiers lacking either a
#'   TRA or TRB row are omitted because no complete paired TCR can be created.
#'   When `TRUE`, identifiers containing only TRA or only TRB rows are retained,
#'   with the missing chain represented by `NA` values. The default is `FALSE`.
#'
#' @param metadata_col A single character string giving the name of the output
#'   metadata column containing the raw TRA and TRB row counts. Values use the
#'   machine-readable format `"TRA=<n>;TRB=<n>"`. The default is
#'   `"chain_counts"`.
#'
#' @param tcr_id_col A single character string giving the name of the generated
#'   paired-TCR identifier column. The default is `"tcr_id"`.
#'
#' @return A tibble containing one row per candidate TRA-TRB combination.
#'   The output contains:
#'
#'   \itemize{
#'     \item the original identifier column;
#'     \item a generated TCR identifier;
#'     \item a chain-count metadata column;
#'     \item all TRA-derived columns prefixed with `"TRA_"`;
#'     \item all TRB-derived columns prefixed with `"TRB_"`.
#'   }
#'
#' @details
#' Pairing is performed with a many-to-many join on the identifier column.
#' Therefore, if an identifier has \eqn{a} TRA rows and \eqn{b} TRB rows, the
#' output contains \eqn{a \times b} rows for that identifier.
#'
#' Identifier groups are returned according to the order in which each
#' identifier first appears in `airr_data`; identifiers are not sorted
#' alphabetically.
#'
#' The chain counts describe the number of raw TRA and TRB input rows, not the
#' number of unique biological sequences. For example, two identical TRA rows
#' and one TRB row produce:
#'
#' \preformatted{
#' chain_counts = "TRA=2;TRB=1"
#' }
#'
#' and two output TCR rows.
#'
#' Rows whose locus is neither TRA nor TRB are ignored with a warning.
#'
#' @examples
#' airr_example <- tibble::tribble(
#'   ~cell_id, ~locus, ~v_call,       ~j_call,      ~junction_aa,
#'   "cell_2", "TRA",  "TRAV1-2*01",  "TRAJ33*01",  "CAVRDSNYQLIW",
#'   "cell_2", "TRB",  "TRBV6-5*01",  "TRBJ1-2*01", "CASSYRGGNTIYF",
#'   "cell_1", "TRA",  "TRAV12-2*01", "TRAJ27*01",  "CAVRDTNTNAGKSTF",
#'   "cell_1", "TRA",  "TRAV8-4*01",  "TRAJ42*01",  "CAVSDGGSQGNLIF",
#'   "cell_1", "TRB",  "TRBV7-9*01",  "TRBJ2-5*01", "CASSLGQETQYF"
#' )
#'
#' paired_tcr <- airr_to_paired_tcr(airr_example)
#'
#' paired_tcr |>
#'   dplyr::select(
#'     cell_id,
#'     tcr_id,
#'     chain_counts,
#'     TRA_junction_aa,
#'     TRB_junction_aa
#'   )
#'
#' @export
airr_to_paired_tcr <- function(
    airr_data,
    id_col = "cell_id",
    locus_col = "locus",
    include_unpaired = FALSE,
    metadata_col = "chain_counts",
    tcr_id_col = "tcr_id"
) {

  #-----------------------------
  # 1. Validate arguments
  #-----------------------------

  if (!is.data.frame(airr_data)) {
    stop(
      "`airr_data` must be a data.frame or tibble.",
      call. = FALSE
    )
  }

  character_arguments <- list(
    id_col = id_col,
    locus_col = locus_col,
    metadata_col = metadata_col,
    tcr_id_col = tcr_id_col
  )

  invalid_character_arguments <- names(character_arguments)[
    !vapply(
      character_arguments,
      function(x) {
        is.character(x) &&
          length(x) == 1L &&
          !is.na(x) &&
          nzchar(x)
      },
      logical(1)
    )
  ]

  if (length(invalid_character_arguments) > 0L) {
    stop(
      "The following arguments must each be one non-empty character string: ",
      paste(invalid_character_arguments, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  if (
    !is.logical(include_unpaired) ||
    length(include_unpaired) != 1L ||
    is.na(include_unpaired)
  ) {
    stop(
      "`include_unpaired` must be either TRUE or FALSE.",
      call. = FALSE
    )
  }

  if (id_col == locus_col) {
    stop(
      "`id_col` and `locus_col` must refer to different columns.",
      call. = FALSE
    )
  }

  if (metadata_col == tcr_id_col) {
    stop(
      "`metadata_col` and `tcr_id_col` must have different names.",
      call. = FALSE
    )
  }

  required_cols <- c(id_col, locus_col)
  missing_cols <- setdiff(required_cols, names(airr_data))

  if (length(missing_cols) > 0L) {
    stop(
      "Missing required column(s): ",
      paste(missing_cols, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  generated_cols <- c(metadata_col, tcr_id_col)
  conflicting_cols <- intersect(generated_cols, names(airr_data))

  if (length(conflicting_cols) > 0L) {
    stop(
      "The input already contains the generated output column(s): ",
      paste(conflicting_cols, collapse = ", "),
      ". Choose different `metadata_col` or `tcr_id_col` names.",
      call. = FALSE
    )
  }

  airr_data <- tibble::as_tibble(airr_data)

  # Record the order in which identifiers first appear in the input.
  #
  # Character conversion is used only for matching positions. The identifier
  # column in the returned data retains its original type.
  original_id_order <- unique(
    as.character(airr_data[[id_col]])
  )

  #-----------------------------
  # 2. Validate identifiers
  #-----------------------------

  id_values <- as.character(airr_data[[id_col]])

  missing_id <- is.na(id_values) | trimws(id_values) == ""

  if (any(missing_id)) {
    stop(
      sum(missing_id),
      " row(s) have a missing or blank value in `",
      id_col,
      "`. Every chain must have an identifier for pairing.",
      call. = FALSE
    )
  }

  # Preserve the original identifier values but normalize locus values.
  airr_data[[locus_col]] <- toupper(
    trimws(as.character(airr_data[[locus_col]]))
  )

  #-----------------------------
  # 3. Retain TRA and TRB rows
  #-----------------------------

  observed_loci <- unique(airr_data[[locus_col]])

  unsupported_loci <- setdiff(
    observed_loci[!is.na(observed_loci)],
    c("TRA", "TRB")
  )

  if (length(unsupported_loci) > 0L) {
    warning(
      "Ignoring rows from unsupported loci: ",
      paste(unsupported_loci, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  chain_data <- airr_data |>
    dplyr::filter(
      .data[[locus_col]] %in% c("TRA", "TRB")
    )

  if (nrow(chain_data) == 0L) {
    stop(
      "No TRA or TRB rows were found in `",
      locus_col,
      "`.",
      call. = FALSE
    )
  }

  # No duplicated() or distinct() operation is performed here.
  # Every TRA and TRB input row is retained.

  #-----------------------------
  # 4. Count raw chain rows
  #-----------------------------

  chain_count_table <- chain_data |>
    dplyr::group_by(.data[[id_col]]) |>
    dplyr::summarise(
      n_TRA = sum(.data[[locus_col]] == "TRA"),
      n_TRB = sum(.data[[locus_col]] == "TRB"),
      .groups = "drop"
    )

  chain_count_table[[metadata_col]] <- sprintf(
    "TRA=%d;TRB=%d",
    chain_count_table$n_TRA,
    chain_count_table$n_TRB
  )

  incomplete_ids <- chain_count_table |>
    dplyr::filter(
      n_TRA == 0L | n_TRB == 0L
    )

  if (!include_unpaired && nrow(incomplete_ids) > 0L) {
    warning(
      nrow(incomplete_ids),
      " identifier(s) lacked either a TRA or TRB chain and were omitted. ",
      "Set `include_unpaired = TRUE` to retain them with missing values.",
      call. = FALSE
    )
  }

  output_metadata <- chain_count_table |>
    dplyr::select(
      dplyr::all_of(id_col),
      dplyr::all_of(metadata_col)
    )

  #-----------------------------
  # 5. Prepare TRA records
  #-----------------------------

  tra_data <- chain_data |>
    dplyr::filter(
      .data[[locus_col]] == "TRA"
    ) |>
    dplyr::rename_with(
      function(column_name) {
        paste0("TRA_", column_name)
      },
      .cols = -dplyr::all_of(id_col)
    ) |>
    dplyr::group_by(.data[[id_col]]) |>
    dplyr::mutate(
      .TRA_input_row = dplyr::row_number()
    ) |>
    dplyr::ungroup()

  #-----------------------------
  # 6. Prepare TRB records
  #-----------------------------

  trb_data <- chain_data |>
    dplyr::filter(
      .data[[locus_col]] == "TRB"
    ) |>
    dplyr::rename_with(
      function(column_name) {
        paste0("TRB_", column_name)
      },
      .cols = -dplyr::all_of(id_col)
    ) |>
    dplyr::group_by(.data[[id_col]]) |>
    dplyr::mutate(
      .TRB_input_row = dplyr::row_number()
    ) |>
    dplyr::ungroup()

  #-----------------------------
  # 7. Generate all TRA × TRB pairs
  #-----------------------------

  if (include_unpaired) {
    paired_data <- dplyr::full_join(
      tra_data,
      trb_data,
      by = id_col,
      relationship = "many-to-many"
    )
  } else {
    paired_data <- dplyr::inner_join(
      tra_data,
      trb_data,
      by = id_col,
      relationship = "many-to-many"
    )
  }

  if (nrow(paired_data) == 0L) {
    stop(
      "No paired TCR rows could be generated. No identifier contained both ",
      "a TRA and TRB row. Set `include_unpaired = TRUE` to retain unpaired ",
      "chains.",
      call. = FALSE
    )
  }

  #-----------------------------
  # 8. Restore original ID order
  #    and add TCR IDs
  #-----------------------------

  paired_data <- paired_data |>
    dplyr::left_join(
      output_metadata,
      by = id_col
    ) |>
    dplyr::mutate(
      .original_id_position = match(
        as.character(.data[[id_col]]),
        original_id_order
      )
    ) |>
    dplyr::arrange(
      .original_id_position,
      .TRA_input_row,
      .TRB_input_row
    ) |>
    dplyr::group_by(.data[[id_col]]) |>
    dplyr::mutate(
      .tcr_number = dplyr::row_number()
    ) |>
    dplyr::ungroup()

  maximum_tcr_number <- max(
    paired_data$.tcr_number,
    na.rm = TRUE
  )

  number_width <- max(
    2L,
    nchar(as.character(maximum_tcr_number))
  )

  paired_data[[tcr_id_col]] <- paste0(
    as.character(paired_data[[id_col]]),
    "__TCR",
    sprintf(
      paste0("%0", number_width, "d"),
      paired_data$.tcr_number
    )
  )

  #-----------------------------
  # 9. Arrange output columns
  #-----------------------------

  paired_data |>
    dplyr::select(
      dplyr::all_of(id_col),
      dplyr::all_of(tcr_id_col),
      dplyr::all_of(metadata_col),
      dplyr::starts_with("TRA_"),
      dplyr::starts_with("TRB_"),
      -dplyr::any_of(
        c(
          ".TRA_input_row",
          ".TRB_input_row",
          ".tcr_number",
          ".original_id_position"
        )
      )
    )
}
