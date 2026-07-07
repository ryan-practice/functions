#' Plot a CDR3/J-gene reconstruction safety check
#'
#' Creates a visual safety-check plot for a single reconstructed TCR chain. The
#' plot aligns the sequencing-derived CDR3 amino acid sequence, the final
#' reconstructed chain, and the corresponding J-gene amino acid sequence. The
#' CDR3/J-gene overlap is highlighted in green, and the appended J-derived GxG
#' motif is highlighted in yellow.
#'
#' This function is intended to help confirm that amino acids already present in
#' the sequencing-derived CDR3 are not duplicated when the J-gene suffix is
#' appended. It supports both alpha-chain TRAJ inputs and beta-chain TRBJ inputs;
#' plot labels are inferred automatically from `j_gene_name`.
#'
#' @param cdr3 Character string containing one sequencing-derived CDR3 amino acid
#'   sequence.
#' @param j_gene_name Character string containing the J-gene name, such as
#'   `"TRAJ42"`, `"TRAJ42*01"`, `"TRBJ2-5"`, or `"TRBJ2-5*01"`.
#' @param j_gene Character string containing the amino acid sequence of the
#'   corresponding J gene.
#' @param final_chain Character string containing the final reconstructed alpha
#'   or beta chain amino acid sequence.
#' @param plot_title Optional character string to use as a custom plot title. In
#'   the current plotting code, the `title = plot_title` line in `labs()` is
#'   commented out; uncomment that line to display the title.
#'
#' @return A `ggplot` object showing the aligned CDR3, final reconstructed chain,
#'   and J-gene sequence. Green tiles indicate the CDR3/J overlap. Yellow tiles
#'   indicate the selected appended J-derived GxG motif.
#'
#' @details
#' The function first identifies the longest suffix of `cdr3` that appears in
#' `j_gene`. This overlap is used to align the J-gene row beneath the matching
#' region of the CDR3 row. The portion of the J gene after this overlap is
#' treated as the appended J-gene suffix.
#'
#' The GxG motif highlight is chosen from motifs that occur in the appended
#' J-gene suffix, rather than from motifs that are already present in the CDR3
#' overlap. For J genes with multiple GxG motifs, the function relies on
#' `select_appended_gxg()` to choose the motif most consistent with the appended
#' J-derived framework region.
#'
#' Special TRAJ cases with known non-standard GxG counts are labeled in the plot:
#' TRAJ16 and TRAJ61 are expected to have zero GxG motifs, while TRAJ10, TRAJ25,
#' TRAJ28, TRAJ35, and TRAJ45 are expected to have two GxG motifs.
#'
#' @examples
#' plot_cdr3_j_safety_check(
#'   cdr3 = "CAMREGTGANSKLT",
#'   j_gene_name = "TRAJ42",
#'   j_gene = "XYTGANSKLTFGKGITLSVRP",
#'   final_chain = "CAMREGTGANSKLTFGKGITLSVRP"
#' )
#'
#' plot_cdr3_j_safety_check(
#'   cdr3 = "ASSQSPGGTQY",
#'   j_gene_name = "TRBJ2-5",
#'   j_gene = "QETQYFGPGTRLLVL",
#'   final_chain = "ASSQSPGGTQYFGPGTRLLVL"
#' )
#'
#' @importFrom dplyr bind_rows mutate case_when
#' @importFrom tibble tibble
#' @importFrom stringr str_locate fixed
#' @importFrom ggplot2 ggplot aes geom_tile geom_text geom_label
#'   scale_fill_manual labs theme_minimal theme element_blank element_text
#'
#' @export
plot_cdr3_j_safety_check <- function(
    cdr3,
    j_gene_name,
    j_gene,
    final_chain,
    plot_title = NULL
) {

  zero_gxg_special_cases <- c("TRAJ16", "TRAJ61")
  two_gxg_special_cases <- c("TRAJ10", "TRAJ25", "TRAJ28", "TRAJ35", "TRAJ45")

  j_gene_name_clean <- clean_j_gene_name(j_gene_name)
  j_type <- infer_j_type(j_gene_name)
  chain_label <- infer_chain_label(j_type)

  overlap <- find_longest_cdr3_j_overlap(cdr3, j_gene)

  if (is.na(overlap)) {
    stop("No CDR3 suffix overlap found in J gene.")
  }

  cdr3_overlap_start <- str_locate(cdr3, fixed(overlap))[1, "start"]
  j_overlap_start <- str_locate(j_gene, fixed(overlap))[1, "start"]

  # Align the J gene so its overlap sits under the matching CDR3 overlap.
  j_plot_start <- cdr3_overlap_start - j_overlap_start + 1

  # First J-gene amino acid that gets appended after the overlap.
  append_start_j <- j_overlap_start + nchar(overlap)

  gxg_info <- select_appended_gxg(
    j_gene = j_gene,
    append_start_j = append_start_j
  )

  cdr3_label <- "CDR3 from sequencing"
  final_label <- paste("Final", chain_label, "chain")
  j_label <- paste(j_type, "gene")

  row_order <- c(
    cdr3_label,
    final_label,
    j_label
  )

  seq_to_df <- function(seq_key, seq_name, aa_seq, start_pos = 1) {
    tibble(
      seq_key = seq_key,
      seq_name = seq_name,
      aa = strsplit(aa_seq, "")[[1]],
      seq_index = seq_along(aa),
      plot_pos = start_pos + seq_index - 1
    )
  }

  plot_df <- bind_rows(
    seq_to_df("cdr3", cdr3_label, cdr3, start_pos = 1),
    seq_to_df("final_chain", final_label, final_chain, start_pos = 1),
    seq_to_df("j_gene", j_label, j_gene, start_pos = j_plot_start)
  )

  overlap_positions <- cdr3_overlap_start:(cdr3_overlap_start + nchar(overlap) - 1)

  if (!is.na(gxg_info$motif_start)) {
    gxg_plot_start <- j_plot_start + gxg_info$motif_start - 1
    gxg_plot_end <- j_plot_start + gxg_info$motif_end - 1
    gxg_positions <- gxg_plot_start:gxg_plot_end
  } else {
    gxg_positions <- integer()
  }

  gxg_note <- gxg_info$gxg_status

  if (j_gene_name_clean %in% zero_gxg_special_cases && gxg_info$n_gxg_total == 0) {
    gxg_note <- paste0(j_gene_name_clean, ": 0 GxG motifs")
  }

  if (j_gene_name_clean %in% two_gxg_special_cases && gxg_info$n_gxg_total == 2) {
    gxg_note <- paste0(
      j_gene_name_clean,
      ": 2 GxG motifs; highlighted appended motif ",
      gxg_info$motif
    )
  }

  plot_df <- plot_df %>%
    mutate(
      seq_name = factor(seq_name, levels = rev(row_order)),

      is_overlap = plot_pos %in% overlap_positions,

      # Yellow marks the appended J-derived GxG in the final chain and J-gene rows.
      # It is intentionally not marked in the CDR3 row.
      is_appended_gxg =
        plot_pos %in% gxg_positions &
        seq_key %in% c("final_chain", "j_gene"),

      highlight = case_when(
        is_appended_gxg ~ "Appended GxG motif",
        is_overlap ~ "CDR3/J overlap",
        TRUE ~ "Other"
      )
    )

  label_df <- tibble(
    x = max(plot_df$plot_pos, na.rm = TRUE),
    y = factor(cdr3_label, levels = rev(row_order)),
    label = gxg_note
  )

  if (is.null(plot_title)) {
    plot_title <- paste0(
      "CDR3 / ",
      j_type,
      " overlap and appended GxG safety check"
    )
  }

  appended_j_suffix <- ifelse(
    append_start_j <= nchar(j_gene),
    substr(j_gene, append_start_j, nchar(j_gene)),
    ""
  )

  ggplot(plot_df, aes(x = plot_pos, y = seq_name)) +
    geom_tile(
      aes(fill = highlight),
      color = "grey70",
      width = 0.95,
      height = 0.85
    ) +
    geom_text(aes(label = aa), size = 5) +
    geom_label(
      data = label_df,
      aes(x = x, y = y, label = label),
      inherit.aes = FALSE,
      hjust = 1,
      vjust = -0.8,
      size = 3.8,
      linewidth = 0.2
    ) +
    scale_fill_manual(
      values = c(
        "CDR3/J overlap" = "palegreen3",
        "Appended GxG motif" = "gold",
        "Other" = "white"
      ),
      breaks = c("CDR3/J overlap", "Appended GxG motif", "Other"),
      name = NULL
    ) +
    labs(
      x = "Aligned amino acid position",
      y = NULL,
      #title = plot_title,
      subtitle = paste0(
        "J gene: ", j_gene_name, "\n",
        "overlap: ", overlap, "\n",
        "appended J suffix: ", appended_j_suffix
      )
    ) +
    theme_minimal(base_size = 14) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(size = 10),
      axis.text.y = element_text(size = 12),
      legend.position = "bottom"
    )
}
