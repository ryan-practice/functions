extract_max_per_group <- function(df, group_col, value_col) {
  df %>%
    group_by({{ group_col }}) %>%
    slice_max(order_by = {{ value_col }}, n = 1, with_ties = FALSE) %>%
    ungroup()
}
