split_column <- function(df, col, delimiter) {
  # Check inputs
  stopifnot(is.data.frame(df))
  if (!col %in% names(df)) stop(paste("Column", col, "not found in data frame"))

  # Split the target column
  split_list <- strsplit(as.character(df[[col]]), delimiter, fixed = TRUE)

  # Figure out the maximum number of splits
  max_len <- max(lengths(split_list))

  # Convert to data.frame
  split_df <- as.data.frame(do.call(rbind, lapply(split_list, function(x) {
    length(x) <- max_len
    x
  })), stringsAsFactors = FALSE)

  # Assign new column names
  new_colnames <- paste0(col, "_", seq_len(max_len))
  names(split_df) <- new_colnames

  # Combine with original data (optionally remove the old column)
  df <- cbind(df, split_df)

  return(df)
}
