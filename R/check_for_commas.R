check_for_commas <- function(df) {
  cat(sapply(df, function(x) any(grepl(",", x))))
}
