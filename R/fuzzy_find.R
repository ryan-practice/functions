fuzzy_find <- function(pattern, string, len_len_moving_window, max_string_distance){
  start_numbers <- c(1:(nchar(string)-(len_moving_window-1)))
  stop_numbers <- start_numbers + (len_moving_window-1)
  v_substr <- Vectorize(substr, vectorize.args = c("start", "stop"))
  
  coordinate_substr <- data.frame("start_coord" = start_numbers,
                                   "end_coord" = stop_numbers,
                                   "substrings" = v_substr(string, start = start_numbers, stop = stop_numbers)
  )
  coordinate_substr$stringdist <- stringdist(pattern, coordinate_substr$substrings)
  coordinate_substr <- coordinate_substr[which(coordinate_substr$stringdist <= max_string_distance),]
}