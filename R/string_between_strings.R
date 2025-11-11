string_between_strings <- function(string, STR1, STR2){
  return(stringr::str_match(string, paste0(STR1, "\\s*(.*?)\\s*", STR2))[,2])
}
