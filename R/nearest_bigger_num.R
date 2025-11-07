nearest_bigger_num <- function(num, vec) {
  which(min(vec[num < vec]) == vec)
}
