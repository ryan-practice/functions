nearest_smaller_num <- function(num, vec) {
  which(max(vec[num > vec]) == vec)
}
