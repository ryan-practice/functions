reverse_complement <- function(nt_sequence){
  reverse_sequence <- unlist(lapply(strsplit(nt_sequence, ""), rev))
  a_indices <- which(reverse_sequence == "a")
  A_indices <- which(reverse_sequence == "A")

  t_indices <- which(reverse_sequence == "t")
  T_indices <- which(reverse_sequence == "T")

  g_indices <- which(reverse_sequence == "g")
  G_indices <- which(reverse_sequence == "G")

  c_indices <- which(reverse_sequence == "c")
  C_indices <- which(reverse_sequence == "C")

  reverse_sequence[a_indices] <- "t"
  reverse_sequence[A_indices] <- "T"

  reverse_sequence[t_indices] <- "a"
  reverse_sequence[T_indices] <- "A"

  reverse_sequence[g_indices] <- "c"
  reverse_sequence[G_indices] <- "C"

  reverse_sequence[c_indices] <- "g"
  reverse_sequence[C_indices] <- "G"

  reverse_comp_seq <- paste(reverse_sequence, collapse = "")

}
