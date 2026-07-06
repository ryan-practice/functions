#' Infer chain label from J-gene type for CDR3 safety-check plotting
#'
#' Converts a J-gene type label into the corresponding TCR chain label. This
#' helper is used by the CDR3/J-gene safety-check plot to choose accurate row
#' labels, such as `"Final alpha chain"` for TRAJ inputs and `"Final beta chain"`
#' for TRBJ inputs.
#'
#' @param j_type Character vector of J-gene type labels, usually produced by
#'   `infer_j_type()`. Expected values include `"TRAJ"` and `"TRBJ"`.
#'
#' @return A character vector with one value per input J-gene type. Values are
#'   `"alpha"` for `"TRAJ"`, `"beta"` for `"TRBJ"`, and `"TCR"` when the chain
#'   type cannot be inferred.
#'
#' @examples
#' infer_chain_label("TRAJ")
#' infer_chain_label("TRBJ")
#' infer_chain_label(c("TRAJ", "TRBJ", "J"))
#'
#' @importFrom dplyr case_when
#'
#' @export
infer_chain_label <- function(j_type) {
  case_when(
    j_type == "TRAJ" ~ "alpha",
    j_type == "TRBJ" ~ "beta",
    TRUE ~ "TCR"
  )
}
