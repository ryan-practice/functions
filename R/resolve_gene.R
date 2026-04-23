#' Verify a gene matches our reference syntax
#'
#' @param input_gene A character vector.
#' @param chain "trav", "trbv", "traj", or "trbj"
#' @return A character vector in uppercase.
#' @export
resolve_gene <- function(input_gene, chain){
  if(!(chain %in% c("trav", "trbv", "traj", "trbj"))){
    stop("must be exactly trav, trbv, traj, trbj")
  }
  test_chain <- as.data.table(get(chain))
  if(chain == "trav" | chain == "trbv"){
    col_name <- "Gene"
  }
  if(chain == "traj"){
    col_name <- "TRAJ"
  }
  if(chain == "trbj"){
    col_name <- "TRBJ"
  }
  input_gene <- toupper(input_gene)
  resolved_gene <- "unresolved"
  if(!(input_gene %in% test_chain[,..col_name][[1]])){

    #strip anything after an asterisk
    input_gene <- sub("\\*.*$", "", input_gene)
    #remove internal spaces
    input_gene <- gsub("\\s+", "", input_gene)
    #remove flanking whitespace
    input_gene <- trimws(input_gene)
    #remove anything after a slash

    #match exactly in trav form because I don't want to accidentally grep trav10 when the trav gene is trav1 for example
    if(length(grep(paste0("^", input_gene, "$"), test_chain[,..col_name][[1]])) != 1){
      stop(paste0("Invalid ", chain, " nomenclature in input file"))
    }else{
      resolved_gene <- test_chain[grep(paste0("^", input_gene, "$"), test_chain[,..col_name][[1]]), ..col_name][[1]]
    }
  }
  if(input_gene == "TRAJ24*02"){
    return("TRAJ24*02")
  }
  return(resolved_gene)
}
