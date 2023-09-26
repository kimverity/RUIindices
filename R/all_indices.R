#' Calculates all index values
#' 
#' @param file A string or phylo object; String of file name containing tree or 
#' the tree itself either as a string or a phylo object
#' @param node_abundances A dataframe containing node abundance data for tree;
#' the first column must contain the node labels and second column the corresponding 
#' node abundance data; if no abundance data is known, don't provide an input and
#' code will assign leaves to be equally abundant and internal nodes to have size
#' zero
#' @returns A dictionary containing all index values (for the key, 
#' if index is 1DS key is D1S)
#' 
#' @examples
#' tree <- '(a:1,(b:0.5,c:0.5)e:0.5)d:0;'
#' tree_abundance <- data.frame("Names" = c("e", "d", "b", "c", "a"), 
#'                              "Values" = c(0, 0, 1/3, 1/3, 1/3))
#' all_indices(tree, tree_abundance)
#' 
#' @export
all_indices <- function(file, node_abundances = FALSE){
  # Calculates all indices
  node <- node(file, node_abundances, "D", 1, FALSE)
  star <- long_star(file, node_abundances, "Star", "D", 0, FALSE)
  long <- long_star(file, node_abundances, "Longitudinal", "D", 0, FALSE)
  
  # List of each index value
  values <- list("D0N"= node$D0N,"D1N" = node$D1N,"J1N" = node$J1N, "D0S" = star$D0S,
                 "D1S" = star$D1S, "J1S" = star$J1S, "D0L" = long$D0L, "D1L" = long$D1L,
                 "J1L" = long$J1L)
  return(values)
}