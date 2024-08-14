#' Calculates qD_i or 1J_i for node i
#' 
#' @param tree A 'phylo' object
#' @param node A number; corresponds to number assigned in 'phylo' object
#' @param abundances Dictionary containing node abundances; output of [abundance_phylo()]
#' @param index_letter A string; a string containing a single letter, either 
#' "D" or "J" to specify desired index, if calculating all indices input does 
#' not matter
#' @param q A number; either 0 or 1 specifying desired index; if calculating
#' J only q=1 will be returned, if calculating all indices input does not matter
#' @param individual A boolean; either TRUE if only one index is desired or
#' FALSE if all indices are desired
#' @returns Either; a single number for the specified index, or a vector of all
#' index values (ordered 1D, 1J, 0D)
#' 
#' @examples
#' tree <- ape::read.tree(text='(a:1,(b:0.5,c:0.5)e:0.5)d:0;')
#' tree_abundance <- data.frame("Names" = c("e", "d", "b", "c", "a"),
#'                              "Values" = c(0, 0, 1/3, 1/3, 1/3))
#' abundances <- abundance_phylo(tree, tree_abundance)
#' calculate_DJ_i(tree, 4, abundances, "D", 1, TRUE)
#' 
#' @export
calculate_DJ_i <- function(tree, node, abundances, index_letter, q, individual){
  
  ancestors <- c(TreeTools::ListAncestors(tree$edge[,1], tree$edge[,2], node), node) # List of ancestors of node
  
  S_i_all <- compute_T_i_S_i(tree, node, abundances) # Run function
  T_i <- S_i_all[[1]] # Select value of T_i
  S_i <- S_i_all[[2]] # Select S_i dataframe
  
  if (T_i != 0){ # Has descendant branch with branch length greater than zero
    # Calculate integral value/s
    DJ_i <- sapply(ancestors, calculate_integral, tree=tree, node=node, S_i=S_i,
                  abundances=abundances, index_letter = index_letter, q = q,
                  individual=individual)
    
    # Add all ancestor contributions and normalise
    if (individual == FALSE){
      dj_i <- apply(DJ_i,1,sum)
      D1_i <- (1/T_i) * dj_i[1]
      J1_i <- (1/T_i) *dj_i[2]
      D0_i <- (1/T_i) *dj_i[3]
      v <- c(D1_i, J1_i, D0_i) 
    }else if (individual == TRUE){
      dj_i <- sum(DJ_i)
      v <- (1/T_i) * dj_i
    }
  }else{
    if (individual == FALSE){
      v <- c(0, 0, 0) 
    }else if (individual == TRUE){
      v <- 0
    }
  }
  
  return(v)
}