#' Calculate S_i_a for tree transformed to a star tree
#' 
#' Transforms tree into a star tree by attaching all of the trees branches to 
#' the root node and then calculates S_i_a 
#' 
#' @param tree A 'phylo' object
#' @param node A number; corresponds to number assigned in 'phylo' object
#' @param abundance_data Dictionary containing node abundances; output of 
#' [abundance_phylo()]
#' @returns A list with two elements; first element is a dataframe containing
#' value of S_i_a and corresponding value of x; second element is a dictionary
#' where each element contains every branch abundance present in each region of
#' x
#' 
#' @examples
#' tree <- ape::read.tree(text='(a:1,(b:0.5,c:0.5)e:0.5)d:0;')
#' tree_abundance <- data.frame("Names" = c("e", "d", "b", "c", "a"), 
#'                              "Values" = c(0, 0, 1/3, 1/3, 1/3))
#' calculate_S_i_a_star(tree,4,tree_abundance)
#' 
#' @export
calculate_S_i_a_star <- function(tree, node, abundance_data) {
  
  temp <- function(a,b){return(a[[b]])}
  
  abundances <- abundance_phylo(tree, abundance_data) # Node/branch abundances
  
  df_node_info <- get_descendant_branches(tree, node) # Get descendant branches from ancestor
  # Set all x to be equal to branch length, i.e. all branches are attached to
  # root node
  df_node_info$x <- df_node_info$Branch_length
  df_node_info_ed <- df_node_info # Create node to delete data from
  
  # Sum abundances of all branches in df_node_info_ed, delete the branch/es 
  # corresponding with shortest branch length/x, repeat until all branches
  # have been summed over (i.e. deleted from dataframe)
  Length <- length(unique(df_node_info$x)) # Preassign length
  df_list <- vector(mode = "list", length = Length) # List for dataframes
  abund_list <- vector(mode = "list", length = Length) # Empty dictionary for abundances
  names(abund_list) <- as.character(c(1:Length)) # Change names to characters
  # Sum over each unique branch length
  for (j in 1:Length) {
    if (nrow(df_node_info_ed) != 0){ # Checks if any branches are left
      end_nodes <- df_node_info_ed$End_node # End nodes of current branches
      # Store all branch abundances present for this value of x
      abund_list[[as.character(j)]] <- sapply(as.character(end_nodes), temp, a = abundances)
      
      abundance_sum <- sum(abund_list[[as.character(j)]]) # Sum abundances 
      
      # Append abundance and corresponding value of x
      df_list[j] <- list(data.frame("S_i" = abundance_sum,
                                    "x" = (min(df_node_info_ed$x))))
      # Remove branch/s corresponding to x value
      df_node_info_ed <- df_node_info_ed[-(which(df_node_info_ed$x %==% min(df_node_info_ed$x))),]
    }
  }
  abund_list <- Filter(Negate(is.null), abund_list)
  # Combine dataframes
  DF_S_i <- Reduce(rbind, df_list)
  return(list(DF_S_i, abund_list))
}
