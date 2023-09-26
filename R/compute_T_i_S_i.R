#' Calculates T_i and S_i for a given node i
#' 
#' @param tree A 'phylo' object
#' @param node A number; corresponds to number assigned in 'phylo' object
#' @param abundances Dictionary containing node abundances; output of [abundance_phylo()]
#' @returns A list with 2 elements; first element is a number corresponding to T_i:
#' second element is a dataframe containing value of S_i and corresponding value 
#' of x
#' 
#' @examples
#' tree <- ape::read.tree(text='(a:1,(b:0.5,c:0.5)e:0.5)d:0;')
#' tree_abundance <- data.frame("Names" = c("e", "d", "b", "c", "a"), 
#'                              "Values" = c(0, 0, 1/3, 1/3, 1/3))
#' abundances <- abundance_phylo(tree, tree_abundance)
#' compute_T_i_S_i(tree, 4, abundances)
#' 
#' @export
compute_T_i_S_i <- function(tree, node, abundances) {
  
  # Choose branches that are immediate descendants of node and record data in dataframe
  df_descen_info <- data.frame("Start_node" = tree$edge[tree$edge[, 1] == node, 1],
                               "End_node" = tree$edge[tree$edge[, 1] == node, 2],
                               "Branch_lengths" = tree$edge.length[which(tree$edge[, 1] == node)])
  # x is distance from root node to end of interval
  DF_S_i <- data.frame("S_i" = numeric(), "x" = numeric()) # Empty dataframe
  
  # Sum abundances of all branches in df_descen_info, delete the branch/es 
  # corresponding with shortest branch length, repeat until all branches
  # have been summed over (i.e. deleted from dataframe)
  T_i <- 0 # Initialise T_i sum
  prev_x <- 0 # Keep track of previous value of x
  # Sum over each unique branch length
  for (j in 1:length(unique(df_descen_info[,"Branch_lengths"]))) {
    abundance_sum <- 0 # Initialise sum
    if (nrow(df_descen_info) != 0){ # Checks if any branches are left
      for (i in 1:nrow(df_descen_info)) { # Sum over all branches in current region
        if (!(is.na(df_descen_info["End_node"][i,]))){
          end_node <- df_descen_info["End_node"][i,] # Select end node of current branch
          abundance_sum <- abundance_sum + abundances[[as.character(end_node)]] # Sum abundances
        }
      }
    }
    
    # Record all info and delete any branches "passed"
    if (nrow(df_descen_info) != 0){ # If there are branches left
      # Append abundance and corresponding value of x
      DF_S_i <- rbind(DF_S_i, 
                      data.frame("S_i" = abundance_sum,
                                 "x" = min(df_descen_info["Branch_lengths"])))
      
      T_i <- T_i + (abundance_sum * (min(df_descen_info["Branch_lengths"]) - prev_x)) # Sum T_i
      prev_x <- min(df_descen_info["Branch_lengths"]) # Update previous x
      
      df_descen_info <- df_descen_info[-(which(df_descen_info["Branch_lengths"] == min(df_descen_info["Branch_lengths"]))),] # Remove branch/s corresponding to x value
    }
  }
  
  return(list(T_i, DF_S_i))
}