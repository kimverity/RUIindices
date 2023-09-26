#' Calculate S_i_a for a given node i and ancestor a
#' 
#' @param tree A 'phylo' object
#' @param node A number; corresponds to number assigned in 'phylo' object
#' @param abundances Dictionary containing node abundances; output of [abundance_phylo()]
#' @param curr_ancestor A number; corresponds to ancestor number assigned 
#' in 'phylo' object
#' @param h A number; distance between node and ancestor
#' @param l_i A number; length of longest direct descendant branch of node
#' @returns A list with two elements; first element is a dataframe containing
#' value of S_i_a and corresponding value of x; second element is a dictionary
#' where each element contains every branch abundance present in each region of
#' x
#' 
#' @examples
#' tree <- ape::read.tree(text='(a:1,(b:0.5,c:0.5)e:0.5)d:0;')
#' tree_abundance <- data.frame("Names" = c("e", "d", "b", "c", "a"), 
#'                              "Values" = c(0, 0, 1/3, 1/3, 1/3))
#' abundances <- abundance_phylo(tree, tree_abundance)
#' h <- distance(tree, 4,5)
#' l_i <- max(tree$edge.length[which(tree$edge[,1] == 5)])
#' calculate_S_i_a(tree,5,abundances,4,h,l_i)
#' 
#' @export
calculate_S_i_a <- function(tree, node, abundances, curr_ancestor, h, l_i) {
  
  temp <- function(a,b){return(a[[b]])}
  
  # Select all descendant branches from ancestor
  if (node != curr_ancestor){ # If not considering itself as ancestor  
    df_node_info <- get_descendant_branches(tree, curr_ancestor) # Get descendant branches from ancestor
    df_node_info <- df_node_info[df_node_info$x > h,] # Select branches that pass node
    
    # Want all branch lengths and x values to be measure from node, not ancestor
    # Select index values of nodes that have descendant branches that start before
    # given node, then adjust x and branch lengths accordingly
    curr_index <- which(df_node_info[,"x"] - df_node_info[,"Branch_length"] < h)
    # Correct branch lengths
    df_node_info[curr_index,"Branch_length"] <- df_node_info[curr_index,"x"] - h
    # Correct x
    df_node_info$x <- df_node_info$x - h
  }else{ # Considering itself as ancestor
    df_node_info <- get_descendant_branches(tree, node) # No corrections needed
  }
  
  df_node_info_ed <- df_node_info # Create dataframe to delete values from
  
  # A branch is in the current region iff its x and branch length are equal, 
  # otherwise it doesn't start at the beginning of current region.
  # Sum over all branches whose x and branch length are equal, deleting the 
  # branch/es corresponding with smallest branch length/x, then transform
  # all x and branch lengths by - length of previous region 
  # i.e. shifting then to be measured from the started of new region
  DF_S_i <- data.frame("S_i" = numeric(), "x" = numeric()) # Empty dataframe
  abund_list <- list() # Empty dictionary for abundances
  x <- 0 # Keep track of distance from node
  # Sum over each branch as all could have different branch lengths
  for (j in 1:nrow(df_node_info)) {
    # Only do until reach the end of nodes longest direct descendant branch, l_i
    if (x != l_i){
      list_ab <- c() # Empty list to hold branch abundances on each iteration
      # Row index of branches where x == branch length
      index_curr_branches <- which(df_node_info_ed[,"x"] %==% df_node_info_ed[,"Branch_length"])
      
      if (nrow(df_node_info_ed) != 0){ # If there are branches left to sum over
        # Store all branch abundances present for this value of x
        abund_list[[as.character(j)]] <- sapply(as.character(end_node <- df_node_info_ed[index_curr_branches,"End_node"]),
                                                temp, a = abundances)
        
        abundance_sum <- sum(abund_list[[as.character(j)]]) # Sum abundances
        
        # Append abundance and corresponding value of x
        DF_S_i <- rbind(DF_S_i, data.frame("S_i" = abundance_sum,
                                           "x" = (min(df_node_info_ed[index_curr_branches,]["x"]) + x)))
        
        prev_x <- min(df_node_info_ed[index_curr_branches,]["x"]) # Update previous x
        x <- x + prev_x # Update x
        
        # Indices of rows to be deleted 
        index_to_be_deleted <- which(df_node_info_ed[,"x"] %==% min(df_node_info_ed[index_curr_branches,]["x"]))
        df_node_info_ed$x <- df_node_info_ed$x - prev_x # "Reset" x = 0 level
        # Shorten current branch lengths by previous x
        df_node_info_ed[index_curr_branches,]["Branch_length"] <- df_node_info_ed[index_curr_branches,]["Branch_length"] - prev_x
        # Remove branch/s
        df_node_info_ed <- df_node_info_ed[-index_to_be_deleted,]
        
      }
    }
  }
  
  return(list(DF_S_i, abund_list))
}