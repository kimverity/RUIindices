#' Calculates all descendant branches of a node
#' 
#' Calculates all descendant branches of a node and records distance from node or "x"
#' 
#' @param tree A tree as a 'phylo' object
#' @param node A number of node corresponding to number assigned in 'phylo'
#' object
#' @returns A dataframe containing for every descendant branch, the start and 
#' end corresponding to numbers assigned in 'phylo' object, branch length and 
#' distance from node to end of descendant branch length
#' 
#' @examples
#' tree <- ape::read.tree(text='(a:1,(b:0.5,c:0.5)e:0.5)d:0;')
#' get_descendant_branches(tree,4)
#'  
#' @export
get_descendant_branches <- function(tree, node){
  root <- FALSE
  if (node == (length(tree$tip.label)+1)){ # Is the root node
    root <- TRUE
    branches <- which(tree$edge[,1] == node) # Select row index of direct descendant branches
    if (length(branches) == 0) { # If node is leaf return empty matrix
      return(data.frame("Start_node" = numeric(), "End_node" = numeric(),
                        "Branch_length" = numeric(), "x" = numeric()))}
    else{
      descendants <- c() # Empty array to store descendants
      # Select row index in tree$edge for all branches that descend from node
      for (i in 1:length(branches)){
        descendants <- append(descendants, 
                              which(descendant_edges(edge = branches[i],
                                                    tree$edge[,1], 
                                                    tree$edge[,2])))
      }
    }  
  }else{ # Not root node
    # Calculate descendant branches from branch connecting node to parent
    # then delete this branch after
    branch <- which(tree$edge[,2] == node) # Select row index of parent branch
    if (length(which(node == c(1:length(tree$tip.label))))) { # If node is leaf return empty matrix
      return(data.frame("Start_node" = numeric(), "End_node" = numeric(),
                        "Branch_length" = numeric(), "x" = numeric()))}
    else{
      descendants <- which(descendant_edges(tree$edge[,1], tree$edge[,2], edge = branch))
    }
  }
    sub_tree <- tree$edge[descendants,] # Select subtree from node
    x <- sapply(sub_tree[,2], distance, top_node=node, tree=tree) # Record distance between node and end of each branch
    
    # x is the distance from node to the end of the branch 
    # i.e. x and branch length are only equal for direct descendant branches
    df_branch_info <- data.frame("Start_node" = sub_tree[,1], "End_node" = sub_tree[,2],
                                 "Branch_length" = tree$edge.length[descendants], 
                                 "x" = x) # Store information

  if (root == FALSE){ # Delete parent branch
    df_branch_info <- df_branch_info[-1,]
  }
  return(df_branch_info)
}