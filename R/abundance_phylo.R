#' Calculates the total abundance descending from each branch of a tree
#' 
#' @param tree A 'phylo' object
#' @param abundance_data A dataframe containing node abundance data for tree;
#' the first column must contain the node labels and second column the corresponding 
#' node abundance data; f no abundance data is known, don't input a dataframe and
#' code will assign leaves to be equally abundant and internal nodes to have size
#' zero
#' @returns A library, key is number assigned to node in 'phylo' object
#' 
#' @examples
#' tree <- ape::read.tree(text='(a:1,(b:0.5,c:0.5)e:0.5)d:0;')
#' tree_abundance <- data.frame("Names" = c("e", "d", "b", "c", "a"), 
#'                              "Values" = c(0, 0, 1/3, 1/3, 1/3))
#' abundance_phylo(tree, tree_abundance)
#' 
#' @export
abundance_phylo <- function(tree, abundance_data) {
  
  node_labels <- append(tree$tip.label, tree$node.label) # all node labels
  
  # Initialize a dictionary to store the total abundance for each branch
  total_abundance <- list()
  
  # Define a recursive function to calculate total abundance
  calculate_abundance <- function(node) {
    
    # Initialize the total abundance for the current node
    node_abundance <- abundance_data[which(abundance_data[,1] == node_labels[node]), 2]
    
    # Loop through each child of the current node
    for (child in tree$edge[tree$edge[, 1] == node, 2]) {
      
      # If the child is not a leaf, recursively calculate the total abundance of its subtree
      if (!(child %in% c(1:length(tree$tip.label)))) {
        child_abundance <- calculate_abundance(child)
        
      }
      # If the child is a leaf, use its abundance as the total abundance of its subtree
      else {
        # Store child/leaf abundance
        child_abundance <- abundance_data[which(abundance_data[,1] == node_labels[child]), 2]
        # Total abundance of leaf is its own abundance
        total_abundance[[as.character(child)]] <<- abundance_data[which(abundance_data[,1] == node_labels[child]), 2]
      }
      
      # Add the total abundance of the child subtree to the current node's total abundance
      node_abundance <- node_abundance + child_abundance
    }
    
    # Store the total abundance for the current node
    total_abundance[[as.character(node)]] <<- node_abundance
    
    # Return the total abundance for the current node
    return(node_abundance)
  }
  
  # Start the recursive calculation from the root node
  calculate_abundance(tree$edge[1, 1])
  
  
  # Return the dictionary of total abundances
  return(total_abundance)
}