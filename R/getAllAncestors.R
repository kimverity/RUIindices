#' Calculates all ancestor nodes of a given node including itself
#' 
#' @param tree A 'phylo' object
#' @param node A number; corresponds to number assigned in 'phylo' object
#' @returns Vector of ancestor node numbers corresponding to numbers 
#' assigned in 'phylo' object
#' 
#' @examples
#' tree <- ape::read.tree(text='(a:1,(b:0.5,c:0.5)e:0.5)d:0;')
#' getAllAncestors(tree,1)
#' 
#' @export
getAllAncestors <- function(tree,node){
  getAncestor <- function(tree, node){ # Returns immediate ancestor of given node
    # Each node only has one parent branch, and the parent branch ends in the node, 
    # select the node at the start of this branch to get parent node
    i <- which(tree$edge[, 2] == node)
    return(tree$edge[i, 1])
  }
  
  ancestors <- c(node) # List for ancestor nodes, include node as an ancestor of itself
  root_node <- length(tree$tip.label) + 1 # Number assigned to root node
  anc_node <- root_node + 1 # Assign so while loop will start
  if (node != root_node){ # Root node itself is its only ancestor
    while (anc_node > root_node){ # Run until reach root node
      anc_node <- getAncestor(tree, node) # select parent of node
      ancestors <- append(ancestors, anc_node) # Append parent to list
      node <- anc_node # Move to parent
    }
  }
  return(ancestors)
}