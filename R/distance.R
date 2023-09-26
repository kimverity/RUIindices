#' Calculates the  distance between two nodes
#' 
#' Calculates the distance between two specified nodes in terms of given
#' branch lengths. 
#' 
#' @param tree A 'phylo' object
#' @param top_node A number; corresponds to number assigned in 'phylo' object;
#' must be ancestor of bottom_node
#' @param bottom_node A number; corresponds to number assigned in 'phylo' object;
#' must be a descendant of top_node
#' @returns A number
#' 
#' @examples
#' tree <- ape::read.tree(text='(a:1,(b:0.5,c:0.5)e:0.5)d:0;')
#' distance(tree,4,1)
#' 
#' @export
distance <- function(tree, top_node, bottom_node){
  curr_node <- bottom_node # Select bottom node as starting node
  dist <- 0 # Initialise distance sum
  while (curr_node != top_node) { # Do until reach desired top node
    # Each node only has one parent branch,
    # select that branch, add the distance between node and parent,
    # move to parent node and repeat
    dist <- dist + tree$edge.length[which(tree$edge[,2] == curr_node)] # Sum distance
    curr_node <- tree$edge[which(tree$edge[,2] == curr_node),1] # Move to parent node
  }
  return(dist)
}