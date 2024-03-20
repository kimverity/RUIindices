#' Identify descendant edges
#'
#' Quickly identify edges that are "descended" from edges in a tree.
#' 
#' @param parent Vector containing first column of edge matrix
#' @param child Vector containing second column of edge matrox
#' @param edge Integer specifying the number of the edge whose children are
#' required (see \code{\link[ape:nodelabels]{edgelabels}()}).
#' @param nEdge number of edges (calculated from `length(parent)` if not
#' supplied).
#' @return `DescendantEdges()` returns a logical vector stating whether each
#' edge in turn is the specified edge or one of its descendants.
#' @examples
#' tree <- ape::read.tree(text='(a:1,(b:0.5,c:0.5)e:0.5)d:0;')
#' desc <- descendant_edges(tree$edge[, 1], tree$edge[, 2], edge = 1)
#' @family tree navigation
#' @export
descendant_edges <- function(parent, child, edge,
                            nEdge = length(parent)) {

  ret <- logical(nEdge)
  edgeSister <- match(parent[edge], parent[-edge])
  if ((!(is.na(edgeSister)))&(edgeSister >= edge)) {
    # Added check for linearity in tree
    # edgeSister is really 1 higher than you think, because we knocked out
    # edge "edge" in the match
    ret[edge:edgeSister] <- TRUE
        
    # Return:
    ret
  } else {
    nextEdge <- edge
    revParent <- rev(parent)
    repeat {
      if (revDescendant <- match(child[nextEdge], revParent, nomatch=FALSE)) {
        nextEdge <- 1 + nEdge - revDescendant
      } else break;
    }
    ret[edge:nextEdge] <- TRUE
        
    # Return:
    ret
    }
}