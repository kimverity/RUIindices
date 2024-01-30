#' convert a data frame tree to igraph
#' @param dft data frame tree
#' @return igraph object and node attribute data frame
#' @export
dfToIgraph <- function(dft) {
  id_map <- as.numeric(factor(dft$Identity))
  parent_map <- match(dft$Parent, dft$Identity)
  edge_matrix <- cbind(parent_map, id_map)
  edge_matrix <- edge_matrix[edge_matrix[, 1] != edge_matrix[, 2], ]
  if ("Generation" %in% names(dft) &&
      !is.null(dft$Generation)) {
    edge_lengths <- dft$OriginTime - dft$Generation
  } else if ("OriginTime" %in% names(dft) &&
             !is.null(dft$OriginTime)) {
    edge_lengths <- dft$OriginTime[id_map] - dft$OriginTime[parent_map]
  } else {
    edge_lengths <- rep(1, nrow(edge_matrix))
  }
  g <- igraph::graph_from_edgelist(as.matrix(edge_matrix), directed = TRUE)
  igraph::E(g)$weight <- edge_lengths
  nadf <- data.frame(names = dft$Identity, values = dft$Population)
  return(list(g, nadf))
}

#' convert igraph to Newick
#' @param g igraph object
#' @return Newick string
#' @export
igraphToNewick <- function(g) {
  if (!igraph::is_tree(g)) {
    stop("igraph object is not a tree")
  }
  if (!igraph::is_directed(g)) {
    stop("igraph object is directed")
  }
  buildNewick <- function(node) {
    children <- igraph::neighbors(g, node, mode = "out")
    if (length(children) == 0) {
      return(as.character(node))
    } else {
      childStrings <- sapply(children, buildNewick, USE.NAMES = FALSE)
      return(paste0("(", paste(childStrings, collapse = ","), ")",
                    as.character(node)))
    }
  }
  root <- igraph::V(g)[igraph::degree(g, mode = "in") == 0]
  if (length(root) != 1) {
    stop("The tree does not have a unique root")
  }
  newickString <- buildNewick(root)
  return(paste0(newickString, ";"))
}

#' convert data frame tree to Newick
#' @param dft data frame tree
#' @return Newick string and node attribute data frame
#' @export
dfToNewick <- function(dft) {
  treelist <- dfToIgraph(dft)
  tree <- igraphToNewick(treelist[[1]])
  return(list(tree, treelist[[2]]))
}
