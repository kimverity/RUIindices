#' Reads in a linear tree and outputs it as a phylo object
#' 
#' Reads in a linear tree in newick format, either as a 
#' string or from a file and outputs it as a phylo object. 
#' If no branch lengthd are given they are assigned to be
#' one.
#' 
#' @param  file A string; String of the filename or the tree itself
#' @returns A 'phylo' object
#' @examples 
#' tree_string <- '(((1:1)4:1)3:1)2:;'
#' tree_phylo <- form_linear_phylo(tree_string)
#' 
#' @export 

form_linear_phylo <- function(file){
  
  # Read in tree
  suppressWarnings({
    tree <- try(scan(file = file, what = "", sep = "\n", quiet = TRUE, skip = 0,
                     comment.char = ""), silent = TRUE)
    if ("try-error" %in% class(tree)) { # If it is not a file, must tree as string
      tree <- file
    }
  })
  
  tree <- gsub("\\[[^]]*\\]", "", tree) # Delete any comments

  colon <- grep(":", tree) # check for colons in tree
  branch_lengths <- FALSE

  # Checks for all common tree forms and their variations
  # e.g. whether tree has branch lengths, brackets around
  # root node 
  if (length(colon) == 1) { # Tree has colons/branch lengths
    # Length of string
    str_len <- unlist(gregexpr(";", tree))
    # Check if there is a root node label
    # if no label assign empty one
    if (substring(tree, str_len - 1, str_len - 1) == ")") { # no root node label
      # Check if branch lengths are empty
      if (unlist(gregexpr(":)", tree)) == -1) { # Branch lengths not empty 
        # indices of right-hand brackets
        ind <- unlist(gregexpr("\\)", tree))
        # select node labels
        node_labels <- mapply(substring, ind - 3, ind - 3,  text = tree)
        # select branch lengths
        branch_lengths <- as.numeric(mapply(substring, ind - 1, ind - 1,  text = tree))
      }else { #  Empty
        # indices of right-hand brackets
        ind <- unlist(gregexpr("\\)", tree))
        # select node labels
        node_labels <- mapply(substring, ind - 2, ind - 2,  text = tree)
      }
      # Add root node, has no label so assign empty label
      node_labels <- append(node_labels, "")
    }else { # Root node label
      # Check if branch lengths are empty
      if (unlist(gregexpr(":)", tree)) == -1) { # Branch lengths not empty 
        # indices of right-hand brackets
        ind <- unlist(gregexpr("\\)", tree))
        # select node labels
        node_labels <- mapply(substring, ind - 3, ind - 3,  text = tree)
        # select branch lengths
        branch_lengths <- as.numeric(mapply(substring, ind - 1, ind - 1,  text = tree))
      }else { #  Empty
        ind <- unlist(gregexpr("\\)", tree))
        node_labels <- mapply(substring, ind - 2, ind - 2,  text = tree)
      }
      # Add root node
      node_labels <- append(node_labels, substring(tree, ind[length(ind)] + 1,
                                                   ind[length(ind)] + 1))
    }
  }else { # No branch lengths
    # Length of string
    str_len <- unlist(gregexpr(";", tree))
    # indices of right-hand brackets
    ind <- unlist(gregexpr("\\)", tree))
    # select node labels
    node_labels <- mapply(substring, ind - 1, ind - 1,  text = tree)
    # Check if there is a root node label
    # if not assign and empty label
    if (substring(tree, str_len - 1, str_len - 1) != ")"){ # Root node has label
      node_labels <- append(node_labels, substring(tree, ind[length(ind)] + 1,
                                                   ind[length(ind)] + 1))
    }else{ # No label, assign empty one
      node_labels <- append(node_labels, "")
    }
  }
  
  no_nodes <- length(node_labels) - 1 # Don't count leaf
  
  # Create edge matrix
  edge_matrix <- matrix(0, no_nodes, 2)
  for (i in 1:(no_nodes-1)) {
    edge_matrix[i,] <- c(i + 1, i + 2)
  }
  edge_matrix[no_nodes,] <- c(no_nodes + 1, 1)

  # Create phylo object
  phylotree <- list()
  phylotree$edge <- edge_matrix
  phylotree$tip.label <- node_labels[1]
  phylotree$node.label <- node_labels[-1]
  phylotree$Nnode <- no_nodes
  class(phylotree) <- "phylo"
  attr(phylotree, "order") <- "cladewise"
  
  # Assign branch lengths
  if (length(branch_lengths) == 1) { # No branch lengths, set to all be equal
    phylotree$edge.length <- rep(1, times = no_nodes)
  }else if (length(branch_lengths) > no_nodes) { # root node had branch length 
    phylotree$edge.length <- branch_lengths[-length(branch_lengths)]
  }else {
    phylotree$edge.length <- branch_lengths
  }
  
  return(phylotree)
}
