#' Calculates node-wise mean 
#' 
#' Calculates the node-wise mean for D for q = 0 or 1, or for J for q = 1, or all.
#' The default, changed by specifying index type, q and individual, is to 
#' calculate all index values for mean type.
#' 
#' @param file A string or phylo object; String of file name containing tree or 
#' the tree itself either as a string or a phylo object
#' @param node_abundances A dataframe containing node abundance data for tree;
#' the first column must contain the node labels and second column the corresponding 
#' node abundance data; if no abundance data is known, don't provide an input and
#' code will assign leaves to be equally abundant and internal nodes to have size
#' zero
#' @param index_letter A string; a string containing a single letter, either 
#' "D" or "J" to specify desired index, if calculating all indices input 
#' does not matter
#' @param q A number; either 0 or 1 specifying desired index; if calculating
#' J only q=1 will be returned, if calculating all indices input does not matter
#' @param individual A boolean; either TRUE if only one index is desired or
#' FALSE if all indices are desired
#' @returns Either; single number for the chosen index or a dictionary with each
#' index value (for the key, if index is 1DS key is D1S)
#' 
#' @examples
#' tree <- '(a:1,(b:0.5,c:0.5)e:0.5)d:0;'
#' tree_abundance <- data.frame("Names" = c("e", "d", "b", "c", "a"),
#'                              "Values" = c(0, 0, 1/3, 1/3, 1/3))
#' node(tree, tree_abundance, "D", 1, TRUE)
#' 
#' @export
node <- function(file, node_abundances = FALSE, index_letter = "D", q = 1,
                 individual = FALSE){
  
  tree <- read_convert(file) # Capitalise inputs
  index_letter <- toupper(index_letter) # Capitalise inputs
  
  # Checks if tree has abundance data
  # If it does it calculates branch/node abundance data.
  # If it doesn't, it assign leaves to be equally abundant and internal nodes 
  # to have size zero
  if (is.data.frame(node_abundances)){ # Tree has abundance data
    abundances <- abundance_phylo(tree, node_abundances) # Calculate branch/node abundances
  }else if (!(is.data.frame(node_abundances))){ # Tree doesn't have abundance data
    num_tips <- tree$edge[1,1] - 1 # Number of tips
    tree$tip.label<- as.character(c(1:num_tips)) # Assign node labels
    num_nodes <- tree$Nnode # Number of nodes
    # Assign tip labels
    tree$node.label<- as.character(c((num_tips+1):(num_tips + num_nodes)))
    
    # Create abundance dataframe
    node_abundances <- data.frame("names" = c(tree$node.label, tree$tip.label),
                                  "values" = rep(c(0, (1/num_tips)), times=c(tree$Nnode, num_tips)))
    abundances <- abundance_phylo(tree, node_abundances) # Calculate branch/node abundances
  }
  
  # Calculates h_bar
  T <- 0 # Initialise sum
  for (i in 1:length(tree$edge[,1])){
    T <- T + (abundances[[as.character(tree$edge[i,2])]] * tree$edge.length[i])
  }
  
  # Numbers of every internal node, corresonding to numbers in phylo object
  nodes <- c((length(tree$tip.label)+1):(length(tree$tip.label) + tree$Nnode))
  
  # Calculates node-averaged indice/s
  if (individual == FALSE){# All indices
    # Calculates the indices values for each node
    DJ <- sapply(nodes, calculate_DJ_i, tree=tree, abundances = abundances, 
                 individual = individual, q = q, index_letter = "D")
    # Normalise the indices
    D1N <- (1/T)*sum(DJ[1,])
    J1N <- (1/T)*sum(DJ[2,])
    D0N <- (1/T)*sum(DJ[3,])
    # List of index values
    List <- list("D1N"= exp(D1N),"J1N" = J1N,"D0N" = exp(D0N))
    return(List)
  }else if (individual == TRUE){ # One index
    # Calculates the index value for each node
    DJ <- sapply(nodes, calculate_DJ_i, tree=tree, abundances = abundances, 
                 individual = individual, index_letter = index_letter, q = q)
    # Normalise
    if (index_letter == "J"){ # Index J
      index <- (1/T)*sum(DJ)
    }else if (!(index_letter == "J")){ # Index D
      index <- exp((1/T)*sum(DJ))
    }
    
    return(index)
  }
}