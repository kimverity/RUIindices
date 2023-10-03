#' Calculates longitudinal or star mean 
#'
#' Calculates the longitudinal or star mean for D for q=0 or 1, or J for q =1,
#' or all for one mean type. The default, changed by specifying index letter, q 
#' and individual, is to calculate all index values for mean type.
#' 
#' @param file A string or phylo object; String of file name containing tree or 
#' the tree itself either as a string or a phylo object
#' @param node_abundances A dataframe containing node abundance data for tree;
#' the first column must contain the node labels and second column the corresponding 
#' node abundance data; if no abundance data is known, don't provide an input and
#' code will assign leaves to be equally abundant and internal nodes to have size
#' zero
#' @param mean_type A string; string containing desired mean type, either "star"
#' or "longitudinal"
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
#' long_star(tree, tree_abundance, "Longitudinal", "D", 1, TRUE)
#' 
#' @export
long_star <- function(file, node_abundances = FALSE, mean_type, index_letter = "D", q = 1,
                      individual = FALSE){
  
  index_letter <- toupper(index_letter) # Capitalise input
  mean_type <- toupper(mean_type) # Capitalise input
  
  tree <- read_convert(file) # Convert tree to phylo object
  node <- tree$edge[1,1] # Select root node
  
  # Checks if tree has abundance data
  # If it does it calculates branch/node abundance data.
  # If it doesn't, it assign leaves to be equally abundant and internal nodes 
  # to have size zero
  if (is.data.frame(node_abundances)) { # Tree has abundance data
    abundances <- abundance_phylo(tree, node_abundances) # Calculate branch/node abundances
  }else if (!(is.data.frame(node_abundances))) { # Tree doesn't have abundance data
    num_tips <- tree$edge[1,1] - 1 # Number of tips
    tree$node.label<- as.character(c(1:num_tips)) # Assign node labels
    num_nodes <- tree$Nnode # Number of nodes
    # Assign node labels
    tree$node.label <- as.character(c((num_tips + 1):(num_tips + num_nodes)))
    
    # Create abundance dataframe
    node_abundances <- data.frame("names" = c(tree$node.label, tree$tip.label),
                                  "values" = rep(c(0, (1/num_tips)), times=c(tree$Nnode, num_tips)))
    abundances <- abundance_phylo(tree, node_abundances) # Calculate branch/node abundances
  }
  
  # Selects mean type and runs corresponding function for S_i_a
  if (mean_type == "LONGITUDINAL") {
    S_i_a_res <- calculate_S_i_a(tree, node, abundances, node, 0, sum(tree$edge.length)) # Run function
  }else if (mean_type == "STAR") {
    S_i_a_res <- calculate_S_i_a_star(tree, node, node_abundances) # Run function
  }
  
  # Calculate index value/s
  # For each region of x in df_S_i_a, calculates index/indices using corresponding
  # abundance list, this contains every branch abundance in this region, 
  # index is calculated and summed over every region of x
  df_S_i_a <- S_i_a_res[[1]] # Select dataframe
  abund_list <- S_i_a_res[[2]] # Select abundance list
  # Initialise sum/s
  if (individual == TRUE) { # One index
    h <- 0
  }else if (individual == FALSE) {# All indices
    h1 <- 0
    h0 <- 0
    j1 <- 0
  }
  T_S_sum <- 0 # Initialise sum
  prev_x <- 0 # Keep track of x
  for (k in 1:length(abund_list)) { # Sum across every region of x 
    S <- df_S_i_a[k,"S_i"] # Select S_i_a
    x <- df_S_i_a[k, "x"] - prev_x # Assign interval length
    abund_vec <- unlist(abund_list[k], use.names = FALSE) # Branch abundances present in interval
    T_S_sum <- T_S_sum + (S * x) # H_bar
    #Calculate index values
    if (individual == TRUE) { # One index
      if (index_letter == "J") { # J for q = 1
        if (length(abund_vec) != 1) { # More than one branch in region
          h <- h + (sum(-(abund_vec) * log((abund_vec / S),
                                           base = length(abund_vec))) * x)
        }else if (length(abund_vec) == 1) {# Only one branch in region
          h <- h + (1 * S * x)
        }
      }
      # For index letter "D", need to check which index is wanted i.e. what q
      if (q == 1) {
        h <- h + (sum(-(abund_vec) * log((abund_vec / S))) * x)
      }else if (q == 0) {
        h <- h + (S * x * log(length(abund_vec), base = exp(1)))
      }
    }else if (individual == FALSE) { # All indices
      h1 <- h1 + (sum(-(abund_vec) * log((abund_vec / S))) * x)
      h0 <- h0 + (S * x * log(length(abund_vec), base = exp(1)))
      if (length(abund_vec) != 1) { # More than one branch in region
        j1 <- j1 + (sum(-(abund_vec) * log((abund_vec / S),
                                           base = length(abund_vec))) * x)
      }else if (length(abund_vec) == 1) { # One branch in region
        j1 <- j1 + (1 * S * x)
      }
    }
    prev_x <- df_S_i_a[k, "x"] # Update x
  }
  # Normalise index/indices
  if (individual == TRUE) { # One index
    if (index_letter == "J") { # Index is J
      H <- (h / T_S_sum)
    }else{
      H <- exp((h / T_S_sum))
    }
  }else if (individual == FALSE) { # All indices
    if (mean_type == "STAR") { # Star mean indices
      H <- list("D0S" = exp((h0 / T_S_sum)), "D1S" = exp((h1 / T_S_sum)),
                "J1S" = j1 / T_S_sum)
    }else if (mean_type == "LONGITUDINAL") { # Lomgitudianl mean indices
      H <- list("D0L" = exp((h0 / T_S_sum)), "D1L" = exp((h1 / T_S_sum)),
                "J1L" = j1 / T_S_sum)
    }
  }
  return(H)
}