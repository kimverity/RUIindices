#' Calculates 0D, 1D or 1J or all
#' 
#' Calculates 0D, 1D or 1J or all for given node i, ancestor a
#' 
#' @import fpCompare
#' 
#' @param tree A 'phylo' object
#' @param node A number; corresponds to number assigned in 'phylo' object
#' @param abundances Dictionary containing node abundances; output of [abundance_phylo()]
#' @param curr_ancestor A number; corresponds to ancestor number assigned 
#' in 'phylo' object
#' @param h A number; distance between node and ancestor
#' @param l_i A number; length of longest direct descendant branch of node
#' @param index_letter A string; a string containing a single letter, either 
#' "D" or "J" to specify desired index, if calculating all indices input does 
#' not matter
#' @param q A number; either 0 or 1 specifying desired index; if calculating
#' J only q=1 will be returned, if calculating all indices input does not matter
#' @param individual A boolean; either TRUE if only one index is desired or
#' FALSE if all indices are desired
#' @returns Either; a list of dataframes if individual is FALSE, each dataframe
#' contains an index value and the corresponding value of x, key is index name;
#' if individual is TRUE, returns just one of these dataframes for the specified
#' index
#' 
#' @examples
#' tree <- ape::read.tree(text='(a:1,(b:0.5,c:0.5)e:0.5)d:0;')
#' tree_abundance <- data.frame("Names" = c("e", "d", "b", "c", "a"), 
#'                              "Values" = c(0, 0, 1/3, 1/3, 1/3))
#' abundances <- abundance_phylo(tree, tree_abundance)
#' h <- distance(tree, 4,4)
#' l_i <- max(tree$edge.length[which(tree$edge[,1] == 4)])
#' calculate_DJ_i_a(tree,4,abundances,4,h,l_i, "D", 1, TRUE)
#' 
#' @export
calculate_DJ_i_a <- function(tree, node, abundances, curr_ancestor, h, l_i, 
                             index_letter, q, individual) {
  
  index_letter <- toupper(index_letter) # Capitalise input
  
  S_i_a_res <- calculate_S_i_a(tree,node,abundances, curr_ancestor, h, l_i) # Run function
  df_S_i_a <- S_i_a_res[[1]] # Select dataframe
  abund_list <- S_i_a_res[[2]] # Select abundance list
  
  # Create empty list for dataframes
  df_In_i_a_list <- vector(mode = "list", length = length(df_S_i_a$x))
  
  # Sum over every region of x
  for (k in 1:length(abund_list)){
    S <- df_S_i_a[k,"S_i"] # Select S_i_a
    # Select abundances in current region
    abund_vec <- unlist(abund_list[k], use.names = FALSE)
    
    if (individual == FALSE){ # Calculating all indices
      D0 <- log(length(abund_vec))
      D1 <- sum(-(abund_vec / S) * log((abund_vec / S)))
      # J is defined to be 1 when only one branch is present
      if (length(abund_vec) > 1 ) { # There is more than one branch in region
        J1 <- sum(-(abund_vec / S) * log((abund_vec / S), base=length(abund_vec)))
      }else if (length(abund_vec) == 1) { # There is one branch in region
        J1 <- 1
      }
      # Store values and corresponding x
      df_In_i_a_list[k] <- list(data.frame("D0" = D0, "D1" = D1, "J1" = J1,
                                           "x" = df_S_i_a[k, "x"]))
    }else if (individual == TRUE){ # Calculating one index
      # Checks desired index and calculates value
      if ((index_letter == "D")&(q == 0)) {
        ind_val <- log(length(abund_vec))
      }else if ((index_letter == "D")&(q == 1)) {
        ind_val <- sum(-(abund_vec / S) * log((abund_vec / S)))
      }else if (index_letter == "J") {
        if (length(abund_vec) > 1 ) { # There is more than one branch in region
          ind_val <- sum(-(abund_vec / S) * log((abund_vec / S), base=length(abund_vec)))
        }else if (length(abund_vec) == 1) { # There is one branch in region
          ind_val <- 1
        }
      }
      # Store value and corresponding x
      df_In_i_a_list[k] <- list(data.frame("In_i_a" = ind_val, "x" = df_S_i_a[k, "x"]))
    }
  }
  # Return either a dictionary of index values or single index
  if (individual == FALSE) {
    # Combine dataframes
    df_In_i_a <- Reduce(rbind, df_In_i_a_list)
    return(list("1DN" = df_In_i_a[,c("D1", "x")], "1JN" = df_In_i_a[,c("J1", "x")],
                "0DN" = df_In_i_a[,c("D0", "x")]))
  }else if (individual == TRUE){
    # Combine dataframes
    df_In_i_a <- Reduce(rbind, df_In_i_a_list)
    return(df_In_i_a)
  }
}