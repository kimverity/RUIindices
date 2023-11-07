#' Calculates integral for each ancestor
#' 
#' @param tree A 'phylo' object
#' @param node A number; corresponds to number assigned in 'phylo' object
#' @param curr_ancestor A number; corresponds to ancestor number assigned 
#' in 'phylo' object
#' @param S_i A dataframe; dataframe output from [compute_T_i_S_i()]
#' @param abundances Dictionary containing node abundances; output of [abundance_phylo()]
#' @param index_letter A string; a string containing a single letter, either 
#' "D" or "J" to specify desired index, if calculating all indices input does 
#' not matter
#' @param q A number; either 0 or 1 specifying desired index; if calculating
#' J only q=1 will be returned, if calculating all indices input does not matter
#' @param individual A boolean; either TRUE if only one index is desired or
#' FALSE if all indices are desired
#' @returns Either; a single number for the desired index, or a list of all indices
#' 
#' @examples
#' tree <- ape::read.tree(text='(a:1,(b:0.5,c:0.5)e:0.5)d:0;')
#' tree_abundance <- data.frame("Names" = c("e", "d", "b", "c", "a"), 
#'                              "Values" = c(0, 0, 1/3, 1/3, 1/3))
#' abundances <- abundance_phylo(tree, tree_abundance)
#' S_i <- compute_T_i_S_i(tree, 4, abundances)[[2]]
#' calculate_integral(tree, 4, 4, S_i, abundances, "D", 1, FALSE)
#' 
#' @export
calculate_integral <- function(tree, node, curr_ancestor, S_i, abundances, index_letter,
                               q, individual){
  
  # Set l_i, longest direct descendant branch of node
  if (length(S_i["x"] > 1)){
    l_i <- max(S_i["x"]) # Select longest branch length of a immediate descendant
  }else{
    l_i <- S_i["x"] # If only one value, select that
  }
  
  # If node is a leaf/has size zero
  if (l_i == 0){
    if (individual == TRUE){
      return(0)
    }else{
      return(c(0,0,0))
    }
  }
  
  h <- distance(tree,curr_ancestor,node) # Distance between current node and ancestor
  
  # Assign distance to parent
  if (curr_ancestor == (length(tree$tip.label) + 1)){ # If ancestor is the root node
    d_parent <- l_i  # By definition its Inf, but integral is only nonzero to l_i
  }else{
    d_parent <- tree$edge.length[which(tree$edge[,2] == curr_ancestor)] + h
  }
  
  # Ancestor integral
  # Calculates integral as sums of areas
  int_2 <- 0 # Initialise sum
  prev_x <- h # Keep track of x, h as integral is from h
  if (l_i <= h){ # Nodes furthest away child is closer than distance to ancestor
    int_2 <- 0
    if (individual == TRUE){
      return(0)
    }else{
      return(c(0,0,0))
    }
  }else{
    current_row <- which(S_i[,"x"] > h)[1] # Start sum over first x that "reaches" ancestor
    while (current_row <= length(S_i[,"x"])) { # Sum over all rows of S_i
      x_s <- S_i[current_row,2] # Current value of x for S_i
      value_s <- S_i[current_row,1] # Current value of S_i for given x
      if (x_s > d_parent){ # Region reaches past ancestor's parent
        int_2 <- int_2 + (value_s * (d_parent - prev_x)) # Integral
        prev_x <- x_s # Update previous x
        break # Exit while loop
      }else { # x_s < d_parent
        int_2 <- int_2 + (value_s * (x_s - prev_x)) # Sum integral
        prev_x <- x_s # Update previous x
        current_row <- current_row + 1 # Move to next row
      }
    }
  }
  
  # Function for index integral
  # Calculates integral as sums of areas
  integral <- function(df){
    prev_x <- 0 # Keep track of x
    current_row_s <- 1 # Start sum over rows of S_i
    current_row_i <- 1 # Start sum over rows of index value
    int_1 <- 0 # Initialise integral sum
    
    # qD/1J integral
    # Sum over all rows of S_i
    while ((current_row_s <= length(S_i[,"x"]))&(current_row_i <= length(df$x))){
      x_s <- S_i[current_row_s,2] # Current value of x for S_i
      x_i <- df[current_row_i, 2] # Current value of x for qD/1J
      value_s <- S_i[current_row_s,1] # Current value of S_i for given x
      value_i <- df[current_row_i,1] # Current value of qD/1J for given x
      if (x_s < x_i){
        int_1 <- int_1 + (value_s * value_i * (x_s - prev_x)) # Sum integral
        prev_x <- x_s # Update previous x
        current_row_s <- current_row_s + 1 # Move to next row
      }else if (isTRUE(all.equal(x_s, x_i))){
        int_1 <- int_1 + (value_s * value_i * (x_s - prev_x)) # Sum integral
        prev_x <- x_s # Update previous x
        current_row_s <- current_row_s + 1 # Move to next row
        current_row_i <- current_row_i + 1 # Move to next row
      }else{ # x_i > x_s
        if (x_i < l_i){
          int_1 <- int_1 + (value_s * value_i * (x_i - prev_x)) # Sum integral
          prev_x <- x_i # Update previous x
        }else{ # x_i > l_i
          int_1 <- int_1 + (value_s * value_i * (l_i - prev_x)) # Sum integral
          prev_x <- x_i # Update previous x
        }
        current_row_i <- current_row_i + 1 # Move to next row
      }
    }
    return(int_1) 
  }
  
  if (int_2 != 0){ # Ancestor integral is not zero, calculate other integral/s
    Sum_i_a <- calculate_DJ_i_a(tree, node, abundances, curr_ancestor, h, l_i, 
                                index_letter, q, individual)
    
    if (individual == TRUE){
      int <- integral(Sum_i_a)*int_2
    }else{
      int <- sapply(Sum_i_a, integral)*int_2
    }
  }else if (int_2 == 0){ # Ancestor integral is zero, return zero
    if (individual == TRUE){
      int <- 0
    }else{
      int <- c(0,0,0)
    }
  }
  
  return(int)
}