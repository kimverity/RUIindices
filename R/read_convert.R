#' Reads in tree and converts to 'phylo' object
#' 
#' Reads in tree from a file or tree can be inputted directly (e.g. as string),
#' and converts to 'phylo' object. If no branch lengths are given, they are 
#' assigned to be one.
#'
#' @param file A string or phylo object; String of file name containing tree or 
#' the tree itself either as a string or a phylo object
#' @returns A 'phylo' object
#' @examples
#' tree <- '(a:1,(b:0.5,c:0.5)e:0.5)d:0;'
#' read_convert(tree)
#' 
#' @export
read_convert <- function(file){
  
  # Checks the format of file/tree and comverts to phylo objext
  suppressWarnings({
    tree <- try(ape::read.tree(file), silent = TRUE) # Newick format
    if ("try-error" %in% class(tree)){
      tree <- try(ape::read.nexus(file), silent = TRUE) # Nexus format
      if ("try-error" %in% class(tree)){
        tree <- try(ape::read.tree(text=file), silent = TRUE) # String in newick format
        if ("try-error" %in% class(tree)){
          if (inherits(file,"phylo")){ # Already phylo object
            tree <- file
          }else if(!(inherits(file,"phylo"))){ # file is none of the above
            return(print("Tree must be in Newick or NEXUS format, or be a phylo object.")) 
          }
        }
      }
    }
  })  
  
  # If no branch lengths, assign all to be one
  if (length(tree$edge.length) == 0){
    tree$edge.length <- rep(1, times = length(tree$edge[,1]))
  }
  
  return(tree)
}