#' Reads in data frame tree (see RUtreebalance repo) and converts to phylo
#' object
#'
#' The df-tree is a data frame with columns: Identity, Parent, Population, and
#' optionally a column called Generation (from demon simulations) which is used
#' to set the branch lengths. If Generation is not included, the branch lengths
#' are all set to 1.
#'
#' @param df_tree A data frame with columns: Identity, Parent, Population, and
#' optionally Generation.
#' @returns A phylo object and a data frame with node abundances.
#' @examples
#' df_tree <- data.frame(Identity=c(0,1,2,3),
#'                    Parent=c(0,0,1,1),
#'                    Population=c(1,1,2,2),
#'                    Generation=c(0,1,3,5))
#' tree <- df_tree_to_phylo(df_tree)
#'
#' @export
df_tree_to_phylo <- function(df_tree) {
    # Check that the df_tree has the correct columns
    if (!all(c("Identity", "Parent", "Population") %in% colnames(df_tree))) {
        stop("df_tree must have columns Identity, Parent & Population")
    }
    na_dataframe <- data.frame("Names" = as.character(df_tree$Identity),
                               "Values" = df_tree$Population)
    phylo_tree <- list()
    phylo_tree$edge <- as.matrix(df_tree[, c("Parent", "Identity")])
    phylo_tree$tip.label <- as.character(setdiff(df_tree$Identity,
                                                 df_tree$Parent))
    phylo_tree$Nnode <- length(unique(df$Parent)) - 1
    edge_lengths <- numeric(nrow(df_tree))
    for (i in 1:nrow(df_tree)) {
        if (df_tree$Identity[i] == 0) {
            edge_lengths[i] <- 0 # root node
        }
        else {
            parent_gen <- df_tree$Generation[df_tree$Identity == df$Parent[i]]
            edge_lengths[i] <- df_tree$Generation[i] - parent_gen
        }
    }
    phylo_tree$edge.length <- edge_lengths
    internal_node_labels <- df$internal_label[df$Identity %in% phylo_tree$edge[,1]]
    phylo_tree$node.label <- as.character(internal_node_labels)
    phylo_tree$node.label <- paste("Node", seq_len(phylo_tree$Nnode))
    return(list(phylo_tree, na_dataframe))
}
