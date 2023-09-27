# RUIindices
 
Robust, universal and interpretable tree indices.

Preprint: (https://doi.org/10.1101/2023.07.17.549219)

The file (or tree) must be in one of the following forms: either a file containing a single tree in Newick or nexus format, or the tree, as a string in Newick format or as a `phylo` object.

Node abundances must be inputted as a dataframe with two columns. The first column must contain the node labels as given by the inputted tree, the second column must contain the corresponding node abundances. If node abundances are omitted, internal nodes will be assigned size zero and the leaves will be assigned to be equally abundant.

If the given tree does not have branch lengths then they will be assigned to be equal with a length of 1.

Example:
```
tree <- ape::read.tree(text='(a:1,(b:0.5,c:0.5)e:0.5)d:0;')
na_dataframe <- data.frame("Names" = c("e", "d", "b", "c", "a"), "Values" = c(0, 0, 1/3, 1/3, 1/3))
all_indices(tree, na_dataframe)
```
