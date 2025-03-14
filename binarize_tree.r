#!/usr/bin/env Rscript

# Load required package
library(ape)

# Read the tree from the Newick file
tree <- read.tree("y_tree.nwk")

# Root the tree using mouse_chrY as the outgroup and resolve the root
#tree <- root(tree, outgroup = "mouse_chrY", resolve.root = TRUE)

# Convert the tree to strictly binary if needed
binary_tree <- multi2di(tree)

# Save the fixed, binary tree
write.tree(binary_tree, file = "binarized_primate_tree.nwk")

cat("Binarized tree saved as binarized_primate_tree.nwk\n")



