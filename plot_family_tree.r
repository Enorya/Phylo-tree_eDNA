#! /data/leuven/340/vsc34088/miniconda3/envs/R-unesco-trees/bin/Rscript

# Set the arguments
args = commandArgs(trailingOnly=TRUE)

# Set working directory
setwd(args[1])

# Load necessary libraries
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggtree))
suppressPackageStartupMessages(library(svglite))
suppressPackageStartupMessages(library(ape))

# Read tree and root it
tree <- read.tree(args[2])
tree2 <- root(tree, outgroup = "Petromyzon_marinus", edgelabel = TRUE)
info <- read.table(args[3], header=TRUE)

# Plot the tree
tree2$node.label = as.numeric(tree2$node.label)
tree2$node.label = round(tree2$node.label * 100)
num <- (length(tree2[["tip.label"]])/10)*2
family_name = paste("rooted tree for",args[4],"family")
plot <- ggtree(tree2) + theme_tree2() + theme(plot.title = element_text(hjust = 0.5, size = 35)) + xlim(0,3) + ggtitle(paste0(family_name)) + geom_nodelab(size = 4, color="darkgreen", vjust=-0.30, hjust = 1.3)
plot2 <- plot %<+% info + geom_tiplab(aes(color=asv), size = 6) + theme(legend.title = element_text(size=40), plot.title = element_text(size = 40, face = "bold"), legend.text = element_text(size=30), legend.key.size = unit(2, 'cm')) + scale_color_manual(values = c("#9c763e","#9e57ad"))
plot2
ggsave(args[5], width = 30, height = num, limitsize = FALSE)
