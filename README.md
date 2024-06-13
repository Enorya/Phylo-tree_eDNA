# Phylo-tree_eDNA
Script to create phylogenetic trees based on taxonomic assignments of Environmental DNA (eDNA) data

## Description
This repository contains 2 different scripts to plot phylogenetic trees based on taxonomic assignments of Environmental DNA (eDNA) data.

There is 2 different methods to generate a phylogenetic tree.

### Location based method
This method gives the possibility to plot all your assigned ASVs of one specific sampling site along with some reference sequences downloaded from the NCBI nucleotide database.

If you made multiple sampling at the same location, with this method, it's possible to color the ASVs depending on how many samplings the species is retrieved. 


### Family based method
This method gives the possibility to plot all the assigned ASVs of one specific taxonomic family accross all sampling sites along with the available reference sequences of the NCBI nucleotide database from this specific family.

With this method, the ASVs are colored depending if they are ASVs or coming from the NCBI database. Also, the names of the ASVs sequences in the tree contain the location where they come from.

## Installation
To be able to use the 2 scripts you need to create 2 new conda environments you need to first clone this repository and then use:
```
conda env create -f unesco-trees.yml
conda env create -f R-unesco-trees.yml
```
