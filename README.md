# Phylo-tree_eDNA
Scripts to create phylogenetic trees based on taxonomic assignments of Environmental DNA (eDNA) data in order to **detect possible miss-assigned** ASVs or **refine the taxonomic assignment** of some ASVs.

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
git clone https://github.com/Enorya/Phylo-tree_eDNA.git
cd Phylo-tree_eDNA/
conda env create -f tree_creation.yml
conda env create -f R_tree_creation.yml
```

## Usage
In order to use the following scripts you can use similar commands as the following ones:
```
location_tree_creation.sh yourSamplingLocationName
family_tree_creation.sh familyTaxonomicName
```

# To do

- Give example file for taxonomy table and sequence table
- Add information in the README to tell to users which amplicon are available by default
- Check if possible to adapt the scripts to have other kind of inputs as the ones used to develop the scripts
- Integrate a check to avoid downloading NCBI data if they are already available (downloaded less than a month ago)
- Idea for a logo: water drop with some DNA strand in it and coming out from it and transforming to a phylogenetic tree
