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
In order to use the TaxonKit dependency you need to download taxonomy information from the NCBI database:
```
wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
tar -xf taxdump.tar.gz
```
In order to download the reference sequences from NCBI you need to prepare a file containing all the accesion numbers and description of the sequence. To do so you can use:
```
mkdir nt_db
cd nt_db
wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nt.gz
gunzip nt.gz
grep "^>" nt.fsa > all_nt_db_acc.txt
```
> [!NOTE]
> A similar file is already provided in this repository but it might be good to redo this step if the file is too old.
>
> (last update: 19 June 2024)

## Usage
In order to use the following scripts you can use similar commands as the following ones:
```
location_tree_creation.sh yourSamplingLocationName
family_tree_creation.sh familyTaxonomicName 'taxonomy_table.tsv' 'sequence_table.tsv' amplicon_name forward_primer reverse_primer
```

## Example use for Family tree
### If you want to use an amplicon from the default list called `list_primers.tsv`:
```
family_tree_creation.sh Dasyatidae '../test/tax_table.tsv' '../test/seq_table.tsv' 12SMifish
```
The available amplicons are:

- Mifish-U from Miya et al. 2015 (171bp) = 12SMifish
- MiMammal-U from Ushio et al. 2017 (171bp) = 12SMimammal
- Teleo from Valentini et al. 2016 (63bp) = 12STeleo
- Leray-COI from Leray et al. 2013 (313bp) = COI
- Vert-16S from Vences et al. 2016 (250bp) = 16SVert

### If you want to use your own amplicon:
```
family_tree_creation.sh Dasyatidae '../test/tax_table.tsv' '../test/seq_table.tsv' 16S ATTCGCCAAGTCAAG GGGTCTCCAAAAGTCGT
```

> [!WARNING]
> The single quotation marks around the table's names are mandatory for proper results, don't forget to put them

# To do

- Give example file for sequence table
- Need to implement cutadapt in conda environment
- Check if possible to adapt the scripts to have other kinds of inputs as the ones used to develop the scripts -> in progress
- Integrate a check to avoid downloading NCBI data if list of accessions less than 1 year old, otherwise download nucleotide db and extract all accession names (possible to implement parameter in the futur to choose if you want to download or not)
- Idea for a logo: water drop with some DNA strand in it and coming out from it and transforming to a phylogenetic tree
