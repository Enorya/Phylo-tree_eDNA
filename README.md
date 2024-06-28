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

<p align="center">
<img src="https://github.com/Enorya/Phylo-tree_eDNA/figures/pipeline_figure_v4.png" alt="family_tree_pipeline" width="500"/>
</p>

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
> (last update: 19 June 2024 with 102960590 sequences)

## Preparation of the files

### Taxonomy table
Your file should look like the `tax_table.tsv` available in the test/ directory
| ASV_name | superkingdom | kingdom | phylum | class | order | family | genus | species | scientificName | taxonRank |
| :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| asv.1 | | Animalia | Chordata | Teleostei | Clupeiformes | Spratelloididae | Jenkinsia | Jenkinsia lamprotaenia | Jenkinsia lamprotaenia | species |
| asv.2 | Eukaryota | | | | | | | | Incertae sedis | kingdom |
| asv.3 | | Animalia | Chordata | Elasmobranchii | Myliobatiformes | | | | Myliobatiformes | order |

### Sequence table
Your file should look like the `seq_table.tsv` available in the test/ directory
| target_gene | pcr_primer_forward | pcr_primer_reverse | DNA_sequence | occurrenceID |
| :---: | :---: | :---: | :---: | :---: |
| 12S | GTCGGTAAAACTCGTGCCAGC | CATAGTGGGGTATCTAATCCCAGTTTG | GTTGGTAAATCTCGTGCCAGCCACCGCGGTCACACGATTAACCCAAGTCAATAGAAGCCGGCGTAAAGAGTGTTTTAGATCACCCCCTCCCCAATAAAGCTAAAACTCACCTGAGTTGTAAAAAACTCCAGTTGACACAAAATAGACTACGAAAGTGGCTTTAACATATCTGAACACACAATAGCTAAGACC | asv.1_S003 |
| 12S | GTCGGTAAAACTCGTGCCAGC | CATAGTGGGGTATCTAATCCCAGTTTG | CACCGCGGTCACACGATTAACCCAAGTCAATAGAAGCCGGCGTAAAGAGTGTTTTAGATCACCCCCTCCCCAATAAAGCTAAAACTCACCTGAGTTGTAAAAAACTCCAGTTGACACAAAATAGACTACGAAAGTGGCTTTAACATATCTGAACACACAATAGCTAAGACC | asv.1_S026 |
| 12S | GTCGGTAAAACTCGTGCCAGC | CATAGTGGGGTATCTAATCCCAGTTTG | GGGTTGGTAAATTTCGTGCCAGCCACCGCGGTCACACGATTAACCCAAGTCAATAGAAGCCGGCGTAAAGAGTGTTTTAGATCACCCCCTCCCCAATAAAGCTAAAACTCACCTGAGTTGTAAAAAACTCCAGTTGACACAAAATAGACTACGAAAGTGGCTTTAACATATCTGAACACACAATAGCTAAGACC | asv.3_S017 |
> [!NOTE]
> Columns with PCR primer sequences are not mandatory in this table.

## Usage of Family tree script
Here is a description of all the parameters of this tool:
```
This script create a phylogenetic tree by taking as input a taxonomic family name, the name of the used amplicon and a taxonomy and sequence tables.

Syntax: script.sh [-t|s|a|f|r|g|l|h]

Options:
-n        name of family taxonomic group (mandatory)
-t        taxonomy table (mandatory)
-s        sequence table (mandatory)
-a        amplicon name (mandatory)
-f        forward primer sequence
-r        reverse primer sequence
-g        gene name (if different than amplicon name)
-l        maximum length of the amplicon (advice: length of the amplicon for outgroup + 100bp)
-h        display this help message.
```

For example, if you want to use an amplicon from the default list called `list_primers.tsv`:
```
family_tree_creation.sh -n Dasyatidae -t '../test/tax_table.tsv' -s '../test/seq_table.tsv' -a 12SMifish -o '../outgroup_sequences/petromyzon_marinus_12SMifish.fa'
```
The available amplicons are:

- Mifish-U from Miya et al. 2015 (171bp) = 12SMifish
- MiMammal-U from Ushio et al. 2017 (171bp) = 12SMimammal
- Teleo from Valentini et al. 2016 (63bp) = Teleo
- Leray-COI from Leray et al. 2013 (313bp) = COI
- Vert-16S from Vences et al. 2016 (250bp) = 16S


If you want to use your own amplicon:
```
family_tree_creation.sh -n Dasyatidae -t '../test/tax_table.tsv' -s '../test/seq_table.tsv' -a 16SVert -f AGTCCCGAAATATAAT -r GCTGTTGTGCCCGAAG -g 16S -l 350 -o '../outgroup_sequences/petromyzon_marinus_12SMifish.fa'
```

If you want to use a different maximum length of the amplicon than the one listed in `list_primers.tsv`:
```
family_tree_creation.sh -n Dasyatidae -t '../test/tax_table.tsv' -s '../test/seq_table.tsv' -a CO1 -l 380 -o outgroup_sequences/petromyzon_marinus_16S.fa
```
> [!TIP]
> Each parameter listed in the help message can be modified, even when using default amplicons, but if you want to change the primers' sequences, you need to change both, otherwise the default sequences will be used.

> [!WARNING]
> The single quotation marks around the table's names are mandatory for proper results, don't forget to put them!

## Usage of Location tree script
> [!CAUTION]
> This script is not finalised for the moment, please take this into consideration if you wish to use it.

# To do

- Idea for a logo: water drop with some DNA strand in it and coming out from it and transforming to a phylogenetic tree
