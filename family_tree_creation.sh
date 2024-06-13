#! /bin/bash

# Move to working directory
cd /lustre1/project/stg_00026/enora/UNESCO/tree_creation

# Load conda environment
source ~/.bashrc
conda activate tree-creation

family=$1
tax_table=$2 # ../pacman-pipeline-results/eDNAexpeditions_batch1_samples/runs/${site}_12SMifish/05-dwca/Full_tax_table_with_lsids.tsv
amplicon=$3 # 12SMifish
primer_forward=$4 # GTCGGTAAAACTCGTGCCAGC
primer_reverse=$5 # CAAACTGGGATTAGATACCCCACTATG

##### Add step to check if tax_table given (if not stop the script execution) #####


##### Need to add check step for amplicon to see if present in database and if not check if primer sequences provided, if not stop script #####


# Checking primer set for cutadapt
if [ "$primer_forward" == "" ] && [ "$primer_reverse" == "" ]
then
        echo -e "You didn't provide primer sequences so we will use the ones present in our primer file."
        primer_forward=`grep "$amplicon" primer_file.tsv | cut -f2`
        primer_reverse=`grep "$amplicon" primer_file.tsv | cut -f3`
elif ["$primer_forward" == "" ] || [ "$primer_reverse" == "" ]
then
        echo -e "You only provided one of the primer sequence so it will not be possible to correctly trim the NCBI reference sequences. Hence$
        primer_forward=`grep "$amplicon" primer_file.tsv | cut -f2`
        primer_reverse=`grep "$amplicon" primer_file.tsv | cut -f3`
else
        echo -e "The following forward ${primer_forward} and reverse ${primer_reverse} sequences will be used to trim the NCBI reference seque$
fi

echo -e "Starting tree creation for the $family family identified as the $family_pacman in the pacman pipeline assignation\n"

mkdir $family

# Check if genus_final.fa exist for this family
FILE=${family}_${amplicon}_genus_final.fa
if [ -f "$FILE" ]
then
        echo -e "$FILE exists, skipping creation steps of this file.\n"
        cd $family
        cp ../${family}_${amplicon}_genus_final.fa ${family}_${amplicon}_genus_final.fa
else
        echo -e "$FILE does not exist, starting creation of the file.\n----------\n"
        # Retrieve genus names of the chosen family
        ## Get taxon ID of the family
        echo -e "Retrieving accession numbers of each sequences of interest in each genus of the $family family...\n----------\n"
        id=`echo $family | taxonkit name2taxid --data-dir . | cut -f2`
        if [ "$id" == "" ]
        then
                echo "This family does not exist in the NCBI taxonomy database, checking the corresponding name in the NCBI database..."
                id=`awk -F'\t' -v fam="$family" ' $6 == fam { print $7 }' *_12SMifish/identified_ASVs.txt | sort | uniq| head -n1 | taxonkit n$
                ranklist=`echo $id | taxonkit lineage -R --data-dir . | cut -f2-3 | sed 's/\t/\n/g' | tac`
                column=`echo $ranklist | head -n1 | tr ';' '\n' | grep -n "^family$" | cut -d':' -f1`
                family_NCBI=`echo $ranklist | sed 's/no rank/norank/g' | sed 's/cellular organisms/cellularOrganisms/g' | sed 's/ /\n/g' | cut$
                echo "The corresponding name for the $family family in the NCBI database is $family_NCBI"
                id=`echo $family_NCBI | taxonkit name2taxid --data-dir . | cut -f2`
        fi

        ## Extract all genus names of this family
        taxonkit list --ids $id -n -r --data-dir . | grep "\[genus\]" | cut -d']' -f2 | sed 's/ //g' | sort | uniq > ${family}_${amplicon}_gen$

        ## Move to family directory
        cd $family

        # Retrieve all NCBI sequences containing "mitochondrial" or "mitochondrion" and "12S" or "complete genome"
        ## Retrieve all the accession numbers
##### Adapt this part to be able to use other amplicons #####
	while read line
        do
                genus=`echo $line`
                grep "$genus" ../all_nt_db_acc.txt | grep -E 'mitochondrion|mitochondrial' > ${genus}_mito.txt # retrieve all sequences contai$
                grep -E 'complete genome' ${genus}_mito.txt | sed 's/ >/\n>/g' | cut -d'.' -f1 | cut -d'>' -f2 >> ${family}_acc-num.txt # from$
                grep -E '12S' ${genus}_mito.txt | sed 's/ >/\n>/g' | cut -d'.' -f1 | cut -d'>' -f2 >> ${family}_acc-num.txt #from these sequen$
                rm ${genus}_mito.txt
        done < ../${family}_${amplicon}_genus.txt

        ## Download the sequences retrieved in the accession file
        echo -e "Downloading all the sequences of the retrieved accession numbers...\n----------\n"
        sort ${family}_acc-num.txt | uniq > ${family}_acc-num_nodup.txt
        mkdir ncbi_acc_files
        while read line
        do
                acc=`echo $line`
                ncbi-acc-download --format fasta --out ncbi_acc_files/${acc}.fasta $acc
        done < ${family}_acc-num_nodup.txt
        cat ncbi_acc_files/*.fasta > ${family}_${amplicon}_genus.fasta

        ## Cleaning temporary files
	rm -r ncbi_acc_files/
        rm ${family}_acc-num.txt

        ## Put all sequences on one line each
        awk '/^>/ { if(NR>1) print "";  printf("%s\n",$0); next; } { printf("%s",$0);}  END {printf("\n");}' < ${family}_${amplicon}_genus.fas$

        ## Clean
        rm ${family}_${amplicon}_genus.fasta
        mv ${family}_${amplicon}_genus_1l.fa ${family}_${amplicon}_genus.fa

        # Retrieve sequence of interest
        ## Activate cutadapt environment
        echo -e "Cutting the retrieved sequences with Cutadapt to reduce the sequences to only the fragment of interest (here Mifish locus)...$
        mamba activate cutadapt

        ## Detect forward primer
        cutadapt -g $primer_forward -o ${family}_${amplicon}_genus_short1.fa -e 4 -j 1 ${family}_${amplicon}_genus.fa > ${family}_${amplicon}_$

        ## Detect reverse primer
        cutadapt -a $primer_reverse -o ${family}_${amplicon}_genus_short2.fa -e 4 -j 1 ${family}_${amplicon}_genus_short1.fa > ${family}_${amp$

        ## Clean
        rm ${family}_${amplicon}_genus_short1.fa ${family}_${amplicon}_genus.fa
        mv ${family}_${amplicon}_genus_short2.fa ${family}_${amplicon}_genus_cut.fa

        # Refine sequence file to make it clean
        ## Activate unesco-trees environment
        echo -e "Refining the sequences and changing sequence's names...\n----------\n"
        mamba activate unesco-trees

        ## Remove duplicated sequences
        seqkit rmdup -s < ${family}_${amplicon}_genus_cut.fa > ${family}_${amplicon}_genus_cut_nodup.fa

        ## Put all sequences on one line each
        awk '/^>/ { if(NR>1) print "";  printf("%s\n",$0); next; } { printf("%s",$0);}  END {printf("\n");}' < ${family}_${amplicon}_genus_cut$

        ## Remove special characters from names and empty sequences
        sed -i 's/://g' ${family}_${amplicon}_genus_cut_nodup_1l.fa | sed 's/(//g' | sed 's/)//g' | sed 's/-//g'
	sed -i '/^[[:space:]]*$/d' ${family}_${amplicon}_genus_cut_nodup_1l.fa

        ## Remove “PREDICTED” words to be able to retrieve full species names
        sed -i 's/PREDICTED //g' ${family}_${amplicon}_genus_cut_nodup_1l.fa

        rm ${family}_${amplicon}_genus_cut.fa ${family}_${amplicon}_genus_cut_nodup.fa
        mv ${family}_${amplicon}_genus_cut_nodup_1l.fa ${family}_${amplicon}_genus_cut.fa

        # Rename the sequences
        while read line
        do
                grep_res=`echo $line | grep -v "^>"`
                if [ -z $grep_res ]
                then
                        name=`echo $line | cut -d' ' -f2-3`
                        grep_name=`grep " ${name} " ${family}_${amplicon}_genus_cut.fa | wc -l`
                        if [ $grep_name -eq 1 ]
                        then
                                sed -i "s|$line|>${name}|g" ${family}_${amplicon}_genus_cut.fa
			else
                                last=`grep "${name}_" ${family}_${amplicon}_genus_cut.fa | tail -n1 | cut -d'_' -f2`
                                num=$(( ${last} + 1 ))
                                new_name=${name}_${num}
                                sed -i "s|$line|>${new_name}|g" ${family}_${amplicon}_genus_cut.fa
                        fi
                fi
        done < ${family}_${amplicon}_genus_cut.fa

        ## Remove empty fasta records
        awk 'BEGIN {RS = ">" ; FS = "\n" ; ORS = ""} $2 {print ">"$0}' ${family}_${amplicon}_genus_cut.fa > ../${family}_${amplicon}_genus_fin$

        cp ../${family}_${amplicon}_genus_final.fa ${family}_${amplicon}_genus_final.fa
fi
echo -e "\nFile containing NCBI sequences for the family $family is created.\n\n----------\n----------\n"

# Add all ASV assigned to this family
## Get all the ASV identified in this family for each site
echo -e "Starting to retrieve ASVs assigned till at least the $family family...\n----------\n"
while read line
do
        site=`echo $line`
        loc_asv="../${site}_12SMifish/${family}_asv.txt"
        if [ -f "$loc_asv" ]
        then
                echo -e "\n"
        else
                echo -e "File containing all identified ASVs till at least family level for $site is missing, starting to create it...\n\n----$
                awk -F'\t' ' $13 == "genus" || $13 == "species" || $13 == "family" { print $0 }' ../pacman-pipeline-results/eDNAexpeditions_ba$
        fi
        awk -F'\t' -v f=$family ' $6 == f { print $1 }' ../${site}_12SMifish/identified_ASVs.txt > ../${site}_12SMifish/${family}_asv.txt
        while read line
        do
                asv=`echo $line`
                species=`grep -P "^$asv\t" ../pacman-pipeline-results/eDNAexpeditions_batch1_samples/runs/${site}_12SMifish/05-dwca/Full_tax_t$
                echo ">"$asv"_"$species"_${site}" >> ../${site}_12SMifish/${family}_asv.fa
                grep "${asv}_" ../pacman-pipeline-results/eDNAexpeditions_batch1_samples/runs/${site}_12SMifish/05-dwca/DNA_extension_table.ts$
        done < ../${site}_12SMifish/${family}_asv.txt
done < ../all_sites.txt

## Put all the ASVs in the same file
cat ../*_12SMifish/${family}_asv.fa > all_${family}_asv.fa

## Add the ASVs to the final fasta file for the alignment
cat all_${family}_asv.fa >> ${family}_genus_final.fa

# Add the outgroup
echo -e "Adding the outgroup...\n----------\n"
cat ../U11880_cut2.fa >> ${family}_genus_final.fa

# Remove the sequences longer than 200 nucleotides
echo -e "Refining the final fasta file containing ASVs and sequences from NCBI database...\n"
seqkit seq -M 280 ${family}_genus_final.fa | awk '/^>/ { if(NR>1) print "";  printf("%s\n",$0); next; } { printf("%s",$0);}  END {printf("\n")$

# Remove duplicated sequences
seqkit rmdup -s < ${family}_genus_final_short.fa > ${family}_genus_final_clean.fa

# Align the sequences with Mafft
echo -e "\nAligning with MAFFT...\n"
mafft --thread -1 --preservecase --adjustdirection ${family}_genus_final_clean.fa > ${family}.aln.fa 2> ${family}_mafft.log

## Replace spaces by underscores
sed -i 's/ /_/g' ${family}.aln.fa

# Curate the alignment with Gblocks
echo -e "Curating the alignment with Gblocks...\n"
Gblocks ${family}.aln.fa -t=d -b4=5 -b5=h -e=.gb > ${family}_gblock.log 2>&1

# Reconstruct phylogenetic tree with FastTree
echo -e "Constructing the phylogenetic tree with FastTree...\n"
FastTreeMP -gtr -nt ${family}.aln.fa.gb > ${family}.nwk 2>${family}_FastTree.log

## Remove "_R_"
sed -i 's/_R_//g' ${family}.nwk

# Create information file to differentiate ASVs and sequences from NCBI database
## Retrieve all names in the tree
echo -e "Creating information file to color differently ASVs and sequences from NCBI database...\n----------\n"
grep "^>" ${family}.aln.fa.gb | sed 's/>//g' | sed 's/_R_//g' > all_labels.txt

## For each name determine if it's an ASV or not and put the information in a file
while read line
do
        name=`echo $line`
        grep_res=`echo $name | grep "asv."`
        if [ -z $grep_res ]
        then
                echo -e "${name}\tdb_name" >> ${family}_info.tsv
        else
                echo -e "${name}\tasv" >> ${family}_info.tsv
        fi
done < all_labels.txt

## Add column's names
sed -i '1itaxa\tasv' ${family}_info.tsv

# Plot the tree
echo -e "Plotting the tree with ggtree package in R...\n"
## Activate environment
mamba activate R-tree-creation

## Run R script
path=`pwd`; Rscript ../plot_tree.r ${path} ${family}.nwk ${family}_info.tsv ${family} ${family}.pdf

## Clean temporary files
rm Rplots.pdf

# Check if final file produced
ENDFILE=${family}.pdf
if [ -f "$ENDFILE" ]
then
        echo -e "\nTree creation script worked, $ENDFILE is available in ${family}/ folder!"
else
        echo -e "\nA problem occurred! $ENDFILE was not created. Please look at the log file to understand what happend."
fi

