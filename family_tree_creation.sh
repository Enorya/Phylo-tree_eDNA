#! /bin/bash

##############
#### Help ####
##############
Help()
{
	# Display Help
	echo "This script create a phylogenetic tree by taking as input a taxonomic family name, the name of the used amplicon and a taxonomy and sequence tables."
	echo
	echo "Syntax: script.sh [-t|s|a|f|r|g|l|o|h]"
	echo
	echo "Options:"
	echo "-n	name of family taxonomic group (mandatory)"
	echo "-t	taxonomy table (mandatory)"
	echo "-s	sequence table (mandatory)"
	echo "-a	amplicon name (mandatory)"
	echo "-f	forward primer sequence"
	echo "-r	reverse primer sequence"
	echo "-g	gene name (if different than amplicon name)"
	echo "-l	maximum length of the amplicon (advice: length of the amplicon for outgroup + 100bp)"
	echo "-o	outgroup fasta file already formatted"
	echo "-h	display this help message."
}

######################
#### Main program ####
######################

###############################
#### Process input options ####
###############################

while getopts hn:t:s:a:f:r:g:l:o: option
do
	case $option in
		(h) # display Help
			Help
			exit;;
		(n) # enter family taxonomic name
			family=$OPTARG;;
		(t) # enter taxonomy table name
			tax_table=$OPTARG;;
		(s) # enter sequence table name
			seq_table=$OPTARG;;
		(a) # enter amplicon name
			amplicon=$OPTARG;;
		(f) # forward primer sequence
			primer_forward=$OPTARG;;
		(r) # reverse primer sequence
			primer_reverse=$OPTARG;;
		(g) # gene name
			gene_name=$OPTARG;;
		(l) # maximum length of the amplicon
			length_amp=$OPTARG;;
		(o) # outgroup
			outgroup=$OPTARG;;
		(\?) # Invalid option
			echo -e "\nError: Invalid option\n\n"
			Help
			exit;;
	esac
done

# Load conda environment
source ~/.bashrc
conda activate tree-creation

echo -e "Checking the parameters given by the user to be sure that all needed informations are provided.\n"

# Checking if family name available

if [ "$family" == "" ]
then
        echo -e "Family name not detected. Please provide a family taxonomy group name."
        exit 1
else
        echo -e "Family name detected: $family"
fi


# Checking if taxonomy table available

if [ "$tax_table" == "" ]
then
	echo -e "Taxonomy table not detected. Please provide a file containing taxonomy information of the ASVs."
	exit 1
else
	echo -e "Taxonomy table detected: $tax_table"
fi

# Checking if sequence table available

if [ "$seq_table" == "" ]
then
        echo -e "Sequence table not detected. Please provide a file containing the sequence of each ASV."
	exit 1
else
        echo -e "Sequence table detected: $seq_table"
fi

# Checking if amplicon name provided and if primer sequences available if primer sequence not in list

if [ "$amplicon" == "" ]
then
	echo -e "No amplicon name detected. Please provide the name of the amplicon you used to start the analysis."
	exit 1
else
	if [ "$primer_forward" == "" ] && [ "$primer_reverse" == "" ]
	then
		grep_amp=`grep -o -P "^$amplicon\t" list_primers.tsv`
		if [ -z $grep_amp ]
		then
			echo -e "You didn't provide the forward and reverse primers of your amplicon and the amplicon name you provided is not in the default list. Please provide the primer sequences or use an amplicon name from the list."
			exit 1
		else
			echo -e "You didn't provide the forward and reverse primers of your amplicon but the amplicon name you provided is in the default list. Default primers from the list will be used for the trimming of the reference sequences.\n"
			primer_forward=`grep "$amplicon" list_primers.tsv | cut -f2`
			primer_reverse=`grep "$amplicon" list_primers.tsv | cut -f3`
			if [ "$gene_name" == "" ]
			then
				echo -e "You didn't provide the real gene name but the amplicon name you provided is in the default list. Default gene name from the list will be used to retrieve the reference sequences.\n"
				gene_name=`grep "$amplicon" list_primers.tsv | cut -f4`
			else
				echo -e "Gene name detected. The following gene name will be used to retrieve the reference sequences: $gene_name"
			fi
			if [ "$length_amp" == "" ]
                        then
                                echo -e "You didn't provide the maximum length of the amplicon but the amplicon name you provided is in the default list. Default maximum length from the list will be used to remove too long sequences before alignment.\n"
                                length_amp=`grep "$amplicon" list_primers.tsv | cut -f5`
                        else
                                echo -e "Maximum length of the amplicon detected. The following length will be used to remove too long sequences before the alignment: $length_amp"
                        fi
			if [ "$outgroup" == "" ]
                        then
                                echo -e "You didn't provide the fasta file of the outgroup but the amplicon name you provided is in the default list. Default outgroup from the list will be used.\n"
				outgroup=./outgroup_sequences/petromyzon_marinus_${amplicon}.fa
                        else
                                echo -e "Fasta file of the outgroup detected. The following outgroup will be used: $outgroup"
                        fi
		fi
	elif ["$primer_forward" == "" ] || [ "$primer_reverse" == "" ]
	then
		grep_amp=`grep -o -P "^$amplicon\t" list_primers.tsv`
                if [ -z $grep_amp ]
                then
                        echo -e "You only provided one primer for your amplicon and the amplicon name you provided is not in the default list. Please provide both forward and reverse primer sequences or use an amplicon name from the list."
                        exit 1
                else
                        echo -e "You only provided one primer for your amplicon but the amplicon name you provided is in the default list. Default primers from the list will be used for the trimming of the reference sequences.\n"
                        primer_forward=`grep "$amplicon" list_primers.tsv | cut -f2`
                        primer_reverse=`grep "$amplicon" list_primers.tsv | cut -f3`
			if [ "$gene_name" == "" ]
                        then
                                echo -e "You didn't provide the real gene name but the amplicon name you provided is in the default list. Default gene name from the list will be used to retrieve the reference sequences.\n"
                                gene_name=`grep "$amplicon" list_primers.tsv | cut -f4`
                        else
                                echo -e "Gene name detected. The following gene name will be used to retrieve the reference sequences: $gene_name"
			fi
			if [ "$length_amp" == "" ]
			then
				echo -e "You didn't provide the maximum length of the amplicon but the amplicon name you provided is in the default list. Default maximum length from the list will be used to remove too long sequences before alignment.\n"
				length_amp=`grep "$amplicon" list_primers.tsv | cut -f5`
			else
				echo -e "Maximum length of the amplicon detected. The following length will be used to remove too long sequences before the alignment: $length_amp"
			fi
			if [ "$outgroup" == "" ]
                        then
                                echo -e "You didn't provide the fasta file of the outgroup but the amplicon name you provided is in the default list. Default outgroup from the list will be used.\n"
                                outgroup=./outgroup_sequences/petromyzon_marinus_${amplicon}.fa
                        else
                                echo -e "Fasta file of the outgroup detected. The following outgroup will be used: $outgroup"
                        fi
                fi
	else
		echo -e "The following forward ${primer_forward} and reverse ${primer_reverse} sequences will be used to trim the NCBI reference sequences.\n"
		if [ "$gene_name" == "" ]
		then
			if [ -z $grep_amp ]
			then
				echo -e "You didn't provide the real gene name and the amplicon name you provided is not in the default list. Please provide a gene name to be able to retrieve the reference sequences.\n"
				exit 1
			else
                                echo -e "You didn't provide the real gene name but the amplicon name you provided is in the default list. Default gene name from the list will be used to retrieve the reference sequences.\n"
                                gene_name=`grep "$amplicon" list_primers.tsv | cut -f4`
			fi
		else
			echo -e "Gene name detected. The following gene name will be used to retrieve the reference sequences: $gene_name\n"
		fi
		if [ "$length_amp" == "" ]
		then
			if [ -z $grep_amp ]
			then
				echo -e "You didn't provide the maximum length of the amplicon and the amplicon name you provided is not in the default list. Please provide a maximum length of the amplicon to be able to remove too long sequences before alignment.\n"
				exit 1
			else
				echo -e "You didn't provide the maximum length of the amplicon but the amplicon name you provided is in the default list. Default maximum length from the list will be used to remove too long sequences before alignment.\n"
				length_amp=`grep "$amplicon" list_primers.tsv | cut -f5`
			fi
		else
			echo -e "Maximum length of the amplicon detected. The following length will be used to remove too long sequences before the alignment: $length_amp\n"
		fi
		if [ "$outgroup" == "" ]
                then
                        if [ -z $grep_amp ]
                        then
                                echo -e "You didn't provide the fasta file of the outgroup and the amplicon name you provided is not in the default list. Please provide the outgroup file to be able to plot your tree.\n"
                                exit 1
                        else
                                echo -e "You didn't provide the fasta file of the outgroup but the amplicon name you provided is in the default list. Default outgroup from the list will be used.\n"
                                outgroup=./outgroup_sequences/petromyzon_marinus_${amplicon}.fa
                        fi
                else
                        echo -e "Fasta file of the outgroup detected. The following outgroup will be used: $outgroup"
                fi
	fi
fi

echo -e "Starting tree creation for the $family family identified as the $family_pacman in the pacman pipeline assignation\n"

mkdir ${family}_${amplicon}

# Check if genus_final.fa exist for this family
FILE=${family}_${amplicon}_genus_final.fa
if [ -f "$FILE" ]
then
	old=`stat -c %y -- $FILE | cut -d' ' -f1 | sed 's/-//g'`
	today=`date -I | sed 's/-//g'`
	let DIFF=($(date +%s -d $today)-$(date +%s -d $old))/86400
	if [ $DIFF -lt 30 ]
	then
		echo -e "$FILE exists, skipping creation steps of this file.\n"
		cd ${family}_${amplicon}
		cp ../${family}_${amplicon}_genus_final.fa ${family}_${amplicon}_genus_final.fa
	elif [ $DIFF -ge 30 ]
	then
		echo -e "$FILE exists but is too old, going to delete $FILE. Please relaunch the analysis after this.\n"
		rm $FILE
		exit 1
	fi
else
	echo -e "$FILE does not exist, starting creation of the file.\n----------\n"
	# Retrieve genus names of the chosen family
	## Get taxon ID of the family
	echo -e "Retrieving accession numbers of each sequences of interest in each genus of the $family family...\n----------\n"
	id=`echo $family | taxonkit name2taxid --data-dir . | cut -f2`
	if [ "$id" == "" ]
	then
		echo "This family does not exist in the NCBI taxonomy database, checking the corresponding name in the NCBI database..."
		all_taxo=`echo ${tax_table} | sed 's/${site}/'*'/g'`
		all_taxo=`echo $all_taxo | sed 's/${amplicon}/'$amplicon'/g'`
		for file in `ls $all_taxo`
		do
			cat $file >> identified_ASVs_${amplicon}.tsv
		done
		id=`awk -F'\t' -v fam="$family" ' $6 == fam { print $7 }' identified_ASVs_${amplicon}.tsv | sort | uniq | head -n1 | taxonkit name2taxid --data-dir . | cut -f2 | head -n1`
		ranklist=`echo $id | taxonkit lineage -R --data-dir . | cut -f2-3 | sed 's/\t/\n/g' | tac`
		column=`echo $ranklist | head -n1 | tr ';' '\n' | grep -n "^family$" | cut -d':' -f1`
		family_NCBI=`echo $ranklist | sed 's/no rank/norank/g' | sed 's/cellular organisms/cellularOrganisms/g' | sed 's/ /\n/g' | cut -d';' -f$column | tail -n1`
		echo "The corresponding name for the $family family in the NCBI database is $family_NCBI"
		id=`echo $family_NCBI | taxonkit name2taxid --data-dir . | cut -f2`
	fi

	## Extract all genus names of this family
	taxonkit list --ids $id -n -r --data-dir . | grep "\[genus\]" | cut -d']' -f2 | sed 's/ //g' | sort | uniq > ${family}_genus.txt

	## Move to family directory
	cd ${family}_${amplicon}

	# Retrieve all NCBI sequences containing "mitochondrial" or "mitochondrion" and "12S" or "complete genome"
	## Retrieve all the accession numbers
	while read line
	do
        	genus=`echo $line`
		grep "$genus" ../all_nt_db_acc.txt | grep -E 'mitochondrion|mitochondrial' > ${genus}_mito.txt # retrieve all sequences containing "mitochondrial" or "mitochondrion" in their name
		grep -E 'complete genome' ${genus}_mito.txt | sed 's/ >/\n>/g' | cut -d'.' -f1 | cut -d'>' -f2 >> ${family}_${amplicon}_acc-num.txt # from these sequences add in the final file the one containing "complete genome"
		grep -E "$gene_name" ${genus}_mito.txt | sed 's/ >/\n>/g' | cut -d'.' -f1 | cut -d'>' -f2 >> ${family}_${amplicon}_acc-num.txt #from these sequences add in the final file the one containing "12S"
		rm ${genus}_mito.txt
	done < ../${family}_genus.txt

	## Download the sequences retrieved in the accession file
	echo -e "Downloading all the sequences of the retrieved accession numbers...\n----------\n"
	sort ${family}_${amplicon}_acc-num.txt | uniq > ${family}_${amplicon}_acc-num_nodup.txt
	mkdir ncbi_acc_files
	while read line
	do
		acc=`echo $line`
		ncbi-acc-download --format fasta --out ncbi_acc_files/${acc}.fasta $acc
	done < ${family}_${amplicon}_acc-num_nodup.txt
	cat ncbi_acc_files/*.fasta > ${family}_${amplicon}_genus.fasta

	## Cleaning temporary files
	rm -r ncbi_acc_files/
	rm ${family}_${amplicon}_acc-num.txt

	## Put all sequences on one line each
	awk '/^>/ { if(NR>1) print "";  printf("%s\n",$0); next; } { printf("%s",$0);}  END {printf("\n");}' < ${family}_${amplicon}_genus.fasta > ${family}_${amplicon}_genus_1l.fa

	## Clean
	rm ${family}_${amplicon}_genus.fasta
	mv ${family}_${amplicon}_genus_1l.fa ${family}_${amplicon}_genus.fa

	# Retrieve sequence of interest
	echo -e "Cutting the retrieved sequences with Cutadapt to reduce the sequences to only the fragment of interest...\n----------\n"

        ## Detect forward primer
        cutadapt -g $primer_forward -o ${family}_${amplicon}_genus_short1.fa -e 4 -j 1 ${family}_${amplicon}_genus.fa > ${family}_${amplicon}_cutadapt1.log 2>&1

        ## Detect reverse primer
        cutadapt -a $primer_reverse -o ${family}_${amplicon}_genus_short2.fa -e 4 -j 1 ${family}_${amplicon}_genus_short1.fa > ${family}_${amplicon}_cutadapt2.log 2>&1

        ## Clean
        rm ${family}_${amplicon}_genus_short1.fa ${family}_${amplicon}_genus.fa
        mv ${family}_${amplicon}_genus_short2.fa ${family}_${amplicon}_genus_cut.fa

        # Refine sequence file to make it clean
	echo -e "Refining the sequences and changing sequence's names...\n----------\n"

	## Remove duplicated sequences
	seqkit rmdup -s < ${family}_${amplicon}_genus_cut.fa > ${family}_${amplicon}_genus_cut_nodup.fa

	## Put all sequences on one line each
	awk '/^>/ { if(NR>1) print "";  printf("%s\n",$0); next; } { printf("%s",$0);}  END {printf("\n");}' < ${family}_${amplicon}_genus_cut_nodup.fa > ${family}_${amplicon}_genus_cut_nodup_1l.fa

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
	awk 'BEGIN {RS = ">" ; FS = "\n" ; ORS = ""} $2 {print ">"$0}' ${family}_${amplicon}_genus_cut.fa > ../${family}_${amplicon}_genus_final.fa

	cp ../${family}_${amplicon}_genus_final.fa ${family}_${amplicon}_genus_final.fa
fi
echo -e "\nFile containing NCBI sequences for the family $family is created.\n\n----------\n----------\n"

# Add all ASV assigned to this family
## Get all the ASV identified in this family for each site
echo -e "Starting to retrieve ASVs assigned till at least the $family family...\n----------\n"
while read line
do
	site=`echo $line`
	loc_asv="../${site}_${amplicon}/identified_ASVs.txt"
	if [ -f "$loc_asv" ]
	then
		echo -e "File containing all identified ASVs till at least family level for $site is already created, checking if file with ASVs of the family is available...\n\n----------\n"
	else
		echo -e "File containing all identified ASVs till at least family level for $site is missing, starting to create it...\n\n----------\n"
		mkdir ../${site}_${amplicon}/
		taxo=`echo "../$tax_table" | sed 's/${site}/'$site'/g'`
		taxo=`echo $taxo | sed 's/${amplicon}/'$amplicon'/g'`
		colTax=`awk -v RS='\t' '/taxonRank/{print NR; exit}' "$taxo"`
		awk -F'\t' -v colTax=$colTax ' $colTax == "genus" || $colTax == "species" || $colTax == "family" { print $0 }' "$taxo" > ../${site}_${amplicon}/identified_ASVs.txt
	fi
	fam_asv="../${site}_${amplicon}/${family}_{amplicon}_asv.fa"
	if [ -f "$fam_asv" ]
	then
		echo -e "File containing all ASVs till at least family level for $family family of $site is available, retrieving the sequences from all sites...\n\n----------\n"
	else
		echo -e "File containing all ASVs till at least family level for $family family of $site is missing, starting to create it...\n\n----------\n"
		taxo=`echo "../${tax_table}" | sed 's/${site}/'$site'/g'`
                taxo=`echo $taxo | sed 's/${amplicon}/'$amplicon'/g'`
		colFam=`awk -v RS='\t' '/family/{print NR; exit}' "$taxo"`
		awk -F'\t' -v f=$family -v name=$colFam ' $name == f { print $1 }' ../${site}_${amplicon}/identified_ASVs.txt > ../${site}_${amplicon}/${family}_${amplicon}_asv.txt
		while read line
		do
			asv=`echo $line`
			colSpe=`awk -v RS='\t' '/scientificName/{print NR; exit}' "$taxo"`
			species=`grep -P "^$asv\t" $taxo | cut -f$colSpe | sed 's/ /_/g'`
			echo ">"$asv"_"$species"_${site}" >> ../${site}_${amplicon}/${family}_${amplicon}_asv.fa
			seq=`echo "../${seq_table}" | sed 's/${site}/'$site'/g'`
			seq=`echo $seq | sed 's/${amplicon}/'$amplicon'/g'`
			colSeq=`awk -v RS='\t' '/DNA_sequence/{print NR; exit}' "$seq"`
			grep "${asv}_" $seq | head -n 1 | cut -f$colSeq >> ../${site}_${amplicon}/${family}_${amplicon}_asv.fa
		done < ../${site}_${amplicon}/${family}_${amplicon}_asv.txt
	fi
done < ../all_sites.txt

## Put all the ASVs in the same file
cat ../*_${amplicon}/${family}_${amplicon}_asv.fa > all_${family}_${amplicon}_asv.fa

## Add the ASVs to the final fasta file for the alignment
cat all_${family}_${amplicon}_asv.fa >> ${family}_${amplicon}_genus_final.fa

# Add the outgroup
echo -e "Adding the outgroup...\n----------\n"
cat ../$outgroup >> ${family}_${amplicon}_genus_final.fa

# Remove the sequences longer than 200 nucleotides
echo -e "Refining the final fasta file containing ASVs and sequences from NCBI database...\n"
seqkit seq -M $length_amp ${family}_${amplicon}_genus_final.fa | awk '/^>/ { if(NR>1) print "";  printf("%s\n",$0); next; } { printf("%s",$0);}  END {printf("\n");}' > ${family}_${amplicon}_genus_final_short.fa

# Align the sequences with Mafft
echo -e "\nAligning with MAFFT...\n"
mafft --thread -1 --preservecase --adjustdirection ${family}_${amplicon}_genus_final_short.fa > ${family}_${amplicon}.aln.fa 2> ${family}_${amplicon}_mafft.log

## Replace spaces by underscores
sed -i 's/ /_/g' ${family}_${amplicon}.aln.fa

# Curate the alignment with Gblocks
echo -e "Curating the alignment with Gblocks...\n"
Gblocks ${family}_${amplicon}.aln.fa -t=d -b4=5 -b5=h -e=.gb > ${family}_${amplicon}_gblock.log 2>&1

# Reconstruct phylogenetic tree with FastTree
echo -e "Constructing the phylogenetic tree with FastTree...\n"
FastTreeMP -gtr -nt ${family}_${amplicon}.aln.fa.gb > ${family}_${amplicon}.nwk 2>${family}_${amplicon}_FastTree.log

## Remove "_R_"
sed -i 's/_R_//g' ${family}_${amplicon}.nwk

# Create information file to differentiate ASVs and sequences from NCBI database
## Retrieve all names in the tree
echo -e "Creating information file to color differently ASVs and sequences from NCBI database...\n----------\n"
grep "^>" ${family}_${amplicon}.aln.fa.gb | sed 's/>//g' | sed 's/_R_//g' > all_labels.txt

## For each name determine if it's an ASV or not and put the information in a file
while read line
do
	name=`echo $line`
	grep_res=`echo $name | grep "asv."`
	if [ -z $grep_res ]
	then
		echo -e "${name}\tdb_name" >> ${family}_${amplicon}_info.tsv
	else
		echo -e "${name}\tasv" >> ${family}_${amplicon}_info.tsv
	fi
done < all_labels.txt

## Add column's names
sed -i '1itaxa\tasv' ${family}_${amplicon}_info.tsv

# Plot the tree
echo -e "Plotting the tree with ggtree package in R...\n"
## Activate environment
mamba activate R-tree-creation

## Run R script
path=`pwd`; Rscript ../plot_family_tree.r ${path} ${family}_${amplicon}.nwk ${family}_${amplicon}_info.tsv ${family} ${family}_${amplicon}.pdf

## Clean temporary files
rm Rplots.pdf

# Check if final file produced
ENDFILE=${family}_${amplicon}.pdf
if [ -f "$ENDFILE" ]
then
	echo -e "\nTree creation script worked, $ENDFILE is available in ${family}_${amplicon}/ folder!\n"
else
	echo -e "\nA problem occurred! $ENDFILE was not created. Please look at the log file to understand what happend."
fi

