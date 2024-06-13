#! /bin/bash

# Move to working directory
cd /lustre1/project/stg_00026/enora/UNESCO/tree_creation

# Load conda environment
source ~/.bashrc
conda activate unesco-trees

# Choice of the location
location=$1
echo -e "\n\nCreation of a phylogenetic tree for the location ${location} starting...\n\n----------\n"

# Create working directory for this location and move into it
mkdir $location
cd $location

# Create file containing all ASV sequences of this sample
## Get column numbers of interesting fields
DNA_col=`awk -v RS='\t' '/DNA_sequence/{print NR; exit}' ../../results_pacman_pip/eDNAexpeditions_batch2_samples/runs/${location}_12SMifish/05-dwca/Full_tax_table_with_lsids.tsv`
spe_col=`awk -v RS='\t' '/scientificName/{print NR; exit}' ../../results_pacman_pip/eDNAexpeditions_batch2_samples/runs/${location}_12SMifish/05-dwca/Full_tax_table_with_lsids.tsv`
rank_col=`awk -v RS='\t' '/taxonRank/{print NR; exit}' ../../results_pacman_pip/eDNAexpeditions_batch2_samples/runs/${location}_12SMifish/05-dwca/Full_tax_table_with_lsids.tsv`
DNA_col=$(( $DNA_col + 1 ))
spe_col=$(( $spe_col + 1 ))
rank_col=$(( $rank_col + 1 ))

## Retrieve ASV sequences
echo -e "Starting to retrieve all the ASV sequences for this location...\n\n----------\n"
grep "^asv" ../../results_pacman_pip/eDNAexpeditions_batch2_samples/runs/${location}_12SMifish/05-dwca/Full_tax_table_with_lsids.tsv | sed 's/\t/|/g' | while read line; do asv=`echo $line | cut -d'|' -f1`; species=`echo $line | cut -d'|' -f${spe_col} | sed 's/ /_/g'`; echo -e ">${asv}_${species}" >> ${location}_ASV.fa; seq=`echo $line | cut -d'|' -f${DNA_col}`; echo $seq >> ${location}_ASV.fa; done

## Remove duplicated sequences
seqkit rmdup -s < ${location}_ASV.fa > ${location}_ASV_nodup.fa

## Transform sequences to get DNA on one line for each ASV
awk '/^>/ { if(NR>1) print "";  printf("%s\n",$0); next; } { printf("%s",$0);}  END {printf("\n");}' < ${location}_ASV_nodup.fa > ${location}_ASV_nodup_correct.fa

# Remove the ASVs identified as human
while read line
do
	grep_res=`echo $line | grep -v "asv." | head -n1 | cut -d' ' -f1`
	if [ -z $grep_res ] # if lines contains asv. do the following
	then
		species=`echo $line | cut -d'_' -f2-` # get species name
		name=`echo $line | cut -d'>' -f2` # get the complete sequence name
		if [[ $species == "Homo_sapiens" || $species == "Incertae_sedis" ]] # check if species is human or assigned incertain
		then
			num=`grep -n "$name" ${location}_ASV_nodup_correct.fa | cut -d":" -f1` # if it is retrieve line number of the asv sequence
			other=$(( $num + 1 ))
			sed -i ''$num','$other'd' ${location}_ASV_nodup_correct.fa # and remove these lines
		fi
	fi
done < ${location}_ASV_nodup_correct.fa

## Remove temporary files
rm ${location}_ASV.fa ${location}_ASV_nodup.fa
mv ${location}_ASV_nodup_correct.fa ${location}_ASV.fa

# Retrieve NCBI sequences to add known species to the tree
echo -e "Starting to retrieve the NCBI sequences corresponding to the genus or species identified with blast and vsearch...\n\n----------\n"
## Retrieve the accession numbers
touch temp_spe_file.txt
while read line
do
	grep_res=`echo $line | grep -v "asv."`
	if [ -z $grep_res ]
	then
		asv=`echo $line | cut -d'_' -f1 | cut -d'>' -f2`
		rank=`grep -P "^$asv\t" ../../results_pacman_pip/eDNAexpeditions_batch2_samples/runs/${location}_12SMifish/05-dwca/Full_tax_table_with_lsids.tsv | cut -f${rank_col}`
		genus=`echo $line | cut -d'_' -f2`
		grep_gen=`grep "$genus" temp_spe_file.txt | head -n1 | cut -d' ' -f1`
		if [ -z $grep_gen ]
		then
			if [ "$rank" == "genus" ]
			then
				grep "$genus " ../all_nt_db_acc.txt | grep -E 'genome|12S' | grep "mitochon" | sort -u -t' ' -k3 | cut -d' ' -f1 | cut -d'>' -f2 >> ${location}_accessions.txt
				echo $genus >> temp_spe_file.txt
			elif [ "$rank" == "species" ]
			then
				# spe=`echo $line | cut -d'_' -f2- | sed 's/_/ /g'`
				grep "$genus " ../all_nt_db_acc.txt | grep -E 'genome|12S' | grep "mitochon" | sort -u -t' ' -k3 | cut -d' ' -f1 | cut -d'>' -f2 >> ${location}_accessions.txt
				echo $genus >> temp_spe_file.txt
			fi
		fi
	fi
done < ${location}_ASV.fa

## Remove temporary files
# rm temp_spe_file.txt

## Download sequences of these accession numbers
echo -e "Downloading the sequences from NCBI...\n\n----------\n"
sort ${location}_accessions.txt | uniq > ${location}_acc_nodup.txt
mkdir ncbi_acc_files
while read line
do
	acc=`echo $line`
	ncbi-acc-download --format fasta --out ncbi_acc_files/${acc}.fasta $acc
done < ${location}_acc_nodup.txt
cat ncbi_acc_files/*.fasta > ${location}_organisms.fasta

## Cleaning temporary files
# rm -r ncbi_acc_files/
# rm ${location}_accessions.txt

awk '/^>/ { if(NR>1) print "";  printf("%s\n",$0); next; } { printf("%s",$0);}  END {printf("\n");}' < ${location}_organisms.fasta > ${location}_organisms_1l.fa

## Clean
rm ${location}_organisms.fasta
mv ${location}_organisms_1l.fa ${location}_organisms.fa

## Rename the sequences to keep only the organism genus and species
sed -i 's/PREDICTED: /PREDICTED:/g' ${location}_organisms.fa
while read line
do
        grep_res=`echo $line | grep -v "^>"`
        if [ -z $grep_res ]
        then
                name=`echo $line | cut -d' ' -f2-3`
                sed -i "s|$line|>${name}|g" ${location}_organisms.fa > sed_output.txt 2>&1
        fi
done < ${location}_organisms.fa

# Put all sequences (organisms and asv) in one file and remove Human sequences
echo -e "Putting all the sequences in one file and removing some unwanted taxa like Hominidae...\n\n----------\n"
## Add asv sequences to the same file
cat ${location}_ASV.fa >> ${location}_organisms_final.fa
cat ${location}_organisms.fa >> ${location}_organisms_final.fa

## Remove all unwanted sequences
seqkit grep -rvif ../removed_organisms.txt ${location}_organisms_final.fa > ${location}_final.fa
awk '/^>/ { if(NR>1) print "";  printf("%s\n",$0); next; } { printf("%s",$0);}  END {printf("\n");}' < ${location}_final.fa > ${location}_final_correct.fa

## Retrieve all unwanted sequences
# grep -A1 -f ../removed_organisms.txt ${location}_organisms_final.fa > ${location}_bad-species.fa

## Remove the unwanted sequences
# grep -v -f ${location}_bad-species.fa ${location}_organisms_final.fa > ${location}_final.fa

## Remove all empty sequences
# awk 'BEGIN {RS = ">" ; FS = "\n" ; ORS = ""} $2 {print ">"$0}' ${location}_final.fa > ${location}_final_correct.fa

## Remove temporary files
# rm ${location}_bad-species.fa
rm ${location}_final.fa
mv ${location}_final_correct.fa ${location}_final.fa

# Apply cutadapt on all the sequences
echo -e "Applying cutadapt to the sequences in order to keep only the fragment of interest...\n\n----------\n"
## activate the conda environment
conda activate cutadapt

## remove forward primer
cutadapt -g GTCGGTAAAACTCGTGCCAGC -o ${location}_final-short.fa -e 4 -j 8 ${location}_final.fa > ${location}_cutadapt1.log 2>&1

## remove reverse primer
cutadapt -a CAAACTGGGATTAGATACCCCACTATG -o ${location}_final-short2.fa -e 4 -j 8 ${location}_final-short.fa > ${location}_cutadapt2.log 2>&1

## Remove temporary files
rm ${location}_final-short.fa

## re-activate unesco-trees conda environment
conda activate unesco-trees

# remove too long sequences which were not cutted with cutadapt (no amplicon in seq?)
seqkit seq -M 250 ${location}_final-short2.fa | awk '/^>/ { if(NR>1) print "";  printf("%s\n",$0); next; } { printf("%s",$0);}  END {printf("\n");}' > ${location}_final-short2_correct.fa

# remove duplicated organisms
echo -e "Removing duplicate sequences and duplicated names...\n\n----------\n"
seqkit rmdup -n < ${location}_final-short2_correct.fa > ${location}_final-short_nodup.fa

## Put the sequences on one line each
awk '/^>/ { if(NR>1) print "";  printf("%s\n",$0); next; } { printf("%s",$0);}  END {printf("\n");}' < ${location}_final-short_nodup.fa > ${location}_final-short2_nodup_correct.fa

## Remove temporary files
rm ${location}_final-short2.fa ${location}_final-short_nodup.fa
mv ${location}_final-short2_nodup_correct.fa ${location}_final-short.fa

# Add the outgroup to the file
echo -e "Adding the outgroup...\n\n----------\n"
cat ../U11880_cut2.fa >> ${location}_final-short.fa

# Create a basic tree (without bootstraps)
## Align the sequences with Mafft
echo -e "Aligning the sequences with Mafft...\n\n----------\n"
mafft --thread -1 --preservecase --adjustdirection ${location}_final-short.fa > ${location}.aln.fa 2> ${location}_mafft.log

## Replace spaces by underscores in the organisms names
sed -i 's/ /_/g' ${location}.aln.fa

## Removed special characters in the names
sed -i 's/://g' ${location}.aln.fa
sed -i 's/_R_asv./asv./g' ${location}.aln.fa
sed -i 's/_R_//g' ${location}.aln.fa

## Curate the alignments with Gblocks
echo -e "Refining the alignment with Gblocks...\n\n----------\n"
Gblocks ${location}.aln.fa -t=d -b4=5 -b5=h -e=.gb > ${location}_gblock.log 2>&1

## Create the tree with Fasttree
echo -e "Creating the phylogenetic tree with Fasttree...\n\n----------\n"
FastTreeMP -gtr -nt ${location}.aln.fa.gb > ${location}.nwk 2>${location}_FastTree.log

# Create a file containing real names of the organisms and a shorter version (necessary for fseqboot because use Phylip format)
#grep "^>" ${sample}.aln.fa.gb | sed 's/>//g' > ${sample}_real_names.txt
#while read line
#do
#	genus=`echo $line | cut -d'_' -f1`
#	gen=${genus:0:3}
#	species=`echo $line | cut -d'_' -f2`
#	spe=${species:0:6}
#	echo $gen"_"$spe >> ${sample}_short_names.txt
#done < ${sample}_real_names.txt
#paste ${sample}_real_names.txt ${sample}_short_names.txt > ${sample}_names.txt

## Remove temporary files
#rm ${sample}_real_names.txt ${sample}_real_species_names.txt ${sample}_short_names.txt ${sample}_short_names_sort.txt ${sample}_real_names_sort.txt

# Resample the tree and create a final tree with bootstrap values ($repeats repeats)
## Change the names to the short versions
#while read line
#do
#	name=`echo $line | cut -d' ' -f1`
#	short=`echo $line | cut -d' ' -f2`
#	sed -i "s|$name|$short|g" ${sample}.aln.fa.gb
#done < ${sample}_names.txt

## Apply seqboot to resample the tree
#fseqboot -reps ${repeats} -sequence ${sample}.aln.fa.gb -outfile ${sample}_${repeats}b.fseqboot

## Run Fasttree to get all the trees
#FastTreeMP -gtr -nt -n ${repeats} ${sample}_${repeats}b.fseqboot > ${sample}_${repeats}b.nwk

## Change the names back to their long versions
#while read line
#do
#	name=`echo $line | cut -d' ' -f1`
#	short=`echo $line | cut -d' ' -f2`
#	sed -i "s|$short|$name|g" ${sample}_${repeats}b.nwk
#done < ${sample}_names.txt

## Run script to get the final tree with bootstrap values
#perl -I ../../mafft_align/ ../../mafft_align/CompareToBootstrap.pl -tree ${sample}.nwk -boot ${sample}_${repeats}b.nwk > ${sample}_${repeats}b-final.nwk

# Create the file to retrieve the number of sampling site for each ASV
echo -e "Creating the table file with quantitative information on the ASV (number of site detecting the ASV)...\n\n----------\n"
## Retrieve all the sequence names
grep "^>" ${location}.aln.fa.gb | sed 's/>//g' > ${location}_quant_names.tsv

## For each ASV retrieve the relative abundance
while read line
do
	asv=`echo $line | cut -d'_' -f1`
	sites=`grep "${asv}_EE" ../../results_pacman_pip/eDNAexpeditions_batch2_samples/runs/${location}_12SMifish/05-dwca/Occurrence_table.tsv | wc -l`
	echo "$sites" >> ${location}_quant_num.tsv
done < ${location}_quant_names.tsv

# Assemble names and percentage in the same file and give headers to the table
paste ${location}_quant_names.tsv ${location}_quant_num.tsv > ${location}_quant.tsv
sed -i '1itaxa\tquantSites' ${location}_quant.tsv

# Remove temporary files
rm ${location}_quant_names.tsv ${location}_quant_num.tsv

# Produce a file with the tree visualization
echo -e "Creating the tree image...\n\n----------\n"
## Change to conda environment for tree construction
mamba activate R-unesco-trees

## call the R script
path=`pwd`
Rscript ../plot_tree_location.r ${path} ${location}.nwk ${location}_quant.tsv ${location} ${location}_root.pdf

## Remove temporary files
rm Rplots.pdf

echo -e "Analysis finished! You can find the tree image in the folder ${location}/ under the name ${location}_root.png"
