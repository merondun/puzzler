#!/bin/bash

#SBATCH --time=48:00:00   
#SBATCH --nodes=1  
#SBATCH --ntasks-per-node=32
#SBATCH --mem=64Gb
#SBATCH --partition=ceres

# use n=32 and mem = 64 typically 
SINGULARITY_TMPDIR=${APPTAINER_TMPDIR}

if [ -z "$1" ]; then
    echo "Error: sample ID positinonal argument is required."
    exit 1
fi

module load miniconda
source activate puzz
module load apptainer

t=32
MEM=64
SAMPLE=$1
# SET THIS!!! 
SPECIES="Guava"

WD=/project/coffea_pangenome/${SPECIES}/Assemblies/20250101_JustinAssemblies
SAMPLE_FILE=/project/coffea_pangenome/${SPECIES}/Assemblies/20250101_JustinAssemblies/samples.csv
PLOIDY=$(awk -F',' -v sample="$SAMPLE" '$1 == sample {print $2}' ${SAMPLE_FILE})
NCHRS=$(awk -F',' -v sample="$SAMPLE" '$1 == sample {print $3}' ${SAMPLE_FILE})
HOM_COV=$(awk -F',' -v sample="$SAMPLE" '$1 == sample {print $4}' ${SAMPLE_FILE})
HIFI=$(awk -F',' -v sample="$SAMPLE" '$1 == sample {print $5}' ${SAMPLE_FILE})
HIC_R1=$(awk -F',' -v sample="$SAMPLE" '$1 == sample {print $6}' ${SAMPLE_FILE})
HIC_R2=$(awk -F',' -v sample="$SAMPLE" '$1 == sample {print $7}' ${SAMPLE_FILE})
REFERENCE=$(awk -F',' -v sample="$SAMPLE" '$1 == sample {print $8}' ${SAMPLE_FILE})
PUZZLER="apptainer exec /project/coffea_pangenome/Software/Merondun/apptainers/puzzler_v1.1.sif"

cat << "EOF"

=======================================================================
__________ ____ _______________________.____     _____________________ 
\______   \    |   \____    /\____    /|    |    \_   _____/\______   \
 |     ___/    |   / /     /   /     / |    |     |    __)_  |       _/
 |    |   |    |  / /     /_  /     /_ |    |___  |        \ |    |   \
 |____|   |______/ /_______ \/_______ \|_______ \/_______  / |____|_  /
                           \/        \/        \/        \/         \/ 
=======================================================================

EOF

echo -e "=======================================================================\nParameters for sample: ${SAMPLE} \nPLOIDY: ${PLOIDY} \nHOM_COV: ${HOM_COV}\nHIFI: ${HIFI}\nHIC_R1: ${HIC_R1}\nHIC_R2: ${HIC_R2}\nREFERENCE: ${REFERENCE}\n=======================================================================\n"

mkdir -p ${WD}/joint_scaffold ${WD}/logs/contact_maps
	
cd ${WD}/${SAMPLE}/01_scaffolding

if [ ! -s ${SAMPLE}_JBAT.review.assembly ]; then
	echo -e "\e[43m~~~~ Post curation assembly file doesn't exist for ${SAMPLE} ~~~~\e[0m" 
else
	if [ ! -s map.txt ]; then
		echo -e "\e[43m~~~~ Extracting final assembly and mapping to reference for ${SAMPLE} ~~~~\e[0m"
		JUICER_PATH=~/symlinks/software/HapHiC/utils
		${JUICER_PATH}/juicer post -o haphic-refsort-post_JBAT ${SAMPLE}_JBAT.review.assembly haphic-refsort_JBAT.liftover.agp all.purged.fa 2> juicer.post.log

		# Renaming: to reference chromosomes 
		$PUZZLER minimap2 -x asm20 ${REFERENCE} haphic-refsort-post_JBAT.FINAL.fa --secondary=no -t ${t} -o 02_orienting/asmpost_to_paf.paf 2> minimap.postjuicer.log
		samtools faidx haphic-refsort-post_JBAT.FINAL.fa
		map_chromosomes --paf 02_orienting/asmpost_to_paf.paf --fai haphic-refsort-post_JBAT.FINAL.fa.fai --out map.txt
	else
		echo -e "\e[42m~~~~ Skipping assembly and mapping for ${SAMPLE}, already exists ~~~~\e[0m"
	fi 

	if [ ! -s haphic_renamed.fa ] && [ -s map.txt ]; then
		echo -e "\e[43m~~~~ Renaming chromosomes for ${SAMPLE} ~~~~\e[0m"

		# This file has $haphic_scafID \t $ref_chr_ID \t $strand
		awk '{OFS="\t"}{print $1, $2, $5, $3, $4}' map.txt > 02_orienting/hap_newID_map.txt 

		# Exclude those chromsomes from map, and then add e.g. h3tg0120320 for the remaining scaffolds 
		awk '{print $1}' 02_orienting/hap_newID_map.txt > 02_orienting/exclude_chr_scafIDs.txt
		grep -vwf 02_orienting/exclude_chr_scafIDs.txt haphic-refsort-post_JBAT.FINAL.fa.fai | awk '{OFS="\t"}{print $1, $1, "+"}' >> 02_orienting/hap_newID_map.txt

		# Must be the same number of scaffolds
		if [ $(cat 02_orienting/hap_newID_map.txt | wc -l) -eq $(cat haphic-refsort-post_JBAT.FINAL.fa.fai | wc -l) ]; then
			echo "Same number of scaffolds, proceeding!"
		else
			echo "Not same number, stop "
			#exit 1
		fi # exit scaffold check 

		### INSPECT! 
		#cat ${WD}/${SAMPLE}/01_scaffolding/02_orienting/hap_newID_map.txt

		# Grab the reverse strands to reverse complement them 
		awk '{if ($3 == "-") print $1}' 02_orienting/hap_newID_map.txt > 02_orienting/reverse_list.txt
		awk '{if ($3 == "+") print $1}' 02_orienting/hap_newID_map.txt > 02_orienting/positive_list.txt
		seqtk subseq haphic-refsort-post_JBAT.FINAL.fa 02_orienting/reverse_list.txt > 02_orienting/to_revercomp.fa
		seqtk subseq haphic-refsort-post_JBAT.FINAL.fa 02_orienting/positive_list.txt > 02_orienting/positive.fa
		seqkit seq --line-width 0 -t DNA -v -r -p 02_orienting/to_revercomp.fa > 02_orienting/revcomp.fa
		cat 02_orienting/positive.fa 02_orienting/revcomp.fa > 02_orienting/orient.fa

		# Detect if there are duplicate scaffolds
		awk '{print $2}' 02_orienting/hap_newID_map.txt | sort | uniq -d > 02_orienting/duplicates.txt

		##### THIS WILL DEAL WITH DUPLICATES, IN CASE THERE ARE 2 SCAFFOLDS CORRESPONDING TO A CHR! eg. Chr1A  
		if [ -s 02_orienting/duplicates.txt ]; then
			
			echo "Multiple scaffolds corresponding to a single Chr for ${SAMPLE}, INSPECT!" > ${WD}/logs/${SAMPLE}.haphic-duplicates.log
			echo -e "\e[41m~~~~ Multiple scaffolds corresponding to a single Chr for ${SAMPLE}, INSPECT!  ~~~~\e[0m"
			exit 1
			
		else 

			echo "Single scaffolds corresponding to a single Chr for ${SAMPLE}"

			# Rename chrs and sort 
			cut -f2 02_orienting/hap_newID_map.txt | sort -V > 02_orienting/chr_sort.list
			seqkit replace --line-width 0 -p "(.*)" -r "{kv}" -k 02_orienting/hap_newID_map.txt 02_orienting/orient.fa > haphic_renamed_unord.fa
			samtools faidx haphic_renamed_unord.fa
			xargs samtools faidx haphic_renamed_unord.fa < 02_orienting/chr_sort.list > haphic_renamed.fa
			rm haphic_renamed_unord.fa haphic_renamed_unord.fa.fai 

		fi  # exit duplicates check

	else
		echo -e "\e[42m~~~~ Skipping renaming chromosomes for ${SAMPLE}, already exists  ~~~~\e[0m"
	fi # exit renaming check 

	cd ${WD}/${SAMPLE}/01_scaffolding
	if [ ! -s haphic_renamed.filtered.bam ] && [ -s haphic_renamed.fa ]; then

		echo -e "\e[43m~~~~ Creating bam for ${SAMPLE} ~~~~\e[0m"

		$PUZZLER bwa index haphic_renamed.fa 2> alignment.indexing.final.log
		$PUZZLER samtools faidx haphic_renamed.fa

		# Align Hi-C reads
		$PUZZLER bwa mem -5SP -t ${t} haphic_renamed.fa ${HIC_R1} ${HIC_R2} 2> alignment.final.log | \
			$PUZZLER samblaster | \
			$PUZZLER samtools view - -@ ${t} -S -h -b -F 3340 | \
			$PUZZLER filter_bam - 1 --nm 3 --threads ${t} --remove-dup | \
			$PUZZLER samtools view - -b -@ ${t} -o haphic_renamed.filtered.bam

	else
		echo -e "\e[42m~~~~ Skipping alignment for final map ${SAMPLE}, already exists ~~~~\e[0m"
	fi # exit final hic mapping check 

	if [ ! -f ${WD}/logs/contact_maps/${SAMPLE}.pdf ]; then

		echo -e "\e[43m~~~~ Creating final contact map for ${SAMPLE} ~~~~\e[0m"
		~/symlinks/software/HapHiC/utils/mock_agp_file.py haphic_renamed.fa > haphic_renamed.agp
		$PUZZLER haphic plot --threads ${t} haphic_renamed.agp haphic_renamed.filtered.bam --bin_size 1000 --min_len 2

		# Copy over the file with the correct pangenome-spec naming
		cp haphic_renamed.fa ${WD}/joint_scaffold/${SAMPLE}.fa
		$PUZZLER samtools faidx ${WD}/joint_scaffold/${SAMPLE}.fa
		cp contact_map.pdf ${WD}/logs/contact_maps/${SAMPLE}.pdf
	else
		echo -e "\e[42m~~~~ Skipping final contact map generating for ${SAMPLE}, already exists ~~~~\e[0m"
	fi # exit final hic plotting check 

fi # exit entire loop file assembly check 