#!/bin/bash

#SBATCH --time=4:00:00   
#SBATCH --nodes=1  
#SBATCH --ntasks-per-node=32
#SBATCH --mem=128Gb
#SBATCH --partition=short

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
SAMPLE=$1

WD=/project/coffea_pangenome/Artocarpus/Assemblies/20250101_JustinAssemblies
SAMPLE_FILE="/project/coffea_pangenome/Artocarpus/Assemblies/20250101_JustinAssemblies/samples.csv"
PLOIDY=$(awk -F',' -v sample="$SAMPLE" '$1 == sample {print $2}' ${SAMPLE_FILE})
NCHRS=$(awk -F',' -v sample="$SAMPLE" '$1 == sample {print $3}' ${SAMPLE_FILE})
HOM_COV=$(awk -F',' -v sample="$SAMPLE" '$1 == sample {print $4}' ${SAMPLE_FILE})
HIFI=$(awk -F',' -v sample="$SAMPLE" '$1 == sample {print $5}' ${SAMPLE_FILE})
HIC_R1=$(awk -F',' -v sample="$SAMPLE" '$1 == sample {print $6}' ${SAMPLE_FILE})
HIC_R2=$(awk -F',' -v sample="$SAMPLE" '$1 == sample {print $7}' ${SAMPLE_FILE})
REFERENCE=$(awk -F',' -v sample="$SAMPLE" '$1 == sample {print $8}' ${SAMPLE_FILE})
PUZZLER="apptainer exec /project/coffea_pangenome/Software/Merondun/apptainers/puzzler_v1.1.sif"
TOTAL_CHRS=$((NCHRS * PLOIDY))

mkdir -p ${WD}/joint_scaffold ${WD}/logs/contact_maps

for IT in pri hap; do 
			
	cd ${WD}/${SAMPLE}/02_${IT}HapHiC

	if [ ! -s ${SAMPLE}.${IT}-MQ1_JBAT.review.assembly ]; then
		echo "~~~~ Post curation assembly file doesn't exist for ${SAMPLE} ${IT} ~~~~" 
	else
		mkdir -p 02_orienting
		if [ ! -s map.txt ]; then
			echo "~~~~ Extracting final ${IT} assembly and mapping to reference for ${SAMPLE} ~~~~"
			JUICER_PATH=~/symlinks/software/HapHiC/utils
			${JUICER_PATH}/juicer post -o haphic-refsort-post_JBAT ${SAMPLE}.${IT}-MQ1_JBAT.review.assembly haphic-refsort-MQ1_JBAT.liftover.agp all.purged.fa

			# Renaming: to reference chromosomes 
			$PUZZLER minimap2 -x asm20 ${REFERENCE} haphic-refsort-post_JBAT.FINAL.fa --secondary=no -t ${t} -o 02_orienting/asmpost_to_paf.paf
			samtools faidx haphic-refsort-post_JBAT.FINAL.fa
			map_chromosomes --paf 02_orienting/asmpost_to_paf.paf --fai haphic-refsort-post_JBAT.FINAL.fa.fai --out map.txt
		else
			echo "~~~~ Skipping assembly and mapping for ${SAMPLE} ${IT} ~~~~"
		fi 

		if [ ! -s haphic_renamed.fa ] && [ -s map.txt ]; then
			echo "~~~~ Renaming chromosomes for ${SAMPLE} ${IT} ~~~~"

			# This file has $haphic_scafID \t $ref_chr_ID \t $strand
			awk '{OFS="\t"}{print $1, $2, $5, $3, $4}' map.txt > 02_orienting/hap_ref_map.txt 
			
			if [ "${IT}" = "hap" ]; then

				# Initialize a counter for each chromosome
				declare -A chr_counter

				# Assign each 'chr' entry an unique haplotype ID 
				awk -v sample="$SAMPLE" '
				BEGIN { OFS="\t" }
				{
					chr = $2
					if (!(chr in chr_counter)) {
						chr_counter[chr] = 1
					} else {
						chr_counter[chr]++
					}
					hap = chr_counter[chr]
					print $1, sample "#" hap "#" chr, $3
				}' 02_orienting/hap_ref_map.txt > 02_orienting/hap_newID_map.txt

			else
				cp 02_orienting/hap_ref_map.txt 02_orienting/hap_newID_map.txt
			fi  # exit haplotype check for chr renaming

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

			# Grab the reverse strands to reverse complement them 
			awk '{if ($3 == "-") print $1}' 02_orienting/hap_newID_map.txt > 02_orienting/reverse_list.txt
			awk '{if ($3 == "+") print $1}' 02_orienting/hap_newID_map.txt > 02_orienting/positive_list.txt
			seqtk subseq haphic-refsort-post_JBAT.FINAL.fa 02_orienting/reverse_list.txt > 02_orienting/to_revercomp.fa
			seqtk subseq haphic-refsort-post_JBAT.FINAL.fa 02_orienting/positive_list.txt > 02_orienting/positive.fa
			seqkit seq --line-width 0 -t DNA -v -r -p 02_orienting/to_revercomp.fa > 02_orienting/revcomp.fa
			cat 02_orienting/positive.fa 02_orienting/revcomp.fa > 02_orienting/orient.fa

			# Detect if there are duplicate scaffolds
			awk '{print $2}' 02_orienting/hap_newID_map.txt | sort | uniq -d > 02_orienting/duplicates.txt

			##### THIS WILL DEAL WITH DUPLICATES, IN CASE THERE ARE 2 SCAFFOLDS CORRESPONDING TO A CHR FROM THE SAME HAPLOTYPE! 
			if [ -s 02_orienting/duplicates.txt ]; then
				
				echo "Multiple scaffolds corresponding to a single Chr for ${SAMPLE} ${IT}, INSPECT!" > ${WD}/logs/${SAMPLE}.${IT}.haphic-duplicates.log

			else 

				echo "Single scaffolds corresponding to a single Chr for ${SAMPLE} ${IT}"

				# Rename chrs and sort 
				seqkit replace --line-width 0 -p "(.*)" -r "{kv}" -k 02_orienting/hap_newID_map.txt 02_orienting/orient.fa | \
					seqkit sort --line-width 0 -n > haphic_renamed.fa

				# For final contact if hap assembly, also sort by CHR for easy chr comparison, header is: SAMPLE # HAP # CHR, so swap CHR and HAP temporarily for name sorting
				if [ "${IT}" = "hap" ]; then
					# Switch headers 
					cat 02_orienting/hap_newID_map.txt | \
						tr '#' '\t' | \
						awk '{OFS="\t"}{ if ($4 ~ /Chr/) {print $2"#"$3"#"$4,$2"#"$4"#"$3,$5} else{print $$0}}' > 02_orienting/hap_samphapchr_map.txt
					seqkit replace --line-width 0 -p "(.*)" -r "{kv}" -k 02_orienting/hap_samphapchr_map.txt haphic_renamed.fa | \
							seqkit sort --line-width 0 -n > pg_renamed.fa
				else
					ln -s haphic_renamed.fa pg_renamed.fa
				fi # exit haplotype check 
			fi  # exit duplicates check

		else
			echo "~~~~ Skipping renaming chromosomes for ${SAMPLE} ${IT}, already exists  ~~~~"
		fi # exit renaming check 

		cd ${WD}/${SAMPLE}/02_${IT}HapHiC
		if [ ! -s pg_renamed.filtered.bam ] && [ -s pg_renamed.fa ]; then

			echo "~~~~ Creating bam for ${SAMPLE} ${IT} ~~~~"

			$PUZZLER bwa index pg_renamed.fa 2> alignment.indexing.finalpg.log
			$PUZZLER samtools faidx pg_renamed.fa

			# Align Hi-C reads
			$PUZZLER bwa mem -5SP -t ${t} pg_renamed.fa ${HIC_R1} ${HIC_R2} | \
				$PUZZLER samblaster | \
				$PUZZLER samtools view - -@ ${t} -S -h -b -F 3340 | \
				$PUZZLER filter_bam - 1 --nm 3 --threads ${t} --remove-dup | \
				$PUZZLER samtools view - -b -@ ${t} -o pg_renamed.filtered.bam 2> alignment.finalpg.log

		else
			echo "~~~~ Bam already created for ${SAMPLE} ${IT} ~~~~"
		fi # exit final hic mapping check 

		if [ ! -f ${WD}/logs/contact_maps/${SAMPLE}.${IT}.pdf ]; then

			echo "~~~~ Creating final contact map for ${SAMPLE} ${IT} ~~~~"
			~/symlinks/software/HapHiC/utils/mock_agp_file.py pg_renamed.fa > pg_renamed.agp
			$PUZZLER haphic plot --threads ${t} pg_renamed.agp pg_renamed.filtered.bam --bin_size 1000 --min_len 2

			# Copy over the file with the correct pangenome-spec naming
			cp haphic_renamed.fa ${WD}/joint_scaffold/${SAMPLE}.${IT}.fa
			$PUZZLER samtools faidx ${WD}/joint_scaffold/${SAMPLE}.${IT}.fa
			cp contact_map.pdf ${WD}/logs/contact_maps/${SAMPLE}.${IT}.pdf
		else
			echo "~~~~ Final contact map already generated for ${SAMPLE} ${IT} ~~~~"
		fi # exit final hic plotting check 

	fi # exit entire loop file assembly check 
done 