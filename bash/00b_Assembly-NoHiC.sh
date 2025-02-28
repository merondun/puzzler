#!/bin/bash

#SBATCH --time=24:00:00   
#SBATCH --nodes=1  
#SBATCH --ntasks-per-node=32
#SBATCH --mem=256Gb
#SBATCH --partition=ceres

# use n=32 and mem = 112 typically 
SINGULARITY_TMPDIR=${APPTAINER_TMPDIR}

module load miniconda
source activate puzz
module load apptainer

if [ -z "$1" ]; then
    echo "Error: sample ID positinonal argument is required."
    exit 1
fi

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

mkdir -p ${WD}/${SAMPLE} ${WD}/logs/juicer ${WD}/logs/haphic
cd ${WD}/${SAMPLE}

##### Generate Primary Assembly #####, add --hom-cov ${HOM_COV} \ if necessary! 
if [ ! -s ${SAMPLE}.p_ctg.gfa ]; then
	${PUZZLER} hifiasm --n-hap ${PLOIDY} --hom-cov ${HOM_COV} --primary -t ${t} -o ${SAMPLE} ${HIFI} 2> ${WD}/logs/${SAMPLE}.hifiasm.log
else
	echo "~~~~ Skipping hifiasm for ${SAMPLE}, already exists ~~~~"
fi 

IT="pri"
mkdir -p ${WD}/${SAMPLE}/02_${IT}HapHiC
cd ${WD}/${SAMPLE}/02_${IT}HapHiC

##### PURGE #####
if [ ! -s all.purged.fa ]; then

	echo "~~~~ Starting Purge_Dups for ${SAMPLE}, primary assembly ~~~~"

	#Ensure fastas have proper contig haplotype names 
	awk '/^S/{print ">"$2;print $3}' ${WD}/${SAMPLE}/${SAMPLE}.p_ctg.gfa > pri.init.fa
	${PUZZLER} split_fa pri.init.fa > pri.split.fa
	${PUZZLER} minimap2 -t ${t} -xasm5 -DP pri.split.fa pri.split.fa | gzip -c > pri.split.self.paf.gz
	${PUZZLER} purge_dups -M1000 -E1000 pri.split.self.paf.gz > pri.dups.bed
	${PUZZLER} get_seqs pri.dups.bed pri.init.fa
	mv purged.fa all.purged.fa
	rm *split.fa *self.paf.gz *dups.bed *init*

else
	echo "~~~~ Skipping purge for ${SAMPLE}, ${IT} ~~~~"
fi

cd ${WD}/${SAMPLE}/02_${IT}HapHiC

##### MAP SCAFFOLDS TO REFERENCE #####
if [ ! -s map.txt ]; then 
	echo "~~~~ Mapping reference chromosomes for ${SAMPLE} ~~~~"
	$PUZZLER minimap2 -x asm20 ${REFERENCE} all.purged.fa --secondary=no -t ${t} -o asm_ref.paf 2> asm_ref.log
	${PUZZLER} samtools faidx all.purged.fa
	map_chromosomes --paf asm_ref.paf --fai all.purged.fa.fai --out map.txt
else
	echo "~~~~ Skipping Mapping reference chromosomes for ${SAMPLE} ~~~~"
fi 

##### RENAME SCAFFOLDS TO REFERENCE CHRS #####
if [ ! -s haphic_renamed.fa ]; then
	mkdir -p 02_orienting
	echo "~~~~ Renaming chromosomes for ${SAMPLE} ~~~~"

	# This file has $haphic_scafID \t $ref_chr_ID \t $strand
	awk '{OFS="\t"}{print $1, $2, $5, $3, $4}' map.txt > 02_orienting/hap_newID_map.txt

	# Exclude those chromsomes from map, and then add e.g. h3tg0120320 for the remaining scaffolds 
	awk '{print $1}' 02_orienting/hap_newID_map.txt > 02_orienting/exclude_chr_scafIDs.txt
	grep -vwf 02_orienting/exclude_chr_scafIDs.txt all.purged.fa.fai | awk '{OFS="\t"}{print $1, $1, "+"}' >> 02_orienting/hap_newID_map.txt

	# Must be the same number of scaffolds
	if [ $(cat 02_orienting/hap_newID_map.txt | wc -l) -eq $(cat all.purged.fa.fai | wc -l) ]; then
		echo "Same number of scaffolds, proceeding!"
	else
		echo "Not same number, stop "
		#exit 1
	fi

	# Grab the reverse strands to reverse complement them 
	awk '{if ($3 == "-") print $1}' 02_orienting/hap_newID_map.txt > 02_orienting/reverse_list.txt
	awk '{if ($3 == "+") print $1}' 02_orienting/hap_newID_map.txt > 02_orienting/positive_list.txt
	seqtk subseq all.purged.fa 02_orienting/reverse_list.txt > 02_orienting/to_revercomp.fa
	seqtk subseq all.purged.fa 02_orienting/positive_list.txt > 02_orienting/positive.fa
	seqkit seq --line-width 0 -t DNA -v -r -p 02_orienting/to_revercomp.fa > 02_orienting/revcomp.fa
	cat 02_orienting/positive.fa 02_orienting/revcomp.fa > 02_orienting/orient.fa

	# Detect if there are duplicate scaffolds
	awk '{print $2}' 02_orienting/hap_newID_map.txt | sort | uniq -d > 02_orienting/duplicates.txt

	##### THIS WILL DEAL WITH DUPLICATES, IN CASE THERE ARE 2 SCAFFOLDS CORRESPONDING TO A CHR FROM THE SAME HAPLOTYPE! 
	if [ -s 02_orienting/duplicates.txt ]; then
		
		echo "Multiple scaffolds corresponding to a single Chr for ${SAMPLE}, INSPECT!" > ${WD}/logs/${SAMPLE}.${IT}.duplicates.log
	else 

		echo "Single scaffolds corresponding to a single Chr for ${SAMPLE}"

		# Rename chrs and sort 
		seqkit replace --line-width 0 -p "(.*)" -r "{kv}" -k 02_orienting/hap_newID_map.txt 02_orienting/orient.fa | \
			seqkit sort --line-width 0 -n > haphic_renamed.fa

	fi

else
	echo "~~~~ Skipping renaming chromosomes for ${SAMPLE}, already exists  ~~~~"
fi