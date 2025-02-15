#!/bin/bash

#SBATCH --time=24:00:00   
#SBATCH --nodes=1  
#SBATCH --ntasks-per-node=32
#SBATCH --mem=128Gb
#SBATCH --partition=short

# use n=32 and mem = 112 typically 
SINGULARITY_TMPDIR=${APPTAINER_TMPDIR}

module load miniconda
source activate puzz
module load apptainer

t=32
MEM=128
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

mkdir -p ${WD}/${SAMPLE} ${WD}/logs/chrjuicer ${WD}/logs/chrjuicer/chrcontact_maps ${WD}/logs/haphic

# It's worthwhile to examine divergence between haplotypes
cd ${WD}/${SAMPLE}/02_hapHapHiC
if [ ! -s chr.divergence.txt ]; then

	echo "~~~~ Estimating Haplotype divergence for ${SAMPLE} ~~~~"
	mkdir -p 00_hapchrs 01_chrassemblies 
	for HAP in $(seq 1 ${PLOIDY}); do 

		#$PUZZLER minimap2 -x asm20 ${REFERENCE} hap${HAP}.purged.fa --secondary=no -t ${t} -o 00_hapchrs/hap${HAP}_ref.paf 2> 00_hapchrs/hap${HAP}.log
		samtools faidx hap${HAP}.purged.fa
		map_chromosomes --paf 00_hapchrs/hap${HAP}_ref.paf --fai hap${HAP}.purged.fa.fai --out hap${HAP}_chrmap.txt

	done 

	# Now loop through chromosomes, extract the chr fasta and estimate divergence 
	for CHRNUM in $(seq 1 ${NCHRS}); do 
		for HAP in $(seq 1 ${PLOIDY}); do 

			CHR=$(printf "Chr%02d" ${CHRNUM})
			echo "Extracting ${CHR} for ${SAMPLE} hap${HAP}"
			scaffold=$(grep -w ${CHR} hap${HAP}_chrmap.txt | sort -k 4 | head -n1 | awk '{print $1}')
			samtools faidx hap${HAP}.purged.fa ${scaffold} | sed "s/^>.*/>${SAMPLE}#${HAP}#${CHR}/" > 00_hapchrs/hap${HAP}_${CHR}.fa 

		done

		cat 00_hapchrs/*${CHR}.fa > 01_chrassemblies/${CHR}.fa 
		${PUZZLER} samtools faidx 01_chrassemblies/${CHR}.fa

		# Calculate divergence, format 
		files=$(find 00_hapchrs -name "*_${CHR}.fa" | tr '\n' ' ')
		mash sketch -o 00_hapchrs/${CHR}.divergence $files
		mash dist 00_hapchrs/${CHR}.divergence.msh 00_hapchrs/${CHR}.divergence.msh | \
			sed "s@00_hapchrs/hap@@g" | sed 's/_/\t/g' | sed 's/.fa//g' | sed 's@\./@@g' | \
			awk '$1 != $3' | sort -k3 | awk '!seen[$1,$3]++ && !seen[$3,$1]++' | awk -v s=${SAMPLE} '{OFS="\t"}{print s, $2, $1, $3, $5, $7}'> 00_hapchrs/${CHR}_divergence.out

	done 

	cat 00_hapchrs/*divergence.out > chr.divergence.txt
	cp chr.divergence.txt ${WD}/logs/divergence_haps/${SAMPLE}.chr.divergence.out 
	awk '{OFS="\t"}{print $1, $2}' 01_chrassemblies/${CHR}.fa.fai > ${WD}/logs/divergence_haps/${SAMPLE}.chrlengths.txt

else 
	echo "~~~~ Skipping haplotype divergence estimation for ${SAMPLE} ~~~~"
fi

# Map reads to whole genome, and then subset.. 

for CHRNUM in $(seq 1 ${NCHRS}); do 
	CHR=$(printf "Chr%02d" ${CHRNUM})

	##### ALIGN HIC: HAPLOTYPES #####
	if [ ! -s 01_chrassemblies/${CHR}.bam ]; then
			echo "~~~~ Creating HiC bam for ${SAMPLE} and ${CHR} ~~~~"
			# Index reference
			${PUZZLER} bwa index 01_chrassemblies/${CHR}.fa
			${PUZZLER} mock_agp_file.py 01_chrassemblies/${CHR}.fa > 01_chrassemblies/${CHR}.agp

			# Align Hi-C reads
			${PUZZLER} bwa mem -5SP -t ${t} 01_chrassemblies/${CHR}.fa ${HIC_R1} ${HIC_R2} | \
				${PUZZLER} samblaster | ${PUZZLER} samtools view - -@ ${t} -S -h -b -F 3340 -o 01_chrassemblies/${CHR}.bam
	else 
		echo "~~~~ Skipping mapping HiC ~~~~"
	fi

	if [ ! -s 01_chrassemblies/${CHR}_JBAT.hic ]; then
			echo "~~~~ Creating juicer file for ${SAMPLE} and ${CHR} ~~~~"
			# Filter alignments
			${PUZZLER} filter_bam 01_chrassemblies/${CHR}.bam 1 --nm 3 --threads ${t} | \
				${PUZZLER} samtools view - -b -@ ${t} -o 01_chrassemblies/${CHR}.filtered.bam

			JUICER_PATH=~/symlinks/software/HapHiC/utils
			${JUICER_PATH}/juicer pre \
				-a -q 1 \
				-o 01_chrassemblies/${CHR}_JBAT \
				01_chrassemblies/${CHR}.filtered.bam \
				01_chrassemblies/${CHR}.agp \
				01_chrassemblies/${CHR}.fa.fai \
				> 01_chrassemblies/${CHR}_JBAT.log 2>&1

			# Extract chromosome sizes and create HiC file
			grep PRE_C_SIZE 01_chrassemblies/${CHR}_JBAT.log | \
				awk '{print $2" "$3}' > 01_chrassemblies/${CHR}_JBAT.chrom.sizes
			java -Xmx${MEM}G -jar ${JUICER_PATH}/juicer_tools.1.9.9_jcuda.0.8.jar pre \
				-r 5000000,4000000,3000000,2000000,1000000,500000,250000,100000,50000 \
				01_chrassemblies/${CHR}_JBAT.txt \
				01_chrassemblies/${CHR}_JBAT.hic \
				01_chrassemblies/${CHR}_JBAT.chrom.sizes

		cp 01_chrassemblies/${CHR}_JBAT.hic ${WD}/logs/chrjuicer/${SAMPLE}_${CHR}.hic
		cp 01_chrassemblies/${CHR}_JBAT.assembly ${WD}/logs/chrjuicer/${SAMPLE}_${CHR}.assembly

		$PUZZLER haphic plot --threads ${t} 01_chrassemblies/${CHR}.agp 01_chrassemblies/${CHR}.filtered.bam --bin_size 1000 --min_len 2
		cp contact_map.pdf ${WD}/logs/chrjuicer/chrcontact_maps/${SAMPLE}.${CHR}.pdf

	else 
		echo "~~~~ Skipping juicer file creation for ${SAMPLE} and ${CHR} ~~~~"
	fi

done
