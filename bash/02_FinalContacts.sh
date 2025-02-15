#!/bin/bash

#SBATCH --time=36:00:00   
#SBATCH --nodes=1  
#SBATCH --ntasks-per-node=32
#SBATCH --mem=64Gb
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
	if [ ! -s pg_renamed.filtered.bam ] && [ -s pg_renamed.fa ]; then

		echo "~~~~ Creating bam for ${IT} and ${SAMPLE} ~~~~"

		$PUZZLER bwa index pg_renamed.fa
		$PUZZLER samtools faidx pg_renamed.fa

		# Align Hi-C reads
		$PUZZLER bwa mem -5SP -t ${t} pg_renamed.fa ${HIC_R1} ${HIC_R2} | \
			$PUZZLER samblaster | \
			$PUZZLER samtools view - -@ ${t} -S -h -b -F 3340 | \
			$PUZZLER filter_bam - 1 --nm 3 --threads ${t} --remove-dup | \
			$PUZZLER samtools view - -b -@ ${t} -o pg_renamed.filtered.bam

	else
		echo "~~~~ Bam already created for ${IT} and ${SAMPLE} ~~~~"
	fi

	if [ ! -f ${WD}/logs/contact_maps/${SAMPLE}.${IT}.pdf ]; then

		~/symlinks/software/HapHiC/utils/mock_agp_file.py pg_renamed.fa > pg_renamed.agp
		$PUZZLER haphic plot --threads ${t} pg_renamed.agp pg_renamed.filtered.bam --bin_size 1000 --min_len 2

		# Copy over the file with the correct pangenome-spec naming
		cp haphic_renamed.fa ${WD}/joint_scaffold/${SAMPLE}.${IT}.fa
		$PUZZLER samtools faidx ${WD}/joint_scaffold/${SAMPLE}.${IT}.fa
		cp contact_map.pdf ${WD}/logs/contact_maps/${SAMPLE}.${IT}.pdf
	else
		echo "~~~~ Final contact map already generated ~~~~"
	fi

done