#!/bin/bash

#SBATCH --time=24:00:00   
#SBATCH --nodes=1  
#SBATCH --ntasks-per-node=32
#SBATCH --mem=256Gb
#SBATCH --partition=short

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
MEM=256
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

mkdir -p ${WD}/${SAMPLE} ${WD}/logs/juicer ${WD}/logs/haphic ${WD}/logs/divergence_haps
cd ${WD}/${SAMPLE}

##### Generate Primary Assembly #####, add --hom-cov ${HOM_COV} \ if necessary! 
if [ ! -s ${SAMPLE}.hic.hap1.p_ctg.gfa ]; then
	${PUZZLER} hifiasm --n-hap ${PLOIDY} --hom-cov ${HOM_COV} --primary -t ${t} -o ${SAMPLE} --h1 ${HIC_R1} --h2 ${HIC_R2} ${HIFI} 2> ${WD}/logs/${SAMPLE}.hifiasm.log
else
	echo "~~~~ Skipping hifiasm for ${SAMPLE}, already exists ~~~~"
fi 

##### Loop through both primary assembly pipeline and haplotype-joint pipeline #####

for IT in pri; do 

	mkdir -p ${WD}/${SAMPLE}/02_${IT}HapHiC
	cd ${WD}/${SAMPLE}/02_${IT}HapHiC

	##### ALIGN HIC #####
	cd ${WD}/${SAMPLE}/02_${IT}HapHiC
	if [ ! -s filtered.MQ${MQ}.bam ]; then
		echo "~~~~ Mapping HiC reads to ${SAMPLE} ${IT} assembly ~~~~" 
		# Align Hi-C reads
		${PUZZLER} filter_bam filtered.MQ1.bam 1 --nm 3 --threads ${t} | \
			${PUZZLER} samtools view - -b -@ ${t} -o filtered.MQ${MQ}.bam

	else
		echo "~~~~ Skipping alignment for ${SAMPLE} ${IT} ~~~~"
	fi

	##### HAPHIC ##### 
	cd ${WD}/${SAMPLE}/02_${IT}HapHiC
	if [ ! -s ${WD}/${SAMPLE}/02_${IT}HapHiC/01_haphicMQ${MQ}/04.build/scaffolds.fa ]; then

		echo "~~~~ Running HapHiC for ${SAMPLE} ${IT} ~~~~" 
		mkdir -p ${WD}/${SAMPLE}/02_${IT}HapHiC/01_haphicMQ${MQ}
		rm -rf ${WD}/${SAMPLE}/02_${IT}HapHiC/01_haphicMQ${MQ}/*

		cd ${WD}/${SAMPLE}/02_${IT}HapHiC/01_haphicMQ${MQ}

		if [ "${IT}" = "hap" ]; then
			NUM_CHRS=${TOTAL_CHRS}
		else
			NUM_CHRS=${NCHRS}
		fi

		${PUZZLER} haphic pipeline ../all.purged.fa ../filtered.MQ${MQ}.bam ${NUM_CHRS} --remove_allelic_links ${PLOIDY} --correct_nrounds 2 --max_inflation 20.0 --threads ${t} --processes ${t} 2> ${WD}/logs/haphic/${SAMPLE}.${IT}.haphic.MQ${MQ}.log

		cp ${WD}/${SAMPLE}/02_${IT}HapHiC/01_haphicMQ${MQ}/04.build/scaffolds.fa ${WD}/${SAMPLE}/02_${IT}HapHiC/haphic.MQ${MQ}.fa
		cp ${WD}/${SAMPLE}/02_${IT}HapHiC/01_haphicMQ${MQ}/01.cluster/HapHiC_cluster.log ${WD}/${SAMPLE}/02_${IT}HapHiC/haphic.MQ${MQ}.log

	else 
		echo "~~~~ Skipping HapHiC for ${SAMPLE} ${IT} ~~~~"
	fi

	##### JUICER #####
	cd ${WD}/${SAMPLE}/02_${IT}HapHiC
	mkdir -p 02_orienting
	if [ ! -s ${WD}/logs/juicer/${SAMPLE}.${IT}-MQ${MQ}_JBAT.hic ]; then

		echo "~~~~ Creating .hic file for juicebox for ${SAMPLE} ${IT} ~~~~"
		# # Orient to jackfruit chromosomes 
		if [ ! -s 02_orienting/asm_to_ref.paf ]; then
			$PUZZLER minimap2 -x asm20 ${REFERENCE} all.purged.fa --secondary=no -t ${t} -o 02_orienting/asm_to_ref.paf
		else
			echo "Skipping initial alignment"
		fi

		$PUZZLER haphic refsort ${WD}/${SAMPLE}/02_${IT}HapHiC/01_haphicMQ${MQ}/04.build/scaffolds.raw.agp 02_orienting/asm_to_ref.paf > refsortMQ${MQ}.agp 2> refsortMQ${MQ}.log
		$PUZZLER samtools faidx haphic.MQ${MQ}.fa
		$PUZZLER samtools faidx all.purged.fa

		JUICER_PATH=~/symlinks/software/HapHiC/utils
		${JUICER_PATH}/juicer pre \
			-a -q 1 \
			-o haphic-refsort-MQ${MQ}_JBAT \
			filtered.MQ${MQ}.bam \
			refsortMQ${MQ}.agp \
			all.purged.fa.fai \
			> haphic-refsort-MQ${MQ}_JBAT.log 2>&1

		# Extract chromosome sizes and create HiC file
		grep PRE_C_SIZE haphic-refsort-MQ${MQ}_JBAT.log | \
			awk '{print $2" "$3}' > chrom.sizesMQ${MQ}
		java -Xmx${MEM}G -jar ${JUICER_PATH}/juicer_tools.1.9.9_jcuda.0.8.jar pre \
			-r 5000000,4000000,3000000,2000000,1500000,1000000,750000,500000,250000,100000,50000 \
			haphic-refsort-MQ${MQ}_JBAT.txt \
			haphic-refsort-MQ${MQ}_JBAT.hic \
			chrom.sizesMQ${MQ}
			
		cp haphic-refsort-MQ${MQ}_JBAT.hic ${WD}/logs/juicer/${SAMPLE}.${IT}-MQ${MQ}_JBAT.hic
		cp haphic-refsort-MQ${MQ}_JBAT.assembly ${WD}/logs/juicer/${SAMPLE}.${IT}-MQ${MQ}_JBAT.assembly
	
	else
		echo "~~~~ No HiC juicer files needed for ${SAMPLE} ${IT} ~~~~"
	fi 

done # terminate IT loop

#### MANUAL CURATION, THEN THROW THE ${SAMPLE}_JBAT.${IT}.review.assembly BACK INTO ${WD}/${SAMPLE}/02_${IT}HapHiC
