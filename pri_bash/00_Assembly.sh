#!/bin/bash

#SBATCH --time=48:00:00   
#SBATCH --nodes=1  
#SBATCH --ntasks-per-node=32
#SBATCH --mem=256Gb
#SBATCH --partition=ceres

# use n=32 and mem = 256 typically 
SINGULARITY_TMPDIR=${APPTAINER_TMPDIR}

t=32
MEM=256

# Default values
SAMPLE=""
MAP_FILE=""

# Parse command-line arguments
while [[ "$#" -gt 0 ]]; do
    case "$1" in
        --sample) SAMPLE="$2"; shift ;;
        --map_file) MAP_FILE="$2"; shift ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
    shift
done

# Validate arguments
if [ -z "$SAMPLE" ]; then
    echo "Error: --sample argument is required."
    exit 1
fi

if [ -z "$MAP_FILE" ]; then
    echo "Error: --map_file argument is required."
    exit 1
fi

module load miniconda
source activate puzz
module load apptainer

# Read CSV line matching sample and assign fields
IFS=',' read -r _ SIF_PATH WD PLOIDY NUM_CHRS HOM_COV HIFI HIC_R1 HIC_R2 REFERENCE < <(awk -F',' -v sample="$SAMPLE" '$1 == sample {print $0}' "${MAP_FILE}")

PUZZLER="apptainer exec ${SIF_PATH}"

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

echo -e "=======================================================================\nParameters for sample: ${SAMPLE} \nCONTAINER: ${SIF_PATH} \nWD: ${WD} \nPLOIDY: ${PLOIDY} \nHOM_COV: ${HOM_COV}\nHIFI: ${HIFI}\nHIC_R1: ${HIC_R1}\nHIC_R2: ${HIC_R2}\nREFERENCE: ${REFERENCE}\n=======================================================================\n"

mkdir -p ${WD}/${SAMPLE} ${WD}/logs/juicer ${WD}/logs/haphic
cd ${WD}/${SAMPLE}

##### Generate Primary Assembly #####
if [ ! -s ${SAMPLE}.hic.p_ctg.gfa ]; then
	echo -e "\e[43m~~~~ Starting hifiasm assembly for ${SAMPLE} ~~~~\e[0m"
	${PUZZLER} hifiasm --n-hap ${PLOIDY} --hom-cov ${HOM_COV} --primary -t ${t} -o ${SAMPLE} --h1 ${HIC_R1} --h2 ${HIC_R2} ${HIFI} 2> ${WD}/logs/${SAMPLE}.hifiasm.log
else
	echo -e "\e[42m~~~~ Skipping hifiasm for ${SAMPLE}, already exists ~~~~\e[0m"
fi 

##### Loop through both primary assembly pipeline and haplotype-joint pipeline #####


mkdir -p ${WD}/${SAMPLE}/01_scaffolding
cd ${WD}/${SAMPLE}/01_scaffolding

##### PURGE #####
if [ ! -s all.purged.fa ]; then

	echo -e "\e[43m~~~~ Starting Purge_Dups for ${SAMPLE}, primary assembly ~~~~\e[0m"

	#Ensure fastas have proper contig haplotype names 
	awk '/^S/{print ">"$2;print $3}' ${WD}/${SAMPLE}/${SAMPLE}.hic.p_ctg.gfa > pri.init.fa
	${PUZZLER} split_fa pri.init.fa > pri.split.fa
	${PUZZLER} minimap2 -t ${t} -xasm5 -DP pri.split.fa pri.split.fa 2> minimap.purge.log | gzip -c > pri.split.self.paf.gz
	${PUZZLER} purge_dups -M1000 -E1000 pri.split.self.paf.gz > pri.dups.bed 2> pri.purge.log
	${PUZZLER} get_seqs pri.dups.bed pri.init.fa 2> pri.getseqs.log
	mv purged.fa all.purged.fa
	rm *split.fa *self.paf.gz *dups.bed *init*

else
	echo -e "\e[42m~~~~ Skipping purge for ${SAMPLE}, already exists ~~~~\e[0m"
fi

##### ALIGN HIC #####
cd ${WD}/${SAMPLE}/01_scaffolding
if [ ! -s filtered.MQ1.bam ]; then
	echo -e "\e[43m~~~~ Mapping HiC reads to ${SAMPLE} pri assembly ~~~~\e[0m" 
	# Index reference
	${PUZZLER} bwa index all.purged.fa 2> alignment.indexing.hic.log

	# Align Hi-C reads
	${PUZZLER} bwa mem -5SP -t ${t} all.purged.fa ${HIC_R1} ${HIC_R2} 2> alignment.firsthic.log | \
		${PUZZLER} samblaster 2> samblaster.firsthic.log | ${PUZZLER} samtools view - -@ ${t} -S -h -b -F 3340 | \
		${PUZZLER} filter_bam - 1 --nm 3 --threads ${t} | \
		${PUZZLER} samtools view - -b -@ ${t} -o filtered.MQ1.bam

else
	echo -e "\e[42m~~~~ Skipping HiC alignment for ${SAMPLE}, already exists ~~~~\e[0m"
fi

##### HAPHIC #####
cd ${WD}/${SAMPLE}/01_scaffolding
if [ ! -s ${WD}/${SAMPLE}/01_scaffolding/haphic/04.build/scaffolds.fa ]; then

	echo -e "\e[43m~~~~ Running HapHiC for ${SAMPLE}  ~~~~\e[0m" 
	mkdir -p ${WD}/${SAMPLE}/01_scaffolding/haphic
	rm -rf ${WD}/${SAMPLE}/01_scaffolding/haphic/*

	cd ${WD}/${SAMPLE}/01_scaffolding/haphic

	${PUZZLER} haphic pipeline ../all.purged.fa ../filtered.MQ1.bam ${NUM_CHRS} --remove_allelic_links ${PLOIDY} --correct_nrounds 2 --max_inflation 20.0 --threads ${t} --processes ${t} 2> ${WD}/logs/haphic/${SAMPLE}.haphic.log

	cp ${WD}/${SAMPLE}/01_scaffolding/haphic/04.build/scaffolds.fa ${WD}/${SAMPLE}/01_scaffolding/haphic.fa
	cp ${WD}/${SAMPLE}/01_scaffolding/haphic/01.cluster/HapHiC_cluster.log ${WD}/${SAMPLE}/01_scaffolding/haphic.cluster.log

else 
	echo -e "\e[42m~~~~ Skipping HapHiC for ${SAMPLE}, already exists ~~~~\e[0m"
fi

##### JUICER #####
cd ${WD}/${SAMPLE}/01_scaffolding
mkdir -p 02_orienting
if [ ! -s ${WD}/logs/juicer/${SAMPLE}_JBAT.hic ]; then

	echo -e "\e[43m~~~~ Creating .hic file for juicebox for ${SAMPLE}  ~~~~\e[0m"
	# # Orient to jackfruit chromosomes 
	if [ ! -s 02_orienting/asm_to_ref.paf ]; then
		${PUZZLER} minimap2 minimap2 -x asm20 ${REFERENCE} all.purged.fa --secondary=no -t ${t} -o 02_orienting/asm_to_ref.paf 2> minimap.ref.log
	else
		echo -e "\e[42m~~~~ Skipping initial alignment between reference and draft ~~~~\e[0m]"
	fi

	${PUZZLER} minimap2 haphic refsort ${WD}/${SAMPLE}/01_scaffolding/haphic/04.build/scaffolds.raw.agp 02_orienting/asm_to_ref.paf > refsort.agp 2> refsort.log
	${PUZZLER} minimap2 samtools faidx haphic.fa
	${PUZZLER} minimap2 samtools faidx all.purged.fa

	${PUZZLER} juicer pre \
		-a -q 1 \
		-o haphic-refsort_JBAT \
		filtered.MQ1.bam \
		refsort.agp \
		all.purged.fa.fai \
		> haphic-refsort_JBAT.log 2>&1

	# Extract chromosome sizes and create HiC file
	grep PRE_C_SIZE haphic-refsort_JBAT.log | \
		awk '{print $2" "$3}' > chrom.sizes
	${PUZZLER} java -Xmx${MEM}G -jar /opt/HapHiC/utils/juicer_tools.1.9.9_jcuda.0.8.jar pre \
		-r 5000000,4000000,3000000,2000000,1500000,1000000,750000,500000,250000,100000,50000 \
		haphic-refsort_JBAT.txt \
		haphic-refsort_JBAT.hic \
		chrom.sizes 2> juicer_pre.log
		
	cp haphic-refsort_JBAT.hic ${WD}/logs/juicer/${SAMPLE}_JBAT.hic
	cp haphic-refsort_JBAT.assembly ${WD}/logs/juicer/${SAMPLE}_JBAT.assembly
	rm haphic-refsort_JBAT.txt 
else
	echo -e "\e[42m~~~~ Skipping juicer HiC file creation for ${SAMPLE}, already exists ~~~~\e[0m"
fi 

#### MANUAL CURATION, THEN THROW THE ${SAMPLE}_JBAT..review.assembly BACK INTO ${WD}/${SAMPLE}/01_scaffolding
