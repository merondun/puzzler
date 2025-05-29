#!/bin/bash

#SBATCH --time=14-00:00:00   
#SBATCH --nodes=1  
#SBATCH --ntasks-per-node=16
#SBATCH --mem=72Gb
#SBATCH --partition=ceres

SINGULARITY_TMPDIR=${APPTAINER_TMPDIR}

module load miniconda
source activate puzz

if [ -z "$1" ]; then
    echo "Error: sample ID positional argument is required."
    exit 1
fi

t=24
MEM=72
SAMPLE=$1

WD=/project/coffea_pangenome/Artocarpus/Assemblies/20250101_JustinAssemblies
SAMPLE_FILE="/project/coffea_pangenome/Artocarpus/Assemblies/20250101_JustinAssemblies/samples.csv"
PUZZLER="apptainer exec /project/coffea_pangenome/Software/Merondun/apptainers/puzzler_v1.1.sif"
PLOIDY=$(awk -F',' -v sample="$SAMPLE" '$1 == sample {print $2}' ${SAMPLE_FILE})
HIFI=$(awk -F',' -v sample="$SAMPLE" '$1 == sample {print $5}' ${SAMPLE_FILE})
GS=$(awk -F',' -v sample="$SAMPLE" '$1 == sample {print $17}' ${SAMPLE_FILE})
GENOME=/project/coffea_pangenome/Artocarpus/Assemblies/20250101_JustinAssemblies/joint_scaffold/HART001.pri.fa
WD=/project/coffea_pangenome/Artocarpus/Assemblies/20250101_JustinAssemblies/panphase

cd ${WD}

mkdir -p bams vcfs subsamp

# Alignment 
if [ ! -s subsamp/${SAMPLE}.10X.fastq.gz ]; then
        TARGET=$(awk -v cov=10 -v gs=${GS} -v p=${PLOIDY} 'BEGIN { printf "%.0f", cov * gs * 1000000 * p } ')
        echo -e "\e[43m~~~~ Subsampling ${TARGET} bases for ${SAMPLE} ~~~~\e[0m"
        bbduk.sh in=${HIFI} out=subsamp/${SAMPLE}.10X.fastq.gz maxbasesout=${TARGET}
else
        echo -e "\e[42m~~~~ Skipping Alignment for ${SAMPLE} ~~~~\e[0m"
fi 

# Alignment 
if [ ! -s bams/${SAMPLE}.hifi.sorted.bam ]; then
        echo -e "\e[43m~~~~ Starting Alignment ${SAMPLE} ~~~~\e[0m"
        minimap2 -t ${t} -ax map-hifi ${GENOME} subsamp/${SAMPLE}.10X.fastq.gz 2> logs/${SAMPLE}.minimap.log | samtools sort --threads ${t} -o bams/${SAMPLE}.hifi.sorted.bam
        samtools index --threads ${t} bams/${SAMPLE}.hifi.sorted.bam
else
        echo -e "\e[42m~~~~ Skipping Alignment for ${SAMPLE} ~~~~\e[0m"
fi 

for CHR in Chr01; do 

        # Call SNP variants 
        MIN_ALT_COUNT=$(( (PLOIDY + 1 ) / 2 ))
        MIN_ALT_FRAC=$( awk -v p=${PLOIDY} 'BEGIN { printf "%.2f", 1 / (p + 1) }' )
        if [ ! -s vcfs/${CHR}_${SAMPLE}.SNP.gvcf.gz ]; then
                echo -e "\e[43m~~~~ Starting SNPs for ${SAMPLE} ${CHR} ~~~~\e[0m"
                freebayes-parallel <(fasta_generate_regions.py ${GENOME}.fai 100000 | grep ${CHR}) ${t} --report-monomorphic --ploidy ${PLOIDY} --skip-coverage 30 --min-coverage 5 --min-alternate-count ${MIN_ALT_COUNT} --min-alternate-fraction ${MIN_ALT_FRAC} -f ${GENOME} ${WD}/bams/${SAMPLE}.hifi.sorted.bam > vcfs/${CHR}_${SAMPLE}.SNP.gvcf
                bcftools view --exclude-types indels vcfs/${CHR}_${SAMPLE}.SNP.gvcf -Oz -o vcfs/${CHR}_${SAMPLE}.SNP.gvcf.gz
                tabix vcfs/${CHR}_${SAMPLE}.SNP.gvcf.gz
        else
                echo -e "\e[42m~~~~ Skipping SNPs for ${SAMPLE} ${CHR} ~~~~\e[0m"
        fi 

        # Phase with whatshap polyphase
        if [ ! -s ${WD}/vcfs/${CHR}_${SAMPLE}.SNP.phased.named.gvcf.gz ]; then
                echo -e "\e[43m~~~~ Phasing ${SAMPLE} ${CHR} ~~~~\e[0m"
                whatshap polyphase -t ${t} \
                        --ploidy $PLOIDY \
                        --distrust-genotypes \
                        --ignore-read-groups \
                        --reference ${GENOME} \
                        --chromosome ${CHR} \
                        -o ${WD}/vcfs/${CHR}_${SAMPLE}.SNP.phased.gvcf.gz \
                        ${WD}/vcfs/${CHR}_${SAMPLE}.SNP.gvcf.gz bams/${SAMPLE}.hifi.sorted.bam
               bcftools reheader -s <(echo "${SAMPLE}") -o ${WD}/vcfs/${CHR}_${SAMPLE}.SNP.phased.named.gvcf.gz ${WD}/vcfs/${CHR}_${SAMPLE}.SNP.phased.gvcf.gz
               bcftools index ${WD}/vcfs/${CHR}_${SAMPLE}.SNP.phased.named.gvcf.gz

        else
                echo -e "\e[42m~~~~ Skipping phasing for ${SAMPLE} ${CHR} ~~~~\e[0m"
        fi 

done 
