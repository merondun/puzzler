#!/bin/bash

#SBATCH --time=100:00:00   
#SBATCH --nodes=1  
#SBATCH --ntasks-per-node=8 
#SBATCH --mem=60Gb
#SBATCH --partition=medium  

if [ -z "$1" ]; then
    echo "Error: sample ID positional argument is required."
    exit 1
fi

t=8
MEM=60
SAMPLE=$1
SPECIES="Artocarpus"

SAMPLE_FILE=/project/coffea_pangenome/${SPECIES}/Assemblies/20250101_JustinAssemblies/samples.csv
HIFI=$(awk -F',' -v sample="$SAMPLE" '$1 == sample {print $5}' ${SAMPLE_FILE})
WD=/project/coffea_pangenome/${SPECIES}/Assemblies/20250101_JustinAssemblies/organelle/reference_based
GENOME=/project/coffea_pangenome/${SPECIES}/Assemblies/20250101_JustinAssemblies/organelle/reference_based/HART001.1.fasta
BAMDIR=/project/coffea_pangenome/${SPECIES}/Assemblies/20250101_JustinAssemblies/organelle/reference_based/bams
READDIR=/project/coffea_pangenome/${SPECIES}/Assemblies/20250101_JustinAssemblies/organelle/reference_based/5gb_reads
SNPDIR=/project/coffea_pangenome/${SPECIES}/Assemblies/20250101_JustinAssemblies/organelle/reference_based/vcfs

module load apptainer
source activate assembly

echo "Mapping SAMPLE: ${SAMPLE}"

# Subset reads so that we have roughly equal inputs 
bbduk.sh in=${HIFI} out=${READDIR}/${SAMPLE}.5gb.fastq.gz maxbasesout=5000000000

# Align, with more relaxed parameters since jackfruit is ~20 MY diverged
minimap2 -ax map-hifi -t ${t} -R @RG\\tID:${SAMPLE}\\tPL:PACBIO\\tLB:${SAMPLE}\\tSM:${SAMPLE} ${GENOME} ${READDIR}/${SAMPLE}.5gb.fastq.gz 2> ${BAMDIR}/${SAMPLE}.minimap.log | \
    samtools view -F 4 -bS - | \
    samtools sort -@ ${t} -o ${BAMDIR}/${SAMPLE}.sorted.bam
samtools index ${BAMDIR}/${SAMPLE}.sorted.bam

# Call variants 
DEEPVAR="/project/coffea_pangenome/Software/Merondun/apptainers/deepvariant_1.5.0.sif"
apptainer exec ${DEEPVAR} run_deepvariant \
    --model_type PACBIO \
    --ref ${GENOME} \
    --reads ${BAMDIR}/${SAMPLE}.sorted.bam \
    --output_vcf ${SNPDIR}/${SAMPLE}.pt.vcf.gz \
    --output_gvcf ${SNPDIR}/${SAMPLE}.pt.gvcf.gz \
    --num_shards ${t}
