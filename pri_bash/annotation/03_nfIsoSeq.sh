#!/bin/bash

#SBATCH --time=8-00:00:00   
#SBATCH --nodes=1  
#SBATCH --cpus-per-task=12
#SBATCH --mem=64Gb
#SBATCH --partition=ceres

module load miniconda
source activate isoseq_ann
module load apptainer

# Submit as e.g. cat Runs.list | xargs -I {} sbatch -J {} 02_RefineCluster.sh {} 
if [ -z "$1" ]; then
    echo "Error: Sample positional argument is required."
    exit 1
fi

SAMPLE=$1
OUTDIR=/project/coffea_pangenome/Artocarpus/Assemblies/20250101_JustinAssemblies/annotation/isoseq
GENOME=/project/coffea_pangenome/Artocarpus/Assemblies/20250101_JustinAssemblies/joint_scaffold/EarlGrey/softmasked/${SAMPLE}.pri.softmasked.fasta

nextflow run nf-core/isoseq \
   --input /project/coffea_pangenome/Artocarpus/RawData/transfer/Breadfruit_BamFiles/skera/${SAMPLE}_samplesheet.csv \
   --outdir ${OUTDIR} \
   --fasta ${GENOME} \
   --primers Barcodes.fa \
   --aligner 'minimap2' \
   --entrypoint 'map' \
   -profile apptainer
