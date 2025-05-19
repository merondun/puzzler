#!/bin/bash

#SBATCH --time=48:00:00   
#SBATCH --nodes=1  
#SBATCH --cpus-per-task=10
#SBATCH --mem=64Gb
#SBATCH --partition=ceres

# Submit as e.g. cat Runs.list | xargs -I {} sbatch -J {} 02_RefineCluster.sh {} 
if [ -z "$1" ]; then
    echo "Error: File prefix positional argument is required."
    exit 1
fi

module load miniconda
source activate isoseq

# trials with m84125_250429_221012_s3.hifi_reads.bcM0001.skera.IsoSeqX_bc01_5p--IsoSeqX_3p
FILE=$1
echo -e "\e[43m~~~~ Refining & clustering file: ${FILE} ~~~~\e[0m"
isoseq refine --require-polya ${FILE}.bam Barcodes.fa ${FILE}.Barcodes.flnc.bam
