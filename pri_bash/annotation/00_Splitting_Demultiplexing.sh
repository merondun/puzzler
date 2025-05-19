#!/bin/bash

#SBATCH --time=48:00:00   
#SBATCH --nodes=1  
#SBATCH --cpus-per-task=10
#SBATCH --partition=ceres

module load miniconda
source activate isoseq

if [ -z "$1" ]; then
    echo "Error: Library prefix file positional argument is required."
    exit 1
fi

LIBRARY=$1

echo -e "\e[43m~~~~ Demultiplexing isoseq file: ${LIBRARY} ~~~~\e[0m"
skera split ${LIBRARY}.bam mas16_primers.fasta ${LIBRARY}.skera.bam
lima ${LIBRARY}.skera.bam Barcodes.fa skera/${LIBRARY}.skera.bam --isoseq --overwrite-biosample-names
