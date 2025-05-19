#!/bin/bash

#SBATCH --time=300:00:00   
#SBATCH --nodes=1  
#SBATCH --ntasks-per-node=30
#SBATCH --mem=512Gb
#SBATCH --partition=ceres

SINGULARITY_TMPDIR=${APPTAINER_TMPDIR}

module load miniconda
source activate getorganelle

if [ -z "$1" ]; then
    echo "Error: sample ID positional argument is required."
    exit 1
fi

t=30
MEM=512
SAMPLE=$1
SPECIES="Guava"

WD=/project/coffea_pangenome/${SPECIES}/Assemblies/20250101_JustinAssemblies
SAMPLE_FILE=/project/coffea_pangenome/${SPECIES}/Assemblies/20250101_JustinAssemblies/samples.csv
HIFI=$(awk -F',' -v sample="$SAMPLE" '$1 == sample {print $5}' ${SAMPLE_FILE})

mkdir -p ${WD}/organelle
cd ${WD}/organelle

# Run this first 
#get_organelle_config.py --add embplant_pt,embplant_mt

get_organelle_from_reads.py -u ${HIFI} -o ${SAMPLE}_2e5 --continue -t ${t} --max-reads 2e5 -R 30 -k 21,45,65,85,105 -F embplant_pt
