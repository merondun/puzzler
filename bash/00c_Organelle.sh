#!/bin/bash

#SBATCH --time=300:00:00   
#SBATCH --nodes=1  
#SBATCH --ntasks-per-node=20
#SBATCH --mem=128Gb
#SBATCH --partition=ceres

SINGULARITY_TMPDIR=${APPTAINER_TMPDIR}

module load miniconda
source activate getorganelle

PMAT="/project/coffea_pangenome/Artocarpus/Assemblies/20250101_JustinAssemblies/organelle/PMAT/bin/PMAT"

if [ -z "$1" ]; then
    echo "Error: sample ID positional argument is required."
    exit 1
fi

t=20
MEM=128
SAMPLE=$1

WD=/project/coffea_pangenome/Artocarpus/Assemblies/20250101_JustinAssemblies
SAMPLE_FILE="/project/coffea_pangenome/Artocarpus/Assemblies/20250101_JustinAssemblies/samples.csv"
HIFI=$(awk -F',' -v sample="$SAMPLE" '$1 == sample {print $5}' ${SAMPLE_FILE})
FC=$(awk -F',' -v sample="$SAMPLE" '$1 == sample {print $12}' ${SAMPLE_FILE} | awk '{print int($1+0.5)}')

mkdir -p ${WD}/organelle/pmat
cd ${WD}/organelle/pmat

${PMAT} autoMito -i ${HIFI} -o ${SAMPLE} -st hifi -g ${FC}m -m -cpu ${t} --type pt