#!/bin/bash

#SBATCH --time=336:00:00   
#SBATCH --nodes=1  
#SBATCH --ntasks-per-node=20
#SBATCH --mem=112Gb
#SBATCH --partition=ceres

#module load miniconda
#source activate earlgrey

SAMPLE=$1
#WD=/project/coffea_pangenome/Artocarpus/Assemblies/20250101_JustinAssemblies/finalized_assemblies
WD=/project/coffea_pangenome/Artocarpus/Assemblies/20250101_JustinAssemblies/joint_scaffold
FASTA="${WD}/${SAMPLE}.pri.fa" 

mkdir -p ${WD}/EarlGrey
cd ${WD}/EarlGrey

samtools faidx $FASTA
earlGrey -d yes -e yes -t 20 -g ${FASTA} -s ${SAMPLE}.pri -o ${WD}/EarlGrey

# HAPLOTYPES
FASTA_HAP="${WD}/${SAMPLE}.hap.fa" 

samtools faidx $FASTA_HAP
earlGrey -d yes -e yes -t 20 -g ${FASTA_HAP} -s ${SAMPLE}.hap -o ${WD}/EarlGrey
