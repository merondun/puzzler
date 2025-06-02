#!/bin/bash

#SBATCH --time=14-00:00:00   
#SBATCH --nodes=1  
#SBATCH --ntasks-per-node=20
#SBATCH --mem=112Gb
#SBATCH --partition=ceres
#SBATCH --account=coffea_pangenome

#module load miniconda
#source activate earlgrey

SAMPLE=$1
WD=/project/coffea_pangenome/Artocarpus/Assemblies/20250101_JustinAssemblies/joint_scaffold

mkdir -p ${WD}/EarlGrey
cd ${WD}/EarlGrey

for IT in pri; do 

    FASTA="${WD}/${SAMPLE}.${IT}.fa" 
    if [ ! -s ${SAMPLE}.${IT}_EarlGrey/${SAMPLE}.${IT}_summaryFiles/${SAMPLE}.${IT}.softmasked.fasta ]  && [ -s ${FASTA} ]; then
        echo -e "\e[43m~~~~ Starting repeat annotation for ${SAMPLE} ${IT} ~~~~\e[0m"
        samtools faidx $FASTA
        earlGrey -d yes -e yes -t 20 -g ${FASTA} -s ${SAMPLE}.${IT} -o ${WD}/EarlGrey

    else
        echo -e "\e[42m~~~~ Skipping repeat annotation for ${SAMPLE} ${IT}, already exists ~~~~\e[0m"
    fi 

done
