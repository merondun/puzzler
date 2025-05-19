#!/bin/bash

#SBATCH --time=48:00:00   
#SBATCH --nodes=1  
#SBATCH --cpus-per-task=20
#SBATCH --mem=256Gb
#SBATCH --partition=ceres

module load miniconda
source activate puzz
module load apptainer

SINGULARITY_TMPDIR=${APPTAINER_TMPDIR}

SAMPLE=$1
SPECIES="Guava"

PUZZLER="apptainer exec /project/coffea_pangenome/Software/Merondun/apptainers/puzzler_v1.1.sif"
WD=/project/coffea_pangenome/${SPECIES}/Assemblies/20250101_JustinAssemblies
OUTDIR=/project/coffea_pangenome/${SPECIES}/Assemblies/20250101_JustinAssemblies/joint_scaffold/summaries
BUSCO_DOWNLOAD_PATH=${OUTDIR}/busco

module load busco5
mkdir -p ${BUSCO_DOWNLOAD_PATH}

if [ ! -d "${BUSCO_DOWNLOAD_PATH}/busco_downloads/lineages/embryophyta_odb10" ]; then
    echo -e "\e[43m~~~~ Downloading BUSCO Dataset ~~~~\e[0m"
    busco --download embryophyta_odb10 --download_path ${BUSCO_DOWNLOAD_PATH}
else
    echo -e "\e[42m~~~~ BUSCO dataset already exists, skipping ~~~~\e[0m"
fi

GENOME=${WD}/joint_scaffold/${SAMPLE}.fa
###### BUSCO ######
if [ ! -s ${OUTDIR}/${SAMPLE}.busco.txt ] && [ -s "${GENOME}" ]; then

    echo -e "\e[43m~~~~ Starting BUSCO Analysis for ${SAMPLE} ~~~~\e[0m"
    mkdir -p ${OUTDIR}/busco/${SAMPLE}
    cd ${OUTDIR}/busco/${SAMPLE}

    busco -i ${GENOME} \
        -l ${BUSCO_DOWNLOAD_PATH}/busco_downloads/lineages/embryophyta_odb10 \
        -m genome \
        -c 20 \
        -o ${SAMPLE} \
        -f

    mv ${SAMPLE}/short_summary.specific.embryophyta_odb10.${SAMPLE}.txt ${OUTDIR}/${SAMPLE}.busco.txt
    cd ${OUTDIR}
    rm -rf ${OUTDIR}/busco/${SAMPLE}

else	
    echo -e "\e[42m~~~~ Skipping BUSCO Analysis for ${SAMPLE}, already exists ~~~~\e[0m"
fi

###### Summary ######
if [ ! -s ${OUTDIR}/${SAMPLE}.summary.txt ]; then
    
    cd ${OUTDIR}
    echo -e "\e[43m~~~~ Summarizing Assembly for ${SAMPLE} ~~~~\e[0m"

    # Get metrics for final assembly

    FINAL_STATS=$(assembly_stats ${GENOME})
    SIZE=$(echo "$FINAL_STATS" | tr ':' '\n' | grep -A 1 'total_bps' | tail -n1 | sed 's/,//g')
    SEQS=$(echo "$FINAL_STATS" | tr ':' '\n' | grep -A 1 'sequence' | tail -n1 | sed 's/,//g')
    CTGS=$(echo "$FINAL_STATS" | tr ':' '\n' | grep -A 1 'sequence' | head -n2 | tail -n1 | sed 's/,//g')
    SCAF_N50=$(echo "$FINAL_STATS" | tr ':' '\n' | grep -A 1 'N50' | tail -n1 | sed 's/,//g')
    CONT_N50=$(echo "$FINAL_STATS" | tr ':' '\n' | grep -A 1 'N50' | head -n2 | tail -n1 | sed 's/,//g')
    GAPS=$((CTGS - SEQS))

    # Parse BUSCO results
    BUSCO_COMPLETE=$(grep "C:" ${OUTDIR}/${SAMPLE}.busco.txt | cut -d'[' -f1 | cut -d':' -f2 | cut -d'%' -f1)
    BUSCO_SINGLE=$(grep "C:" ${OUTDIR}/${SAMPLE}.busco.txt | cut -d'[' -f2 | cut -d'%' -f1 | sed 's/S://g')

    # What percentage of the assembly is in chrs1-11, excluding putative haplotype chrs (ending in 'A')
    samtools faidx ${GENOME}
    grep 'chr' ${GENOME}.fai | awk '$1 !~ /A/' | cut -f1 > ${SAMPLE}.chrs
    samtools faidx ${GENOME} -r ${SAMPLE}.chrs > ${SAMPLE}.chr.fa
    CHR_STATS=$(assembly_stats ${SAMPLE}.chr.fa)
    CHR_SIZE=$(echo "$CHR_STATS" | tr ':' '\n' | grep -A 1 'total_bps' | tail -n1 | sed 's/,//g')
    CHR_PROP=$(perl -e "print sprintf('%.4f', $CHR_SIZE / $SIZE)")
    NUM_CHRS=$(cat ${SAMPLE}.chrs | wc -l)

    # Write summary
    echo -e "ID\tSample\tBUSCO_Complete\tBUSCO_singlecopy\tSizeBP\tWithinChrsBP\tPropWithinChrs\tChrs\tSequences\tContigs\tGaps\tContigN50\tScafN50" > ${SAMPLE}.summary.txt
    echo -e "${SAMPLE}\t${SAMPLE}\t$BUSCO_COMPLETE\t$BUSCO_SINGLE\t$SIZE\t$CHR_SIZE\t$CHR_PROP\t$NUM_CHRS\t$SEQS\t$CTGS\t$GAPS\t$CONT_N50\t$SCAF_N50" >> ${SAMPLE}.summary.txt

else	
    echo -e "\e[42m~~~~ Skipping Assembly Summary for ${SAMPLE}, already exists ~~~~\e[0m"

fi
