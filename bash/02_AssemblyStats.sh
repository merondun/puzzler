#!/bin/bash

#SBATCH --time=48:00:00   
#SBATCH --nodes=1  
#SBATCH --cpus-per-task=20
#SBATCH --mem=256Gb
#SBATCH --partition=ceres

module load miniconda
source activate puzz
module load apptainer

SAMPLE=$1

PUZZLER="apptainer exec /project/coffea_pangenome/Software/Merondun/apptainers/puzzler_v1.1.sif"
WD=/project/coffea_pangenome/Artocarpus/Assemblies/20250101_JustinAssemblies
OUTDIR=/project/coffea_pangenome/Artocarpus/Assemblies/20250101_JustinAssemblies/joint_scaffold/summaries
BUSCO_DOWNLOAD_PATH=${OUTDIR}/busco

module load busco5
mkdir -p ${BUSCO_DOWNLOAD_PATH}

if [ ! -d "${BUSCO_DOWNLOAD_PATH}/busco_downloads/lineages/embryophyta_odb10" ]; then
    echo -e "\e[43m~~~~ Downloading BUSCO Dataset ~~~~\e[0m"
    busco --download embryophyta_odb10 --download_path ${BUSCO_DOWNLOAD_PATH}
else
    echo -e "\e[42m~~~~ BUSCO dataset already exists, skipping ~~~~\e[0m"
fi

for IT in pri hap; do 

    GENOME=${WD}/joint_scaffold/${SAMPLE}.${IT}.fa
    ID="${SAMPLE}.${IT}"
    ###### BUSCO ######
    if [ ! -s ${OUTDIR}/${ID}.busco.txt ] && [ -s "${GENOME}" ]; then
    
        echo -e "\e[43m~~~~ Starting BUSCO Analysis for ${ID} ~~~~\e[0m"
        mkdir -p ${OUTDIR}/busco/${ID}
        cd ${OUTDIR}/busco/${ID}

        busco -i ${GENOME} \
            -l ${BUSCO_DOWNLOAD_PATH}/busco_downloads/lineages/embryophyta_odb10 \
            -m genome \
            -c 20 \
            -o ${ID} \
            -f

        mv ${ID}/short_summary.specific.embryophyta_odb10.${ID}.txt ${OUTDIR}/${ID}.busco.txt
        cd ${OUTDIR}
        rm -rf ${OUTDIR}/busco/${ID}

    else	
		echo -e "\e[42m~~~~ Skipping BUSCO Analysis for ${ID}, already exists ~~~~\e[0m"
	fi

    ###### Summary ######
    if [ ! -s ${OUTDIR}/${ID}.summary.txt ]; then
        
        cd ${OUTDIR}
        echo -e "\e[43m~~~~ Summarizing Assembly for ${ID} ~~~~\e[0m"

        # Get metrics for final assembly

        FINAL_STATS=$(assembly_stats ${GENOME})
        SIZE=$(echo "$FINAL_STATS" | tr ':' '\n' | grep -A 1 'total_bps' | tail -n1 | sed 's/,//g')
        SEQS=$(echo "$FINAL_STATS" | tr ':' '\n' | grep -A 1 'sequence' | tail -n1 | sed 's/,//g')
        CTGS=$(echo "$FINAL_STATS" | tr ':' '\n' | grep -A 1 'sequence' | head -n2 | tail -n1 | sed 's/,//g')
        SCAF_N50=$(echo "$FINAL_STATS" | tr ':' '\n' | grep -A 1 'N50' | tail -n1 | sed 's/,//g')
        CONT_N50=$(echo "$FINAL_STATS" | tr ':' '\n' | grep -A 1 'N50' | head -n2 | tail -n1 | sed 's/,//g')
        GAPS=$((CTGS - SEQS))

        # Parse BUSCO results
        BUSCO_COMPLETE=$(grep "C:" ${OUTDIR}/${ID}.busco.txt | cut -d'[' -f1 | cut -d':' -f2 | cut -d'%' -f1)
        BUSCO_SINGLE=$(grep "C:" ${OUTDIR}/${ID}.busco.txt | cut -d'[' -f2 | cut -d'%' -f1 | sed 's/S://g')

        # What percentage of the assembly is in chrs1-28? 
        samtools faidx ${GENOME}
        grep 'Chr' ${GENOME}.fai | cut -f1 > ${ID}.chrs
        samtools faidx ${GENOME} -r ${ID}.chrs > ${ID}.chr.fa
        CHR_STATS=$(assembly_stats ${ID}.chr.fa)
        CHR_SIZE=$(echo "$CHR_STATS" | tr ':' '\n' | grep -A 1 'total_bps' | tail -n1 | sed 's/,//g')
        CHR_PROP=$(perl -e "print sprintf('%.4f', $CHR_SIZE / $SIZE)")
        NUM_CHRS=$(cat ${ID}.chrs | wc -l)

        # Write summary
        echo -e "ID\tSample\tAssembly\tBUSCO_Complete\tBUSCO_singlecopy\tSizeBP\tWithinChrsBP\tPropWithinChrs\tChrs\tSequences\tContigs\tGaps\tContigN50\tScafN50" > ${ID}.summary.txt
        echo -e "${ID}\t${SAMPLE}\t${IT}\t$BUSCO_COMPLETE\t$BUSCO_SINGLE\t$SIZE\t$CHR_SIZE\t$CHR_PROP\t$NUM_CHRS\t$SEQS\t$CTGS\t$GAPS\t$CONT_N50\t$SCAF_N50" >> ${ID}.summary.txt

    else	
        echo -e "\e[42m~~~~ Skipping Assembly Summary for ${ID}, already exists ~~~~\e[0m"

    fi

done