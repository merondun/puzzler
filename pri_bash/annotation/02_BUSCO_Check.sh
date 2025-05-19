#!/bin/bash

#SBATCH --time=48:00:00   
#SBATCH --nodes=1  
#SBATCH --cpus-per-task=20
#SBATCH --mem=1024Gb
#SBATCH --partition=ceres

module load miniconda
source activate puzz
module load apptainer

if [ -z "$1" ]; then
    echo "Error: File prefix positional argument is required."
    exit 1
fi

# Submit as m84125_250429_221012_s3.hifi_reads.bcM0001.skera.IsoSeqX_bc01_5p--IsoSeqX_3p.Barcodes.flnc
FILE=$1
t=20

WD=/project/coffea_pangenome/Artocarpus/RawData/transfer/Breadfruit_BamFiles/skera
OUTDIR=/project/coffea_pangenome/Artocarpus/RawData/transfer/Breadfruit_BamFiles/skera/busco
BUSCO_DOWNLOAD_PATH=/project/coffea_pangenome/Artocarpus/RawData/transfer/IsoSeqRawData/init

module load busco5
mkdir -p ${BUSCO_DOWNLOAD_PATH}

if [ ! -d "${BUSCO_DOWNLOAD_PATH}/busco_downloads/lineages/embryophyta_odb10" ]; then
    echo -e "\e[43m~~~~ Downloading BUSCO Dataset ~~~~\e[0m"
    busco --download embryophyta_odb10 --download_path ${BUSCO_DOWNLOAD_PATH}
else
    echo -e "\e[42m~~~~ BUSCO dataset already exists, skipping ~~~~\e[0m"
fi

#samtools fasta ${FILE}.bam > ${FILE}.fasta
#sed -i 's@/@_@g' ${FILE}.fasta
GENOME=$(realpath ${FILE}.fasta)

###### BUSCO ######
if [ ! -s ${OUTDIR}/${FILE}.busco.txt ] && [ -s ${GENOME} ]; then

    echo -e "\e[43m~~~~ Starting BUSCO Analysis for ${SAMPLE} ~~~~\e[0m"
    mkdir -p ${OUTDIR}/${FILE}
    cd ${OUTDIR}/${FILE}

    busco -i ${GENOME} \
        -l ${BUSCO_DOWNLOAD_PATH}/busco_downloads/lineages/embryophyta_odb10 \
        -m genome \
        -c ${t} \
        -o ${FILE} \
        --metaeuk_rerun_parameters="-s=2,--remove-tmp-files=1" \
        --metaeuk_parameters="-s=2,--remove-tmp-files=1" \
        -f

    mv ${FILE}/short_summary.specific.embryophyta_odb10.${FILE}.txt ${OUTDIR}/${FILE}.busco.txt
    cd ${OUTDIR}
    rm -rf ${OUTDIR}/${FILE}

else
        echo -e "\e[42m~~~~ Skipping BUSCO Analysis for ${FILE}, already exists ~~~~\e[0m"
fi
