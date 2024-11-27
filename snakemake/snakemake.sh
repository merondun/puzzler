#!/bin/bash
#SBATCH --time=48:00:00   
#SBATCH --nodes=1  
#SBATCH --cpus-per-task=64
#SBATCH --mem=370Gb
#SBATCH --partition=short

module load miniconda
source activate snakemake
module load busco5

cat << "EOF"
=======================================================================
__________ ____ _______________________.____     _____________________ 
\______   \    |   \____    /\____    /|    |    \_   _____/\______   \
 |     ___/    |   / /     /   /     / |    |     |    __)_  |       _/
 |    |   |    |  / /     /_  /     /_ |    |___  |        \ |    |   \
 |____|   |______/ /_______ \/_______ \|_______ \/_______  / |____|_  /
                           \/        \/        \/        \/         \/ 
=======================================================================

EOF

# Check arguments
if [ "$#" -ne 1 ]; then
    echo "Error: Incorrect number of arguments"
    echo "Usage: sbatch $0 <sample_name>"
    echo "Example: sbatch $0 HART001"
    exit 1
fi

SAMPLE=$1

WORKDIR="/project/coffea_pangenome/Artocarpus/Assemblies/20241115_JustinAssemblies"
CONTAINER="/project/coffea_pangenome/Software/Merondun/apptainers/puzzler.sif"

# Check if sample exists in CSV
if ! grep -q "^${SAMPLE}," samples.csv; then
    echo "Error: Sample ${SAMPLE} not found in samples.csv"
    exit 1
fi

# Get ploidy and hifi path from csv for this sample
PLOIDY=$(awk -F',' -v sample="$SAMPLE" '$1 == sample {print $2}' samples.csv)
HIFI=$(awk -F',' -v sample="$SAMPLE" '$1 == sample {print $3}' samples.csv)
HIC_R1=$(awk -F',' -v sample="$SAMPLE" '$1 == sample {print $4}' samples.csv)
HIC_R2=$(awk -F',' -v sample="$SAMPLE" '$1 == sample {print $5}' samples.csv)

# Check if HiC fields are empty or contain only whitespace
if [ -z "${HIC_R1// }" ] || [ -z "${HIC_R2// }" ]; then
    HIC_R1=""
    HIC_R2=""
fi

# Create a directory for each sample
mkdir -p sample_runs/${SAMPLE}
cd sample_runs/${SAMPLE}

# Copy necessary files
cp ../../Snakefile .
cp ../../samples.csv .

# Create temporary config
cat > config_${SAMPLE}.yaml << EOF
workdir: "${WORKDIR}"
container: "${CONTAINER}"
reference: "/project/coffea_pangenome/Artocarpus/WholeGenomeAlignments/fastas/ASM2540343.fa"
samples:
    ${SAMPLE}:
        hifi: "${HIFI}"
        ploidy: ${PLOIDY}
        hic_r1: "${HIC_R1}"
        hic_r2: "${HIC_R2}"
EOF

echo "Starting assembly pipeline for ${SAMPLE}"
echo "Ploidy: ${PLOIDY}"
echo "HiFi reads: ${HIFI}"

if [ -n "$HIC_R1" ] && [ -n "$HIC_R2" ]; then
    echo "HiC reads: ${HIC_R1}, ${HIC_R2}"
    echo "Running HiC-enabled pipeline"
else
    echo "No HiC reads provided, running standard pipeline"
fi

snakemake -s Snakefile \
    --configfile config_${SAMPLE}.yaml \
    --use-singularity \
    -j 1 \
    Primary_Assemblies/${SAMPLE}.HifiasmHiFi-HiC-PurgeDups-RagTag.fa
