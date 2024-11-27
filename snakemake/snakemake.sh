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

dry_run=""

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --sample)
            SAMPLE="$2"
            shift 2
            ;;
        --reference)
            REFERENCE="$2"
            shift 2
            ;;
        --wd)
            WORKDIR="$2"
            shift 2
            ;;
        --container)
            CONTAINER="$2"
            shift 2
            ;;
        --dry-run)
            dry_run="--dry-run"
            shift 1
            ;;
        *)
            echo "Error: Unknown argument $1"
            exit 1
            ;;
    esac
done

# Check if all required arguments are provided
if [ -z "$SAMPLE" ] || [ -z "$REFERENCE" ] || [ -z "$WORKDIR" ] || [ -z "$CONTAINER" ]; then
    echo "Error: Missing required arguments"
    echo "Usage: sbatch $0 --sample <sample_name> --reference <reference_path> --wd <working_directory> --container <container_path>"
    echo "Example: sbatch $0 --sample HART001 --reference /path/to/close_species_reference.fasta --wd /path/to/workdir --container /path/to/container.sif"
    exit 1
fi

# Check if sample exists in CSV
if ! grep -q "^${SAMPLE}," samples.csv; then
    echo "Error: Sample ${SAMPLE} not found in samples.csv"
    exit 1
fi

# Get ploidy, hifi, and HiC paths from csv for this sample
PLOIDY=$(awk -F',' -v sample="$SAMPLE" '$1 == sample {print $2}' samples.csv)
HIFI=$(awk -F',' -v sample="$SAMPLE" '$1 == sample {print $3}' samples.csv)
HIC_R1=$(awk -F',' -v sample="$SAMPLE" '$1 == sample {print $4}' samples.csv)
HIC_R2=$(awk -F',' -v sample="$SAMPLE" '$1 == sample {print $5}' samples.csv)

# Verify HiC data is present
if [ -z "${HIC_R1}" ] || [ -z "${HIC_R2}" ]; then
    echo "Error: HiC data is required. Both HiC R1 and R2 must be specified in samples.csv"
    exit 1
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
reference: "${REFERENCE}"
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
echo "HiC reads: ${HIC_R1}, ${HIC_R2}"

snakemake -s Snakefile \
    --configfile config_${SAMPLE}.yaml \
    --use-apptainer \
    ${dry_run} \
    -j 1