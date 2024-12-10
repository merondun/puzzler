#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64Gb
#SBATCH --partition=short

module load miniconda
source activate snakemake
module load apptainer

cat << "EOF"
=========================================================================
  /$$$$$$   /$$$$$$  /$$      /$$  /$$$$$$  /$$   /$$  /$$$$$$  /$$      
 /$$__  $$ /$$__  $$| $$$    /$$$ /$$__  $$| $$  | $$ /$$__  $$| $$      
| $$  \ $$| $$  \__/| $$$$  /$$$$| $$  \ $$| $$  | $$| $$  \ $$| $$      
| $$$$$$$$|  $$$$$$ | $$ $$/$$ $$| $$  | $$| $$  | $$| $$$$$$$$| $$      
| $$__  $$ \____  $$| $$  $$$| $$| $$  | $$| $$  | $$| $$__  $$| $$      
| $$  | $$ /$$  \ $$| $$\  $ | $$| $$/$$ $$| $$  | $$| $$  | $$| $$      
| $$  | $$|  $$$$$$/| $$ \/  | $$|  $$$$$$/|  $$$$$$/| $$  | $$| $$$$$$$$
|__/  |__/ \______/ |__/     |__/ \____ $$$ \______/ |__/  |__/|________/
                                       \__/                              
=========================================================================
EOF

dry_run=""
unlock=""

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --id)
            ID="$2"
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
        --unlock)
            unlock="--unlock"
            shift 1
            ;;
        *)
            echo "Error: Unknown argument $1"
            exit 1
            ;;
    esac
done

# Check if all required arguments are provided
if [ -z "$ID" ] || [ -z "$WORKDIR" ] || [ -z "$CONTAINER" ]; then
    echo "Error: Missing required arguments"
    echo "Usage: sbatch $0 --ID <assembly id> --wd <working_directory> --container <container_path>"
    echo "Example: sbatch $0 --ID HART001_1 --wd /path/to/workdir --container /path/to/container.sif"
    exit 1
fi

# Check if file exists in CSV
if ! grep -q "^${ID}," samples.csv; then
    echo "Error: ID ${ID} not found in samples.csv"
    exit 1
fi

# Get ploidy, hifi, and HiC paths from csv for this sample
FILE=$(awk -F',' -v id="${ID}" '$1 == id {print $2}' samples.csv)
SAMPLE=$(awk -F',' -v id="${ID}" '$1 == id {print $3}' samples.csv)
HAP=$(awk -F',' -v id="${ID}" '$1 == id {print $4}' samples.csv)

# Verify fasta data is present
if [ -z "${FILE}" ]; then
    echo "Error: fasta is required, specify it in samples.csv"
    exit 1
fi

# Create sample directory
mkdir -p sample_runs/${ID}
cd sample_runs/${ID}

# Copy Snakefile
cp ../../Snakefile .

# Create config for this sample
cat > config_${ID}.yaml << EOF
workdir: "${WORKDIR}"
container: "${CONTAINER}"
samples:
    ${ID}:
        file: "${FILE}"
        id: "${ID}"
        sample: "${SAMPLE}"
        hap: "${HAP}"

EOF

# Run snakemake
snakemake -s Snakefile \
    --configfile config_${ID}.yaml \
    --use-apptainer \
    ${dry_run} \
    ${unlock} \
    --cores 16