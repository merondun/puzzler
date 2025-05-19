#!/bin/bash
#SBATCH --time=180:00:00   
#SBATCH --nodes=1  
#SBATCH --ntasks-per-node=20
#SBATCH --mem=64Gb
#SBATCH --partition=ceres

module load apptainer
SINGULARITY_TMPDIR=${APPTAINER_TMPDIR}
export AUGUSTUS_CONFIG_PATH=/project/coffea_pangenome/Software/Merondun/Augustus/config

SAMPLE=$1
masked=/project/coffea_pangenome/Artocarpus/Assemblies/20250101_JustinAssemblies/joint_scaffold/EarlGrey/softmasked
WD=/project/coffea_pangenome/Artocarpus/Assemblies/20250101_JustinAssemblies/annotation/${SAMPLE}

rm -rf ${WD} /project/coffea_pangenome/Software/Merondun/Augustus/config/species/${SAMPLE}
mkdir -p ${WD}

#if [ ! -s chr.divergence.txt ] && [ "${IT}" = "hap" ]; then

apptainer exec \
  ~/symlinks/software/apptainers/braker3.sif braker.pl \
  --AUGUSTUS_CONFIG_PATH=/project/coffea_pangenome/Software/Merondun/Augustus/config \
  --workingdir ${WD} --threads 20 --species=${SAMPLE} --genome=${masked}/${SAMPLE}.pri.softmasked.fasta \
	--rnaseq_sets_ids=SRR5997516,SRR5997517,SRR5997518,SRR5997519,SRR5997520,SRR5997521,SRR5997522,SRR5997523,SRR5997524,SRR5997525,SRR5997526,SRR5997527,SRR5997536,SRR5997537,SRR5997538,SRR5997539 \
	--rnaseq_sets_dirs=/project/coffea_pangenome/Artocarpus/Expression/unzipped
