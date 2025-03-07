#!/bin/bash

#SBATCH --time=48:00:00   
#SBATCH --nodes=1  
#SBATCH --ntasks-per-node=5
#SBATCH --mem=75Gb
#SBATCH --partition=short

# e.g. sbatch 01_Mash.sh /project/coffea_pangenome/Artocarpus/Assemblies/20250101_JustinAssemblies/joint_scaffold
WD=$1
module load miniconda
source activate puzz

cd ${WD}
mash sketch -o AllvAll.msh *pri.fa
mash dist AllvAll.msh AllvAll.msh > AllvAll_pairwise_distances.txt

