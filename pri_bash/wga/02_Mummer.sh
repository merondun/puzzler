#!/bin/bash

#SBATCH --time=48:00:00   
#SBATCH --nodes=1  
#SBATCH --ntasks-per-node=10
#SBATCH --partition=ceres

# Check if the correct number of arguments is provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <first_pair> <second_pair>"
    exit 1
fi

#Submit as cat Pairs.list  | xargs -L 1 sbatch 02_Mummer.sh
module load miniconda
source activate wga

# This comes almost directly from scripts from https://github.com/lurebgi/monkParakeet
# Use a file with $SAMPLE_1 \t $SAMPLE_2 within pairs.list, submit
target=$1  #reference
query=$2 #target
GENOME_DIR=/project/coffea_pangenome/Artocarpus/Comparative_Paper/assemblies

mkdir -p aln
cd aln
align_len=1000 # the cutoff for filtering small alignments, you may need to reduce it for species with large divergence

# generating nucmer alignments
nucmer -l 50 -b 400  -t 10 ${GENOME_DIR}/${target}.pri.chr.fa  ${GENOME_DIR}/${query}.pri.chr.fa -p ${target}-${query}
delta-filter -1 -l 400 ${target}-${query}.delta > ${target}-${query}.delta.filt
show-coords -H -c -l -o -r -T ${target}-${query}.delta.filt > ${target}-${query}.delta.filt.coords

# preparing simple files, this comes straight from Parakeet paper 
cat ${target}-${query}.delta.filt.coords  |  awk '$3<$4 && $5>"'$align_len'" && $7<99' | awk '!a[$12"_"$1]++' | awk '!a[$12"_"$2]++' |  awk '$12==ref && $13==ref2{  if($1-end<100000 && $3-end2<100000 && $3-end2>-50000){end=$2;end2=$4;n++;y=0}else{y++; if(y==1){next;}else{ print ref"__"st"\t"ref"__"end"\t"ref2"__"st2"\t"ref2"__"end2"\t"n"\t+"; st=$1;end=$2;st2=$3;end2=$4;  n=0} }  }$12!=ref || $13!=ref2{y++; if(y==1){next;}else{print ref"__"st"\t"ref"__"end"\t"ref2"__"st2"\t"ref2"__"end2"\t"n"\t+"; n=0;ref=$12; ref2=$13; st=$1;end=$2;st2=$3;end2=$4;}}' | sed 1d | awk '$5>2' > ${target}-${query}.simple
cat ${target}-${query}.delta.filt.coords | awk '$3>$4 && $5>"'$align_len'" && $7<99' | awk '!a[$12"_"$1]++' | awk '!a[$12"_"$2]++' | awk '$5=$4' | awk '$4=$3' | awk '$3=$5' | awk '$12==ref && $13==ref2{  if($1-end<100000 && st2-$4<100000 && st2-$4>-50000){end=$2;st2=$3;n++;y=0}else{y++; if(y==1){next;}else{  print ref"__"st"\t"ref"__"end"\t"ref2"__"st2"\t"ref2"__"end2"\t"n"\t-"; st=$1;end=$2;st2=$3;end2=$4;  n=0}}   } $12!=ref || $13!=ref2{y++; if(y==1){next;}else{ print ref"__"st"\t"ref"__"end"\t"ref2"__"st2"\t"ref2"__"end2"\t"n"\t-"; n=0;ref=$12; ref2=$13; st=$1;end=$2;st2=$3;end2=$4;}} ' | sed 1d | awk '$5>2' >>  ${target}-${query}.simple

# preparing bed files this comes straight from Parakeet paper 
cat ${target}-${query}.delta.filt.coords |  awk '$5>"'$align_len'" && $7<99' | awk '!a[$12"_"$1]++' | awk '!a[$12"_"$2]++' | awk '{print $12"\t"$1"\t"$2"\t"$12"__"$1"\t0\t+"}' > ${target}.bed.common
cat ${target}-${query}.delta.filt.coords |  awk '$5>"'$align_len'" && $7<99' | awk '!a[$12"_"$1]++' | awk '!a[$12"_"$2]++' | awk '$3<$4{print $13"\t"$3"\t"$4"\t"$13"__"$2"\t0\t+"}$3>$4{print $13"\t"$4"\t"$3"\t"$13"__"$3"\t0\t+"}' > ${query}.bed.common

cat ${target}-${query}.simple | sed 's/__/\t/g' | awk '{print $1"\t"$2"\t"$2+1"\t"$1"__"$2"\t0\t+\n"$3"\t"$4"\t"$4+1"\t"$3"__"$4"\t0\t+"}' | cat - ${target}.bed.common > ${target}.bed.${target}-${query}  
cat ${target}-${query}.simple | sed 's/__/\t/g' | awk '{print $5"\t"$6"\t"$6+1"\t"$5"__"$6"\t0\t+\n"$7"\t"$8"\t"$8+1"\t"$7"__"$8"\t0\t+"}' | cat - ${query}.bed.common > ${query}.bed.${target}-${query}
