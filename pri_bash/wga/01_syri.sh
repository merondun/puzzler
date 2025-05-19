#!/bin/bash

#SBATCH --time=48:00:00   
#SBATCH --nodes=1  
#SBATCH --ntasks-per-node=12
#SBATCH --mem=150Gb
#SBATCH --partition=ceres

module load miniconda 

# Check if the correct number of arguments is provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <first_pair> <second_pair>"
    exit 1
fi

target=$1  #reference
query=$2 #target

tfile=$(realpath fastas/${target}.fa)
qfile=$(realpath fastas/${query}.fa)

echo "WGA between ${target} and ${query}"

mkdir -p out work
cd work

if [ ! -s ${target}_${query}.bam ]; then
	echo -e "\e[43m~~~~ Initial Alignment between ${target} and ${query} ~~~~\e[0m"
    source activate syri
    # Using minimap2 for generating alignment. Any other whole genome alignment tool can also be used.
    minimap2 -t 10 -ax asm5 -H -f 100 --rmq=no --secondary=no --eqx ${tfile} ${qfile} | samtools view -@ 10 -bS - | samtools sort -@ 10 -o ${target}_${query}.bam 2> work/${target}_${query}.initialmap.log
    conda deactivate
fi 

if [ ! -s ${target}_${query}fixchr.qry.filtered.fa ]; then
	echo -e "\e[43m~~~~ Inverting inverted chromosomes for ${target} and ${query} ~~~~\e[0m"
    source activate fixchr
    # Re-orient so that chromosomes are largely in syntenic order
    fixchr -c ${target}_${query}.bam -r ${tfile} -q ${qfile} -F B --prefix ${target}_${query} --contig_size 100000
    conda deactivate
fi

if [ ! -s ${target}_${query}.reorient.bam ]; then
	echo -e "\e[43m~~~~ Second Alignment between ${target} and ${query} ~~~~\e[0m"
    source activate syri
    # Re-run minimap2
    minimap2 -t 10 -ax asm5 -H -f 100 --rmq=no --secondary=no --eqx ${target}_${query}fixchr.ref.filtered.fa ${target}_${query}fixchr.qry.filtered.fa | samtools view -@ 10 -bS - | samtools sort -@ 10 -o ${target}_${query}.reorient.bam 2> work/${target}_${query}.reorientmap.log
fi 

if [ ! -s ${target}_${query}syri.summary ]; then

    # Run syri
    echo -e "\e[43m~~~~ Running SYRI between ${target} and ${query} ~~~~\e[0m"
    syri --nc 10 -c ${target}_${query}.reorient.bam -r ${target}_${query}fixchr.ref.filtered.fa -q ${target}_${query}fixchr.qry.filtered.fa -k -F B --prefix ${target}_${query} 2> ${target}_${query}.std.syrilog
    # Migrate
    cp ${target}_${query}*summary ${target}_${query}*vcf ${target}_${query}*out ../out/
else
    echo -e "\e[42m~~~~ Skipping SYRI between ${target} and ${query}, already exists ~~~~\e[0m"
fi