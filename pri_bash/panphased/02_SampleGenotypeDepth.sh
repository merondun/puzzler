#!/bin/bash

#SBATCH --time=2-00:00:00   
#SBATCH --nodes=1  
#SBATCH --ntasks-per-node=8
#SBATCH --mem=36Gb
#SBATCH --partition=ceres

SINGULARITY_TMPDIR=${APPTAINER_TMPDIR}

module load miniconda
source activate panphase

t=8
MEM=36

SAMPLE_FILE="/project/coffea_pangenome/Artocarpus/Assemblies/20250101_JustinAssemblies/samples.csv"
GENOME=/project/coffea_pangenome/Artocarpus/Assemblies/20250101_JustinAssemblies/joint_scaffold/HART001.pri.fa
WD=/project/coffea_pangenome/Artocarpus/Assemblies/20250101_JustinAssemblies/panphase

cd ${WD}

# Merge all the samples for each chromosome 
for CHR in Chr01 Chr02 Chr03 Chr04; do 
	bcftools merge vcfs/${CHR}_*named.vcf.gz -Oz -o vcfs/${CHR}.full.vcf.gz
done 

# Merge all the chromosomes.. 
bcftools concat vcfs/*.full.vcf.gz | \
	bcftools view --min-alleles 2 --max-alleles 2 --types snps | \
    bcftools annotate -x INFO,^FORMAT/GT,^FORMAT/PS -Oz -o vcfs/full.snps.vcf.gz
bcftools index vcfs/full.snps.vcf.gz

# Create a backup merged file
cp vcfs/full.snps.vcf.gz vcfs/full.snps.vcf.gz.bak

# Output bed of SNP pos for mosdepth
bcftools query -f '%CHROM\t%POS0\t%POS\n' vcfs/full.snps.vcf.gz > vcfs/full.snps.bed

# For each sample, assess coverage at that site 
echo '##FORMAT=<ID=PLOIDY,Number=1,Type=Integer,Description="Sample-specific ploidy">' > annotation_header.txt
echo '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Sample-specific allele depth">' >> annotation_header.txt

for SAMPLE in $(cat Samples.list); do 
    PLOIDY=$(awk -F',' -v sample="$SAMPLE" '$1 == sample {print $2}' ${SAMPLE_FILE})
	echo -e "\e[43m~~~~ Extracting genotype coverage for ${SAMPLE}, ploidy: ${PLOIDY}N ~~~~\e[0m"

	# Extract depth for each SNP 
	echo -e "#CHROM\tPOS\tFORMAT/DP\tFORMAT/PLOIDY" > vcfs/${SAMPLE}.tmp.regions.tsv
	mosdepth --by vcfs/full.snps.bed --no-per-base vcfs/${SAMPLE}.tmp ${WD}/bams/${SAMPLE}.hifi.sorted.bam
	zcat vcfs/${SAMPLE}.tmp.regions.bed.gz | awk -v p=${PLOIDY} '{OFS="\t"}{print $1, $3, $4, p}' >> vcfs/${SAMPLE}.tmp.regions.tsv
	bgzip -c vcfs/${SAMPLE}.tmp.regions.tsv > vcfs/${SAMPLE}.tmp.regions.tsv.gz
	tabix -s 1 -b 2 -e 2 vcfs/${SAMPLE}.tmp.regions.tsv.gz

    # Calculate FMT/DP
    bcftools view vcfs/full.snps.vcf.gz --samples ${SAMPLE} -Ou | \
    	bcftools annotate -x INFO,^FORMAT/DP -h annotation_header.txt | \
    	bcftools annotate -a vcfs/${SAMPLE}.tmp.regions.tsv.gz -c CHROM,POS,FORMAT/DP,FORMAT/PLOIDY -o vcfs/${SAMPLE}.annotation.vcf.gz
    bcftools index vcfs/${SAMPLE}.annotation.vcf.gz
    
    # Add FMT/DP for that sample
    bcftools annotate -a vcfs/${SAMPLE}.annotation.vcf.gz -c FORMAT/DP,FORMAT/PLOIDY -Oz -o vcfs/full.snps.ann.vcf.gz vcfs/full.snps.vcf.gz
    mv vcfs/full.snps.ann.vcf.gz vcfs/full.snps.vcf.gz
    bcftools index -f vcfs/full.snps.vcf.gz
               
done 

# rm vcfs/*tsv* vcfs/*annotation* 