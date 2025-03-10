#!/bin/bash

#SBATCH --time=24:00:00   
#SBATCH --nodes=1  
#SBATCH --ntasks-per-node=64
#SBATCH --mem=512Gb
#SBATCH --partition=short

# use n=32 and mem = 112 typically 
SINGULARITY_TMPDIR=${APPTAINER_TMPDIR}

module load miniconda
source activate puzz
module load apptainer

if [ -z "$1" ]; then
    echo "Error: sample ID positinonal argument is required."
    exit 1
fi

t=64
MEM=512
SAMPLE=$1

WD=/project/coffea_pangenome/Artocarpus/Assemblies/20250101_JustinAssemblies
SAMPLE_FILE="/project/coffea_pangenome/Artocarpus/Assemblies/20250101_JustinAssemblies/samples.csv"
PLOIDY=$(awk -F',' -v sample="$SAMPLE" '$1 == sample {print $2}' ${SAMPLE_FILE})
NCHRS=$(awk -F',' -v sample="$SAMPLE" '$1 == sample {print $3}' ${SAMPLE_FILE})
HOM_COV=$(awk -F',' -v sample="$SAMPLE" '$1 == sample {print $4}' ${SAMPLE_FILE})
HIFI=$(awk -F',' -v sample="$SAMPLE" '$1 == sample {print $5}' ${SAMPLE_FILE})
HIC_R1=$(awk -F',' -v sample="$SAMPLE" '$1 == sample {print $6}' ${SAMPLE_FILE})
HIC_R2=$(awk -F',' -v sample="$SAMPLE" '$1 == sample {print $7}' ${SAMPLE_FILE})
REFERENCE=$(awk -F',' -v sample="$SAMPLE" '$1 == sample {print $8}' ${SAMPLE_FILE})
PUZZLER="apptainer exec /project/coffea_pangenome/Software/Merondun/apptainers/puzzler_v1.1.sif"

mkdir -p ${WD}/${SAMPLE} ${WD}/logs/juicer ${WD}/logs/haphic
cd ${WD}/${SAMPLE}

#rm -rf 02_HapHiC 01_purgedups 01_maphapstogether 02_priHapHiC

##### Generate Primary Assembly #####, add --hom-cov ${HOM_COV} \ if necessary! 
if [ ! -s ${SAMPLE}.hic.hap1.p_ctg.gfa ]; then
	${PUZZLER} hifiasm --n-hap ${PLOIDY} --hom-cov ${HOM_COV} --primary -t ${t} -o ${SAMPLE} --h1 ${HIC_R1} --h2 ${HIC_R2} ${HIFI} 2> ${WD}/logs/${SAMPLE}.hifiasm.log
else
	echo "~~~~ Skipping hifiasm for ${SAMPLE}, already exists ~~~~"
fi 

##### Loop through both primary assembly pipeline and haplotype-joint pipeline #####

for IT in pri hap; do 

	mkdir -p ${WD}/${SAMPLE}/02_${IT}HapHiC
	cd ${WD}/${SAMPLE}/02_${IT}HapHiC

	##### PURGE #####
	if [ ! -s all.purged.fa ]; then

		if [ "${IT}" = "hap" ]; then

			### Purge and collect gfas: Haplotypes!
			for HAP in $(seq 1 ${PLOIDY}); do 

				if [ ! -s hap${HAP}.purged.fa ]; then
					
				echo "~~~~ Starting Purge_Dups for ${SAMPLE}, Haplotype ${HAP} ~~~~"

				# Ensure fastas have proper contig haplotype names 
				awk '/^S/{print ">"$2;print $3}' ${WD}/${SAMPLE}/${SAMPLE}.hic.hap${HAP}.p_ctg.gfa > hap${HAP}.init.fa
				split_fa hap${HAP}.init.fa > hap${HAP}.split.fa
				minimap2 -t ${t} -xasm5 -DP hap${HAP}.split.fa hap${HAP}.split.fa | gzip -c > hap${HAP}.split.self.paf.gz
				purge_dups -M1000 -E1000 hap${HAP}.split.self.paf.gz > hap${HAP}.dups.bed
				get_seqs hap${HAP}.dups.bed hap${HAP}.init.fa
				# Ensure haplotype IDs are correct 
				sed "s/h1tg/h${HAP}tg/g" purged.fa > hap${HAP}.purged.fa
				mv hap.fa hap${HAP}.hap.fa

				# Also, ensure that the .gfa files have correct haplotype IDs, which they don't as of current hifiasm version! 
				sed "s/h1tg/h${HAP}tg/g" ${WD}/${SAMPLE}/${SAMPLE}.hic.hap${HAP}.p_ctg.gfa > hap${HAP}.gfa

				else	
					echo "~~~~ Skipping ${SAMPLE}, Haplotype ${HAP} as it already exists ~~~~"
					continue
				fi

			done 

			cat hap*purged.fa > all.purged.fa
			rm *split.fa *self.paf.gz *dups.bed *init*

		else 

			echo "~~~~ Starting Purge_Dups for ${SAMPLE}, primary assembly ~~~~"

			#Ensure fastas have proper contig haplotype names 
			awk '/^S/{print ">"$2;print $3}' ${WD}/${SAMPLE}/${SAMPLE}.hic.p_ctg.gfa > pri.init.fa
			split_fa pri.init.fa > pri.split.fa
			minimap2 -t ${t} -xasm5 -DP pri.split.fa pri.split.fa | gzip -c > pri.split.self.paf.gz
			purge_dups -M1000 -E1000 pri.split.self.paf.gz > pri.dups.bed
			get_seqs pri.dups.bed pri.init.fa
			mv purged.fa all.purged.fa
			rm *split.fa *self.paf.gz *dups.bed *init*

		fi

	else
		echo "~~~~ Skipping purge for ${SAMPLE}, ${IT} ~~~~"
	fi

done # Exit hap pri iteration fork 
 
##### ALIGN HIC for Primary Assembly #####
cd ${WD}/${SAMPLE}/02_priHapHiC
if [ ! -s filtered.MQ1.bam ]; then
	echo "~~~~ Mapping HiC reads to ${SAMPLE} assembly ~~~~" 
	# Index reference
	${PUZZLER} bwa index all.purged.fa

	# Align Hi-C reads
	${PUZZLER} bwa mem -5SP -t ${t} all.purged.fa ${HIC_R1} ${HIC_R2} | \
		${PUZZLER} samblaster | \
		${PUZZLER} samtools view - -@ ${t} -S -h -b -F 3340 | \
		${PUZZLER} filter_bam - 1 --nm 3 --threads ${t} | \
		${PUZZLER} samtools view - -b -@ ${t} -o filtered.MQ1.bam

else
	echo "~~~~ Skipping HiC alignment for ${SAMPLE} ~~~~"
fi

##### HAPHIC for Primary Assembly ##### 
if [ ! -s ${WD}/${SAMPLE}/02_priHapHiC/01_haphicMQ1/04.build/scaffolds.fa ]; then

	echo "~~~~ Running HapHiC for ${SAMPLE} ~~~~" 
	mkdir -p ${WD}/${SAMPLE}/02_priHapHiC/01_haphicMQ1${M}
	rm -rf ${WD}/${SAMPLE}/02_priHapHiC/01_haphicMQ1/*

	cd ${WD}/${SAMPLE}/02_priHapHiC/01_haphicMQ1

	${PUZZLER} haphic pipeline ../all.purged.fa ../filtered.MQ1.bam ${NCHRS} --remove_allelic_links ${PLOIDY} --correct_nrounds 2 --max_inflation 20.0 --threads ${t} --processes ${t} 2> ${WD}/logs/haphic/${SAMPLE}.pri.haphic.MQ1.log

	cp ${WD}/${SAMPLE}/02_${IT}HapHiC/01_haphicMQ1/04.build/scaffolds.fa ${WD}/${SAMPLE}/02_priHapHiC/haphic.MQ1.fa
	cp ${WD}/${SAMPLE}/02_${IT}HapHiC/01_haphicMQ1/01.cluster/HapHiC_cluster.log ${WD}/${SAMPLE}/02_priHapHiC/haphic.MQ1.log

else 
	echo "~~~~ Skipping HapHiC for ${SAMPLE} ~~~~"
fi

##### JUICER #####
cd ${WD}/${SAMPLE}/02_priHapHiC
mkdir -p 02_orienting

if [ ! -s ${WD}/logs/juicer/${SAMPLE}.pri-MQ1_JBAT.hic ]; then

	echo "~~~~ Creating .hic file for juicebox for ${SAMPLE} ~~~~"
	# # Orient to jackfruit chromosomes 
	if [ ! -s 02_orienting/asm_to_paf.paf ]; then
		$PUZZLER minimap2 -x asm20 ${REFERENCE} all.purged.fa --secondary=no -t ${t} -o 02_orienting/asm_to_paf.paf
	else
		echo "Skipping initial alignment"
	fi

	$PUZZLER haphic refsort ${WD}/${SAMPLE}/02_priHapHiC/01_haphicMQ1/04.build/scaffolds.raw.agp 02_orienting/asm_to_paf.paf > refsortMQ1.agp 2> refsortMQ1.log
	$PUZZLER samtools faidx haphic.MQ1.fa
	$PUZZLER samtools faidx all.purged.fa

	JUICER_PATH=~/symlinks/software/HapHiC/utils
	${JUICER_PATH}/juicer pre \
		-a -q 1 \
		-o haphic-refsort-MQ1_JBAT \
		filtered.MQ1.bam \
		refsortMQ1.agp \
		all.purged.fa.fai \
		> haphic-refsort-MQ1_JBAT.log 2>&1

	# Extract chromosome sizes and create HiC file
	grep PRE_C_SIZE haphic-refsort-MQ1_JBAT.log | \
		awk '{print $2" "$3}' > chrom.sizesMQ1
	java -Xmx${MEM}G -jar ${JUICER_PATH}/juicer_tools.1.9.9_jcuda.0.8.jar pre \
		-r 5000000,4000000,3000000,2000000,1500000,1000000,500000,250000,100000,50000 \
		haphic-refsort-MQ1_JBAT.txt \
		haphic-refsort-MQ1_JBAT.hic \
		chrom.sizesMQ1
		
	cp haphic-refsort-MQ1_JBAT.hic ${WD}/logs/juicer/${SAMPLE}.pri-MQ1_JBAT.hic
	cp haphic-refsort-MQ1_JBAT.assembly ${WD}/logs/juicer/${SAMPLE}.pri-MQ1_JBAT.assembly

else
	echo "~~~~ No HiC juicer files needed for ${SAMPLE} ${IT} and MQ1 ~~~~"
fi 

#### MANUAL CURATION, THEN THROW THE ${SAMPLE}_JBAT.pri.review.assembly BACK INTO ${WD}/${SAMPLE}/02_priHapHiC


##### HAPLOTYPES ######

mkdir -p ${WD}/logs/chrjuicer ${WD}/logs/chrjuicer/chrcontact_maps

# It's worthwhile to examine divergence between haplotypes
cd ${WD}/${SAMPLE}/02_hapHapHiC
if [ ! -s chr.divergence.txt ]; then

	echo "~~~~ Estimating Haplotype divergence for ${SAMPLE} ~~~~"
	mkdir -p 00_hapchrs 01_chrassemblies 
	for HAP in $(seq 1 ${PLOIDY}); do 

		$PUZZLER minimap2 -x asm20 ${REFERENCE} hap${HAP}.purged.fa --secondary=no -t ${t} -o 00_hapchrs/hap${HAP}_ref.paf 2> 00_hapchrs/hap${HAP}.log
		samtools faidx hap${HAP}.purged.fa
		map_chromosomes --paf 00_hapchrs/hap${HAP}_ref.paf --fai hap${HAP}.purged.fa.fai --out hap${HAP}_chrmap.txt

	done 

	# Now loop through chromosomes, extract the chr fasta and estimate divergence 
	for CHRNUM in $(seq 1 ${NCHRS}); do 
		for HAP in $(seq 1 ${PLOIDY}); do 

			CHR=$(printf "Chr%02d" ${CHRNUM})
			echo "Extracting ${CHR} for ${SAMPLE} hap${HAP}"
			scaffold=$(grep -w ${CHR} hap${HAP}_chrmap.txt | sort -k 4 | head -n1 | awk '{print $1}')
			strand=$(grep -w ${CHR} hap${HAP}_chrmap.txt | sort -k 4 | head -n1 | awk '{print $5}')

			if [ "$strand" == "-" ]; then
				echo "~~~~ ${SAMPLE} hap${HAP} ${CHR} is negative sense, revcomp ~~~~"
				samtools faidx hap${HAP}.purged.fa ${scaffold} | \
				seqkit seq --line-width 0 -t DNA -v -r -p | \
				sed "s/^>.*/>${SAMPLE}#${HAP}#${CHR}/" > 00_hapchrs/hap${HAP}_${CHR}.fa 
			else
				echo "~~~~ ${SAMPLE} hap${HAP} ${CHR} is positive sense, maintain ~~~~"
				samtools faidx hap${HAP}.purged.fa ${scaffold} | \
				sed "s/^>.*/>${SAMPLE}#${HAP}#${CHR}/" > 00_hapchrs/hap${HAP}_${CHR}.fa 
			fi

		done

		cat 00_hapchrs/*${CHR}.fa > 01_chrassemblies/${CHR}.fa 
		${PUZZLER} samtools faidx 01_chrassemblies/${CHR}.fa

		# Calculate divergence, format 
		files=$(find 00_hapchrs -name "*_${CHR}.fa" | tr '\n' ' ')
		mash sketch -o 00_hapchrs/${CHR}.divergence $files
		mash dist 00_hapchrs/${CHR}.divergence.msh 00_hapchrs/${CHR}.divergence.msh | \
			sed "s@00_hapchrs/hap@@g" | sed 's/_/\t/g' | sed 's/.fa//g' | sed 's@\./@@g' | \
			awk '$1 != $3' | sort -k3 | awk '!seen[$1,$3]++ && !seen[$3,$1]++' | awk -v s=${SAMPLE} '{OFS="\t"}{print s, $2, $1, $3, $5, $7}'> 00_hapchrs/${CHR}_divergence.out

	done 

	cat 00_hapchrs/*divergence.out > chr.divergence.txt
	cp chr.divergence.txt ${WD}/logs/divergence_haps/${SAMPLE}.chr.divergence.out 
	awk '{OFS="\t"}{print $1, $2}' 01_chrassemblies/${CHR}.fa.fai > ${WD}/logs/divergence_haps/${SAMPLE}.chrlengths.txt

else 
	echo "~~~~ Skipping haplotype divergence estimation for ${SAMPLE} ~~~~"
fi

cd ${WD}/${SAMPLE}/02_hapHapHiC
# Map reads to whole genome, and then subset.. 
if [ ! -s gw.bam ]; then

	echo "~~~~ Creating HiC bam for ${SAMPLE} ~~~~"
	# Index reference
	cat 01_chrassemblies/Chr*fa > allhaps_chr.fa
	${PUZZLER} samtools faidx allhaps_chr.fa
	${PUZZLER} bwa index allhaps_chr.fa
	${PUZZLER} mock_agp_file.py allhaps_chr.fa > allhaps_chr.agp

	# Align Hi-C reads, filter for MQ 1 
	${PUZZLER} bwa mem -5SP -t ${t} allhaps_chr.fa ${HIC_R1} ${HIC_R2} | \
		${PUZZLER} samblaster | \
		${PUZZLER} samtools view - -@ ${t} -S -h -b -F 3340 | \
		${PUZZLER} filter_bam - 1 --nm 3 --threads ${t} | \
		${PUZZLER} samtools view - -b -@ ${t} -o allhaps_chr.filtered.bam
else 
	echo "~~~~ Skipping mapping HiC ~~~~"
fi

# JUICER
if [ ! -s allhaps_chr_JBAT.hic ]; then
	echo "~~~~ Creating juicer file for ${SAMPLE} ~~~~"

	JUICER_PATH=~/symlinks/software/HapHiC/utils
	${JUICER_PATH}/juicer pre \
		-a -q 1 \
		-o allhaps_chr_JBAT \
		allhaps_chr.filtered.bam \
		allhaps_chr.agp \
		allhaps_chr.fa.fai \
		> allhaps_chr_JBAT.log 2>&1

	# Extract chromosome sizes and create HiC file
	grep PRE_C_SIZE allhaps_chr_JBAT.log | \
		awk '{print $2" "$3}' > allhaps_chr_JBAT.chrom.sizes
	java -Xmx${MEM}G -jar ${JUICER_PATH}/juicer_tools.1.9.9_jcuda.0.8.jar pre \
		-r 5000000,4000000,3000000,2000000,1500000,1000000,500000,250000,100000,50000 \
		allhaps_chr_JBAT.txt \
		allhaps_chr_JBAT.hic \
		allhaps_chr_JBAT.chrom.sizes

	cp allhaps_chr_JBAT.hic ${WD}/logs/chrjuicer/${SAMPLE}.hic
	cp allhaps_chr_JBAT.assembly ${WD}/logs/chrjuicer/${SAMPLE}.assembly

	$PUZZLER haphic plot --threads ${t} allhaps_chr.agp allhaps_chr.filtered.bam --bin_size 1000 --min_len 2
	mv contact_map.pdf ${WD}/logs/chrjuicer/chrcontact_maps/${SAMPLE}.allhapchr.pdf

else 
	echo "~~~~ Skipping juicer file creation for ${SAMPLE} ~~~~"
fi
