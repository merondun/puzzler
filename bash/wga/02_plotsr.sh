cd out

# Run this first to generate .cl files 
# for i in $(ls *fai | sed 's/.pri.chr.fa.fai//g'); do awk '{print $1,$2}' ${i}*fai > ../out/${i}.cl; done
plotsr --genomes ../full_genomes.txt --itx -S 0.7 -o ~/symlinks/bf/mero_scratch/figures/20250423_BF_Pangenome_SYRI-Chr10.png -W7 -H 10 -f6 -s 100000 \
	--sr HART001_1_HART001_2syri.out \
	--sr HART001_2_HART030_1syri.out \
	--sr HART030_1_HART030_2syri.out \
	--sr HART030_2_HART050_1syri.out \
	--sr HART050_1_HART050_2syri.out \
	--sr HART050_2_HART069_1syri.out \
	--sr HART069_1_HART069_2syri.out \
	--sr HART069_2_HART032_1syri.out \
	--sr HART032_1_HART032_2syri.out \
	--sr HART032_2_HART032_3syri.out \
	--sr HART032_3_HART033_1syri.out \
	--sr HART033_1_HART033_2syri.out \
	--sr HART033_2_HART033_3syri.out \
	--sr HART033_3_HART038_1syri.out \
	--sr HART038_1_HART038_2syri.out \
	--sr HART038_2_HART038_3syri.out \
	--sr HART038_3_HART046_1syri.out \
	--sr HART046_1_HART046_2syri.out \
	--sr HART046_2_HART046_3syri.out \
	--sr HART046_3_HART049_1syri.out \
	--sr HART049_1_HART049_2syri.out \
	--sr HART049_2_HART049_3syri.out \
	--sr HART049_3_HART053_1syri.out \
	--sr HART053_1_HART053_2syri.out \
	--sr HART053_2_HART053_3syri.out \
	--sr HART053_3_H6_1syri.out \
	--sr H6_1_H6_2syri.out \
	--sr H6_2_H6_3syri.out \
	--sr H6_3_H6_4syri.out
