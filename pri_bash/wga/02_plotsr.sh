cd out

# Run this first to generate .cl files 
# for i in $(ls *fai | sed 's/.pri.chr.fa.fai//g'); do awk '{print $1,$2}' ${i}*fai > ../out/${i}.cl; done
plotsr --genomes ../full_genomes.txt --itx -S 0.7 -o ../20250312_WGA-100KB.pdf -W7 -H 10 -f6 -s 100000 \
	--sr HART001_HART030syri.out \
	--sr HART030_HART050syri.out \
	--sr HART050_HART069syri.out \
	--sr HART069_HART032syri.out \
	--sr HART032_HART033syri.out \
	--sr HART033_HART038syri.out \
	--sr HART038_HART046syri.out \
	--sr HART046_HART049syri.out \
	--sr HART049_HART053syri.out \
	--sr HART053_H6syri.out \
	--sr H6_HART067syri.out \
	--sr HART067_HART063syri.out \
	--sr HART063_HART061syri.out \
	--sr HART061_HART062syri.out \
	--sr HART062_HART058syri.out \
	--sr HART058_HART027syri.out \
	--sr HART027_HART060syri.out \
	--sr HART060_N97_50syri.out