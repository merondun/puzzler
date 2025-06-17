# Needlecast fungus Genome Assembly with **`puzzler`**

[Needlecast fungus](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_964106605.1/) has a small genome (25 Mb) and 14 chromosomes, so it's good for a quick tutorial. 

## File Preparation

### Reads

First, grab the [HiFi](https://www.ncbi.nlm.nih.gov/sra/ERX12134039[accn]) and [HiC](https://www.ncbi.nlm.nih.gov/sra/ERX12138301[accn]) data from the SRA. Let's only use 100K HiFi reads and 10M HiC paired reads so it runs quickly.

```bash
#!/bin/bash

#SBATCH --time=4-00:00:00   
#SBATCH --nodes=1  
#SBATCH --ntasks-per-node=1
#SBATCH --partition=ceres

WD=/90daydata/coffea_pangenome/puzzler_trials

mkdir Rhizosphaera_kalkhoffii
cd Rhizosphaera_kalkhoffii
module load sratoolkit/3.2.0

# FUNGUS
#HIFI ERR12760830	   HIC ERR12765203	
prefetch ERR12760830	
fasterq-dump --split-files ERR12760830
head -n 200000 ERR12760830.fastq | gzip -c > ../concat_reads/Fungus.HiFi.fastq.gz

prefetch ERR12765203 --max-size 0	
fasterq-dump --split-files ERR12765203	
head -n 1000000 ERR12765203_1.fastq | gzip -c > ../concat_reads/Fungus.HiC.R1.fastq.gz
head -n 1000000 ERR12765203_2.fastq | gzip -c > ../concat_reads/Fungus.HiC.R2.fastq.gz
```

### [Optional] Simplify Reference Chr Names

I will name the chromosomes according to the published reference assembly:

```bash
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/964/106/605/GCA_964106605.1_gdRhiKalk1.hap1.1/GCA_964106605.1_gdRhiKalk1.hap1.1_genomic.fna.gz
gunzip GCA_964106605.1_gdRhiKalk1.hap1.1_genomic.fna.gz
```

I don't want the re-named chromosomes to be this cryptic, so I will strip all the text between `>` and `chromosome :`. **Note that this is preference only**. Otherwise, your chromosomes would be named e.g.  `OZ066564.1 Rhizosphaera kalkhoffii genome assembly, chromosome: 1` 

```
grep '>' GCA_964106605.1_gdRhiKalk1.hap1.1_genomic.fna
>OZ066564.1 Rhizosphaera kalkhoffii genome assembly, chromosome: 1
>OZ066565.1 Rhizosphaera kalkhoffii genome assembly, chromosome: 2
>OZ066566.1 Rhizosphaera kalkhoffii genome assembly, chromosome: 3
>OZ066567.1 Rhizosphaera kalkhoffii genome assembly, chromosome: 4
>OZ066568.1 Rhizosphaera kalkhoffii genome assembly, chromosome: 5
>OZ066569.1 Rhizosphaera kalkhoffii genome assembly, chromosome: 6
>OZ066570.1 Rhizosphaera kalkhoffii genome assembly, chromosome: 7
>OZ066571.1 Rhizosphaera kalkhoffii genome assembly, chromosome: 8
>OZ066572.1 Rhizosphaera kalkhoffii genome assembly, chromosome: 9
>OZ066573.1 Rhizosphaera kalkhoffii genome assembly, chromosome: 10
>OZ066574.1 Rhizosphaera kalkhoffii genome assembly, chromosome: 11
>OZ066575.1 Rhizosphaera kalkhoffii genome assembly, chromosome: 12
>OZ066576.1 Rhizosphaera kalkhoffii genome assembly, chromosome: 13
>OZ066577.1 Rhizosphaera kalkhoffii genome assembly, chromosome: 14
>OZ066578.1 Rhizosphaera kalkhoffii genome assembly, organelle: mitochondrion
```

After stripping:

```
sed -i 's/>.*chromosome: />chr/g' GCA_964106605.1_gdRhiKalk1.hap1.1_genomic.fna

grep '>' GCA_964106605.1_gdRhiKalk1.hap1.1_genomic.fna
>chr1
>chr2
>chr3
>chr4
>chr5
>chr6
>chr7
>chr8
>chr9
>chr10
>chr11
>chr12
>chr13
>chr14
>OZ066578.1 Rhizosphaera kalkhoffii genome assembly, organelle: mitochondrion
```



## [Optional] Determine homozygous coverage & genome size with genomescope2

Using HiFi reads, we can determine homozygous coverage peaks, which also gives a sense of HiFi data quality. We ideally want enough HiFi data that we have a clean peak that is easily separated from the errors. 

```bash
module load apptainer
SAMPLE=Fungus
WD=/90daydata/coffea_pangenome/puzzler_trials/assemblies/kmers
HIFI=/90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/${SAMPLE}.HiFi.fastq.gz
SIF_PATH=/home/justin.merondun/apptainer/puzzler_v1.5.sif
PUZZLER="apptainer exec ${SIF_PATH}"
K=31

$PUZZLER FastK -v -t1 -k${K} -M70 -T3 -NFastK_Table_${SAMPLE} ${HIFI}
$PUZZLER Histex -G FastK_Table_${SAMPLE} > ${SAMPLE}.hist
$PUZZLER genomescope2 --input ${SAMPLE}.hist --output . --kmer_length ${K} --name_prefix ${SAMPLE}
rm FastK*
```

Will output this plot:

![genomescope2](/examples/figs/Fungus_linear_plot.png)

The left-peak corresponds to around 18, so we would set `--hom_cov` to 32 if it was diploid. 

It can sometimes be difficult to disentangle haploid and diploid samples using genomescope2, and we don't know the ploidy of this species. 

Instead of specifying `--hom_cov` directly, I will just let use `hifiasm` determine appropriate levels, and check the `.log`.

## `puzzler`

### Step 1: Draft Assembly

:bulb: The most important part is to prepare the map file `samples.tsv` with these columns, **in this specific order.**

:heavy_exclamation_mark:***Use full paths***!! 



* **sample:** Sample ID, all assembly work will be saved in `$WD/Fungus`.
* **runtime:** Either "apptainer", "singularity", or "conda". If other runtime, ensure `$runtime exec puzzler.sif` works.
* **container:** Path to the apptainer `.sif`. If all software available on path, simply write 'conda'. 
* **wd:** Path to working directory to store all files.
* **hifi:** Path to HiFi reads.
* **hic_r1:** Path to HiC R1.
* **hic_r2:** Path to HiC R2.
* **chromosomes:** Number of chromosomes, or best guess. The pipeline will attempt +/- 4 your estimate if unknown.

***OPTIONAL columns*** 
*Specify "NA" and the script will skip respective components.*

* **reference:** Path to related species genome for chromosome naming. Scaffolds will be renamed to the closest syntenic chromosome **using their scaffold naming convention**.
* **hom_cov:** Homozygous peak coverage.
* **blob_database:** Directory to save all blobtools databases.
* **busco_lineage:** Busco odb10 version lineage.
* **busco_database:** Directory to save busco dbs.



| sample | runtime   | container                                        | wd                                                    | hifi                                                         | hic_r1                                                       | hic_r2                                                       | num_chrs | reference                                                    | hom_cov | blob_database                                             | busco_lineage | busco_database                                             |
| ------ | --------- | ------------------------------------------------ | ----------------------------------------------------- | ------------------------------------------------------------ | ------------------------------------------------------------ | ------------------------------------------------------------ | -------- | ------------------------------------------------------------ | ------- | --------------------------------------------------------- | ------------- | ---------------------------------------------------------- |
| Fungus | apptainer | /home/justin.merondun/apptainer/puzzler_v1.7.sif | /90daydata/coffea_pangenome/puzzler_trials/assemblies | /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Fungus.HiFi.fastq.gz | /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Fungus.HiC.R1.fastq.gz | /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Fungus.HiC.R2.fastq.gz | 14       | /90daydata/coffea_pangenome/puzzler_trials/raw_data/references/GCA_964106605.1_gdRhiKalk1.hap1.1_genomic.fna | NA      | /90daydata/coffea_pangenome/puzzler_trials/blob_downloads | fungi_odb10   | /90daydata/coffea_pangenome/puzzler_trials/busco_downloads |

I set `hom_cov` to "NA", so it won't specify this for `hifiasm`, instead letting it select this level internally. 

Because the assembly and input files (400K HiFi reads; 310mb file size, 2M paired HiC reads; 373mb each) are small, I simply run this on a compute node with `--threads 4`and `--mem 32`  Gb of memory. 

```bash
=======================================================================
__________ ____ _______________________.____     _____________________
\______   \    |   \____    /\____    /|    |    \_   _____/\______   \
 |     ___/    |   / /     /   /     / |    |     |    __)_  |       _/
 |    |   |    |  / /     /_  /     /_ |    |___  |        \ |    |   \
 |____|   |______/ /_______ \/_______ \|_______ \/_______  / |____|_  /
                           \/        \/        \/        \/         \/
=======================================================================

=======================================================================
Parameters for sample: Fungus
RUNTIME: apptainer
CONTAINER: /home/justin.merondun/apptainer/puzzler_v1.7.sif
WD: /project/90daydata/coffea_pangenome/puzzler_trials/assemblies
HIFI: /project/90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Fungus.HiFi.fastq.gz
HIC_R1: /project/90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Fungus.HiC.R1.fastq.gz
HIC_R2: /project/90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Fungus.HiC.R2.fastq.gz
NUMBER CHRS: 14
REFERENCE: /project/90daydata/coffea_pangenome/puzzler_trials/raw_data/references/GCA_964106605.1_gdRhiKalk1.hap1.1_genomic.fna
HOM_COV: NA
BLOB_DB: /project/90daydata/coffea_pangenome/puzzler_trials/blob_downloads
BUSCO_LINEAGE: fungi_odb10
BUSCO_DB: /project/90daydata/coffea_pangenome/puzzler_trials/busco_downloads
Cores Requested: 4
Cores Available: 4
RAM Requested: 32
Memory Available: 327.3 GB
=======================================================================

~~~~ Assembling genome for Fungus ~~~~
~~~~ Starting hifiasm assembly for Fungus ~~~~
```

:alarm_clock: This is a good time to check that all of your paths and parameters look appropriate. 

 We can then check the output `Fungus/01_hifiasm/Fungus.hifiasm.log` file to ensure an appropriate `--hom_cov` was selected:

```bash
 1: ****************************************************************************************************> 905745
 2: ****************************************************************************************************> 63784
 3: ********************************** 12624
 4: ************** 5273
 5: ******** 2871
 6: **** 1643
 7: *** 1091
 8: ** 842
 9: ** 677
10: * 383
11: * 306
12: * 298
13: * 413
14: * 443
15: ** 577
16: ** 904
17: *** 1281
18: ***** 1752
19: ****** 2310
20: ********* 3375
21: ************ 4347
22: **************** 5891
23: ********************** 8012
24: **************************** 10458
25: *********************************** 12942
26: ***************************************** 14919
27: ************************************************** 18171
28: ********************************************************* 20927
29: *************************************************************** 22993
30: ********************************************************************** 25613
31: ******************************************************************************* 28887
32: ************************************************************************************* 31331
33: ******************************************************************************************* 33257
34: ************************************************************************************************ 35302
35: **************************************************************************************************** 36699
36: ************************************************************************************************** 36086
37: *********************************************************************************************** 34944
38: ********************************************************************************************** 34616
39: *************************************************************************************** 31801
40: ************************************************************************************ 30938
41: ******************************************************************************* 28950
42: ************************************************************************** 27087
43: ****************************************************************** 24120
44: ********************************************************** 21313
45: *********************************************** 17428
46: **************************************** 14650
47: ********************************** 12535
48: ***************************** 10638
49: ************************* 9016
50: ******************** 7335
51: *************** 5624
52: ************ 4336
53: ********* 3417
54: ******* 2709
55: ***** 1883
56: **** 1415
57: *** 1273
58: *** 981
59: ** 821
60: ** 719
61: ** 606
62: * 521
63: * 476
64: * 458
65: * 450
66: * 507
67: * 450
68: * 467
69: * 401
70: * 325
71: * 401
72: * 327
73: * 408
74: * 312
75: * 314
76: * 324
77: * 373
78: * 313
79: * 294
80: * 248
81: * 271
82: * 342
83: * 228
84: * 289
85: * 349
86: * 338
87: * 318
88: * 311
89: * 299
90: * 236
91: * 261
92: * 222
93: * 230
94: * 215
95: * 208
[M::ha_hist_line]  rest: ************************ 8842
[M::ha_analyze_count] left: none
[M::ha_analyze_count] right: none
[M::ha_pt_gen] peak_hom: 35; peak_het: -1
```

`hifiasm` automatically selected 35, while `genomescope2` estimated left peak at 18.1 (if assuming diploid, `--hom_cov` =36).

Puzzler continues up until the juicebox curation stage. The output indicates that puzzler was unable to successfully run [HapHiC](https://github.com/zengxiaofei/HapHiC) with the specified number of chromosomes (from `samples.tsv`). Puzzler will automatically atttempt HapHiC with +/- 4 chromosomes, to see if it will successfully run with those chromosomes numbers instead. If that fails, puzzler will instead run [YAHS](https://github.com/c-zhou/yahs), which in my experience has always run successfully. 

Of the many genomes I have run so far, Fungus is the only genome not successfully scaffolded with HapHiC, so maybe it has something to do with the small size of chromosomes. 

```bash
=======================================================================
__________ ____ _______________________.____     _____________________ 
\______   \    |   \____    /\____    /|    |    \_   _____/\______   \
 |     ___/    |   / /     /   /     / |    |     |    __)_  |       _/
 |    |   |    |  / /     /_  /     /_ |    |___  |        \ |    |   \
 |____|   |______/ /_______ \/_______ \|_______ \/_______  / |____|_  /
                           \/        \/        \/        \/         \/ 
=======================================================================

=======================================================================
Parameters for sample: Fungus 
CONTAINER: /home/justin.merondun/apptainer/puzzler_v1.6.sif 
WD: /90daydata/coffea_pangenome/puzzler_trials/assemblies 
PLOIDY: 1 
NUMBER CHRS: 14
HIFI: /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Fungus.HiFi.fastq.gz
HIC_R1: /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Fungus.HiC.R1.fastq.gz
HIC_R2: /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Fungus.HiC.R2.fastq.gz
REFERENCE: /90daydata/coffea_pangenome/puzzler_trials/raw_data/references/GCA_964106605.1_gdRhiKalk1.hap1.1_genomic.fna
HOM_COV: NA
RUNTIME: apptainer
=======================================================================

~~~~ Starting hifiasm assembly for Fungus ~~~~
~~~~ Starting Purge_Dups for Fungus ~~~~
~~~~ Mapping HiC reads to Fungus draft ~~~~
~~~~ Running HapHiC for Fungus  ~~~~
~~~~ HapHiC for Fungus with 14 chrs failed, trying: 10 ~~~~
~~~~ HapHiC for Fungus with 14 chrs failed, trying: 11 ~~~~
~~~~ HapHiC for Fungus with 14 chrs failed, trying: 12 ~~~~
~~~~ HapHiC for Fungus with 14 chrs failed, trying: 13 ~~~~
~~~~ HapHiC for Fungus with 14 chrs failed, trying: 14 ~~~~
~~~~ HapHiC for Fungus with 14 chrs failed, trying: 15 ~~~~
~~~~ HapHiC for Fungus with 14 chrs failed, trying: 16 ~~~~
~~~~ HapHiC for Fungus with 14 chrs failed, trying: 17 ~~~~
~~~~ HapHiC for Fungus with 14 chrs failed, trying: 18 ~~~~
~~~~ HapHiC for Fungus failed, scaffolding with YAHS instead ~~~~
~~~~ Creating .hic file for juicebox for Fungus  ~~~~
~~~~ Post curation assembly file doesn't exist for Fungus - Run Juicebox! ~~~~
```

In any case, the pipeline will also run `haphic refsort`, provided a related species is provided, which will align your draft scaffolds in order of the reference genome for juicebox curation. **This is not reference-guided assembly**, and will not rename or otherwise modify the scaffolds. Please see the [HapHiC](https://github.com/zengxiaofei/HapHiC) `refsort` description for more details. 

If you did not provide a `$REFERENCE` file, it will still output the `.hic` and `.assembly` files, just without running `haphic refsort`. 

### Step 2: Juicebox Curation

Now we are ready to load the `$WD/juicer_files/Fungus_JBAT.hic` and `$WD/juicer_files/Fungus_JBAT.assembly` file into juicebox. 

Perform our edits. Pull in the `.hic` file, import the `.assembly` file with `Assembly > Import Map Assembly`. There's very minimal edits on this genome, I just merge two blocks which originate from one chromosome: 

![Fungus_Edits_Juicebox]()

Afterwards, save the file with `Assembly > Export Assembly`, and save the  `Fungus_JBAT.review.assembly` to the directory `$WD/juicer_files`, retaining the default file name! 

### Step 3: Chromosome Renaming, Remapping HiC

Simply submit `puzzler` with the same parameters as before, and the script will assign chromosome names, output a final HiC contact map on the final assembly, and output assembly statistics:

```bash
puzzler --sample Fungus --map samples.tsv 

=======================================================================
__________ ____ _______________________.____     _____________________ 
\______   \    |   \____    /\____    /|    |    \_   _____/\______   \
 |     ___/    |   / /     /   /     / |    |     |    __)_  |       _/
 |    |   |    |  / /     /_  /     /_ |    |___  |        \ |    |   \
 |____|   |______/ /_______ \/_______ \|_______ \/_______  / |____|_  /
                           \/        \/        \/        \/         \/ 
=======================================================================

=======================================================================
Parameters for sample: Fungus 
CONTAINER: /home/justin.merondun/apptainer/puzzler_v1.6.sif 
WD: /90daydata/coffea_pangenome/puzzler_trials/assemblies 
PLOIDY: 1 
NUMBER CHRS: 14
HIFI: /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Fungus.HiFi.fastq.gz
HIC_R1: /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Fungus.HiC.R1.fastq.gz
HIC_R2: /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Fungus.HiC.R2.fastq.gz
REFERENCE: /90daydata/coffea_pangenome/puzzler_trials/raw_data/references/GCA_964106605.1_gdRhiKalk1.hap1.1_genomic.fna
HOM_COV: 
RUNTIME: apptainer
=======================================================================

~~~~ Skipping hifiasm for Fungus: /90daydata/coffea_pangenome/puzzler_trials/assemblies/Fungus/01_hifiasm/asm.hic.p_ctg.gfa exists ~~~~
~~~~ Skipping purge for Fungus: /90daydata/coffea_pangenome/puzzler_trials/assemblies/Fungus/02_purge_dups/p_ctg.purged.fa exists ~~~~
~~~~ Skipping HiC alignment for Fungus: /90daydata/coffea_pangenome/puzzler_trials/assemblies/Fungus/03_haphic/filtered.bam exists ~~~~
~~~~ Skipping HapHiC for Fungus: /90daydata/coffea_pangenome/puzzler_trials/assemblies/Fungus/03_haphic/haphic/04.build/scaffolds.fa exists ~~~~
~~~~ Skipping juicer HiC file creation for Fungus: /90daydata/coffea_pangenome/puzzler_trials/assemblies/primary_asm/juicer_files/Fungus_JBAT.hic exists ~~~~
~~~~ Skipping draft-reference mapping for Fungus: /90daydata/coffea_pangenome/puzzler_trials/assemblies/Fungus/05_postjuicebox/map.txt exists ~~~~
~~~~ Renaming chromosomes for Fungus ~~~~
~~~~ Scaffold sanity check passed for renaming, proceeding! ~~~~
~~~~ Single scaffolds corresponding to a single Chr for Fungus ~~~~
~~~~ Creating final HiC bam for Fungus ~~~~
```

Which outputs a final HiC contact map at `$WD/primary_asm/stats/$SAMPLE.pdf`: 

![final hic contacts](/examples/figs/Fungus_Final.png)

### Step 3: If No `$REFERENCE` 

:bomb: If you **DO NOT** have a somewhat-related species to assign chromosomes (note that we only really need weak alignments, I have done this with ~60 MY diverged species), then you can be finished after juicebox curation. 

After Juicebox curation, just run this, which will extract the final assembly:

```bash
SAMPLE=Fungus
WD=/90daydata/coffea_pangenome/puzzler_trials/assemblies
PUZZLER="apptainer exec /home/justin.merondun/apptainer/puzzler_v1.6.sif"

cd ${WD}/${SAMPLE}/05_postjuicebox
${PUZZLER} juicer post \
    -o haphic-refsort-post_JBAT \
    ${WD}/primary_asm/juicer_files/${SAMPLE}_JBAT.review.assembly \
    ${WD}/${SAMPLE}/04_juicer/haphic-refsort_JBAT.liftover.agp \
    ${WD}/${SAMPLE}/02_purge_dups/p_ctg.purged.fa 2> ${SAMPLE}.juicer.post.log
```

Which will output: `haphic-refsort-post_JBAT.FINAL.fa`. Your final assembly. 