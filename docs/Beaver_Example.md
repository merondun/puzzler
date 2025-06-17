# Beaver Genome Assembly with **`puzzler`**

[Beaver](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_047511655.1/) has a rather large genome, 19 autosomes and assembled X for female. 

## File Preparation

### Reads

First, grab the HiFi and HiC data from GenomeArk. 

```bash
#!/bin/bash

#SBATCH --time=2-00:00:00   
#SBATCH --nodes=1  
#SBATCH --ntasks-per-node=1
#SBATCH --partition=ceres

WD=/90daydata/coffea_pangenome/puzzler_trials

mkdir Castor_canadensis
cd Castor_canadensis
wget https://genomeark.s3.amazonaws.com/species/Castor_canadensis/mCasCan1/genomic_data/pacbio_hifi/m84091_240920_000814_s4.hifi_reads.bc1018.fastq.gz
wget https://genomeark.s3.amazonaws.com/species/Castor_canadensis/mCasCan1/genomic_data/pacbio_hifi/m84091_240927_164339_s1.hifi_reads.bc1018.fastq.gz
wget https://genomeark.s3.amazonaws.com/species/Castor_canadensis/mCasCan1/genomic_data/arima/mCasCan1_L1_S5_R1_001.fastq.gz
wget https://genomeark.s3.amazonaws.com/species/Castor_canadensis/mCasCan1/genomic_data/arima/mCasCan1_L1_S5_R2_001.fastq.gz
```

### [Optional] Simplify Reference Chr Names

I will name the chromosomes according to the published reference assembly:

```bash
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/047/511/655/GCF_047511655.1_mCasCan1.hap1v2/GCF_047511655.1_mCasCan1.hap1v2_genomic.fna.gz
gunzip GCF_047511655.1_mCasCan1.hap1v2_genomic.fna.gz
```

I don't want the re-named chromosomes to be this cryptic, so I will strip all the text between `>` and `chromosome :`.

After stripping:

```
# Strip text before and after 'chromosome N'
sed -i 's/>.*chromosome />chr/g' GCF_047511655.1_mCasCan1.hap1v2_genomic.fna
sed -i 's/,.*//g' GCF_047511655.1_mCasCan1.hap1v2_genomic.fna

grep '>' GCF_047511655.1_mCasCan1.hap1v2_genomic.fna
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
```

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



| sample | runtime   | container                                        | wd                                                    | hifi                                                         | hic_r1                                                       | hic_r2                                                       | num_chrs | reference                                                    | hom_cov | blob_database                                             | busco_lineage  | busco_database                                             |
| ------ | --------- | ------------------------------------------------ | ----------------------------------------------------- | ------------------------------------------------------------ | ------------------------------------------------------------ | ------------------------------------------------------------ | -------- | ------------------------------------------------------------ | ------- | --------------------------------------------------------- | -------------- | ---------------------------------------------------------- |
| Beaver | apptainer | /home/justin.merondun/apptainer/puzzler_v1.7.sif | /90daydata/coffea_pangenome/puzzler_trials/assemblies | /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Beaver.HiFi.fastq.gz | /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Beaver.HiC.R1.fastq.gz | /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Beaver.HiC.R2.fastq.gz | 20       | /90daydata/coffea_pangenome/puzzler_trials/raw_data/references/GCF_047511655.1_mCasCan1.hap1v2_genomic.fna | NA      | /90daydata/coffea_pangenome/puzzler_trials/blob_downloads | mammalia_odb10 | /90daydata/coffea_pangenome/puzzler_trials/busco_downloads |

I set `hom_cov` to "NA", so it won't specify this for `hifiasm`, instead letting it select this level internally. 

```bash
sbatch -J asm_Beaver puzzler --sample Beaver --map samples.tsv
...
cat slurm-15604116.out

=======================================================================
__________ ____ _______________________.____     _____________________ 
\______   \    |   \____    /\____    /|    |    \_   _____/\______   \
 |     ___/    |   / /     /   /     / |    |     |    __)_  |       _/
 |    |   |    |  / /     /_  /     /_ |    |___  |        \ |    |   \
 |____|   |______/ /_______ \/_______ \|_______ \/_______  / |____|_  /
                           \/        \/        \/        \/         \/ 
=======================================================================

=======================================================================
Parameters for sample: Beaver 
RUNTIME: apptainer
CONTAINER: /home/justin.merondun/apptainer/puzzler_v1.7.sif 
WD: /90daydata/coffea_pangenome/puzzler_trials/assemblies 
HIFI: /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Beaver.HiFi.fastq.gz
HIC_R1: /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Beaver.HiC.R1.fastq.gz
HIC_R2: /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Beaver.HiC.R2.fastq.gz
NUMBER CHRS: 20
REFERENCE: /90daydata/coffea_pangenome/puzzler_trials/raw_data/references/GCF_047511655.1_mCasCan1.hap1v2_genomic.fna
HOM_COV: NA
BLOB_DB: /90daydata/coffea_pangenome/puzzler_trials/blob_downloads
BUSCO_LINEAGE: mammalia_odb10
BUSCO_DB: /90daydata/coffea_pangenome/puzzler_trials/busco_downloads
Cores Requested: 64
Cores Available: 64
RAM Requested: 512
Memory Available: 512.0 GB
=======================================================================

~~~~ Starting hifiasm assembly for Beaver ~~~~
~~~~ Starting Purge_Dups for Beaver ~~~~
~~~~ Mapping HiC reads to Beaver draft ~~~~
~~~~ Running HapHiC for Beaver  ~~~~
~~~~ Creating .hic file for juicebox for Beaver, with reference alignment  ~~~~
~~~~ Post curation assembly file doesn't exist for Beaver - Run Juicebox! ~~~~
```

:alarm_clock: This is a good time to check that all of your paths and parameters look appropriate. 

 We can then check the output `Beaver/01_hifiasm/Beaver.hifiasm.log` file to ensure an appropriate `--hom_cov` was selected:

```bash
  1: **************************** 476915
  2: ***** 86664
  3: *** 55124
  4: ** 35861
  5: ** 26443
  6: * 21344
  7: * 18989
  8: * 18888
  9: * 20222
 10: * 23206
 11: ** 28461
 12: ** 36350
 13: *** 46244
 14: *** 57888
 15: **** 72688
 16: ***** 91997
 17: ******* 113428
 18: ******** 137842
 19: ********* 161804
 20: *********** 191264
 21: ************* 226008
 22: *************** 259766
 23: ***************** 298037
 24: ******************** 341507
 25: ********************** 387685
 26: ************************* 431145
 27: **************************** 477356
 28: ****************************** 525676
 29: ********************************* 573187
 30: ************************************ 614176
 31: ************************************** 652660
 32: **************************************** 684777
 33: ***************************************** 711899
 34: ******************************************* 735799
 35: ******************************************** 753431
 36: ******************************************** 764911
 37: ********************************************* 769853
 38: ******************************************** 766664
 39: ******************************************** 761691
 40: ******************************************** 760897
 41: ******************************************** 763545
 42: ******************************************** 763706
 43: ******************************************** 765086
 44: ********************************************* 770508
 45: ********************************************* 778444
 46: ********************************************** 799132
 47: *********************************************** 815358
 48: ************************************************* 837651
 49: ************************************************** 866006
 50: ***************************************************** 907093
 51: ******************************************************* 947103
 52: ********************************************************** 998091
 53: ************************************************************* 1045608
 54: **************************************************************** 1105535
 55: ******************************************************************** 1164970
 56: *********************************************************************** 1218954
 57: ************************************************************************** 1271892
 58: ****************************************************************************** 1337123
 59: ********************************************************************************* 1396242
 60: ************************************************************************************ 1455921
 61: **************************************************************************************** 1517596
 62: ******************************************************************************************* 1564817
 63: ********************************************************************************************* 1611250
 64: ************************************************************************************************ 1650964
 65: ************************************************************************************************** 1684947
 66: *************************************************************************************************** 1707885
 67: **************************************************************************************************** 1721196
 68: **************************************************************************************************** 1724779
 69: **************************************************************************************************** 1722429
 70: *************************************************************************************************** 1710780
 71: ************************************************************************************************** 1695993
 72: ************************************************************************************************* 1666265
 73: ********************************************************************************************** 1629205
 74: ******************************************************************************************** 1578837
 75: **************************************************************************************** 1522875
 76: ************************************************************************************* 1462132
 77: ********************************************************************************* 1393915
 78: ***************************************************************************** 1320491
 79: ************************************************************************* 1251732
 80: ******************************************************************** 1165715
 81: *************************************************************** 1086583
 82: ********************************************************** 1006928
 83: ****************************************************** 926857
 84: ************************************************* 850094
 85: ********************************************* 772947
 86: ***************************************** 704671
 87: ************************************* 635036
 88: ********************************* 567381
 89: ***************************** 505765
 90: ************************** 447555
 91: *********************** 394188
 92: ******************** 345559
 93: ****************** 303634
 94: *************** 265141
 95: ************* 231366
 96: ************ 201452
 97: ********** 171954
 98: ********* 148171
 99: ******* 128251
100: ****** 111167
101: ****** 98421
102: ***** 86603
103: **** 74068
104: **** 66567
105: *** 59049
106: *** 51351
107: *** 45144
108: ** 39650
109: ** 36851
110: ** 34180
111: ** 32281
112: ** 30613
113: ** 29733
114: ** 27750
115: ** 26319
116: ** 26536
117: * 25376
118: * 25385
119: * 25121
120: * 24975
121: * 24862
122: * 24327
123: * 24077
124: * 24168
125: * 24281
126: * 24406
127: * 24542
128: * 24800
129: * 24625
130: * 24782
131: * 24794
132: * 24593
133: * 24566
134: * 24570
135: * 24419
136: * 24713
137: * 24674
138: * 24768
139: * 24567
140: * 24294
141: * 24101
142: * 23821
143: * 23184
144: * 23007
145: * 22352
146: * 21840
147: * 21744
148: * 21037
149: * 20786
150: * 19753
151: * 19562
152: * 19277
153: * 18707
154: * 17750
155: * 16942
156: * 16240
157: * 16264
158: * 15460
159: * 14955
160: * 14414
161: * 13999
162: * 13442
163: * 13380
164: * 12926
165: * 12598
166: * 11881
167: * 11939
168: * 11788
169: * 11136
170: * 10851
171: * 10504
172: * 10203
173: * 9858
174: * 9795
175: * 9400
176: * 9013
177: * 8850
[M::ha_hist_line]  rest: ***************************************** 710416
[M::ha_analyze_count] left: none
[M::ha_analyze_count] right: none
[M::ha_pt_gen] peak_hom: 68; peak_het: -1
```

`hifiasm` automatically selected 68 (`peak_hom: 68`), which is right on the diploid peak. 

In any case, the pipeline will also run `haphic refsort`, provided a related species is provided, which will align your draft scaffolds in order of the reference genome for juicebox curation. **This is not reference-guided assembly**, and will not rename or otherwise modify the scaffolds. Please see the [HapHiC](https://github.com/zengxiaofei/HapHiC) `refsort` description for more details. 

If you did not provide a `$REFERENCE` file, it will still output the `.hic` and `.assembly` files, just without running `haphic refsort`. 

### Step 2: Juicebox Curation

Now we are ready to load the `$WD/juicer_files/Beaver_JBAT.hic` and `$WD/juicer_files/Beaver_JBAT.assembly` file into juicebox. 

Perform our edits. Pull in the `.hic` file, import the `.assembly` file with `Assembly > Import Map Assembly`. There's very minimal edits on this genome, I just merge two blocks which originate from one chromosome: 

![Fungus_Edits_Juicebox]()

Afterwards, save the file with `Assembly > Export Assembly`, and save the  `Beaver_JBAT.review.assembly` to the directory `$WD/juicer_files`, retaining the default file name! 

### Step 3: Chromosome Renaming, QC

Simply submit `puzzler` with the same parameters as before, and the script will assign chromosome names, output a final HiC contact map on the final assembly, and output assembly statistics:

```bash
puzzler --sample Beaver --map samples.tsv 

cat slurm-15672866.out

=======================================================================
__________ ____ _______________________.____     _____________________ 
\______   \    |   \____    /\____    /|    |    \_   _____/\______   \
 |     ___/    |   / /     /   /     / |    |     |    __)_  |       _/
 |    |   |    |  / /     /_  /     /_ |    |___  |        \ |    |   \
 |____|   |______/ /_______ \/_______ \|_______ \/_______  / |____|_  /
                           \/        \/        \/        \/         \/ 
=======================================================================

=======================================================================
Parameters for sample: Beaver 
RUNTIME: apptainer
CONTAINER: /home/justin.merondun/apptainer/puzzler_v1.7.sif 
WD: /90daydata/coffea_pangenome/puzzler_trials/assemblies 
HIFI: /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Beaver.HiFi.fastq.gz
HIC_R1: /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Beaver.HiC.R1.fastq.gz
HIC_R2: /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Beaver.HiC.R2.fastq.gz
NUMBER CHRS: 20
REFERENCE: /90daydata/coffea_pangenome/puzzler_trials/raw_data/references/GCF_047511655.1_mCasCan1.hap1v2_genomic.fna
HOM_COV: NA
BLOB_DB: /90daydata/coffea_pangenome/puzzler_trials/blob_downloads
BUSCO_LINEAGE: mammalia_odb10
BUSCO_DB: /90daydata/coffea_pangenome/puzzler_trials/busco_downloads
Cores Requested: 64
Cores Available: 64
RAM Requested: 512
Memory Available: 512.0 GB
=======================================================================

~~~~ Assembling genome for Beaver ~~~~
~~~~ Skipping hifiasm for Beaver: /90daydata/coffea_pangenome/puzzler_trials/assemblies/Beaver/01_hifiasm/asm.hic.p_ctg.gfa exists ~~~~
~~~~ Skipping purge for Beaver: /90daydata/coffea_pangenome/puzzler_trials/assemblies/Beaver/02_purge_dups/p_ctg.purged.fa exists ~~~~
~~~~ Skipping HiC alignment for Beaver: /90daydata/coffea_pangenome/puzzler_trials/assemblies/Beaver/03_haphic/filtered.bam exists ~~~~
~~~~ Skipping HapHiC for Beaver: /90daydata/coffea_pangenome/puzzler_trials/assemblies/Beaver/03_haphic/haphic/04.build/scaffolds.fa exists ~~~~
~~~~ Skipping juicer HiC file creation for Beaver: /90daydata/coffea_pangenome/puzzler_trials/assemblies/juicer_files/Beaver_JBAT.hic exists ~~~~
~~~~ Extracting post-curation assembly and mapping to reference for Beaver ~~~~
~~~~ Genome size is 2GB, running default minimap2 command ~~~~
~~~~ Renaming chromosomes for Beaver ~~~~
~~~~ Scaffold sanity check passed for renaming, proceeding! ~~~~
~~~~ Multiple scaffolds corresponding to single Chr for Beaver, Renaming them e.g. Chr1A, Chr1B.. ~~~~
~~~~ Assessing genome quality for Beaver ~~~~
~~~~ Creating final HiC bam for Beaver ~~~~
~~~~ Creating final contact map for Beaver ~~~~
~~~~ Starting HiFi Alignment for Beaver ~~~~
~~~~ Running BUSCO for Beaver using lineage: /90daydata/coffea_pangenome/puzzler_trials/assemblies/mammalia_odb10 ~~~~
```

There is one output message in red font here: **"Multiple scaffolds corresponding to single Chr for Beaver, Renaming them e.g. Chr1A, Chr1B.."**. This message means that there are multiple scaffolds in our draft which correspond to a single chromosome in the reference. This will mainly happen when there's a smaller section of a chromosome which isn't scaffolded with the main sequence. 

`map_chromosomes` extracts synteny from the minimap2 `.paf` file, and assigns likely chromosomes based on the highest percent match against the reference. It will assign the largest scaffold in your assembly as the main chromosome (e.g. Chr1), and then all the smaller scaffolds which also have synteny to Chr1 as Chr1A, Chr1B, etc. 

`puzzler` continues and maps HiC and HiFi reads, outputs the final HiC contact map, runs BUSCO according to the specified `$BUSCO_LINEAGE` (mammalia_odb10), and runs yak to assess base quality. 



Which outputs a final HiC contact map at `$WD/primary_asm/stats/$SAMPLE.pdf`: 

![final hic contacts](/examples/figs/Beaver_Final.png)


