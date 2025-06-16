# Needlecast fungus Genome Assembly with **`puzzler`**

[Needlecast fungus](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_964106605.1/) has a small genome (25 Mb) and 14 chromosomes, so it's good for a quick tutorial. 

## File Preparation

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



## Determine homozygous coverage with genomescope2

Using HiFi reads, determine homozygous coverage peaks:

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

The left-peak corresponds to around 14, so we would set `--hom_cov` to 28. Instead, I will just let use hifiasm determine appropriate levels, and check the `.log`.

## Draft Assembly

The most important part is to prepare the map file `samples.tsv` with these columns, in this specific order. 

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
      |

I left the `hom_cov` column empty, so it won't specify this for `hifiasm`, instead letting hifiasm select this level internally. 

Because the assembly and input files are small, I simply run this on a compute node with 4 cores and 36 Gb of memory:

```bash
puzzler --sample Fungus --map samples.tsv --threads 4 --mem 36

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
CONTAINER: /home/justin.merondun/apptainer/puzzler_v1.5.sif 
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

~~~~ Starting hifiasm assembly for Fungus ~~~~
~~~~ Starting Purge_Dups for Fungus ~~~~
~~~~ Skipping purge for Fungus: no duplicates found! ~~~~
~~~~ Mapping HiC reads to Fungus draft ~~~~
~~~~ Running HapHiC for Fungus  ~~~~
```

:alarm_clock: This is a good time to check that all of your paths and parameters look appropriate. 

 We can then check the output `Fungus/01_hifiasm/Fungus.hifiasm.log` file to ensure an appropriate `--hom_cov` was selected:

```bash
[M::ha_hist_line]     1: ****************************************************************************************************> 593244
[M::ha_hist_line]     2: ****************************************************************************************************> 40550
[M::ha_hist_line]     3: ************************************************ 12704
[M::ha_hist_line]     4: ************************************** 9994
[M::ha_hist_line]     5: *********************************** 9165
[M::ha_hist_line]     6: ****************************************** 10998
[M::ha_hist_line]     7: ******************************************* 11210
[M::ha_hist_line]     8: ******************************************* 11324
[M::ha_hist_line]     9: *************************************************** 13249
[M::ha_hist_line]    10: ************************************************************ 15787
[M::ha_hist_line]    11: ************************************************************* 15949
[M::ha_hist_line]    12: ***************************************************************** 16971
[M::ha_hist_line]    13: *********************************************************************** 18676
[M::ha_hist_line]    14: ******************************************************************************* 20572
[M::ha_hist_line]    15: *********************************************************************************** 21659
[M::ha_hist_line]    16: **************************************************************************************** 23154
[M::ha_hist_line]    17: ******************************************************************************************** 24120
[M::ha_hist_line]    18: ********************************************************************************************** 24568
[M::ha_hist_line]    19: ************************************************************************************************* 25442
[M::ha_hist_line]    20: **************************************************************************************************** 26141
[M::ha_hist_line]    21: **************************************************************************************************** 26206
[M::ha_hist_line]    22: ************************************************************************************************ 25169
[M::ha_hist_line]    23: ********************************************************************************************** 24528
[M::ha_hist_line]    24: ************************************************************************************************* 25467
[M::ha_hist_line]    25: ********************************************************************************************** 24699
[M::ha_hist_line]    26: ******************************************************************************************* 23770
[M::ha_hist_line]    27: ************************************************************************************** 22573
[M::ha_hist_line]    28: ******************************************************************************* 20825
[M::ha_hist_line]    29: ******************************************************************************* 20710
[M::ha_hist_line]    30: ************************************************************************* 19117
[M::ha_hist_line]    31: *********************************************************************** 18523
[M::ha_hist_line]    32: ******************************************************************* 17560
[M::ha_hist_line]    33: ******************************************************************* 17552
[M::ha_hist_line]    34: ********************************************************* 15052
[M::ha_hist_line]    35: ************************************************** 13078
[M::ha_hist_line]    36: ************************************************ 12449
[M::ha_hist_line]    37: ***************************************** 10692
[M::ha_hist_line]    38: *************************************** 10336
[M::ha_hist_line]    39: *********************************** 9066
[M::ha_hist_line]    40: ****************************** 7739
[M::ha_hist_line]    41: ************************ 6334
[M::ha_hist_line]    42: ******************* 5061
[M::ha_hist_line]    43: ************** 3797
[M::ha_hist_line]    44: ************ 3082
[M::ha_hist_line]    45: ********** 2711
[M::ha_hist_line]    46: ******* 1795
[M::ha_hist_line]    47: ***** 1200
[M::ha_hist_line]    48: ***** 1342
[M::ha_hist_line]    49: **** 952
[M::ha_hist_line]    50: *** 686
[M::ha_hist_line]    51: ** 541
[M::ha_hist_line]    52: * 393
[M::ha_hist_line]    53: * 266
[M::ha_hist_line]    54: * 238
[M::ha_hist_line]    55: * 163
[M::ha_hist_line]    56: * 138
[M::ha_hist_line]  rest: ******************* 4852
[M::ha_analyze_count] left: none
[M::ha_analyze_count] right: none
[M::ha_pt_gen] peak_hom: 21; peak_het: -1
```

`hifiasm` automatically selected 21, instead of `genomescope2` estimate of 28, which is actually a bigger gap than I've seen with larger genomes. 

We successfully got to the `haphic pipeline` step, but now we ran into an issue because haphic could not identify our specified number of chromosomes (n=14).

```bash
#!/bin/bash

#SBATCH --time=48:00:00   
#SBATCH --nodes=1  
#SBATCH --ntasks-per-node=10
#SBATCH --mem=56Gb
#SBATCH --partition=ceres




if [ ! -s ${WD}/${SAMPLE}/02_${IT}HapHiC/01_haphicMQ1/04.build/scaffolds.fa ]; then
rm -rf ${WD}/${SAMPLE}/02_${IT}HapHiC/01_haphicMQ1/*

    # Attempt NUM_CHRS values from 2 to 40
    for CHR_ATTEMPT in {2..40}; do

        echo -e "\e[43m~~~~ Running HapHiC for ${SAMPLE} ${IT} ${CHR_ATTEMPT} ~~~~\e[0m"
        cd ${WD}/${SAMPLE}/02_${IT}HapHiC/01_haphicMQ1

        ${PUZZLER} haphic pipeline ../all.purged.fa ../filtered.MQ1.bam ${CHR_ATTEMPT} --remove_allelic_links ${PLOIDY} --correct_nrounds 2 --max_inflation 20.0 --threads ${t} --processes ${t} 2> ${WD}/${SAMPLE}/02_${IT}HapHiC/${SAMPLE}.${IT}.${CHR_ATTEMPT}.haphic.LOOP.log

        # Check if scaffolds.fa exists
        if [ -s ${WD}/${SAMPLE}/02_${IT}HapHiC/01_haphicMQ1/04.build/scaffolds.fa ]; then
            echo -e "\e[32mScaffolds.fa created successfully with NUM_CHRS=${CHR_ATTEMPT}\e[0m"
            break
        else
            echo -e "\e[31mFailed with NUM_CHRS=${CHR_ATTEMPT}, retrying...\e[0m"
            rm -rf ${WD}/${SAMPLE}/02_${IT}HapHiC/01_haphicMQ1/*
        fi
    done
fi
```





Otherwise, now it's load the `$WD/primary_asm/juicer_files/Fungus_JBAT.hic` and `$WD/primary_asm/juicer_files/Fungus_JBAT.assembly` file into juicebox. 





