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

### Simplify Reference Chr Names

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

Now, all of the scaffolds matching these chromosomes based on synteny will just be labled 'chr1', 'chr2'. 

## [Optional] Determine homozygous coverage & genome size with genomescope2

Using HiFi reads, we can determine homozygous coverage peaks, which also gives a sense of HiFi data quality. We ideally want enough HiFi data that we have a clean peak that is easily separated from the errors. 

```bash
module load apptainer
SAMPLE=Fungus
WD=/90daydata/coffea_pangenome/puzzler_trials/assemblies/kmers
HIFI=/90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/${SAMPLE}.HiFi.fastq.gz
SIF_PATH=/home/justin.merondun/apptainer/puzzler_v1.8.sif
# Note that depending on your architecture, you might need to bind all the paths. Run puzzler once with your samples.tsv to learn the full Puzzler command, e.g. apptainer exec --bind /90daydata/coffea_pangenome/puzzler_trials:/90daydata/coffea_pangenome/puzzler_trials --bind /90daydata/coffea_pangenome/puzzler_trials/assemblies:/90daydata/coffea_pangenome/puzzler_trials/assemblies --bind /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads:/90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads --bind /home/justin.merondun/apptainer:/home/justin.merondun/apptainer /home/justin.merondun/apptainer/puzzler_v1.8.sif FastK

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

## `puzzler` Workflow

I will analyze 4 analytical forks for the Fungus sample:

1) With `juicer` manual curation.
2) Same with (1), except without a `$REFERENCE` relative species.
3) Without `jucier` manual curation. 
4) Same with (3), except without a `$REFERENCE` relative species.

This will give us 4 assemblies to compare, and show what the pipeline can do. 

:stop_sign: Analyses (1) and (2) will stop at the manual curation step, giving us the `.hic` and `.assembly` file within `$WD/juicer_files` to load into juicebox. 

:rabbit: Analyses (3) and (4) go all the way to the end, one-command genome assembly and QC. **This is not recommended, as manual curation is essential for checking misassemblies!** 

_____

### Inputs

**Puzzler script:**

First and foremost, I ran this fresh on a different HPC, so I downloaded the puzzler script and added it to my path:

```bash
wget https://raw.githubusercontent.com/merondun/puzzler/main/bin/puzzler
chmod +x puzzler
INSTALL_PATH=$(dirname "$(realpath puzzler)")
export PATH="$PATH:$INSTALL_PATH"
```

Then, change the SLURM scheduler so it can be submitted directly.

```
vim puzzler
########### EDIT THIS BLOCK WITH SLURM & APPTAINER/SINGULARITY SETTINGS ############
#SBATCH --time=2-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mem=64Gb
#SBATCH --partition=atlas
#SBATCH --account=coffea_pangenome

module load apptainer &> /dev/null || true
#module load singularity &> /dev/null || true
SINGULARITY_TMPDIR=$APPTAINER_TMPDIR
########### EDIT THIS BLOCK WITH SLURM & APPTAINER/SINGULARITY SETTINGS ############
:wq
```

**Puzzler container:**

I can simply grab the container with all the software with:

```bash
module load apptainer
apptainer pull --arch amd64 library://merondun/default/puzzler:v1.8

ls -lhtr ~/apptainer/
total 3.4G
-rwxr-x--- 1 justin.merondun justin.merondun 3.4G Jun 16 14:41 puzzler_v1.8.sif
-rwxr-x--- 1 justin.merondun justin.merondun  49K Jun 17 16:40 puzzler
```

Then I just navigate to my `$WD` and make sure it's still on my path:

```bash
cd /project/90daydata/coffea_pangenome/puzzler_trials/assemblies/fungus_tutorial
puzzler -h
Usage: puzzler [OPTIONS]

Options:
  -s, --sample SAMPLE   Sample name (required)
  -m, --map FILE        Path to .tsv/.csv map file (required)
  --threads t           Number of threads (optional; default 64)
  --mem MEM             Memory allocation (optional; default 512)
  --no_juice            Skip juicer file creation entirely (optional; not recommended!)
  -v, --version         Show version and exit
  -h, --help            Show help and exit

  Required --map Structure:
  The provided map file (e.g., samples.txt) must contain the following columns in this order:
  RUNTIME CONTAINER WD HIFI HIC_R1 HIC_R2 NUM_CHRS REFERENCE HOM_COV BLOB_DB BUSCO_LINEAGE BUSCO_DB
  For optional columns (REFERENCE - BUSCO_DB), write NA if undesired.
```

Good to go. 

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

***OPTIONAL columns***: *Specify "NA" if not needed! The script will skip respective components if "NA".*

* **reference:** Path to related species genome for chromosome naming. Scaffolds will be renamed to the closest syntenic chromosome **using their scaffold naming convention**.
* **hom_cov:** Homozygous peak coverage.
* **blob_database:** Directory to save all blobtools databases.
* **busco_lineage:** Busco odb10 version lineage.
* **busco_database:** Directory to save busco dbs.

For this example with 4 forks, my file will look like this, where the only difference is `$SAMPLE` (where files are stored) and `$REFERENCE`. 

| sample              | runtime   | container                                        | wd                                                           | hifi                                                         | hic_r1                                                       | hic_r2                                                       | num_chrs | reference                                                    | hom_cov | blob_database | busco_lineage | busco_database                                             |
| ------------------- | --------- | ------------------------------------------------ | ------------------------------------------------------------ | ------------------------------------------------------------ | ------------------------------------------------------------ | ------------------------------------------------------------ | -------- | ------------------------------------------------------------ | ------- | ------------- | ------------- | ---------------------------------------------------------- |
| Fungus_Ref          | apptainer | /home/justin.merondun/apptainer/puzzler_v1.8.sif | /90daydata/coffea_pangenome/puzzler_trials/assemblies/fungus_tutorial | /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Fungus.HiFi.fastq.gz | /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Fungus.HiC.R1.fastq.gz | /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Fungus.HiC.R2.fastq.gz | 14       | /90daydata/coffea_pangenome/puzzler_trials/raw_data/references/GCA_964106605.1_gdRhiKalk1.hap1.1_genomic.fna | NA      | NA            | fungi_odb10   | /90daydata/coffea_pangenome/puzzler_trials/busco_downloads |
| Fungus_NoRef        | apptainer | /home/justin.merondun/apptainer/puzzler_v1.8.sif | /90daydata/coffea_pangenome/puzzler_trials/assemblies/fungus_tutorial | /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Fungus.HiFi.fastq.gz | /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Fungus.HiC.R1.fastq.gz | /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Fungus.HiC.R2.fastq.gz | 14       | NA                                                           | NA      | NA            | fungi_odb10   | /90daydata/coffea_pangenome/puzzler_trials/busco_downloads |
| Fungus_NoJuice      | apptainer | /home/justin.merondun/apptainer/puzzler_v1.8.sif | /90daydata/coffea_pangenome/puzzler_trials/assemblies/fungus_tutorial | /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Fungus.HiFi.fastq.gz | /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Fungus.HiC.R1.fastq.gz | /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Fungus.HiC.R2.fastq.gz | 14       | /90daydata/coffea_pangenome/puzzler_trials/raw_data/references/GCA_964106605.1_gdRhiKalk1.hap1.1_genomic.fna | NA      | NA            | fungi_odb10   | /90daydata/coffea_pangenome/puzzler_trials/busco_downloads |
| Fungus_NoRefNoJuice | apptainer | /home/justin.merondun/apptainer/puzzler_v1.8.sif | /90daydata/coffea_pangenome/puzzler_trials/assemblies/fungus_tutorial | /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Fungus.HiFi.fastq.gz | /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Fungus.HiC.R1.fastq.gz | /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Fungus.HiC.R2.fastq.gz | 14       | NA                                                           | NA      | NA            | fungi_odb10   | /90daydata/coffea_pangenome/puzzler_trials/busco_downloads |

### Recommended: Juicer & Reference

#### Step 1: Draft Assembly

Submit: `sbatch -J Fungus_Ref puzzler -s Fungus_Ref -m fungus_example.tsv --threads 16 --mem 64` 

```bash
cat slurm-15858752.out

=======================================================================
__________ ____ _______________________.____     _____________________ 
\______   \    |   \____    /\____    /|    |    \_   _____/\______   \
 |     ___/    |   / /     /   /     / |    |     |    __)_  |       _/
 |    |   |    |  / /     /_  /     /_ |    |___  |        \ |    |   \
 |____|   |______/ /_______ \/_______ \|_______ \/_______  / |____|_  /
                           \/        \/        \/        \/         \/ 
=======================================================================

=======================================================================
Parameters for sample: Fungus_Ref 
RUNTIME: apptainer
CONTAINER: /home/justin.merondun/apptainer/puzzler_v1.8.sif 
WD: /90daydata/coffea_pangenome/puzzler_trials/assemblies/fungus_tutorial 
HIFI: /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Fungus.HiFi.fastq.gz
HIC_R1: /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Fungus.HiC.R1.fastq.gz
HIC_R2: /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Fungus.HiC.R2.fastq.gz
NUMBER CHRS: 14
REFERENCE: /90daydata/coffea_pangenome/puzzler_trials/raw_data/references/GCA_964106605.1_gdRhiKalk1.hap1.1_genomic.fna
HOM_COV: NA
BLOB_DB: NA
BUSCO_LINEAGE: fungi_odb10
BUSCO_DB: /90daydata/coffea_pangenome/puzzler_trials/busco_downloads
Cores Requested: 16
Cores Available: 16
RAM Requested: 64
Memory Available: 369.7 GB
PUZZLER command: apptainer exec --bind /90daydata/coffea_pangenome/puzzler_trials:/90daydata/coffea_pangenome/puzzler_trials --bind /90daydata/coffea_pangenome/puzzler_trials/assemblies:/90daydata/coffea_pangenome/puzzler_trials/assemblies --bind /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads:/90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads --bind /90daydata/coffea_pangenome/puzzler_trials/raw_data/references:/90daydata/coffea_pangenome/puzzler_trials/raw_data/references --bind /home/justin.merondun/apptainer:/home/justin.merondun/apptainer /home/justin.merondun/apptainer/puzzler_v1.8.sif
=======================================================================

~~~ [+0m] Checking software availability, this will take about 1 minute ~~~
~~~ [+1m] Software check complete ~~~
~~~ [+1m] Assembling genome for Fungus_Ref ~~~
~~~ [+1m] Running juicer, script will stop after .hic files created ~~~
~~~ [+1m] Starting hifiasm assembly for Fungus_Ref ~~~
~~~ [+6m] Starting Purge_Dups for Fungus_Ref ~~~
~~~ [+7m] Mapping HiC reads to Fungus_Ref draft ~~~
~~~ [+9m] Running HapHiC for Fungus_Ref  ~~~

❌ Command failed in in /90daydata/coffea_pangenome/puzzler_trials/assemblies/fungus_tutorial/Fungus_Ref/03_haphic/haphic: "${PUZZLER} haphic pipeline ${WD}/${SAMPLE}/02_purge_dups/p_ctg.purged.fa ../filtered.bam ${NUM_CHRS} --correct_nrounds 2 --max_inflation 10.0 --threads ${t} --processes ${t} 2> ../${SAMPLE}.haphic.log" (line 289)

~~~ [+11m] HapHiC for Fungus_Ref with 14 chrs failed, trying: 10 ~~~

❌ Command failed in in /90daydata/coffea_pangenome/puzzler_trials/assemblies/fungus_tutorial/Fungus_Ref/03_haphic/haphic: "${PUZZLER} haphic pipeline ${WD}/${SAMPLE}/02_purge_dups/p_ctg.purged.fa ../filtered.bam ${CHR_ATTEMPT} --correct_nrounds 2 --max_inflation 10.0 --threads ${t} --processes ${t} 2> ../${SAMPLE}.haphic_N_chrs${CHR_ATTEMPT}.log" (line 312)

~~~ [+13m] HapHiC for Fungus_Ref with 14 chrs failed, trying: 11 ~~~

❌ Command failed in in /90daydata/coffea_pangenome/puzzler_trials/assemblies/fungus_tutorial/Fungus_Ref/03_haphic/haphic: "${PUZZLER} haphic pipeline ${WD}/${SAMPLE}/02_purge_dups/p_ctg.purged.fa ../filtered.bam ${CHR_ATTEMPT} --correct_nrounds 2 --max_inflation 10.0 --threads ${t} --processes ${t} 2> ../${SAMPLE}.haphic_N_chrs${CHR_ATTEMPT}.log" (line 312)

~~~ [+15m] HapHiC for Fungus_Ref with 14 chrs failed, trying: 12 ~~~

❌ Command failed in in /90daydata/coffea_pangenome/puzzler_trials/assemblies/fungus_tutorial/Fungus_Ref/03_haphic/haphic: "${PUZZLER} haphic pipeline ${WD}/${SAMPLE}/02_purge_dups/p_ctg.purged.fa ../filtered.bam ${CHR_ATTEMPT} --correct_nrounds 2 --max_inflation 10.0 --threads ${t} --processes ${t} 2> ../${SAMPLE}.haphic_N_chrs${CHR_ATTEMPT}.log" (line 312)

~~~ [+18m] HapHiC for Fungus_Ref with 14 chrs failed, trying: 13 ~~~

❌ Command failed in in /90daydata/coffea_pangenome/puzzler_trials/assemblies/fungus_tutorial/Fungus_Ref/03_haphic/haphic: "${PUZZLER} haphic pipeline ${WD}/${SAMPLE}/02_purge_dups/p_ctg.purged.fa ../filtered.bam ${CHR_ATTEMPT} --correct_nrounds 2 --max_inflation 10.0 --threads ${t} --processes ${t} 2> ../${SAMPLE}.haphic_N_chrs${CHR_ATTEMPT}.log" (line 312)

~~~ [+20m] HapHiC for Fungus_Ref with 14 chrs failed, trying: 14 ~~~

❌ Command failed in in /90daydata/coffea_pangenome/puzzler_trials/assemblies/fungus_tutorial/Fungus_Ref/03_haphic/haphic: "${PUZZLER} haphic pipeline ${WD}/${SAMPLE}/02_purge_dups/p_ctg.purged.fa ../filtered.bam ${CHR_ATTEMPT} --correct_nrounds 2 --max_inflation 10.0 --threads ${t} --processes ${t} 2> ../${SAMPLE}.haphic_N_chrs${CHR_ATTEMPT}.log" (line 312)

~~~ [+23m] HapHiC for Fungus_Ref with 14 chrs failed, trying: 15 ~~~

❌ Command failed in in /90daydata/coffea_pangenome/puzzler_trials/assemblies/fungus_tutorial/Fungus_Ref/03_haphic/haphic: "${PUZZLER} haphic pipeline ${WD}/${SAMPLE}/02_purge_dups/p_ctg.purged.fa ../filtered.bam ${CHR_ATTEMPT} --correct_nrounds 2 --max_inflation 10.0 --threads ${t} --processes ${t} 2> ../${SAMPLE}.haphic_N_chrs${CHR_ATTEMPT}.log" (line 312)

~~~ [+25m] HapHiC for Fungus_Ref with 14 chrs failed, trying: 16 ~~~

❌ Command failed in in /90daydata/coffea_pangenome/puzzler_trials/assemblies/fungus_tutorial/Fungus_Ref/03_haphic/haphic: "${PUZZLER} haphic pipeline ${WD}/${SAMPLE}/02_purge_dups/p_ctg.purged.fa ../filtered.bam ${CHR_ATTEMPT} --correct_nrounds 2 --max_inflation 10.0 --threads ${t} --processes ${t} 2> ../${SAMPLE}.haphic_N_chrs${CHR_ATTEMPT}.log" (line 312)

~~~ [+28m] HapHiC for Fungus_Ref with 14 chrs failed, trying: 17 ~~~

❌ Command failed in in /90daydata/coffea_pangenome/puzzler_trials/assemblies/fungus_tutorial/Fungus_Ref/03_haphic/haphic: "${PUZZLER} haphic pipeline ${WD}/${SAMPLE}/02_purge_dups/p_ctg.purged.fa ../filtered.bam ${CHR_ATTEMPT} --correct_nrounds 2 --max_inflation 10.0 --threads ${t} --processes ${t} 2> ../${SAMPLE}.haphic_N_chrs${CHR_ATTEMPT}.log" (line 312)

~~~ [+30m] HapHiC for Fungus_Ref with 14 chrs failed, trying: 18 ~~~

❌ Command failed in in /90daydata/coffea_pangenome/puzzler_trials/assemblies/fungus_tutorial/Fungus_Ref/03_haphic/haphic: "${PUZZLER} haphic pipeline ${WD}/${SAMPLE}/02_purge_dups/p_ctg.purged.fa ../filtered.bam ${CHR_ATTEMPT} --correct_nrounds 2 --max_inflation 10.0 --threads ${t} --processes ${t} 2> ../${SAMPLE}.haphic_N_chrs${CHR_ATTEMPT}.log" (line 312)

~~~ [+32m] HapHiC for Fungus_Ref failed, scaffolding with YAHS instead ~~~
~~~ [+33m] Creating .hic file for juicebox for Fungus_Ref ~~~
~~~ [+34m] Post curation assembly file missing for Fungus_Ref: Run Juicebox & place in /90daydata/coffea_pangenome/puzzler_trials/assemblies/fungus_tutorial/juicer_files/Fungus_Ref_JBAT.review.assembly ~~~
```

:alarm_clock: This is a good time to check that all of your paths and parameters look appropriate. 

 We can then check the output `Fungus_Ref/01_hifiasm/Fungus_Ref.hifiasm.log` file to ensure an appropriate `--hom_cov` was selected:

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

:exclamation: **❌ Command failed:**

The output indicates that puzzler was unable to successfully run [HapHiC](https://github.com/zengxiaofei/HapHiC) with the specified number of chromosomes (from `samples.tsv`). Puzzler will automatically atttempt HapHiC with +/- 4 chromosomes, to see if it will successfully run with those chromosomes numbers instead. If that fails, puzzler will instead run [YAHS](https://github.com/c-zhou/yahs), which in my experience has always run successfully.

As long as you see:

`~~~~ HapHiC for Fungus_Ref failed, scaffolding with YAHS instead ~~~~` 

Afterwards, then everything is fine. 

Of the many genomes I have run so far, Fungus is the only genome not successfully scaffolded with HapHiC, so maybe it has something to do with the small size of chromosomes. 

**Runtime:**

It completed the assembly in 34 minutes.

``` 
Cores per node: 16
CPU Utilized: 02:25:14
CPU Efficiency: 26.58% of 09:06:24 core-walltime
Job Wall-clock time: 00:34:09
Memory Utilized: 17.50 GB
Memory Efficiency: 27.34% of 64.00 GB (64.00 GB/node)
```

In any case, the pipeline will also run `haphic refsort`, provided a related species is provided, which will align your draft scaffolds in order of the reference genome for juicebox curation. **This is not reference-guided assembly**, and will not rename or otherwise modify the scaffolds. Please see the [HapHiC](https://github.com/zengxiaofei/HapHiC) `refsort` description for more details. 

If you did not provide a `$REFERENCE` file, it will still output the `.hic` and `.assembly` files, just without running `haphic refsort`. 

#### Step 2: Juicebox Curation

Now we are ready to load the `$WD/juicer_files/Fungus_Ref_JBAT.hic` and `$WD/juicer_files/Fungus_Ref_JBAT.assembly` file into juicebox. 

Perform our edits. Pull in the `.hic` file, import the `.assembly` file with `Assembly > Import Map Assembly`. There's very minimal edits on this genome, I just merge two blocks which originate from one chromosome: 

To watch the video, please download the [Fungus_Ref.mp4](https://zenodo.org/records/15693026/files/Fungus_Ref.mp4?download=1) file from the manuscript's [Zenodo link](https://doi.org/10.5281/zenodo.15693025). 

Afterwards, save the file with `Assembly > Export Assembly`, and save the  `Fungus_Ref_JBAT.review.assembly` to the directory `$WD/juicer_files`, retaining the default file name! 

#### Step 3: Chromosome Renaming, Remapping HiC

Simply submit `puzzler` with the same parameters as before, and the script will assign chromosome names, output a final HiC contact map on the final assembly, and output assembly statistics. **I skip blobtools, because that will take several hours to blast against the nt database...**. 

We simply resubmit the exact same command, and the script will take up where it left off: 

`sbatch -J asm_Fungus_Ref puzzler -s Fungus_Ref -m samples.tsv --threads 16 --mem 64` 

```bash
cat slurm-15866418.out

=======================================================================
__________ ____ _______________________.____     _____________________ 
\______   \    |   \____    /\____    /|    |    \_   _____/\______   \
 |     ___/    |   / /     /   /     / |    |     |    __)_  |       _/
 |    |   |    |  / /     /_  /     /_ |    |___  |        \ |    |   \
 |____|   |______/ /_______ \/_______ \|_______ \/_______  / |____|_  /
                           \/        \/        \/        \/         \/ 
=======================================================================

=======================================================================
Parameters for sample: Fungus_Ref 
RUNTIME: apptainer
CONTAINER: /home/justin.merondun/apptainer/puzzler_v1.8.sif 
WD: /90daydata/coffea_pangenome/puzzler_trials/assemblies/fungus_tutorial 
HIFI: /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Fungus.HiFi.fastq.gz
HIC_R1: /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Fungus.HiC.R1.fastq.gz
HIC_R2: /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Fungus.HiC.R2.fastq.gz
NUMBER CHRS: 14
REFERENCE: /90daydata/coffea_pangenome/puzzler_trials/raw_data/references/GCA_964106605.1_gdRhiKalk1.hap1.1_genomic.fna
HOM_COV: NA
BLOB_DB: NA
BUSCO_LINEAGE: fungi_odb10
BUSCO_DB: /90daydata/coffea_pangenome/puzzler_trials/busco_downloads
Cores Requested: 16
Cores Available: 16
RAM Requested: 64
Memory Available: 287.3 GB
PUZZLER command: apptainer exec --bind /90daydata/coffea_pangenome/puzzler_trials:/90daydata/coffea_pangenome/puzzler_trials --bind /90daydata/coffea_pangenome/puzzler_trials/assemblies:/90daydata/coffea_pangenome/puzzler_trials/assemblies --bind /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads:/90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads --bind /90daydata/coffea_pangenome/puzzler_trials/raw_data/references:/90daydata/coffea_pangenome/puzzler_trials/raw_data/references --bind /home/justin.merondun/apptainer:/home/justin.merondun/apptainer /home/justin.merondun/apptainer/puzzler_v1.8.sif
=======================================================================

~~~ [+0m] Checking software availability, this will take about 1 minute ~~~
~~~ [+0m] Software check complete ~~~
~~~ [+0m] Assembling genome for Fungus_Ref ~~~
~~~ [+0m] Running juicer, script will stop after .hic files created ~~~
~~~ [+0m] Skipping hifiasm for Fungus_Ref: /90daydata/coffea_pangenome/puzzler_trials/assemblies/fungus_tutorial/Fungus_Ref/01_hifiasm/asm.hic.p_ctg.gfa exists ~~~
~~~ [+0m] Skipping purge for Fungus_Ref: /90daydata/coffea_pangenome/puzzler_trials/assemblies/fungus_tutorial/Fungus_Ref/02_purge_dups/p_ctg.purged.fa exists ~~~
~~~ [+0m] Skipping HiC alignment for Fungus_Ref: /90daydata/coffea_pangenome/puzzler_trials/assemblies/fungus_tutorial/Fungus_Ref/03_haphic/filtered.bam exists ~~~
~~~ [+0m] Skipping HapHiC for Fungus_Ref: /90daydata/coffea_pangenome/puzzler_trials/assemblies/fungus_tutorial/Fungus_Ref/03_haphic/haphic/04.build/scaffolds.fa exists ~~~
~~~ [+0m] Skipping juicer HiC file creation for Fungus_Ref: /90daydata/coffea_pangenome/puzzler_trials/assemblies/fungus_tutorial/juicer_files/Fungus_Ref_JBAT.hic exists ~~~
~~~ [+0m] Extracting post-curation assembly and mapping to reference for Fungus_Ref ~~~
~~~ [+1m] Renaming chromosomes for Fungus_Ref ~~~
~~~ [+1m] Scaffold sanity check passed for renaming, proceeding! ~~~
~~~ [+1m] Single scaffolds corresponding to a single Chr for Fungus_Ref ~~~
~~~ [+2m] Assessing genome quality for Fungus_Ref ~~~
~~~ [+2m] Creating final HiC bam for Fungus_Ref ~~~
~~~ [+4m] Creating final contact map for Fungus_Ref, plotting named chromosomes ~~~
~~~ [+4m] Running YAK on Fungus_Ref ~~~
~~~ [+5m] BUSCO lineage dataset already exists, skipping ~~~
~~~ [+5m] Running BUSCO for Fungus_Ref using lineage: fungi_odb10 ~~~
~~~ [+17m] Skipping blobtools for Fungus_Ref, not desired ~~~
~~~ [+17m] Summarizing Assembly for Fungus_Ref ~~~
~~~ [+18m] Your final assembly for Fungus_Ref is: /90daydata/coffea_pangenome/puzzler_trials/assemblies/fungus_tutorial/primary_asm/Fungus_Ref.fa ~~~
~~~ [+18m] Your final assembly stats for Fungus_Ref are in: /90daydata/coffea_pangenome/puzzler_trials/assemblies/fungus_tutorial/primary_asm/stats/Fungus_Ref.summary.txt ~~~
```

Which only ran for about 18 minutes. 

```
Cores per node: 16
CPU Utilized: 02:48:53
CPU Efficiency: 57.26% of 04:54:56 core-walltime
Job Wall-clock time: 00:18:26
Memory Utilized: 52.10 GB
Memory Efficiency: 81.41% of 64.00 GB (64.00 GB/node)
```

The script skips all the previously run steps, because we have the `.complete` files in the `$WD/Fungus_Ref` dir, as well as any already-completed step-specific files (e.g. `$WD/juicer_files/$SAMPLE_JBAT.review.assembly` for post-curation assembly and renaming). 

```
drwxr-s---. 2 justin.merondun proj-coffea_pangenome 4096 Jun 26 15:21 01_hifiasm
drwxr-s---. 2 justin.merondun proj-coffea_pangenome 4096 Jun 26 15:47 02_purge_dups
drwxr-s---. 2 justin.merondun proj-coffea_pangenome 4096 Jun 26 15:47 03_haphic
drwxr-s---. 2 justin.merondun proj-coffea_pangenome 4096 Jun 26 15:48 04_juicer
drwxr-s---. 2 justin.merondun proj-coffea_pangenome 4096 Jun 26 16:04 05_postjuicebox
drwxr-s---. 2 justin.merondun proj-coffea_pangenome 4096 Jun 26 16:07 06_realign_hic_hifi
drwxr-s---. 2 justin.merondun proj-coffea_pangenome 4096 Jun 26 16:19 07_busco_yak_blob
-rw-r-----. 1 justin.merondun proj-coffea_pangenome    0 Jun 26 15:23 align_hic.complete
-rw-r-----. 1 justin.merondun proj-coffea_pangenome    0 Jun 26 16:20 busco.complete
-rw-r-----. 1 justin.merondun proj-coffea_pangenome    0 Jun 26 15:21 hifiasm.complete
-rw-r-----. 1 justin.merondun proj-coffea_pangenome    0 Jun 26 15:48 juicer.complete
-rw-r-----. 1 justin.merondun proj-coffea_pangenome    0 Jun 26 15:22 purge_dups.complete
-rw-r-----. 1 justin.merondun proj-coffea_pangenome    0 Jun 26 16:06 qc_align_hic.complete
-rw-r-----. 1 justin.merondun proj-coffea_pangenome    0 Jun 26 15:47 scaffolding.complete
-rw-r-----. 1 justin.merondun proj-coffea_pangenome    0 Jun 26 16:08 yak.complete
```

Which outputs a final HiC contact map at `$WD/primary_asm/stats/$SAMPLE.pdf`: 

![final hic contacts](/examples/figs/Fungus_Final.png)

And the final assembly stats here:

```
cat primary_asm/stats/Fungus_Ref.summary.txt
Sample  SizeBP  Sequences       Contigs Gaps    ContigN50       ScafN50 WithinChrsBP    PropWithinChrs  BUSCO_Complete  BUSCO_singlecopy        YAK_CV  YAK_QV
Fungus_Ref       24817076        22      23     1        1776619         1806824         24324650       0.9802  99.1    99.1    0.999   69.904
```



### Juicer - No Reference

Submit: `sbatch -J Fungus_NoRef puzzler -s Fungus_NoRef -m fungus_example.tsv --threads 16 --mem 64` 

```
cat slurm-15858948.out

=======================================================================
__________ ____ _______________________.____     _____________________ 
\______   \    |   \____    /\____    /|    |    \_   _____/\______   \
 |     ___/    |   / /     /   /     / |    |     |    __)_  |       _/
 |    |   |    |  / /     /_  /     /_ |    |___  |        \ |    |   \
 |____|   |______/ /_______ \/_______ \|_______ \/_______  / |____|_  /
                           \/        \/        \/        \/         \/ 
=======================================================================

=======================================================================
Parameters for sample: Fungus_NoRef 
RUNTIME: apptainer
CONTAINER: /home/justin.merondun/apptainer/puzzler_v1.8.sif 
WD: /90daydata/coffea_pangenome/puzzler_trials/assemblies/fungus_tutorial 
HIFI: /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Fungus.HiFi.fastq.gz
HIC_R1: /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Fungus.HiC.R1.fastq.gz
HIC_R2: /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Fungus.HiC.R2.fastq.gz
NUMBER CHRS: 14
REFERENCE: NA
HOM_COV: NA
BLOB_DB: NA
BUSCO_LINEAGE: fungi_odb10
BUSCO_DB: /90daydata/coffea_pangenome/puzzler_trials/busco_downloads
Cores Requested: 16
Cores Available: 16
RAM Requested: 64
Memory Available: 363.8 GB
PUZZLER command: apptainer exec --bind /90daydata/coffea_pangenome/puzzler_trials:/90daydata/coffea_pangenome/puzzler_trials --bind /90daydata/coffea_pangenome/puzzler_trials/assemblies:/90daydata/coffea_pangenome/puzzler_trials/assemblies --bind /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads:/90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads --bind /home/justin.merondun/apptainer:/home/justin.merondun/apptainer /home/justin.merondun/apptainer/puzzler_v1.8.sif
=======================================================================

~~~ [+0m] Checking software availability, this will take about 1 minute ~~~
~~~ [+0m] Software check complete ~~~
~~~ [+0m] Assembling genome for Fungus_NoRef ~~~
~~~ [+0m] Running juicer, script will stop after .hic files created ~~~
~~~ [+0m] Starting hifiasm assembly for Fungus_NoRef ~~~
~~~ [+6m] Starting Purge_Dups for Fungus_NoRef ~~~
~~~ [+6m] Mapping HiC reads to Fungus_NoRef draft ~~~
~~~ [+8m] Running HapHiC for Fungus_NoRef  ~~~

❌ Command failed in in /90daydata/coffea_pangenome/puzzler_trials/assemblies/fungus_tutorial/Fungus_NoRef/03_haphic/haphic: "${PUZZLER} haphic pipeline ${WD}/${SAMPLE}/02_purge_dups/p_ctg.purged.fa ../filtered.bam ${NUM_CHRS} --correct_nrounds 2 --max_inflation 10.0 --threads ${t} --processes ${t} 2> ../${SAMPLE}.haphic.log" (line 289)

~~~ [+10m] HapHiC for Fungus_NoRef with 14 chrs failed, trying: 10 ~~~

❌ Command failed in in /90daydata/coffea_pangenome/puzzler_trials/assemblies/fungus_tutorial/Fungus_NoRef/03_haphic/haphic: "${PUZZLER} haphic pipeline ${WD}/${SAMPLE}/02_purge_dups/p_ctg.purged.fa ../filtered.bam ${CHR_ATTEMPT} --correct_nrounds 2 --max_inflation 10.0 --threads ${t} --processes ${t} 2> ../${SAMPLE}.haphic_N_chrs${CHR_ATTEMPT}.log" (line 312)

~~~ [+12m] HapHiC for Fungus_NoRef with 14 chrs failed, trying: 11 ~~~

❌ Command failed in in /90daydata/coffea_pangenome/puzzler_trials/assemblies/fungus_tutorial/Fungus_NoRef/03_haphic/haphic: "${PUZZLER} haphic pipeline ${WD}/${SAMPLE}/02_purge_dups/p_ctg.purged.fa ../filtered.bam ${CHR_ATTEMPT} --correct_nrounds 2 --max_inflation 10.0 --threads ${t} --processes ${t} 2> ../${SAMPLE}.haphic_N_chrs${CHR_ATTEMPT}.log" (line 312)

~~~ [+15m] HapHiC for Fungus_NoRef with 14 chrs failed, trying: 12 ~~~

❌ Command failed in in /90daydata/coffea_pangenome/puzzler_trials/assemblies/fungus_tutorial/Fungus_NoRef/03_haphic/haphic: "${PUZZLER} haphic pipeline ${WD}/${SAMPLE}/02_purge_dups/p_ctg.purged.fa ../filtered.bam ${CHR_ATTEMPT} --correct_nrounds 2 --max_inflation 10.0 --threads ${t} --processes ${t} 2> ../${SAMPLE}.haphic_N_chrs${CHR_ATTEMPT}.log" (line 312)

~~~ [+17m] HapHiC for Fungus_NoRef with 14 chrs failed, trying: 13 ~~~

❌ Command failed in in /90daydata/coffea_pangenome/puzzler_trials/assemblies/fungus_tutorial/Fungus_NoRef/03_haphic/haphic: "${PUZZLER} haphic pipeline ${WD}/${SAMPLE}/02_purge_dups/p_ctg.purged.fa ../filtered.bam ${CHR_ATTEMPT} --correct_nrounds 2 --max_inflation 10.0 --threads ${t} --processes ${t} 2> ../${SAMPLE}.haphic_N_chrs${CHR_ATTEMPT}.log" (line 312)

~~~ [+20m] HapHiC for Fungus_NoRef with 14 chrs failed, trying: 14 ~~~

❌ Command failed in in /90daydata/coffea_pangenome/puzzler_trials/assemblies/fungus_tutorial/Fungus_NoRef/03_haphic/haphic: "${PUZZLER} haphic pipeline ${WD}/${SAMPLE}/02_purge_dups/p_ctg.purged.fa ../filtered.bam ${CHR_ATTEMPT} --correct_nrounds 2 --max_inflation 10.0 --threads ${t} --processes ${t} 2> ../${SAMPLE}.haphic_N_chrs${CHR_ATTEMPT}.log" (line 312)

~~~ [+22m] HapHiC for Fungus_NoRef with 14 chrs failed, trying: 15 ~~~

❌ Command failed in in /90daydata/coffea_pangenome/puzzler_trials/assemblies/fungus_tutorial/Fungus_NoRef/03_haphic/haphic: "${PUZZLER} haphic pipeline ${WD}/${SAMPLE}/02_purge_dups/p_ctg.purged.fa ../filtered.bam ${CHR_ATTEMPT} --correct_nrounds 2 --max_inflation 10.0 --threads ${t} --processes ${t} 2> ../${SAMPLE}.haphic_N_chrs${CHR_ATTEMPT}.log" (line 312)

~~~ [+24m] HapHiC for Fungus_NoRef with 14 chrs failed, trying: 16 ~~~

❌ Command failed in in /90daydata/coffea_pangenome/puzzler_trials/assemblies/fungus_tutorial/Fungus_NoRef/03_haphic/haphic: "${PUZZLER} haphic pipeline ${WD}/${SAMPLE}/02_purge_dups/p_ctg.purged.fa ../filtered.bam ${CHR_ATTEMPT} --correct_nrounds 2 --max_inflation 10.0 --threads ${t} --processes ${t} 2> ../${SAMPLE}.haphic_N_chrs${CHR_ATTEMPT}.log" (line 312)

~~~ [+27m] HapHiC for Fungus_NoRef with 14 chrs failed, trying: 17 ~~~

❌ Command failed in in /90daydata/coffea_pangenome/puzzler_trials/assemblies/fungus_tutorial/Fungus_NoRef/03_haphic/haphic: "${PUZZLER} haphic pipeline ${WD}/${SAMPLE}/02_purge_dups/p_ctg.purged.fa ../filtered.bam ${CHR_ATTEMPT} --correct_nrounds 2 --max_inflation 10.0 --threads ${t} --processes ${t} 2> ../${SAMPLE}.haphic_N_chrs${CHR_ATTEMPT}.log" (line 312)

~~~ [+29m] HapHiC for Fungus_NoRef with 14 chrs failed, trying: 18 ~~~

❌ Command failed in in /90daydata/coffea_pangenome/puzzler_trials/assemblies/fungus_tutorial/Fungus_NoRef/03_haphic/haphic: "${PUZZLER} haphic pipeline ${WD}/${SAMPLE}/02_purge_dups/p_ctg.purged.fa ../filtered.bam ${CHR_ATTEMPT} --correct_nrounds 2 --max_inflation 10.0 --threads ${t} --processes ${t} 2> ../${SAMPLE}.haphic_N_chrs${CHR_ATTEMPT}.log" (line 312)

~~~ [+31m] HapHiC for Fungus_NoRef failed, scaffolding with YAHS instead ~~~
~~~ [+32m] Creating .hic file for juicebox for Fungus_NoRef ~~~
~~~ [+32m] Post curation assembly file missing for Fungus_NoRef: Run Juicebox & place in /90daydata/coffea_pangenome/puzzler_trials/assemblies/fungus_tutorial/juicer_files/Fungus_NoRef_JBAT.review.assembly ~~~
```

Ran in 32 minutes, and of course we didn't need nearly that many resources:

```bash
Cores per node: 16
CPU Utilized: 02:32:09
CPU Efficiency: 29.04% of 08:44:00 core-walltime
Job Wall-clock time: 00:32:45
Memory Utilized: 17.49 GB
Memory Efficiency: 27.33% of 64.00 GB (64.00 GB/node)
```

Perform juicebox edits....

To watch the video, please download the [Fungus_NoRef.mp4](https://zenodo.org/records/15693026/files/Fungus_NoRef.mp4?download=1) file from the manuscript's [Zenodo link](https://doi.org/10.5281/zenodo.15693025). 

and then resubmit: `sbatch -J Fungus_NoRef puzzler -s Fungus_NoRef -m fungus_example.tsv --threads 16 --mem 64` 

```
cat slurm-15866422.out

=======================================================================
__________ ____ _______________________.____     _____________________ 
\______   \    |   \____    /\____    /|    |    \_   _____/\______   \
 |     ___/    |   / /     /   /     / |    |     |    __)_  |       _/
 |    |   |    |  / /     /_  /     /_ |    |___  |        \ |    |   \
 |____|   |______/ /_______ \/_______ \|_______ \/_______  / |____|_  /
                           \/        \/        \/        \/         \/ 
=======================================================================

=======================================================================
Parameters for sample: Fungus_NoRef 
RUNTIME: apptainer
CONTAINER: /home/justin.merondun/apptainer/puzzler_v1.8.sif 
WD: /90daydata/coffea_pangenome/puzzler_trials/assemblies/fungus_tutorial 
HIFI: /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Fungus.HiFi.fastq.gz
HIC_R1: /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Fungus.HiC.R1.fastq.gz
HIC_R2: /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Fungus.HiC.R2.fastq.gz
NUMBER CHRS: 14
REFERENCE: NA
HOM_COV: NA
BLOB_DB: NA
BUSCO_LINEAGE: fungi_odb10
BUSCO_DB: /90daydata/coffea_pangenome/puzzler_trials/busco_downloads
Cores Requested: 16
Cores Available: 16
RAM Requested: 64
Memory Available: 2159.1 GB
PUZZLER command: apptainer exec --bind /90daydata/coffea_pangenome/puzzler_trials:/90daydata/coffea_pangenome/puzzler_trials --bind /90daydata/coffea_pangenome/puzzler_trials/assemblies:/90daydata/coffea_pangenome/puzzler_trials/assemblies --bind /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads:/90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads --bind /home/justin.merondun/apptainer:/home/justin.merondun/apptainer /home/justin.merondun/apptainer/puzzler_v1.8.sif
=======================================================================

~~~ [+0m] Checking software availability, this will take about 1 minute ~~~
~~~ [+0m] Software check complete ~~~
~~~ [+0m] Assembling genome for Fungus_NoRef ~~~
~~~ [+0m] Running juicer, script will stop after .hic files created ~~~
~~~ [+0m] Skipping hifiasm for Fungus_NoRef: /90daydata/coffea_pangenome/puzzler_trials/assemblies/fungus_tutorial/Fungus_NoRef/01_hifiasm/asm.hic.p_ctg.gfa exists ~~~
~~~ [+0m] Skipping purge for Fungus_NoRef: /90daydata/coffea_pangenome/puzzler_trials/assemblies/fungus_tutorial/Fungus_NoRef/02_purge_dups/p_ctg.purged.fa exists ~~~
~~~ [+0m] Skipping HiC alignment for Fungus_NoRef: /90daydata/coffea_pangenome/puzzler_trials/assemblies/fungus_tutorial/Fungus_NoRef/03_haphic/filtered.bam exists ~~~
~~~ [+0m] Skipping HapHiC for Fungus_NoRef: /90daydata/coffea_pangenome/puzzler_trials/assemblies/fungus_tutorial/Fungus_NoRef/03_haphic/haphic/04.build/scaffolds.fa exists ~~~
~~~ [+0m] Skipping juicer HiC file creation for Fungus_NoRef: /90daydata/coffea_pangenome/puzzler_trials/assemblies/fungus_tutorial/juicer_files/Fungus_NoRef_JBAT.hic exists ~~~
~~~ [+0m] No reference provided for Fungus_NoRef: simply extracting assembly to: /90daydata/coffea_pangenome/puzzler_trials/assemblies/fungus_tutorial/primary_asm/Fungus_NoRef.fa ~~~
~~~ [+0m] No reference provided for Fungus_NoRef: no chromosome re-naming, so final assembly already: /90daydata/coffea_pangenome/puzzler_trials/assemblies/fungus_tutorial/primary_asm/Fungus_NoRef.fa ~~~
~~~ [+0m] Assessing genome quality for Fungus_NoRef ~~~
~~~ [+0m] Creating final HiC bam for Fungus_NoRef ~~~
~~~ [+2m] Creating final contact map for Fungus_NoRef, no ref specified, plotting scaffolds > 2mb ~~~
~~~ [+3m] Running YAK on Fungus_NoRef ~~~
~~~ [+4m] BUSCO lineage dataset already exists, skipping ~~~
~~~ [+4m] Running BUSCO for Fungus_NoRef using lineage: fungi_odb10 ~~~
~~~ [+15m] Skipping blobtools for Fungus_NoRef, not desired ~~~
~~~ [+15m] Summarizing Assembly for Fungus_NoRef ~~~
~~~ [+15m] Your final assembly for Fungus_NoRef is: /90daydata/coffea_pangenome/puzzler_trials/assemblies/fungus_tutorial/primary_asm/Fungus_NoRef.fa ~~~
~~~ [+15m] Your final assembly stats for Fungus_NoRef are in: /90daydata/coffea_pangenome/puzzler_trials/assemblies/fungus_tutorial/primary_asm/stats/Fungus_NoRef.summary.txt ~~~
```

With only 16 cores, took around 16 minutes:

```bash
Cores per node: 16
CPU Utilized: 02:13:45
CPU Efficiency: 52.96% of 04:12:32 core-walltime
Job Wall-clock time: 00:15:47
Memory Utilized: 51.13 GB
Memory Efficiency: 79.89% of 64.00 GB (64.00 GB/node)
```

### No Juicer - With Reference

Just add the `--no_juice` flag to skip juicer. Not recommended! 

Submit: `sbatch -J Fungus_NoJuice puzzler -s Fungus_NoJuice -m fungus_example.tsv --no_juice --threads 16 --mem 64` 

```
cat slurm-15866287.out

=======================================================================
__________ ____ _______________________.____     _____________________ 
\______   \    |   \____    /\____    /|    |    \_   _____/\______   \
 |     ___/    |   / /     /   /     / |    |     |    __)_  |       _/
 |    |   |    |  / /     /_  /     /_ |    |___  |        \ |    |   \
 |____|   |______/ /_______ \/_______ \|_______ \/_______  / |____|_  /
                           \/        \/        \/        \/         \/ 
=======================================================================

=======================================================================
Parameters for sample: Fungus_NoJuice 
RUNTIME: apptainer
CONTAINER: /home/justin.merondun/apptainer/puzzler_v1.8.sif 
WD: /90daydata/coffea_pangenome/puzzler_trials/assemblies/fungus_tutorial 
HIFI: /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Fungus.HiFi.fastq.gz
HIC_R1: /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Fungus.HiC.R1.fastq.gz
HIC_R2: /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Fungus.HiC.R2.fastq.gz
NUMBER CHRS: 14
REFERENCE: /90daydata/coffea_pangenome/puzzler_trials/raw_data/references/GCA_964106605.1_gdRhiKalk1.hap1.1_genomic.fna
HOM_COV: NA
BLOB_DB: NA
BUSCO_LINEAGE: fungi_odb10
BUSCO_DB: /90daydata/coffea_pangenome/puzzler_trials/busco_downloads
Cores Requested: 16
Cores Available: 16
RAM Requested: 64
Memory Available: 359.8 GB
PUZZLER command: apptainer exec --bind /90daydata/coffea_pangenome/puzzler_trials:/90daydata/coffea_pangenome/puzzler_trials --bind /90daydata/coffea_pangenome/puzzler_trials/assemblies:/90daydata/coffea_pangenome/puzzler_trials/assemblies --bind /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads:/90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads --bind /90daydata/coffea_pangenome/puzzler_trials/raw_data/references:/90daydata/coffea_pangenome/puzzler_trials/raw_data/references --bind /home/justin.merondun/apptainer:/home/justin.merondun/apptainer /home/justin.merondun/apptainer/puzzler_v1.8.sif
=======================================================================

~~~ [+0m] Checking software availability, this will take about 1 minute ~~~
~~~ [+0m] Software check complete ~~~
~~~ [+0m] Assembling genome for Fungus_NoJuice ~~~
~~~ [+0m] Skipping juicer, no manual curation (not recommended!) ~~~
~~~ [+0m] Starting hifiasm assembly for Fungus_NoJuice ~~~
~~~ [+5m] Starting Purge_Dups for Fungus_NoJuice ~~~
~~~ [+6m] Mapping HiC reads to Fungus_NoJuice draft ~~~
~~~ [+7m] Running HapHiC for Fungus_NoJuice  ~~~

❌ Command failed in in /90daydata/coffea_pangenome/puzzler_trials/assemblies/fungus_tutorial/Fungus_NoJuice/03_haphic/haphic: "${PUZZLER} haphic pipeline ${WD}/${SAMPLE}/02_purge_dups/p_ctg.purged.fa ../filtered.bam ${NUM_CHRS} --correct_nrounds 2 --max_inflation 10.0 --threads ${t} --processes ${t} 2> ../${SAMPLE}.haphic.log" (line 289)

~~~ [+9m] HapHiC for Fungus_NoJuice with 14 chrs failed, trying: 10 ~~~

❌ Command failed in in /90daydata/coffea_pangenome/puzzler_trials/assemblies/fungus_tutorial/Fungus_NoJuice/03_haphic/haphic: "${PUZZLER} haphic pipeline ${WD}/${SAMPLE}/02_purge_dups/p_ctg.purged.fa ../filtered.bam ${CHR_ATTEMPT} --correct_nrounds 2 --max_inflation 10.0 --threads ${t} --processes ${t} 2> ../${SAMPLE}.haphic_N_chrs${CHR_ATTEMPT}.log" (line 312)

~~~ [+12m] HapHiC for Fungus_NoJuice with 14 chrs failed, trying: 11 ~~~

❌ Command failed in in /90daydata/coffea_pangenome/puzzler_trials/assemblies/fungus_tutorial/Fungus_NoJuice/03_haphic/haphic: "${PUZZLER} haphic pipeline ${WD}/${SAMPLE}/02_purge_dups/p_ctg.purged.fa ../filtered.bam ${CHR_ATTEMPT} --correct_nrounds 2 --max_inflation 10.0 --threads ${t} --processes ${t} 2> ../${SAMPLE}.haphic_N_chrs${CHR_ATTEMPT}.log" (line 312)

~~~ [+14m] HapHiC for Fungus_NoJuice with 14 chrs failed, trying: 12 ~~~

❌ Command failed in in /90daydata/coffea_pangenome/puzzler_trials/assemblies/fungus_tutorial/Fungus_NoJuice/03_haphic/haphic: "${PUZZLER} haphic pipeline ${WD}/${SAMPLE}/02_purge_dups/p_ctg.purged.fa ../filtered.bam ${CHR_ATTEMPT} --correct_nrounds 2 --max_inflation 10.0 --threads ${t} --processes ${t} 2> ../${SAMPLE}.haphic_N_chrs${CHR_ATTEMPT}.log" (line 312)

~~~ [+17m] HapHiC for Fungus_NoJuice with 14 chrs failed, trying: 13 ~~~

❌ Command failed in in /90daydata/coffea_pangenome/puzzler_trials/assemblies/fungus_tutorial/Fungus_NoJuice/03_haphic/haphic: "${PUZZLER} haphic pipeline ${WD}/${SAMPLE}/02_purge_dups/p_ctg.purged.fa ../filtered.bam ${CHR_ATTEMPT} --correct_nrounds 2 --max_inflation 10.0 --threads ${t} --processes ${t} 2> ../${SAMPLE}.haphic_N_chrs${CHR_ATTEMPT}.log" (line 312)

~~~ [+19m] HapHiC for Fungus_NoJuice with 14 chrs failed, trying: 14 ~~~

❌ Command failed in in /90daydata/coffea_pangenome/puzzler_trials/assemblies/fungus_tutorial/Fungus_NoJuice/03_haphic/haphic: "${PUZZLER} haphic pipeline ${WD}/${SAMPLE}/02_purge_dups/p_ctg.purged.fa ../filtered.bam ${CHR_ATTEMPT} --correct_nrounds 2 --max_inflation 10.0 --threads ${t} --processes ${t} 2> ../${SAMPLE}.haphic_N_chrs${CHR_ATTEMPT}.log" (line 312)

~~~ [+22m] HapHiC for Fungus_NoJuice with 14 chrs failed, trying: 15 ~~~

❌ Command failed in in /90daydata/coffea_pangenome/puzzler_trials/assemblies/fungus_tutorial/Fungus_NoJuice/03_haphic/haphic: "${PUZZLER} haphic pipeline ${WD}/${SAMPLE}/02_purge_dups/p_ctg.purged.fa ../filtered.bam ${CHR_ATTEMPT} --correct_nrounds 2 --max_inflation 10.0 --threads ${t} --processes ${t} 2> ../${SAMPLE}.haphic_N_chrs${CHR_ATTEMPT}.log" (line 312)

~~~ [+24m] HapHiC for Fungus_NoJuice with 14 chrs failed, trying: 16 ~~~

❌ Command failed in in /90daydata/coffea_pangenome/puzzler_trials/assemblies/fungus_tutorial/Fungus_NoJuice/03_haphic/haphic: "${PUZZLER} haphic pipeline ${WD}/${SAMPLE}/02_purge_dups/p_ctg.purged.fa ../filtered.bam ${CHR_ATTEMPT} --correct_nrounds 2 --max_inflation 10.0 --threads ${t} --processes ${t} 2> ../${SAMPLE}.haphic_N_chrs${CHR_ATTEMPT}.log" (line 312)

~~~ [+26m] HapHiC for Fungus_NoJuice with 14 chrs failed, trying: 17 ~~~

❌ Command failed in in /90daydata/coffea_pangenome/puzzler_trials/assemblies/fungus_tutorial/Fungus_NoJuice/03_haphic/haphic: "${PUZZLER} haphic pipeline ${WD}/${SAMPLE}/02_purge_dups/p_ctg.purged.fa ../filtered.bam ${CHR_ATTEMPT} --correct_nrounds 2 --max_inflation 10.0 --threads ${t} --processes ${t} 2> ../${SAMPLE}.haphic_N_chrs${CHR_ATTEMPT}.log" (line 312)

~~~ [+29m] HapHiC for Fungus_NoJuice with 14 chrs failed, trying: 18 ~~~

❌ Command failed in in /90daydata/coffea_pangenome/puzzler_trials/assemblies/fungus_tutorial/Fungus_NoJuice/03_haphic/haphic: "${PUZZLER} haphic pipeline ${WD}/${SAMPLE}/02_purge_dups/p_ctg.purged.fa ../filtered.bam ${CHR_ATTEMPT} --correct_nrounds 2 --max_inflation 10.0 --threads ${t} --processes ${t} 2> ../${SAMPLE}.haphic_N_chrs${CHR_ATTEMPT}.log" (line 312)

~~~ [+31m] HapHiC for Fungus_NoJuice failed, scaffolding with YAHS instead ~~~
~~~ [+31m] Skipping juicer HiC file creation for Fungus_NoJuice: --no_juice given, not recommended! ~~~
~~~ [+31m] Skipping juicer extraction for Fungus_NoJuice, --no_juice requested! (not recommended!) ~~~
~~~ [+32m] Renaming chromosomes for Fungus_NoJuice ~~~
~~~ [+32m] Scaffold sanity check passed for renaming, proceeding! ~~~
~~~ [+32m] Multiple scaffolds corresponding to single Chr for Fungus_NoJuice, Renaming them e.g. Chr1A, Chr1B.. ~~~
~~~ [+32m] Assessing genome quality for Fungus_NoJuice ~~~
~~~ [+32m] Creating final HiC bam for Fungus_NoJuice ~~~
~~~ [+34m] Creating final contact map for Fungus_NoJuice, plotting named chromosomes ~~~
~~~ [+34m] Running YAK on Fungus_NoJuice ~~~
~~~ [+35m] BUSCO lineage dataset already exists, skipping ~~~
~~~ [+35m] Running BUSCO for Fungus_NoJuice using lineage: fungi_odb10 ~~~
~~~ [+48m] Skipping blobtools for Fungus_NoJuice, not desired ~~~
~~~ [+48m] Summarizing Assembly for Fungus_NoJuice ~~~
~~~ [+48m] Your final assembly for Fungus_NoJuice is: /90daydata/coffea_pangenome/puzzler_trials/assemblies/fungus_tutorial/primary_asm/Fungus_NoJuice.fa ~~~
~~~ [+48m] Your final assembly stats for Fungus_NoJuice are in: /90daydata/coffea_pangenome/puzzler_trials/assemblies/fungus_tutorial/primary_asm/stats/Fungus_NoJuice.summary.txt ~~~
```

Runtime:

```
Cores per node: 16
CPU Utilized: 05:32:55
CPU Efficiency: 42.51% of 13:03:12 core-walltime
Job Wall-clock time: 00:48:57
Memory Utilized: 53.50 GB
Memory Efficiency: 83.60% of 64.00 GB (64.00 GB/node)
```

### No Juicer - No Reference

Just add the `--no_juice` flag to skip juicer. Not recommended! 

Submit: `sbatch -J Fungus_NoRefNoJuice puzzler -s Fungus_NoRefNoJuice -m fungus_example.tsv --no_juice --threads 16 --mem 64` 

```
cat slurm-15858351.out

=======================================================================
__________ ____ _______________________.____     _____________________ 
\______   \    |   \____    /\____    /|    |    \_   _____/\______   \
 |     ___/    |   / /     /   /     / |    |     |    __)_  |       _/
 |    |   |    |  / /     /_  /     /_ |    |___  |        \ |    |   \
 |____|   |______/ /_______ \/_______ \|_______ \/_______  / |____|_  /
                           \/        \/        \/        \/         \/ 
=======================================================================

=======================================================================
Parameters for sample: Fungus_NoRefNoJuice 
RUNTIME: apptainer
CONTAINER: /home/justin.merondun/apptainer/puzzler_v1.8.sif 
WD: /90daydata/coffea_pangenome/puzzler_trials/assemblies/fungus_tutorial 
HIFI: /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Fungus.HiFi.fastq.gz
HIC_R1: /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Fungus.HiC.R1.fastq.gz
HIC_R2: /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Fungus.HiC.R2.fastq.gz
NUMBER CHRS: 14
REFERENCE: NA
HOM_COV: NA
BLOB_DB: NA
BUSCO_LINEAGE: fungi_odb10
BUSCO_DB: /90daydata/coffea_pangenome/puzzler_trials/busco_downloads
Cores Requested: 16
Cores Available: 16
RAM Requested: 64
Memory Available: 317.3 GB
PUZZLER command: apptainer exec --bind /90daydata/coffea_pangenome/puzzler_trials:/90daydata/coffea_pangenome/puzzler_trials --bind /90daydata/coffea_pangenome/puzzler_trials/assemblies:/90daydata/coffea_pangenome/puzzler_trials/assemblies --bind /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads:/90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads --bind /home/justin.merondun/apptainer:/home/justin.merondun/apptainer /home/justin.merondun/apptainer/puzzler_v1.8.sif
=======================================================================

~~~ [+0m] Checking software availability, this will take about 1 minute ~~~
~~~ [+0m] Software check complete ~~~
~~~ [+0m] Assembling genome for Fungus_NoRefNoJuice ~~~
~~~ [+0m] Skipping juicer, no manual curation (not recommended!) ~~~
~~~ [+0m] Starting hifiasm assembly for Fungus_NoRefNoJuice ~~~
~~~ [+6m] Starting Purge_Dups for Fungus_NoRefNoJuice ~~~
~~~ [+7m] Mapping HiC reads to Fungus_NoRefNoJuice draft ~~~
~~~ [+8m] Running HapHiC for Fungus_NoRefNoJuice  ~~~

❌ Command failed in in /90daydata/coffea_pangenome/puzzler_trials/assemblies/fungus_tutorial/Fungus_NoRefNoJuice/03_haphic/haphic: "${PUZZLER} haphic pipeline ${WD}/${SAMPLE}/02_purge_dups/p_ctg.purged.fa ../filtered.bam ${NUM_CHRS} --correct_nrounds 2 --max_inflation 10.0 --threads ${t} --processes ${t} 2> ../${SAMPLE}.haphic.log" (line 289)

~~~ [+11m] HapHiC for Fungus_NoRefNoJuice with 14 chrs failed, trying: 10 ~~~

❌ Command failed in in /90daydata/coffea_pangenome/puzzler_trials/assemblies/fungus_tutorial/Fungus_NoRefNoJuice/03_haphic/haphic: "${PUZZLER} haphic pipeline ${WD}/${SAMPLE}/02_purge_dups/p_ctg.purged.fa ../filtered.bam ${CHR_ATTEMPT} --correct_nrounds 2 --max_inflation 10.0 --threads ${t} --processes ${t} 2> ../${SAMPLE}.haphic_N_chrs${CHR_ATTEMPT}.log" (line 312)

~~~ [+13m] HapHiC for Fungus_NoRefNoJuice with 14 chrs failed, trying: 11 ~~~

❌ Command failed in in /90daydata/coffea_pangenome/puzzler_trials/assemblies/fungus_tutorial/Fungus_NoRefNoJuice/03_haphic/haphic: "${PUZZLER} haphic pipeline ${WD}/${SAMPLE}/02_purge_dups/p_ctg.purged.fa ../filtered.bam ${CHR_ATTEMPT} --correct_nrounds 2 --max_inflation 10.0 --threads ${t} --processes ${t} 2> ../${SAMPLE}.haphic_N_chrs${CHR_ATTEMPT}.log" (line 312)

~~~ [+15m] HapHiC for Fungus_NoRefNoJuice with 14 chrs failed, trying: 12 ~~~

❌ Command failed in in /90daydata/coffea_pangenome/puzzler_trials/assemblies/fungus_tutorial/Fungus_NoRefNoJuice/03_haphic/haphic: "${PUZZLER} haphic pipeline ${WD}/${SAMPLE}/02_purge_dups/p_ctg.purged.fa ../filtered.bam ${CHR_ATTEMPT} --correct_nrounds 2 --max_inflation 10.0 --threads ${t} --processes ${t} 2> ../${SAMPLE}.haphic_N_chrs${CHR_ATTEMPT}.log" (line 312)

~~~ [+18m] HapHiC for Fungus_NoRefNoJuice with 14 chrs failed, trying: 13 ~~~

❌ Command failed in in /90daydata/coffea_pangenome/puzzler_trials/assemblies/fungus_tutorial/Fungus_NoRefNoJuice/03_haphic/haphic: "${PUZZLER} haphic pipeline ${WD}/${SAMPLE}/02_purge_dups/p_ctg.purged.fa ../filtered.bam ${CHR_ATTEMPT} --correct_nrounds 2 --max_inflation 10.0 --threads ${t} --processes ${t} 2> ../${SAMPLE}.haphic_N_chrs${CHR_ATTEMPT}.log" (line 312)

~~~ [+20m] HapHiC for Fungus_NoRefNoJuice with 14 chrs failed, trying: 14 ~~~

❌ Command failed in in /90daydata/coffea_pangenome/puzzler_trials/assemblies/fungus_tutorial/Fungus_NoRefNoJuice/03_haphic/haphic: "${PUZZLER} haphic pipeline ${WD}/${SAMPLE}/02_purge_dups/p_ctg.purged.fa ../filtered.bam ${CHR_ATTEMPT} --correct_nrounds 2 --max_inflation 10.0 --threads ${t} --processes ${t} 2> ../${SAMPLE}.haphic_N_chrs${CHR_ATTEMPT}.log" (line 312)

~~~ [+22m] HapHiC for Fungus_NoRefNoJuice with 14 chrs failed, trying: 15 ~~~

❌ Command failed in in /90daydata/coffea_pangenome/puzzler_trials/assemblies/fungus_tutorial/Fungus_NoRefNoJuice/03_haphic/haphic: "${PUZZLER} haphic pipeline ${WD}/${SAMPLE}/02_purge_dups/p_ctg.purged.fa ../filtered.bam ${CHR_ATTEMPT} --correct_nrounds 2 --max_inflation 10.0 --threads ${t} --processes ${t} 2> ../${SAMPLE}.haphic_N_chrs${CHR_ATTEMPT}.log" (line 312)

~~~ [+25m] HapHiC for Fungus_NoRefNoJuice with 14 chrs failed, trying: 16 ~~~

❌ Command failed in in /90daydata/coffea_pangenome/puzzler_trials/assemblies/fungus_tutorial/Fungus_NoRefNoJuice/03_haphic/haphic: "${PUZZLER} haphic pipeline ${WD}/${SAMPLE}/02_purge_dups/p_ctg.purged.fa ../filtered.bam ${CHR_ATTEMPT} --correct_nrounds 2 --max_inflation 10.0 --threads ${t} --processes ${t} 2> ../${SAMPLE}.haphic_N_chrs${CHR_ATTEMPT}.log" (line 312)

~~~ [+27m] HapHiC for Fungus_NoRefNoJuice with 14 chrs failed, trying: 17 ~~~

❌ Command failed in in /90daydata/coffea_pangenome/puzzler_trials/assemblies/fungus_tutorial/Fungus_NoRefNoJuice/03_haphic/haphic: "${PUZZLER} haphic pipeline ${WD}/${SAMPLE}/02_purge_dups/p_ctg.purged.fa ../filtered.bam ${CHR_ATTEMPT} --correct_nrounds 2 --max_inflation 10.0 --threads ${t} --processes ${t} 2> ../${SAMPLE}.haphic_N_chrs${CHR_ATTEMPT}.log" (line 312)

~~~ [+29m] HapHiC for Fungus_NoRefNoJuice with 14 chrs failed, trying: 18 ~~~

❌ Command failed in in /90daydata/coffea_pangenome/puzzler_trials/assemblies/fungus_tutorial/Fungus_NoRefNoJuice/03_haphic/haphic: "${PUZZLER} haphic pipeline ${WD}/${SAMPLE}/02_purge_dups/p_ctg.purged.fa ../filtered.bam ${CHR_ATTEMPT} --correct_nrounds 2 --max_inflation 10.0 --threads ${t} --processes ${t} 2> ../${SAMPLE}.haphic_N_chrs${CHR_ATTEMPT}.log" (line 312)

~~~ [+31m] HapHiC for Fungus_NoRefNoJuice failed, scaffolding with YAHS instead ~~~
~~~ [+32m] Skipping juicer HiC file creation for Fungus_NoRefNoJuice: --no_juice given, not recommended! ~~~
~~~ [+32m] No reference provided for Fungus_NoRefNoJuice: simply extracting assembly to: /90daydata/coffea_pangenome/puzzler_trials/assemblies/fungus_tutorial/primary_asm/Fungus_NoRefNoJuice.fa ~~~
~~~ [+32m] Skipping juicer extraction for Fungus_NoRefNoJuice, --no_juice requested! (not recommended!) ~~~
~~~ [+32m] No reference provided for Fungus_NoRefNoJuice: no chromosome re-naming, so final assembly already: /90daydata/coffea_pangenome/puzzler_trials/assemblies/fungus_tutorial/primary_asm/Fungus_NoRefNoJuice.fa ~~~
~~~ [+32m] Assessing genome quality for Fungus_NoRefNoJuice ~~~
~~~ [+32m] Creating final HiC bam for Fungus_NoRefNoJuice ~~~
~~~ [+34m] Creating final contact map for Fungus_NoRefNoJuice, no ref specified, plotting scaffolds > 2mb ~~~
~~~ [+34m] Running YAK on Fungus_NoRefNoJuice ~~~
~~~ [+35m] BUSCO lineage dataset already exists, skipping ~~~
~~~ [+35m] Running BUSCO for Fungus_NoRefNoJuice using lineage: fungi_odb10 ~~~
~~~ [+46m] Skipping blobtools for Fungus_NoRefNoJuice, not desired ~~~
~~~ [+46m] Summarizing Assembly for Fungus_NoRefNoJuice ~~~
~~~ [+47m] Your final assembly for Fungus_NoRefNoJuice is: /90daydata/coffea_pangenome/puzzler_trials/assemblies/fungus_tutorial/primary_asm/Fungus_NoRefNoJuice.fa ~~~
~~~ [+47m] Your final assembly stats for Fungus_NoRefNoJuice are in: /90daydata/coffea_pangenome/puzzler_trials/assemblies/fungus_tutorial/primary_asm/stats/Fungus_NoRefNoJuice.summary.txt ~~~
```

Runtime (47 minutes including BUSCO):

```
Cores per node: 16
CPU Utilized: 05:34:13
CPU Efficiency: 44.38% of 12:33:04 core-walltime
Job Wall-clock time: 00:47:04
Memory Utilized: 53.57 GB
Memory Efficiency: 83.70% of 64.00 GB (64.00 GB/node)
```

## Outputs

Within `$WD/primary_asm/stats/`, we have our summary statistics.

The samples which have a `$REFERENCE` have 2 additional columns: `WithinChrsBp` and `PropWithinChrs`, showing how much of the assembly is within named chromosomes:

```bash
awk 'NR == 1 || FNR > 1' Fungus_Ref*summary.txt
Sample  SizeBP  Sequences       Contigs Gaps    ContigN50       ScafN50 WithinChrsBP    PropWithinChrs  BUSCO_Complete  BUSCO_singlecopy        YAK_CV  YAK_QV
Fungus_RefNoJuice        24816976        23      23     0        1776619         1776619         24117944       0.9718  99.1    99.1    0.999   69.904
Fungus_Ref       24817076        22      23     1        1776619         1806824         24324650       0.9802  99.1    99.1    0.999   69.904
```

| Sample            | SizeBP   | Sequences | Contigs | Gaps | ContigN50 | ScafN50 | WithinChrsBP | PropWithinChrs | BUSCO_Complete | BUSCO_singlecopy | YAK_CV | YAK_QV |
| ----------------- | -------- | --------- | ------- | ---- | --------- | ------- | ------------ | -------------- | -------------- | ---------------- | ------ | ------ |
| Fungus_RefNoJuice | 24816976 | 23        | 23      | 0    | 1776619   | 1776619 | 24117944     | 0.9718         | 99.1           | 99.1             | 0.999  | 69.904 |
| Fungus_Ref        | 24817076 | 22        | 23      | 1    | 1776619   | 1806824 | 24324650     | 0.9802         | 99.1           | 99.1             | 0.999  | 69.904 |

While the assemblies with no `$REFERENCE` do not have those columns:

```bash
awk 'NR == 1 || FNR > 1' Fungus_NoRef*summary.txt
Sample  SizeBP  Sequences       Contigs Gaps    ContigN50       ScafN50 BUSCO_Complete  BUSCO_singlecopy        YAK_CV  YAK_QV
Fungus_NoRefNoJuice      24816976        23      23     0        1776619         1776619        99.5    99.5    0.999   69.904
Fungus_NoRef     24817076        22      23     1        1776619         1806824        99.5    99.5    0.999   69.904
```

| Sample              | SizeBP   | Sequences | Contigs | Gaps | ContigN50 | ScafN50 | BUSCO_Complete | BUSCO_singlecopy | YAK_CV | YAK_QV |
| ------------------- | -------- | --------- | ------- | ---- | --------- | ------- | -------------- | ---------------- | ------ | ------ |
| Fungus_NoRefNoJuice | 24816976 | 23        | 23      | 0    | 1776619   | 1776619 | 99.5           | 99.5             | 0.999  | 69.904 |
| Fungus_NoRef        | 24817076 | 22        | 23      | 1    | 1776619   | 1806824 | 99.5           | 99.5             | 0.999  | 69.904 |

Overall, all assemblies look good. 

![fungus_assemblies](/examples/figs/Fungus_Assemblies_HiC.png)

## Alignments

I will perform a quick WGA and create a dotplot between the assemblies. Since the no `$REFERENCE` assemblies do not have named chromosomes, I will just align any sequences > 200 Kb:

```bash
PUZZLER="apptainer exec /home/justin.merondun/apptainer/puzzler_v1.8.sif"

for i in Fungus_Ref Fungus_NoRef Fungus_RefNoJuice Fungus_NoRefNoJuice RhiKalk1.hap1; do 

    echo "Extracting chromosomes for ${i}"
    $PUZZLER samtools faidx ${i}.fa
    #egrep 'chr|Chr' ${i}.fa.fai | egrep -v 'SUPER|unloc' | awk '{print $1}' > $i.tmp 
    awk '$2 > 2e5' $i.fa.fai | awk '{print $1}' > $i.tmp
    $PUZZLER seqtk subseq ${i}.fa $i.tmp > ../../alignments/fastas/${i}.fa

done
```

I will compare each assembly to the published reference:

```bash
for i in Fungus_Ref Fungus_NoRef Fungus_RefNoJuice Fungus_NoRefNoJuice; do 

    echo "Alignment for ${i}"
    ${PUZZLER} minimap2 -t ${t} -x asm20 --secondary=no -c -D --max-chain-skip 100 --max-chain-iter 1000 --frag yes ${i}.fa RhiKalk1.hap1.fa > fungus_pafs/${i}.paf
    ${PUZZLER} paf2dotplot.R fungus_pafs/${i}.paf -r 1e4 -m 1e4 

done 
```

Alignments:

![fungus_alignments](/examples/figs/Fungus_Assemblies_Alignments.png)