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
#HIFI ERR12760830      HIC ERR12765203  
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
grep -qxF "export PATH=\$PATH:$INSTALL_PATH" ~/.bashrc || echo "export PATH=\$PATH:$INSTALL_PATH" >> ~/.bashrc
export PATH="$PATH:$INSTALL_PATH"
```

And then changed the SLURM scheduler so I can submit the script directly. I also modify the default values for threads and memory so I don't need to specify `--t 48 --mem 384` every time to `puzzler`. 

```
vim puzzler
########### EDIT THIS BLOCK WITH SLURM & APPTAINER/SINGULARITY SETTINGS ############
#SBATCH --time=2-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --mem=256Gb
#SBATCH --partition=atlas
#SBATCH --account=coffea_pangenome

module load apptainer &> /dev/null || true
#module load singularity &> /dev/null || true
SINGULARITY_TMPDIR=$APPTAINER_TMPDIR
########### EDIT THIS BLOCK WITH SLURM & APPTAINER/SINGULARITY SETTINGS ############

# Default values
SAMPLE=""
MAP_FILE=""
t=48
MEM=384
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
cd /project/90daydata/coffea_pangenome/puzzler_trials/assemblies
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

| sample              | runtime   | container                                        | wd                                                    | hifi                                                         | hic_r1                                                       | hic_r2                                                       | num_chrs | reference                                                    | hom_cov | blob_database | busco_lineage | busco_database                                             |
| ------------------- | --------- | ------------------------------------------------ | ----------------------------------------------------- | ------------------------------------------------------------ | ------------------------------------------------------------ | ------------------------------------------------------------ | -------- | ------------------------------------------------------------ | ------- | ------------- | ------------- | ---------------------------------------------------------- |
| Fungus_Ref          | apptainer | /home/justin.merondun/apptainer/puzzler_v1.7.sif | /90daydata/coffea_pangenome/puzzler_trials/assemblies | /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Fungus.HiFi.fastq.gz | /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Fungus.HiC.R1.fastq.gz | /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Fungus.HiC.R2.fastq.gz | 14       | /90daydata/coffea_pangenome/puzzler_trials/raw_data/references/GCA_964106605.1_gdRhiKalk1.hap1.1_genomic.fna | NA      | NA            | fungi_odb10   | /90daydata/coffea_pangenome/puzzler_trials/busco_downloads |
| Fungus_NoRef        | apptainer | /home/justin.merondun/apptainer/puzzler_v1.7.sif | /90daydata/coffea_pangenome/puzzler_trials/assemblies | /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Fungus.HiFi.fastq.gz | /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Fungus.HiC.R1.fastq.gz | /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Fungus.HiC.R2.fastq.gz | 14       | NA                                                           | NA      | NA            | fungi_odb10   | /90daydata/coffea_pangenome/puzzler_trials/busco_downloads |
| Fungus_NoJuice      | apptainer | /home/justin.merondun/apptainer/puzzler_v1.7.sif | /90daydata/coffea_pangenome/puzzler_trials/assemblies | /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Fungus.HiFi.fastq.gz | /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Fungus.HiC.R1.fastq.gz | /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Fungus.HiC.R2.fastq.gz | 14       | /90daydata/coffea_pangenome/puzzler_trials/raw_data/references/GCA_964106605.1_gdRhiKalk1.hap1.1_genomic.fna | NA      | NA            | fungi_odb10   | /90daydata/coffea_pangenome/puzzler_trials/busco_downloads |
| Fungus_NoRefNoJuice | apptainer | /home/justin.merondun/apptainer/puzzler_v1.7.sif | /90daydata/coffea_pangenome/puzzler_trials/assemblies | /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Fungus.HiFi.fastq.gz | /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Fungus.HiC.R1.fastq.gz | /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Fungus.HiC.R2.fastq.gz | 14       | NA                                                           | NA      | NA            | fungi_odb10   | /90daydata/coffea_pangenome/puzzler_trials/busco_downloads |

### Recommended: Juicer & Reference

#### Step 1: Draft Assembly

Submit: `sbatch -J asm_Fungus_Ref puzzler -s Fungus_Ref -m samples.tsv --threads 12 --mem 64` 

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
Parameters for sample: Fungus_Ref 
RUNTIME: apptainer
CONTAINER: /home/justin.merondun/apptainer/puzzler_v1.7.sif 
WD: /90daydata/coffea_pangenome/puzzler_trials/assemblies 
HIFI: /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Fungus.HiFi.fastq.gz
HIC_R1: /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Fungus.HiC.R1.fastq.gz
HIC_R2: /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Fungus.HiC.R2.fastq.gz
NUMBER CHRS: 14
REFERENCE: /90daydata/coffea_pangenome/puzzler_trials/raw_data/references/GCA_964106605.1_gdRhiKalk1.hap1.1_genomic.fna
HOM_COV: NA
BLOB_DB: NA
BUSCO_LINEAGE: fungi_odb10
BUSCO_DB: /90daydata/coffea_pangenome/puzzler_trials/busco_downloads
Cores Requested: 64
Cores Available: 64
RAM Requested: 512
Memory Available: 2241.3 GB
=======================================================================

~~~~ Assembling genome for Fungus_Ref ~~~~
~~~~ Running juicer, script will stop after .hic files created ~~~~
~~~~ Starting hifiasm assembly for Fungus_Ref ~~~~
~~~~ Starting Purge_Dups for Fungus_Ref ~~~~
~~~~ Mapping HiC reads to Fungus_Ref draft ~~~~
~~~~ Running HapHiC for Fungus_Ref  ~~~~
~~~~ HapHiC for Fungus_Ref with 14 chrs failed, trying: 10 ~~~~

❌ Command failed in in /90daydata/coffea_pangenome/puzzler_trials/assemblies/Fungus_NoRefNoJuice/03_haphic/haphic: "${PUZZLER} haphic pipeline ${WD}/${SAMPLE}/02_purge_dups/p_ctg.purged.fa ../filtered.bam ${CHR_ATTEMPT} --correct_nrounds 2 --max_inflation 10.0 --threads ${t} --processes ${t} 2> ../${SAMPLE}.haphic_N_chrs${CHR_ATTEMPT}.log" (line 257)

~~~~ HapHiC for Fungus_Ref with 14 chrs failed, trying: 11 ~~~~

❌ Command failed in in /90daydata/coffea_pangenome/puzzler_trials/assemblies/Fungus_NoRefNoJuice/03_haphic/haphic: "${PUZZLER} haphic pipeline ${WD}/${SAMPLE}/02_purge_dups/p_ctg.purged.fa ../filtered.bam ${CHR_ATTEMPT} --correct_nrounds 2 --max_inflation 10.0 --threads ${t} --processes ${t} 2> ../${SAMPLE}.haphic_N_chrs${CHR_ATTEMPT}.log" (line 257)

~~~~ HapHiC for Fungus_Ref with 14 chrs failed, trying: 12 ~~~~

❌ Command failed in in /90daydata/coffea_pangenome/puzzler_trials/assemblies/Fungus_NoRefNoJuice/03_haphic/haphic: "${PUZZLER} haphic pipeline ${WD}/${SAMPLE}/02_purge_dups/p_ctg.purged.fa ../filtered.bam ${CHR_ATTEMPT} --correct_nrounds 2 --max_inflation 10.0 --threads ${t} --processes ${t} 2> ../${SAMPLE}.haphic_N_chrs${CHR_ATTEMPT}.log" (line 257)

~~~~ HapHiC for Fungus_Ref with 14 chrs failed, trying: 13 ~~~~

❌ Command failed in in /90daydata/coffea_pangenome/puzzler_trials/assemblies/Fungus_NoRefNoJuice/03_haphic/haphic: "${PUZZLER} haphic pipeline ${WD}/${SAMPLE}/02_purge_dups/p_ctg.purged.fa ../filtered.bam ${CHR_ATTEMPT} --correct_nrounds 2 --max_inflation 10.0 --threads ${t} --processes ${t} 2> ../${SAMPLE}.haphic_N_chrs${CHR_ATTEMPT}.log" (line 257)

~~~~ HapHiC for Fungus_Ref with 14 chrs failed, trying: 14 ~~~~

❌ Command failed in in /90daydata/coffea_pangenome/puzzler_trials/assemblies/Fungus_NoRefNoJuice/03_haphic/haphic: "${PUZZLER} haphic pipeline ${WD}/${SAMPLE}/02_purge_dups/p_ctg.purged.fa ../filtered.bam ${CHR_ATTEMPT} --correct_nrounds 2 --max_inflation 10.0 --threads ${t} --processes ${t} 2> ../${SAMPLE}.haphic_N_chrs${CHR_ATTEMPT}.log" (line 257)

~~~~ HapHiC for Fungus_Ref with 14 chrs failed, trying: 15 ~~~~

❌ Command failed in in /90daydata/coffea_pangenome/puzzler_trials/assemblies/Fungus_NoRefNoJuice/03_haphic/haphic: "${PUZZLER} haphic pipeline ${WD}/${SAMPLE}/02_purge_dups/p_ctg.purged.fa ../filtered.bam ${CHR_ATTEMPT} --correct_nrounds 2 --max_inflation 10.0 --threads ${t} --processes ${t} 2> ../${SAMPLE}.haphic_N_chrs${CHR_ATTEMPT}.log" (line 257)

~~~~ HapHiC for Fungus_Ref with 14 chrs failed, trying: 16 ~~~~

❌ Command failed in in /90daydata/coffea_pangenome/puzzler_trials/assemblies/Fungus_NoRefNoJuice/03_haphic/haphic: "${PUZZLER} haphic pipeline ${WD}/${SAMPLE}/02_purge_dups/p_ctg.purged.fa ../filtered.bam ${CHR_ATTEMPT} --correct_nrounds 2 --max_inflation 10.0 --threads ${t} --processes ${t} 2> ../${SAMPLE}.haphic_N_chrs${CHR_ATTEMPT}.log" (line 257)

~~~~ HapHiC for Fungus_Ref with 14 chrs failed, trying: 17 ~~~~

❌ Command failed in in /90daydata/coffea_pangenome/puzzler_trials/assemblies/Fungus_NoRefNoJuice/03_haphic/haphic: "${PUZZLER} haphic pipeline ${WD}/${SAMPLE}/02_purge_dups/p_ctg.purged.fa ../filtered.bam ${CHR_ATTEMPT} --correct_nrounds 2 --max_inflation 10.0 --threads ${t} --processes ${t} 2> ../${SAMPLE}.haphic_N_chrs${CHR_ATTEMPT}.log" (line 257)

~~~~ HapHiC for Fungus_Ref with 14 chrs failed, trying: 18 ~~~~

❌ Command failed in in /90daydata/coffea_pangenome/puzzler_trials/assemblies/Fungus_NoRefNoJuice/03_haphic/haphic: "${PUZZLER} haphic pipeline ${WD}/${SAMPLE}/02_purge_dups/p_ctg.purged.fa ../filtered.bam ${CHR_ATTEMPT} --correct_nrounds 2 --max_inflation 10.0 --threads ${t} --processes ${t} 2> ../${SAMPLE}.haphic_N_chrs${CHR_ATTEMPT}.log" (line 257)

~~~~ HapHiC for Fungus_Ref failed, scaffolding with YAHS instead ~~~~
~~~~ Creating .hic file for juicebox for Fungus_Ref ~~~~
~~~~ Genome size is 0GB, running default minimap2 command ~~~~
~~~~ Post curation assembly file missing for Fungus_Ref: Run Juicebox & place in /90daydata/coffea_pangenome/puzzler_trials/assemblies/juicer_files/Fungus_Ref_JBAT.review.assembly ~~~~
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

It completed the assembly in 24 minutes, and of course with this tiny genome we could have dramatically reduced our resources to `--threads 8 --mem 32` 

``` 
Nodes: 1
Cores per node: 64
CPU Utilized: 04:31:37
CPU Efficiency: 17.08% of 1-02:30:24 core-walltime
Job Wall-clock time: 00:24:51
Memory Utilized: 22.98 GB
Memory Efficiency: 4.49% of 512.00 GB
```

In any case, the pipeline will also run `haphic refsort`, provided a related species is provided, which will align your draft scaffolds in order of the reference genome for juicebox curation. **This is not reference-guided assembly**, and will not rename or otherwise modify the scaffolds. Please see the [HapHiC](https://github.com/zengxiaofei/HapHiC) `refsort` description for more details. 

If you did not provide a `$REFERENCE` file, it will still output the `.hic` and `.assembly` files, just without running `haphic refsort`. 

#### Step 2: Juicebox Curation

Now we are ready to load the `$WD/juicer_files/Fungus_Ref_JBAT.hic` and `$WD/juicer_files/Fungus_Ref_JBAT.assembly` file into juicebox. 

Perform our edits. Pull in the `.hic` file, import the `.assembly` file with `Assembly > Import Map Assembly`. There's very minimal edits on this genome, I just merge two blocks which originate from one chromosome: 

![Fungus_Edits_Juicebox]()

Afterwards, save the file with `Assembly > Export Assembly`, and save the  `Fungus_Ref_JBAT.review.assembly` to the directory `$WD/juicer_files`, retaining the default file name! 

#### Step 3: Chromosome Renaming, Remapping HiC

Simply submit `puzzler` with the same parameters as before, and the script will assign chromosome names, output a final HiC contact map on the final assembly, and output assembly statistics. **I skip blobtools, because that will take several hours to blast against the nt database...**. 

We simply resubmit the exact same command, and the script will take up where it left off: 

`sbatch -J asm_Fungus_Ref puzzler -s Fungus_Ref -m samples.tsv` 

```bash
cat slurm-15728608.out

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
CONTAINER: /home/justin.merondun/apptainer/puzzler_v1.7.sif 
WD: /90daydata/coffea_pangenome/puzzler_trials/assemblies 
HIFI: /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Fungus.HiFi.fastq.gz
HIC_R1: /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Fungus.HiC.R1.fastq.gz
HIC_R2: /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Fungus.HiC.R2.fastq.gz
NUMBER CHRS: 14
REFERENCE: /90daydata/coffea_pangenome/puzzler_trials/raw_data/references/GCA_964106605.1_gdRhiKalk1.hap1.1_genomic.fna
HOM_COV: NA
BLOB_DB: NA
BUSCO_LINEAGE: fungi_odb10
BUSCO_DB: /90daydata/coffea_pangenome/puzzler_trials/busco_downloads
Cores Requested: 64
Cores Available: 64
RAM Requested: 512
Memory Available: 2247.0 GB
=======================================================================

~~~~ Assembling genome for Fungus_Ref ~~~~
~~~~ Running juicer, script will stop after .hic files created ~~~~
~~~~ Skipping hifiasm for Fungus_Ref: /90daydata/coffea_pangenome/puzzler_trials/assemblies/Fungus_Ref/01_hifiasm/asm.hic.p_ctg.gfa exists ~~~~
~~~~ Skipping purge for Fungus_Ref: /90daydata/coffea_pangenome/puzzler_trials/assemblies/Fungus_Ref/02_purge_dups/p_ctg.purged.fa exists ~~~~
~~~~ Skipping HiC alignment for Fungus_Ref: /90daydata/coffea_pangenome/puzzler_trials/assemblies/Fungus_Ref/03_haphic/filtered.bam exists ~~~~
~~~~ Skipping HapHiC for Fungus_Ref: /90daydata/coffea_pangenome/puzzler_trials/assemblies/Fungus_Ref/03_haphic/haphic/04.build/scaffolds.fa exists ~~~~
~~~~ Skipping juicer HiC file creation for Fungus_Ref: /90daydata/coffea_pangenome/puzzler_trials/assemblies/juicer_files/Fungus_Ref_JBAT.hic exists ~~~~
~~~~ Extracting post-curation assembly and mapping to reference for Fungus_Ref ~~~~
~~~~ Renaming chromosomes for Fungus_Ref ~~~~
~~~~ Scaffold sanity check passed for renaming, proceeding! ~~~~
~~~~ Single scaffolds corresponding to a single Chr for Fungus_Ref ~~~~
~~~~ Assessing genome quality for Fungus_Ref ~~~~
~~~~ Creating final HiC bam for Fungus_Ref ~~~~
~~~~ Creating final contact map for Fungus_Ref, plotting named chromosomes ~~~~
~~~~ Running YAK on Fungus_Ref ~~~~
~~~~ BUSCO lineage dataset already exists, skipping ~~~~
~~~~ Running BUSCO for Fungus_Ref using lineage: fungi_odb10 ~~~~
~~~~ Skipping blobtools for Fungus_Ref, not desired ~~~~
~~~~ Summarizing Assembly for Fungus_Ref ~~~~
~~~~ Your final assembly for Fungus_Ref is: /90daydata/coffea_pangenome/puzzler_trials/assemblies/primary_asm/Fungus_Ref.fa ~~~~
~~~~ Your final assembly stats for Fungus_Ref are in: /90daydata/coffea_pangenome/puzzler_trials/assemblies/primary_asm/stats/Fungus_Ref.summary.txt ~~~~
```

Which only ran for about 10 minutes (note I left off `--threads 12 and --mem 64`, so this would typically take a big longer) 

```
Cores per node: 64
CPU Utilized: 02:23:10
CPU Efficiency: 27.85% of 08:34:08 core-walltime
Job Wall-clock time: 00:08:02
Memory Utilized: 53.61 GB
Memory Efficiency: 10.47% of 512.00 GB
```

The script skips all the previously run steps, because we have the `.complete` files in the `$WD/Fungus_Ref` dir, as well as any already-completed step-specific files (e.g. `$WD/juicer_files/$SAMPLE_JBAT.review.assembly` for post-curation assembly and renaming). 

```
ls -l $WD/Fungus_Ref
01_hifiasm/
02_purge_dups/
03_haphic/
04_juicer/
05_postjuicebox/
06_realign_hic_hifi/
07_busco_yak_blob/
align_hic.complete
hifiasm.complete
juicer.complete
purge_dups.complete
qc_align_hic.complete
scaffolding.complete
yak.complete
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

Submit: `sbatch -J asm_Fungus_NoRef puzzler -s Fungus_NoRef -m samples.tsv --threads 12 --mem 64` 

```
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
CONTAINER: /home/justin.merondun/apptainer/puzzler_v1.7.sif 
WD: /90daydata/coffea_pangenome/puzzler_trials/assemblies 
HIFI: /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Fungus.HiFi.fastq.gz
HIC_R1: /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Fungus.HiC.R1.fastq.gz
HIC_R2: /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Fungus.HiC.R2.fastq.gz
NUMBER CHRS: 14
REFERENCE: NA
HOM_COV: NA
BLOB_DB: NA
BUSCO_LINEAGE: fungi_odb10
BUSCO_DB: /90daydata/coffea_pangenome/puzzler_trials/busco_downloads
Cores Requested: 64
Cores Available: 64
RAM Requested: 512
Memory Available: 2241.3 GB
=======================================================================

~~~~ Assembling genome for Fungus_NoRef ~~~~
~~~~ Running juicer, script will stop after .hic files created ~~~~
~~~~ Starting hifiasm assembly for Fungus_NoRef ~~~~
~~~~ Starting Purge_Dups for Fungus_NoRef ~~~~
~~~~ Mapping HiC reads to Fungus_NoRef draft ~~~~
~~~~ Running HapHiC for Fungus_NoRef  ~~~~
~~~~ HapHiC for Fungus_NoRef with 14 chrs failed, trying: 10 ~~~~

❌ Command failed in in /90daydata/coffea_pangenome/puzzler_trials/assemblies/Fungus_NoRefNoJuice/03_haphic/haphic: "${PUZZLER} haphic pipeline ${WD}/${SAMPLE}/02_purge_dups/p_ctg.purged.fa ../filtered.bam ${CHR_ATTEMPT} --correct_nrounds 2 --max_inflation 10.0 --threads ${t} --processes ${t} 2> ../${SAMPLE}.haphic_N_chrs${CHR_ATTEMPT}.log" (line 257)

~~~~ HapHiC for Fungus_NoRef with 14 chrs failed, trying: 11 ~~~~

❌ Command failed in in /90daydata/coffea_pangenome/puzzler_trials/assemblies/Fungus_NoRefNoJuice/03_haphic/haphic: "${PUZZLER} haphic pipeline ${WD}/${SAMPLE}/02_purge_dups/p_ctg.purged.fa ../filtered.bam ${CHR_ATTEMPT} --correct_nrounds 2 --max_inflation 10.0 --threads ${t} --processes ${t} 2> ../${SAMPLE}.haphic_N_chrs${CHR_ATTEMPT}.log" (line 257)

~~~~ HapHiC for Fungus_NoRef with 14 chrs failed, trying: 12 ~~~~

❌ Command failed in in /90daydata/coffea_pangenome/puzzler_trials/assemblies/Fungus_NoRefNoJuice/03_haphic/haphic: "${PUZZLER} haphic pipeline ${WD}/${SAMPLE}/02_purge_dups/p_ctg.purged.fa ../filtered.bam ${CHR_ATTEMPT} --correct_nrounds 2 --max_inflation 10.0 --threads ${t} --processes ${t} 2> ../${SAMPLE}.haphic_N_chrs${CHR_ATTEMPT}.log" (line 257)

~~~~ HapHiC for Fungus_NoRef with 14 chrs failed, trying: 13 ~~~~

❌ Command failed in in /90daydata/coffea_pangenome/puzzler_trials/assemblies/Fungus_NoRefNoJuice/03_haphic/haphic: "${PUZZLER} haphic pipeline ${WD}/${SAMPLE}/02_purge_dups/p_ctg.purged.fa ../filtered.bam ${CHR_ATTEMPT} --correct_nrounds 2 --max_inflation 10.0 --threads ${t} --processes ${t} 2> ../${SAMPLE}.haphic_N_chrs${CHR_ATTEMPT}.log" (line 257)

~~~~ HapHiC for Fungus_NoRef with 14 chrs failed, trying: 14 ~~~~

❌ Command failed in in /90daydata/coffea_pangenome/puzzler_trials/assemblies/Fungus_NoRefNoJuice/03_haphic/haphic: "${PUZZLER} haphic pipeline ${WD}/${SAMPLE}/02_purge_dups/p_ctg.purged.fa ../filtered.bam ${CHR_ATTEMPT} --correct_nrounds 2 --max_inflation 10.0 --threads ${t} --processes ${t} 2> ../${SAMPLE}.haphic_N_chrs${CHR_ATTEMPT}.log" (line 257)

~~~~ HapHiC for Fungus_NoRef with 14 chrs failed, trying: 15 ~~~~

❌ Command failed in in /90daydata/coffea_pangenome/puzzler_trials/assemblies/Fungus_NoRefNoJuice/03_haphic/haphic: "${PUZZLER} haphic pipeline ${WD}/${SAMPLE}/02_purge_dups/p_ctg.purged.fa ../filtered.bam ${CHR_ATTEMPT} --correct_nrounds 2 --max_inflation 10.0 --threads ${t} --processes ${t} 2> ../${SAMPLE}.haphic_N_chrs${CHR_ATTEMPT}.log" (line 257)

~~~~ HapHiC for Fungus_NoRef with 14 chrs failed, trying: 16 ~~~~

❌ Command failed in in /90daydata/coffea_pangenome/puzzler_trials/assemblies/Fungus_NoRefNoJuice/03_haphic/haphic: "${PUZZLER} haphic pipeline ${WD}/${SAMPLE}/02_purge_dups/p_ctg.purged.fa ../filtered.bam ${CHR_ATTEMPT} --correct_nrounds 2 --max_inflation 10.0 --threads ${t} --processes ${t} 2> ../${SAMPLE}.haphic_N_chrs${CHR_ATTEMPT}.log" (line 257)

~~~~ HapHiC for Fungus_NoRef with 14 chrs failed, trying: 17 ~~~~

❌ Command failed in in /90daydata/coffea_pangenome/puzzler_trials/assemblies/Fungus_NoRefNoJuice/03_haphic/haphic: "${PUZZLER} haphic pipeline ${WD}/${SAMPLE}/02_purge_dups/p_ctg.purged.fa ../filtered.bam ${CHR_ATTEMPT} --correct_nrounds 2 --max_inflation 10.0 --threads ${t} --processes ${t} 2> ../${SAMPLE}.haphic_N_chrs${CHR_ATTEMPT}.log" (line 257)

~~~~ HapHiC for Fungus_NoRef with 14 chrs failed, trying: 18 ~~~~

❌ Command failed in in /90daydata/coffea_pangenome/puzzler_trials/assemblies/Fungus_NoRefNoJuice/03_haphic/haphic: "${PUZZLER} haphic pipeline ${WD}/${SAMPLE}/02_purge_dups/p_ctg.purged.fa ../filtered.bam ${CHR_ATTEMPT} --correct_nrounds 2 --max_inflation 10.0 --threads ${t} --processes ${t} 2> ../${SAMPLE}.haphic_N_chrs${CHR_ATTEMPT}.log" (line 257)

~~~~ HapHiC for Fungus_NoRef failed, scaffolding with YAHS instead ~~~~
~~~~ Creating .hic file for juicebox for Fungus_NoRef ~~~~
~~~~ Post curation assembly file missing for Fungus_NoRef: Run Juicebox & place in /90daydata/coffea_pangenome/puzzler_trials/assemblies/juicer_files/Fungus_NoRef_JBAT.review.assembly ~~~~
```

Ran in 24 minutes, and of course we didn't need nearly that many resources:

```bash
Nodes: 1
Cores per node: 64
CPU Utilized: 04:34:08
CPU Efficiency: 17.20% of 1-02:33:36 core-walltime
Job Wall-clock time: 00:24:54
Memory Utilized: 19.84 GB
Memory Efficiency: 3.87% of 512.00 GB
```

Perform juicebox edits....

![juicebox]()

and then resubmit: `sbatch -J asm_Fungus_NoRef puzzler -s Fungus_NoRef -m samples.tsv --threads 12 --mem 64` 

```
cat slurm-15728609.out

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
CONTAINER: /home/justin.merondun/apptainer/puzzler_v1.7.sif 
WD: /90daydata/coffea_pangenome/puzzler_trials/assemblies 
HIFI: /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Fungus.HiFi.fastq.gz
HIC_R1: /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Fungus.HiC.R1.fastq.gz
HIC_R2: /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Fungus.HiC.R2.fastq.gz
NUMBER CHRS: 14
REFERENCE: NA
HOM_COV: NA
BLOB_DB: NA
BUSCO_LINEAGE: fungi_odb10
BUSCO_DB: /90daydata/coffea_pangenome/puzzler_trials/busco_downloads
Cores Requested: 12
Cores Available: 12
RAM Requested: 64
Memory Available: 2243.4 GB
=======================================================================

~~~~ Assembling genome for Fungus_NoRef ~~~~
~~~~ Running juicer, script will stop after .hic files created ~~~~
~~~~ Skipping hifiasm for Fungus_NoRef: /90daydata/coffea_pangenome/puzzler_trials/assemblies/Fungus_NoRef/01_hifiasm/asm.hic.p_ctg.gfa exists ~~~~
~~~~ Skipping purge for Fungus_NoRef: /90daydata/coffea_pangenome/puzzler_trials/assemblies/Fungus_NoRef/02_purge_dups/p_ctg.purged.fa exists ~~~~
~~~~ Skipping HiC alignment for Fungus_NoRef: /90daydata/coffea_pangenome/puzzler_trials/assemblies/Fungus_NoRef/03_haphic/filtered.bam exists ~~~~
~~~~ Skipping HapHiC for Fungus_NoRef: /90daydata/coffea_pangenome/puzzler_trials/assemblies/Fungus_NoRef/03_haphic/haphic/04.build/scaffolds.fa exists ~~~~
~~~~ Skipping juicer HiC file creation for Fungus_NoRef: /90daydata/coffea_pangenome/puzzler_trials/assemblies/juicer_files/Fungus_NoRef_JBAT.hic exists ~~~~
~~~~ No reference provided for Fungus_NoRef: simply extracting assembly to: /90daydata/coffea_pangenome/puzzler_trials/assemblies/primary_asm/Fungus_NoRef.fa ~~~~
~~~~ No reference provided for Fungus_NoRef: no chromosome re-naming, so final assembly already: /90daydata/coffea_pangenome/puzzler_trials/assemblies/primary_asm/Fungus_NoRef.fa ~~~~
~~~~ Assessing genome quality for Fungus_NoRef ~~~~
~~~~ Creating final HiC bam for Fungus_NoRef ~~~~
~~~~ Creating final contact map for Fungus_NoRef, no ref specified, plotting scaffolds > 2mb ~~~~
~~~~ Running YAK on Fungus_NoRef ~~~~
~~~~ BUSCO lineage dataset already exists, skipping ~~~~
~~~~ Running BUSCO for Fungus_NoRef using lineage: fungi_odb10 ~~~~
~~~~ Skipping blobtools for Fungus_NoRef, not desired ~~~~
~~~~ Summarizing Assembly for Fungus_NoRef ~~~~
~~~~ Your final assembly for Fungus_NoRef is: /90daydata/coffea_pangenome/puzzler_trials/assemblies/primary_asm/Fungus_NoRef.fa ~~~~
~~~~ Your final assembly stats for Fungus_NoRef are in: /90daydata/coffea_pangenome/puzzler_trials/assemblies/primary_asm/stats/Fungus_NoRef.summary.txt ~~~~
```

With only 12 cores, took around 15 minutes:

```
Nodes: 1
Cores per node: 12
CPU Utilized: 02:12:31
CPU Efficiency: 72.81% of 03:02:00 core-walltime
Job Wall-clock time: 00:15:10
Memory Utilized: 49.71 GB
Memory Efficiency: 77.68% of 64.00 GB
```





### No Juicer - With Reference

Just add the `--no_juice` flag to skip juicer. Not recommended! 

Submit: `sbatch -J asm_Fungus_RefNoJuice puzzler --no_juice -s Fungus_RefNoJuice -m samples.tsv --threads 12 --mem 64` 

```
cat slurm-15728087.out

=======================================================================
__________ ____ _______________________.____     _____________________ 
\______   \    |   \____    /\____    /|    |    \_   _____/\______   \
 |     ___/    |   / /     /   /     / |    |     |    __)_  |       _/
 |    |   |    |  / /     /_  /     /_ |    |___  |        \ |    |   \
 |____|   |______/ /_______ \/_______ \|_______ \/_______  / |____|_  /
                           \/        \/        \/        \/         \/ 
=======================================================================

=======================================================================
Parameters for sample: Fungus_RefNoJuice 
RUNTIME: apptainer
CONTAINER: /home/justin.merondun/apptainer/puzzler_v1.7.sif 
WD: /90daydata/coffea_pangenome/puzzler_trials/assemblies 
HIFI: /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Fungus.HiFi.fastq.gz
HIC_R1: /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Fungus.HiC.R1.fastq.gz
HIC_R2: /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Fungus.HiC.R2.fastq.gz
NUMBER CHRS: 14
REFERENCE: /90daydata/coffea_pangenome/puzzler_trials/raw_data/references/GCA_964106605.1_gdRhiKalk1.hap1.1_genomic.fna
HOM_COV: NA
BLOB_DB: NA
BUSCO_LINEAGE: fungi_odb10
BUSCO_DB: /90daydata/coffea_pangenome/puzzler_trials/busco_downloads
Cores Requested: 12
Cores Available: 12
RAM Requested: 64
Memory Available: 1492.1 GB
=======================================================================

~~~~ Assembling genome for Fungus_RefNoJuice ~~~~
~~~~ Skipping juicer, no manual curation (not recommended!) ~~~~
~~~~ Starting hifiasm assembly for Fungus_RefNoJuice ~~~~
~~~~ Starting Purge_Dups for Fungus_RefNoJuice ~~~~
~~~~ Mapping HiC reads to Fungus_RefNoJuice draft ~~~~
~~~~ Running HapHiC for Fungus_RefNoJuice  ~~~~

❌ Command failed in in /90daydata/coffea_pangenome/puzzler_trials/assemblies/Fungus_RefNoJuice/03_haphic/haphic: "${PUZZLER} haphic pipeline ${WD}/${SAMPLE}/02_purge_dups/p_ctg.purged.fa ../filtered.bam ${NUM_CHRS} --correct_nrounds 2 --max_inflation 20.0 --threads ${t} --processes ${t} 2> ../${SAMPLE}.haphic.log" (line 234)

~~~~ HapHiC for Fungus_RefNoJuice with 14 chrs failed, trying: 10 ~~~~

❌ Command failed in in /90daydata/coffea_pangenome/puzzler_trials/assemblies/Fungus_RefNoJuice/03_haphic/haphic: "${PUZZLER} haphic pipeline ${WD}/${SAMPLE}/02_purge_dups/p_ctg.purged.fa ../filtered.bam ${CHR_ATTEMPT} --correct_nrounds 2 --max_inflation 10.0 --threads ${t} --processes ${t} 2> ../${SAMPLE}.haphic_N_chrs${CHR_ATTEMPT}.log" (line 257)

~~~~ HapHiC for Fungus_RefNoJuice with 14 chrs failed, trying: 11 ~~~~

❌ Command failed in in /90daydata/coffea_pangenome/puzzler_trials/assemblies/Fungus_RefNoJuice/03_haphic/haphic: "${PUZZLER} haphic pipeline ${WD}/${SAMPLE}/02_purge_dups/p_ctg.purged.fa ../filtered.bam ${CHR_ATTEMPT} --correct_nrounds 2 --max_inflation 10.0 --threads ${t} --processes ${t} 2> ../${SAMPLE}.haphic_N_chrs${CHR_ATTEMPT}.log" (line 257)

~~~~ HapHiC for Fungus_RefNoJuice with 14 chrs failed, trying: 12 ~~~~

❌ Command failed in in /90daydata/coffea_pangenome/puzzler_trials/assemblies/Fungus_RefNoJuice/03_haphic/haphic: "${PUZZLER} haphic pipeline ${WD}/${SAMPLE}/02_purge_dups/p_ctg.purged.fa ../filtered.bam ${CHR_ATTEMPT} --correct_nrounds 2 --max_inflation 10.0 --threads ${t} --processes ${t} 2> ../${SAMPLE}.haphic_N_chrs${CHR_ATTEMPT}.log" (line 257)

~~~~ HapHiC for Fungus_RefNoJuice with 14 chrs failed, trying: 13 ~~~~

❌ Command failed in in /90daydata/coffea_pangenome/puzzler_trials/assemblies/Fungus_RefNoJuice/03_haphic/haphic: "${PUZZLER} haphic pipeline ${WD}/${SAMPLE}/02_purge_dups/p_ctg.purged.fa ../filtered.bam ${CHR_ATTEMPT} --correct_nrounds 2 --max_inflation 10.0 --threads ${t} --processes ${t} 2> ../${SAMPLE}.haphic_N_chrs${CHR_ATTEMPT}.log" (line 257)

~~~~ HapHiC for Fungus_RefNoJuice with 14 chrs failed, trying: 14 ~~~~

❌ Command failed in in /90daydata/coffea_pangenome/puzzler_trials/assemblies/Fungus_RefNoJuice/03_haphic/haphic: "${PUZZLER} haphic pipeline ${WD}/${SAMPLE}/02_purge_dups/p_ctg.purged.fa ../filtered.bam ${CHR_ATTEMPT} --correct_nrounds 2 --max_inflation 10.0 --threads ${t} --processes ${t} 2> ../${SAMPLE}.haphic_N_chrs${CHR_ATTEMPT}.log" (line 257)

~~~~ HapHiC for Fungus_RefNoJuice with 14 chrs failed, trying: 15 ~~~~

❌ Command failed in in /90daydata/coffea_pangenome/puzzler_trials/assemblies/Fungus_RefNoJuice/03_haphic/haphic: "${PUZZLER} haphic pipeline ${WD}/${SAMPLE}/02_purge_dups/p_ctg.purged.fa ../filtered.bam ${CHR_ATTEMPT} --correct_nrounds 2 --max_inflation 10.0 --threads ${t} --processes ${t} 2> ../${SAMPLE}.haphic_N_chrs${CHR_ATTEMPT}.log" (line 257)

~~~~ HapHiC for Fungus_RefNoJuice with 14 chrs failed, trying: 16 ~~~~

❌ Command failed in in /90daydata/coffea_pangenome/puzzler_trials/assemblies/Fungus_RefNoJuice/03_haphic/haphic: "${PUZZLER} haphic pipeline ${WD}/${SAMPLE}/02_purge_dups/p_ctg.purged.fa ../filtered.bam ${CHR_ATTEMPT} --correct_nrounds 2 --max_inflation 10.0 --threads ${t} --processes ${t} 2> ../${SAMPLE}.haphic_N_chrs${CHR_ATTEMPT}.log" (line 257)

~~~~ HapHiC for Fungus_RefNoJuice with 14 chrs failed, trying: 17 ~~~~

❌ Command failed in in /90daydata/coffea_pangenome/puzzler_trials/assemblies/Fungus_RefNoJuice/03_haphic/haphic: "${PUZZLER} haphic pipeline ${WD}/${SAMPLE}/02_purge_dups/p_ctg.purged.fa ../filtered.bam ${CHR_ATTEMPT} --correct_nrounds 2 --max_inflation 10.0 --threads ${t} --processes ${t} 2> ../${SAMPLE}.haphic_N_chrs${CHR_ATTEMPT}.log" (line 257)

~~~~ HapHiC for Fungus_RefNoJuice with 14 chrs failed, trying: 18 ~~~~

❌ Command failed in in /90daydata/coffea_pangenome/puzzler_trials/assemblies/Fungus_RefNoJuice/03_haphic/haphic: "${PUZZLER} haphic pipeline ${WD}/${SAMPLE}/02_purge_dups/p_ctg.purged.fa ../filtered.bam ${CHR_ATTEMPT} --correct_nrounds 2 --max_inflation 10.0 --threads ${t} --processes ${t} 2> ../${SAMPLE}.haphic_N_chrs${CHR_ATTEMPT}.log" (line 257)

~~~~ HapHiC for Fungus_RefNoJuice failed, scaffolding with YAHS instead ~~~~
~~~~ Skipping juicer HiC file creation for Fungus_RefNoJuice: --no_juice given, not recommended! ~~~~
~~~~ Skipping juicer extraction for Fungus_RefNoJuice, --no_juice requested! (not recommended!) ~~~~
~~~~ Renaming chromosomes for Fungus_RefNoJuice ~~~~
~~~~ Scaffold sanity check passed for renaming, proceeding! ~~~~
~~~~ Multiple scaffolds corresponding to single Chr for Fungus_RefNoJuice, Renaming them e.g. Chr1A, Chr1B.. ~~~~
~~~~ Assessing genome quality for Fungus_RefNoJuice ~~~~
~~~~ Creating final HiC bam for Fungus_RefNoJuice ~~~~
~~~~ Creating final contact map for Fungus_RefNoJuice, plotting named chromosomes ~~~~
~~~~ Running YAK on Fungus_RefNoJuice ~~~~
~~~~ BUSCO lineage dataset already exists, skipping ~~~~
~~~~ Running BUSCO for Fungus_RefNoJuice using lineage: fungi_odb10 ~~~~
~~~~ Skipping blobtools for Fungus_RefNoJuice, not desired ~~~~
~~~~ Summarizing Assembly for Fungus_RefNoJuice ~~~~
~~~~ Your final assembly for Fungus_RefNoJuice is: /90daydata/coffea_pangenome/puzzler_trials/assemblies/primary_asm/Fungus_RefNoJuice.fa ~~~~
~~~~ Your final assembly stats for Fungus_RefNoJuice are in: /90daydata/coffea_pangenome/puzzler_trials/assemblies/primary_asm/stats/Fungus_RefNoJuice.summary.txt ~~~~
```

Runtime:

```
Nodes: 1
Cores per node: 12
CPU Utilized: 05:27:02
CPU Efficiency: 51.40% of 10:36:12 core-walltime
Job Wall-clock time: 00:53:01
Memory Utilized: 51.74 GB
Memory Efficiency: 80.85% of 64.00 GB
```





### No Juicer - No Reference

Just add the `--no_juice` flag to skip juicer. Not recommended! 

Submit: `sbatch -J asm_Fungus_NoRefNoJuice puzzler --no_juice -s Fungus_NoRefNoJuice -m samples.tsv --threads 12 --mem 64` 

```
cat slurm-15728050.out

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
CONTAINER: /home/justin.merondun/apptainer/puzzler_v1.7.sif 
WD: /90daydata/coffea_pangenome/puzzler_trials/assemblies 
HIFI: /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Fungus.HiFi.fastq.gz
HIC_R1: /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Fungus.HiC.R1.fastq.gz
HIC_R2: /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Fungus.HiC.R2.fastq.gz
NUMBER CHRS: 14
REFERENCE: NA
HOM_COV: NA
BLOB_DB: NA
BUSCO_LINEAGE: fungi_odb10
BUSCO_DB: /90daydata/coffea_pangenome/puzzler_trials/busco_downloads
Cores Requested: 12
Cores Available: 12
RAM Requested: 64
Memory Available: 1494.8 GB
=======================================================================

~~~~ Assembling genome for Fungus_NoRefNoJuice ~~~~
~~~~ Skipping juicer, no manual curation (not recommended!) ~~~~
~~~~ Starting hifiasm assembly for Fungus_NoRefNoJuice ~~~~
~~~~ Starting Purge_Dups for Fungus_NoRefNoJuice ~~~~
~~~~ Mapping HiC reads to Fungus_NoRefNoJuice draft ~~~~
~~~~ Running HapHiC for Fungus_NoRefNoJuice  ~~~~

❌ Command failed in in /90daydata/coffea_pangenome/puzzler_trials/assemblies/Fungus_NoRefNoJuice/03_haphic/haphic: "${PUZZLER} haphic pipeline ${WD}/${SAMPLE}/02_purge_dups/p_ctg.purged.fa ../filtered.bam ${NUM_CHRS} --correct_nrounds 2 --max_inflation 20.0 --threads ${t} --processes ${t} 2> ../${SAMPLE}.haphic.log" (line 234)

~~~~ HapHiC for Fungus_NoRefNoJuice with 14 chrs failed, trying: 10 ~~~~

❌ Command failed in in /90daydata/coffea_pangenome/puzzler_trials/assemblies/Fungus_NoRefNoJuice/03_haphic/haphic: "${PUZZLER} haphic pipeline ${WD}/${SAMPLE}/02_purge_dups/p_ctg.purged.fa ../filtered.bam ${CHR_ATTEMPT} --correct_nrounds 2 --max_inflation 10.0 --threads ${t} --processes ${t} 2> ../${SAMPLE}.haphic_N_chrs${CHR_ATTEMPT}.log" (line 257)

~~~~ HapHiC for Fungus_NoRefNoJuice with 14 chrs failed, trying: 11 ~~~~

❌ Command failed in in /90daydata/coffea_pangenome/puzzler_trials/assemblies/Fungus_NoRefNoJuice/03_haphic/haphic: "${PUZZLER} haphic pipeline ${WD}/${SAMPLE}/02_purge_dups/p_ctg.purged.fa ../filtered.bam ${CHR_ATTEMPT} --correct_nrounds 2 --max_inflation 10.0 --threads ${t} --processes ${t} 2> ../${SAMPLE}.haphic_N_chrs${CHR_ATTEMPT}.log" (line 257)

~~~~ HapHiC for Fungus_NoRefNoJuice with 14 chrs failed, trying: 12 ~~~~

❌ Command failed in in /90daydata/coffea_pangenome/puzzler_trials/assemblies/Fungus_NoRefNoJuice/03_haphic/haphic: "${PUZZLER} haphic pipeline ${WD}/${SAMPLE}/02_purge_dups/p_ctg.purged.fa ../filtered.bam ${CHR_ATTEMPT} --correct_nrounds 2 --max_inflation 10.0 --threads ${t} --processes ${t} 2> ../${SAMPLE}.haphic_N_chrs${CHR_ATTEMPT}.log" (line 257)

~~~~ HapHiC for Fungus_NoRefNoJuice with 14 chrs failed, trying: 13 ~~~~

❌ Command failed in in /90daydata/coffea_pangenome/puzzler_trials/assemblies/Fungus_NoRefNoJuice/03_haphic/haphic: "${PUZZLER} haphic pipeline ${WD}/${SAMPLE}/02_purge_dups/p_ctg.purged.fa ../filtered.bam ${CHR_ATTEMPT} --correct_nrounds 2 --max_inflation 10.0 --threads ${t} --processes ${t} 2> ../${SAMPLE}.haphic_N_chrs${CHR_ATTEMPT}.log" (line 257)

~~~~ HapHiC for Fungus_NoRefNoJuice with 14 chrs failed, trying: 14 ~~~~

❌ Command failed in in /90daydata/coffea_pangenome/puzzler_trials/assemblies/Fungus_NoRefNoJuice/03_haphic/haphic: "${PUZZLER} haphic pipeline ${WD}/${SAMPLE}/02_purge_dups/p_ctg.purged.fa ../filtered.bam ${CHR_ATTEMPT} --correct_nrounds 2 --max_inflation 10.0 --threads ${t} --processes ${t} 2> ../${SAMPLE}.haphic_N_chrs${CHR_ATTEMPT}.log" (line 257)

~~~~ HapHiC for Fungus_NoRefNoJuice with 14 chrs failed, trying: 15 ~~~~

❌ Command failed in in /90daydata/coffea_pangenome/puzzler_trials/assemblies/Fungus_NoRefNoJuice/03_haphic/haphic: "${PUZZLER} haphic pipeline ${WD}/${SAMPLE}/02_purge_dups/p_ctg.purged.fa ../filtered.bam ${CHR_ATTEMPT} --correct_nrounds 2 --max_inflation 10.0 --threads ${t} --processes ${t} 2> ../${SAMPLE}.haphic_N_chrs${CHR_ATTEMPT}.log" (line 257)

~~~~ HapHiC for Fungus_NoRefNoJuice with 14 chrs failed, trying: 16 ~~~~

❌ Command failed in in /90daydata/coffea_pangenome/puzzler_trials/assemblies/Fungus_NoRefNoJuice/03_haphic/haphic: "${PUZZLER} haphic pipeline ${WD}/${SAMPLE}/02_purge_dups/p_ctg.purged.fa ../filtered.bam ${CHR_ATTEMPT} --correct_nrounds 2 --max_inflation 10.0 --threads ${t} --processes ${t} 2> ../${SAMPLE}.haphic_N_chrs${CHR_ATTEMPT}.log" (line 257)

~~~~ HapHiC for Fungus_NoRefNoJuice with 14 chrs failed, trying: 17 ~~~~

❌ Command failed in in /90daydata/coffea_pangenome/puzzler_trials/assemblies/Fungus_NoRefNoJuice/03_haphic/haphic: "${PUZZLER} haphic pipeline ${WD}/${SAMPLE}/02_purge_dups/p_ctg.purged.fa ../filtered.bam ${CHR_ATTEMPT} --correct_nrounds 2 --max_inflation 10.0 --threads ${t} --processes ${t} 2> ../${SAMPLE}.haphic_N_chrs${CHR_ATTEMPT}.log" (line 257)

~~~~ HapHiC for Fungus_NoRefNoJuice with 14 chrs failed, trying: 18 ~~~~

❌ Command failed in in /90daydata/coffea_pangenome/puzzler_trials/assemblies/Fungus_NoRefNoJuice/03_haphic/haphic: "${PUZZLER} haphic pipeline ${WD}/${SAMPLE}/02_purge_dups/p_ctg.purged.fa ../filtered.bam ${CHR_ATTEMPT} --correct_nrounds 2 --max_inflation 10.0 --threads ${t} --processes ${t} 2> ../${SAMPLE}.haphic_N_chrs${CHR_ATTEMPT}.log" (line 257)

~~~~ HapHiC for Fungus_NoRefNoJuice failed, scaffolding with YAHS instead ~~~~
~~~~ Skipping juicer HiC file creation for Fungus_NoRefNoJuice: --no_juice given, not recommended! ~~~~
~~~~ No reference provided for Fungus_NoRefNoJuice: simply extracting assembly to: /90daydata/coffea_pangenome/puzzler_trials/assemblies/primary_asm/Fungus_NoRefNoJuice.fa ~~~~
~~~~ Skipping juicer extraction for Fungus_NoRefNoJuice, --no_juice requested! (not recommended!) ~~~~
~~~~ No reference provided for Fungus_NoRefNoJuice: no chromosome re-naming, so final assembly already: /90daydata/coffea_pangenome/puzzler_trials/assemblies/primary_asm/Fungus_NoRefNoJuice.fa ~~~~
~~~~ Assessing genome quality for Fungus_NoRefNoJuice ~~~~
~~~~ Creating final HiC bam for Fungus_NoRefNoJuice ~~~~
~~~~ Creating final contact map for Fungus_NoRefNoJuice, no ref specified, plotting scaffolds > 2mb ~~~~
~~~~ Running YAK on Fungus_NoRefNoJuice ~~~~
~~~~ BUSCO lineage dataset already exists, skipping ~~~~
~~~~ Running BUSCO for Fungus_NoRefNoJuice using lineage: fungi_odb10 ~~~~
~~~~ Skipping blobtools for Fungus_NoRefNoJuice, not desired ~~~~
~~~~ Summarizing Assembly for Fungus_NoRefNoJuice ~~~~
~~~~ Your final assembly for Fungus_NoRefNoJuice is: /90daydata/coffea_pangenome/puzzler_trials/assemblies/primary_asm/Fungus_NoRefNoJuice.fa ~~~~
~~~~ Your final assembly stats for Fungus_NoRefNoJuice are in: /90daydata/coffea_pangenome/puzzler_trials/assemblies/primary_asm/stats/Fungus_NoRefNoJuice.summary.txt ~~~~
```

Runtime:

```
Nodes: 1
Cores per node: 12
CPU Utilized: 05:20:07
CPU Efficiency: 51.35% of 10:23:24 core-walltime
Job Wall-clock time: 00:51:57
Memory Utilized: 54.41 GB
Memory Efficiency: 85.01% of 64.00 GB
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
PUZZLER="apptainer exec /home/justin.merondun/apptainer/puzzler_v1.7.sif"

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