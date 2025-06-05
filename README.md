![Puzzler](/examples/figs/logo.png)

Simple check-point aware shell script for high-throughput genome assembly.

> What's it for?!

Making many genomes, with a focus on pangenomics. Primarily designed for container-capable SLURM resources. 

<!-- TOC start (generated with https://github.com/derlin/bitdowntoc) -->

- [Installation](#installation)
- [Workflow](#workflow)
- [Quick Start](#quick-start)
- [Details](#details)
- [Outputs](#outputs)
- [FAQ](#faq)
- [Contact](#contact)

<!-- TOC end --> 

<!-- TOC --><a name="installation"></a>
## Installation

Currently, only a containerized release is supported because we shouldn't waste our time with software install. The singularity / apptainer `.sif` can be yanked with: `apptainer pull --arch amd64 library://merondun/default/puzzler:latest` 

The workhorse script `puzzler` is simply a check-point aware bash script. It can be submitted directly with slurm, e.g. `sbatch puzzler --sample Fungus --map samples.tsv`. 

For installation:  

1) Download the container `.sif`
2) Download the `puzzler` shell script
3) Modify your SLURM and apptainer/singularity loading within that `puzzler` script: 

```
wget https://raw.githubusercontent.com/merondun/puzzler/main/bin/puzzler
vim puzzler

########### EDIT THIS BLOCK WITH APPTAINER/SINGULARITY/MODULES ############
#SBATCH --time=10-00:00:00   
#SBATCH --nodes=1  
#SBATCH --ntasks-per-node=64
#SBATCH --mem=512Gb
#SBATCH --partition=ceres
#SBATCH --account=coffea_pangenome

RUNTIME="apptainer"
#RUNTIME="singularity"

module load apptainer
SINGULARITY_TMPDIR=$APPTAINER_TMPDIR
########### EDIT THIS BLOCK WITH APPTAINER/SINGULARITY/MODULES ############
:wq
```

Test: 

```
puzzler -h
Usage: puzzler [OPTIONS]

Options:
  --sample SAMPLE       Sample name (required)
  --map_file FILE       Path to .csv map file (required)
  --threads t           Number of threads (optional)
  --mem MEM             Memory allocation (optional)
  --help                Show this help message and exit
```

The (large, ~2.5Gb) container contains many common tools including:

* [hifiasm v0.25.0-r726](https://github.com/chhylp123/hifiasm)
* [purge_dups v1.2.6](https://github.com/dfguan/purge_dups)
* [minimap2 v2.29-r1283](https://github.com/lh3/minimap2)
* [bwa-mem2 v2.2.1](https://github.com/bwa-mem2/bwa-mem2)
* [samblaster v0.1.26](https://github.com/GregoryFaust/samblaster)
* [samtools v1.21](https://github.com/samtools/samtools)
* [HapHiC v1.0.6](https://github.com/zengxiaofei/HapHiC)
* [juicer v1.2](https://github.com/aidenlab/juicer)
* [merothon v0.4.2](https://github.com/merondun/merothon)
* [seqkit v2.10.0](https://bioinf.shenwei.me/seqkit/)
* [seqtk v1.4-r122](https://github.com/lh3/seqtk)
* [assembly_stats v0.1.2](https://github.com/MikeTrizna/assembly_stats)

1b) As a last resort, you can create a conda environment from the `environment.yml` file in `/apptainer/`, and install [HapHiC](https://github.com/zengxiaofei/HapHiC) and [assembly_stats](https://github.com/MikeTrizna/assembly_stats). J

<!-- TOC --><a name="workflow"></a>
## Workflow

The pipeline creates a primary genome assembly using both HiFi and HiC data. The workflow is: 

:pushpin: **`puzzler`**

1) [Hifiasm](https://github.com/chhylp123/hifiasm) assembly using HiFi + HiC reads.
2) Purging of haplotigs using [purge_dups](https://github.com/dfguan/purge_dups). Typically with sufficient HiFi data for diploid organisms, `hifiasm` adequately purges most coverage-based diplotigs. This `puzzler` implementation is not coverage based - only sequence similarity based, which works well in many species with variable ploidy. 
3) Scaffolding with [HapHiC](https://github.com/zengxiaofei/HapHiC). 
4) Creating [Juicebox](https://github.com/aidenlab/Juicebox) manual curation inputs (`.hic`, `.assembly`). 

:mag: **Manual Curation in Juicebox** 

:pushpin: **`puzzler`**

5) Finalizing assembly, assign chromosome names based on a related reference. 
6) Re-map HiC for final assembly check. 
7) Output assembly statistics 

<!-- TOC --><a name="quick-start"></a>
## Quick Start

**1. Make pipeline map file** 

:boom: *Most important!*

Prepare a `samples.tsv` which outlines all necessary pipeline components. An example template is found in [/examples/samples.tsv](/examples/samples.tsv)

```
sample	container	wd	ploidy	chromosomes	hifi	hic_r1	hic_r2	Reference	hom_cov
Crane	/home/justin.merondun/apptainer/puzzler_v1.5.sif	/90daydata/coffea_pangenome/puzzler_trials/assemblies	2	40	/90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Crane.HiFi.fastq.gz	/90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Crane.HiC.R1.fastq.gz	/90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Crane.HiC.R2.fastq.gz	/90daydata/coffea_pangenome/puzzler_trials/raw_data/references/GCA_964106855.1_bGruGru1.hap1.1_genomic.fna	48
Fish	/home/justin.merondun/apptainer/puzzler_v1.5.sif	/90daydata/coffea_pangenome/puzzler_trials/assemblies	2	24	/90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Fish.HiFi.fastq.gz	/90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Fish.HiC.R1.fastq.gz	/90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Fish.HiC.R2.fastq.gz	/90daydata/coffea_pangenome/puzzler_trials/raw_data/references/GCA_951216825.1_fEleAnt2.1_genomic.fna	42
```

:page_with_curl: Map file descriptions (:exclamation: Use full paths!!!)

* **sample:** Sample ID, all assembly work will be saved in `$WD/$SAMPLE``
* **container:** Path to the apptainer `.sif`
* **wd:** Path to working directory to store all files 
* **ploidy:** Ploidy of organism
* **chromosomes:** Number of chromosomes
* **hifi:** Path to HiFi reads 
* **hic_r1:** Path to HiC R1 
* **hic_r2:** Path to HiC R2
* **reference:** Path to related species genome for chromosome naming. Scaffolds will be renamed to the closest syntenic chromosome using their naming convention. 
* **hom_cov:** Homozygous peak coverage (OPTIONAL; see below) 

***Homozygous peak coverage*** is the homozygous peak covewrage identified from k-mer coverage in the HiFi library. I prefer to quickly run genomescope2, where the `*_linear_plot.png` indicates the left peak with `kcov:`, which you can multiple by ploidy to get `--hom_cov`. 

If you do not run `genomescope2`, I **highly** recommend you inspect the hifiasm log file (in $WD/$SAMPLE/$SAMPLE.hifiasm.log).

For more details about this homozygous coverage and for a full example, see the [Fungus example](/docs/Fungus_Example.md). 


<!-- TOC --><a name="details"></a>
## Details

:fire: `puzzler` is checkpoint-based, so it will skip already-completed tasks based on file existence (files with a size > 0). You can therefore trick the script by creating or copying manual files into the appropriate directories. 

For instance, if you submit the command `puzzler --sample Fungus --map samples.tsv` and you had already completed hifiasm, purge duplicates, and haphic steps, yet miss the juicer outputs, the script will create the juicebox `.hic` files. 

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
Parameters for sample: Fungus
CONTAINER: /home/justin.merondun/apptainer/puzzler_v1.5.sif 
WD: /90daydata/coffea_pangenome/puzzler_trials/assemblies 
PLOIDY: 1 
NUMBER CHRS: 3
HIFI: /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Fungus.HiFi.fastq.gz
HIC_R1: /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Fungus.HiC.R1.fastq.gz
HIC_R2: /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Fungus.HiC.R2.fastq.gz
REFERENCE: /90daydata/coffea_pangenome/puzzler_trials/raw_data/references/GCA_964106605.1_gdRhiKalk1.hap1.1_genomic.fna
HOM_COV: 
RUNTIME: apptainer
=======================================================================

~~~~ Starting hifiasm assembly for Fungus ~~~~
~~~~ Starting Purge_Dups for Fungus ~~~~
~~~~ Mapping HiC reads to Fungus draft ~~~~
```


:one: **Step 1:** `sbatch puzzler --sample Fungus --map samples.tsv`

This script will check for these files, and create them if they do not exist in these directories:

| File                                               | Description                  | If missing, run |
| -------------------------------------------------- | ---------------------------- | --------------- |
| $WD/Fungus/01_hifiasm/asm.hic.p_ctg.gfa           | Primary contig assembly      | hifiasm         |
| $WD/Fungus/02_purge_dups/p_ctg.purged.fa          | Purged contig assembly       | purge_dups      |
| $WD/Fungus/03_haphic/filtered.bam                 | Filtered re-mapped HiC reads | bwa mem         |
| $WD/Fungus/03_haphic/haphic/04.build/scaffolds.fa | HapHiC scaffolded assembly   | haphic pipeline |
| $WD/primary_asm/juicer_files/Fungus_JBAT.hic      | Juicebox input file          | juicer pre      |

<br/>

:mag: **Manual curation:** Juicebox! 



Open Juicebox, and drag the `.hic` file into the window. Import the `.assembly` file using `Assembly > Import Map Assembly`. Make any adjustments if necessary (only about 50% of my genomes need it), and then export the file with `Assembly > Export Assembly`. 

This will create e.g. `${SAMPLE}_JBAT.review.assembly`. Maintain that file name, and copy it to `${WD}/primary_asm/juicer_files/${SAMPLE}_JBAT.review.assembly`. 

<br/>

:one: **Step:** resubmit `sbatch puzzler --sample Fungus --map samples.tsv`

This script will check for these files, and create them if they do not exist in these directories. The script will not start if you haven't added the Juicebox `.review` file. 

| File                                                      |                            | If missing, run |
| --------------------------------------------------------- | -------------------------- | --------------- |
| $WD/primary_asm/juicer_files/Fungus_JBAT.review.assembly | Juicebox output file       | NOTHING!        |
| $WD/Fungus/05_postjuicebox/map.txt                       | Scaffold renaming file map | minimap2        |
| $WD/Fungus/05_postjuicebox/final_asm.fa                  | Renamed scaffold file      | seqkit renaming |
| $WD/Fungus/06_realign_hic_hifi/final_asm.filtered.bam    | Final HiC re-mapped file   | bwa mem         |
| $WD/primary_asm/stats/Fungus.pdf                         | Final contact map          | haphic plot     |
| $WD/primary_asm/stats/Fungus.stats.txt                   | Calculate Assembly Stats   | assembly_stats  |

<br/>

:volcano: **NOTE:** If you don't care about renaming chromosomes to a reference, you can simply run juicer:

```
WD=/90daydata/coffea_pangenome/puzzler_trials/assemblies
SAMPLE=Fungus
SIF_PATH=/home/justin.merondun/apptainer/puzzler_v1.4.sif
PUZZLER="apptainer exec ${SIF_PATH}"

cd ${WD}/${SAMPLE}/05_postjuicebox
${PUZZLER} juicer post -o final_asm.fa \
  ${WD}/primary_asm/juicer_files/${SAMPLE}_JBAT.review.assembly \
  ${WD}/${SAMPLE}/04_juicer/haphic-refsort_JBAT.liftover.agp \
  ${WD}/${SAMPLE}/02_purge_dups/p_ctg.purged.fa
```

You can then re-submit the script `puzzler --sample Fungus --map samples.tsv`, and the script will simply run the final HiC mapping check, outputting the juicebox-curated assembly in `$WD/primary_asm/$SAMPLE.fa` and the final contact map in `WD/primary_asm/stats/$SAMPLE.pdf`. 

<br/>

:bomb: :warning: **If you then encounter this warning during post-curation** `puzzler`: 

`~~~~ Multiple scaffolds corresponding to a single Chr for ${SAMPLE} INSPECT!  ~~~~` 

This means that there are multiple scaffolds/chromosomes in your draft which correspond to a single chromosome in the reference genome. The script will automatically rename the duplicates into e.g. `chr1` `chr1A` `chr1B` according to length. 

You should probably inspect `${WD}/${SAMPLE}/05_postjuicebox/chromosome_naming_map.txt` and to ensure that you are happy with this naming s cheme.

You can also manually modify any of these by changing the value in the second column.

For example, in this file:

```
head ${WD}/${SAMPLE}/05_postjuicebox/chromosome_naming_map.txt
scaffold_1      chr1    +       42.46%  35386063
scaffold_2      chr1  +       41.33%  9847223
scaffold_3      chr11   +       47.76%  28452307
scaffold_4      chr2    +       39.75%  31132341
```

There are two scaffolds with a high similarity (~42%) and a large amount of aligned bases (35 Mb and 9.8 Mb) to reference chr1. 

The script will automatically rename these chr1 and chr1A, but you could also modify this manually and resubmit:

```
head ${WD}/${SAMPLE}/05_postjuicebox/chromosome_naming_map.txt
scaffold_1      chr1    +       42.46%  35386063
scaffold_2      chr1A  +       41.33%  9847223
scaffold_3      chr11   +       47.76%  28452307
scaffold_4      chr2    +       39.75%  31132341
```

If you did any manual edits and want to re-generate `final_asm.fa`, just fill in these fields and re-run: 

```
WD=/90daydata/coffea_pangenome/puzzler_trials/assemblies
SAMPLE=Fungus
SIF_PATH=/home/justin.merondun/apptainer/puzzler_v1.4.sif
PUZZLER="apptainer exec ${SIF_PATH}"

cd ${WD}/${SAMPLE}/05_postjuicebox/

# Sort them according to ID
awk '$2 ~ /^chr/ {print $0}' chromosome_naming_map.txt | sort -k2,2V | awk '{print $2}' > sorted_chr.txt
awk '$2 ~ /^scaffold_/ {print $0}' chromosome_naming_map.txt | sort -k2,2V | awk '{print $2}' >> sorted_chr.txt

# Rename
${PUZZLER} seqkit replace --line-width 0 -p "(.*)" -r "{kv}" -k chromosome_naming_map.txt orient.fa > haphic_renamed_unord.fa 2> seqkit_renaming.log
${PUZZLER} samtools faidx haphic_renamed_unord.fa
xargs ${PUZZLER} samtools faidx haphic_renamed_unord.fa < sorted_chr.txt > final_asm.fa
```

As above, you can re-submit now to generate the contact maps: `puzzler --sample Fungus --map samples.tsv`.

<!-- TOC --><a name="faq"></a>
## FAQ

> Do I need to specify the Hi-C enzyme motifs for HapHiC?

[HapHiC](https://github.com/zengxiaofei/HapHiC) is very good at recognizing fragments without specification, I haven't needed to do this - although I have mainly tested with Arima kits. 

> The script fails on the HapHic step. 

This can be difficult to diagnose, so please inspect the `${WD}/${SAMPLE}/03_haphic/${SAMPLE}.haphic.log`, and refer to [HapHiC](https://github.com/zengxiaofei/HapHiC) about e.g. "chromosomes are grouped together". Typically I see that error when there is insufficient or low quality HiC data. 

<!-- TOC --><a name="outputs"></a>
## Outputs 

Relative to the `$WD` path, the outputs will be: 

* The final assembly will be found in: `${WD}/primary_asm/${SAMPLE}.fa`.

* The final assembly stats will be in: `${WD}/primary_asm/stats/${SAMPLE}.summary.txt`, e.g.: 

* The final contact map for assessment will be in: `${WD}/primary_asm/stats/${SAMPLE}.pdf` and should look something like this: 

![Final_Map](/examples/figs/contact_example.png)


Each directory within `$WD/$SAMPLE` corresponds to an assembly step, where all logs for all steps are saved in `.log` files. Please inspect these before opening an issue to see what the problem might be. 


<!-- TOC --><a name="contact"></a>
## Contact

Either make an issue or send an email to Justin at heritabilities [@] gmail.com 

Please ensure you cite the developers of the main  `puzzler`:

* [hifiasm v0.25.0-r726](https://github.com/chhylp123/hifiasm)
* [purge_dups v1.2.6](https://github.com/dfguan/purge_dups)
* [minimap2 v2.29-r1283](https://github.com/lh3/minimap2)
* [bwa-mem2 v2.2.1](https://github.com/bwa-mem2/bwa-mem2)
* [samblaster v0.1.26](https://github.com/GregoryFaust/samblaster)
* [samtools v1.21](https://github.com/samtools/samtools)
* [HapHiC v1.0.6](https://github.com/zengxiaofei/HapHiC)
* [juicer v1.2](https://github.com/aidenlab/juicer)
* [merothon v0.4.2](https://github.com/merondun/merothon)
* [seqkit v2.10.0](https://bioinf.shenwei.me/seqkit/)
* [seqtk v1.4-r122](https://github.com/lh3/seqtk)
* [assembly_stats v0.1.2](https://github.com/MikeTrizna/assembly_stats)
* [genomescope2 v2.0](https://github.com/tbenavi1/genomescope2.0)