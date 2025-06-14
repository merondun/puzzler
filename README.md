![Puzzler](/examples/figs/logo.png)

Simple check-point aware shell script wrapper for high-throughput genome assembly.

> What's it for?!

Scalable genome assembly. As simple and portable as possible: only requires the `.sif` and a table with file paths. Optimized for container-capable SLURM resources. 

<!-- TOC start (generated with https://github.com/derlin/bitdowntoc) -->

- [Installation](#installation)
- [Workflow](#workflow)
- [Quick Start](#quick-start)
- [Details](#details)
- [Visual Workflow](#visual)
- [Outputs](#outputs)
- [FAQ](#faq)
- [Contact](#contact)
- [Citation](#citation)
- [Changelog](#changelog)

<!-- TOC end --> 

<!-- TOC --><a name="installation"></a>
## Installation

Currently, only a containerized release is supported because we shouldn't waste our time with software install and [apptainer](https://apptainer.org/docs/admin/main/installation.html) is obtainable without root. The singularity / apptainer `.sif` can be yanked with: `apptainer pull --arch amd64 library://merondun/default/puzzler:latest` 

The workhorse script `puzzler` is simply a check-point aware bash script. It can be submitted directly with slurm, e.g. `sbatch puzzler --sample Fungus --map samples.tsv`. 

For installation:  

1) Download the container `.sif`
2) Download the `puzzler` shell script, ideally add to path (e.g. `chmod +x /path/puzzler; echo 'export PATH=$PATH:/path/puzzler' >> ~/.bashrc`)
3) Modify your SLURM/HPC settings within the `puzzler` script, if necessary. 

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

module load apptainer
########### EDIT THIS BLOCK WITH APPTAINER/SINGULARITY/MODULES ############
:wq
```

Test: 

```
puzzler -h
Usage: puzzler [OPTIONS]

Options:
  --sample SAMPLE       Sample name (required)
  --map FILE            Path to .tsv/.csv map file (required)
  --threads t           Number of threads (optional; default 64)
  --mem MEM             Memory allocation (optional; default 512)
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

**Conda install**

1b) As a last resort, create conda environment from `/apptainer/environment.yml`, and install [HapHiC](https://github.com/zengxiaofei/HapHiC) and [assembly_stats](https://github.com/MikeTrizna/assembly_stats). If you go this route, please troubleshoot software and inspect the `.logs` before posting issues. 

:exclamation: If you go the conda route, you must make sure that `juicer pre` is available on path (included in HapHic installation), and you will need to modify path to `java -Xmx${MEM}G -jar /opt/HapHiC/utils/juicer_tools.1.9.9_jcuda.0.8.jar` to the correct path within `puzzler`. Simply search for that line in the script and replace the `/opt/` path with the path to your `.jarfile`. 


<!-- TOC --><a name="workflow"></a>
## Workflow

The pipeline creates a collapsed completely *de novo* primary genome assembly using both HiFi and HiC data. The workflow is: 

:pushpin: **`puzzler`**

1) [Hifiasm](https://github.com/chhylp123/hifiasm) assembly using HiFi + HiC reads.
2) Purging of haplotigs using [purge_dups](https://github.com/dfguan/purge_dups). Typically with sufficient HiFi data for diploid organisms, `hifiasm` adequately purges most coverage-based diplotigs. This `puzzler` implementation is not coverage based - only sequence similarity based, which works well in many species with variable ploidy. 
3) Scaffolding with [HapHiC](https://github.com/zengxiaofei/HapHiC). 
4) Creating [Juicebox](https://github.com/aidenlab/Juicebox) manual curation inputs (`.hic`, `.assembly`). 

:mag: **Manual Curation in Juicebox** 

:pushpin: **`puzzler`**

5) Finalizing assembly, assign chromosome names based on a related reference using [merothon](https://github.com/merondun/merothon). 
6) Re-map HiC for final assembly check with [HapHiC](https://github.com/zengxiaofei/HapHiC). 

You may then wish to proceed with `puzzle_quality`, which runs a variety of quality checks on assembly fasta files. 

<!-- TOC --><a name="quick-start"></a>
## Quick Start

`puzzler` only requires two inputs: `--sample`, `--map`. 

**1. Make pipeline map file** 

:boom: *Most important!*

Prepare a `samples.tsv` (either tab or comma separated) which outlines all necessary pipeline components. An example template is found in [/examples/samples.tsv](/examples/samples.tsv)

```
sample	container	wd	ploidy	chromosomes	hifi	hic_r1	hic_r2	Reference	hom_cov
Crane	/home/justin.merondun/apptainer/puzzler_v1.5.sif	/90daydata/coffea_pangenome/puzzler_trials/assemblies	2	40	/90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Crane.HiFi.fastq.gz	/90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Crane.HiC.R1.fastq.gz	/90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Crane.HiC.R2.fastq.gz	/90daydata/coffea_pangenome/puzzler_trials/raw_data/references/GCA_964106855.1_bGruGru1.hap1.1_genomic.fna	48
Fish	/home/justin.merondun/apptainer/puzzler_v1.5.sif	/90daydata/coffea_pangenome/puzzler_trials/assemblies	2	24	/90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Fish.HiFi.fastq.gz	/90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Fish.HiC.R1.fastq.gz	/90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Fish.HiC.R2.fastq.gz	/90daydata/coffea_pangenome/puzzler_trials/raw_data/references/GCA_951216825.1_fEleAnt2.1_genomic.fna	42
```

:page_with_curl: Map file descriptions (:exclamation: Use full paths!!!)

* **sample:** Sample ID, all assembly work will be saved in `$WD/$SAMPLE``
* **container:** Path to the apptainer `.sif`. If using conda, simply write 'conda'. 
* **wd:** Path to working directory to store all files 
* **ploidy:** Ploidy of organism, or best guess. 
* **chromosomes:** Number of chromosomes, or best guess. The pipeline will attempt +/- 4 your estimate in case you don't know.
* **hifi:** Path to HiFi reads 
* **hic_r1:** Path to HiC R1 
* **hic_r2:** Path to HiC R2
* **reference:** Path to related species genome for chromosome naming. Scaffolds will be renamed to the closest syntenic chromosome **using their scaffold naming convention**. 
* **hom_cov:** Homozygous peak coverage (OPTIONAL; see below, otherwise write "NA") 

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


:one: **Step 1:** `[sbatch] puzzler --sample Fungus --map samples.tsv`

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

:two: **Step 2:** resubmit `[sbatch] puzzler --sample Fungus --map samples.tsv`

This script will check for these files, and create them if they do not exist in these directories. The script will not start if you haven't added the Juicebox `.review` file or if you do not have a related reference genome. 

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
SIF_PATH=/home/justin.merondun/apptainer/puzzler_v1.6.sif
PUZZLER="apptainer exec ${SIF_PATH}"

cd ${WD}/${SAMPLE}/05_postjuicebox
${PUZZLER} juicer post -o final_asm.fa \
  ${WD}/primary_asm/juicer_files/${SAMPLE}_JBAT.review.assembly \
  ${WD}/${SAMPLE}/04_juicer/haphic-refsort_JBAT.liftover.agp \
  ${WD}/${SAMPLE}/02_purge_dups/p_ctg.purged.fa
```

If you still want the final HiC map, you can then re-submit the script `puzzler --sample Fungus --map samples.tsv`, and the script will simply run the final HiC mapping check, outputting the juicebox-curated assembly in `$WD/primary_asm/$SAMPLE.fa` and the final contact map in `WD/primary_asm/stats/$SAMPLE.pdf`. 

<br/>

<!-- TOC --><a name="visual"></a>
## Visual Workflow

Below is an exhaustive workflow outline documenting the inputs and outputs for each `puzzler` and `puzzle_quality` step. Only four hands-on steps are required from raw reads to assembly summary statistics. 

![workflow](/examples/figs/Workflow.png)

<!-- TOC --><a name="faq"></a>
## FAQ

> Do I need to specify the Hi-C enzyme motifs for HapHiC?

[HapHiC](https://github.com/zengxiaofei/HapHiC) is very good at recognizing fragments without specification, I haven't needed to do this - although I have mainly tested with Arima kits. 

> The script fails on the HapHic step. 

This can be difficult to diagnose, so please inspect the `${WD}/${SAMPLE}/03_haphic/${SAMPLE}.haphic.log`, and refer to [HapHiC](https://github.com/zengxiaofei/HapHiC) about e.g. "chromosomes are grouped together". Typically I see that error when there is insufficient or low quality HiC data. 

> I can't launch the container, "FATAL:   container creation failed: failed to resolve session directory /var/apptainer/mnt/session: lstat /var/apptainer: no such file or directory"

I see this warning if you try to run apptainer on a login node, which is typically a setting set by IT. You can try modifying some paths IF NECESSARY, but I really recommend testing on a compute node or submitting to a job scheduler. 

Potentially try:

```
export APPTAINER_TMPDIR=/your/tmp/dir
export APPTAINER_CACHEDIR=$HOME/.apptainer_cache
mkdir -p "$APPTAINER_TMPDIR"
apptainer exe $SIF hifiasm
```

> I encountered the warning: Multiple scaffolds corresponding to a single Chr for ${SAMPLE} INSPECT!

This means that there are multiple scaffolds/chromosomes in your draft which correspond to a single chromosome in the reference genome. The script will automatically rename the duplicates into e.g. `chr1` `chr1A` `chr1B` according to length. 

You should probably inspect `${WD}/${SAMPLE}/05_postjuicebox/chromosome_naming_map.txt` and to ensure that you are happy with this naming scheme.

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

If you did any manual edits and want to re-generate `final_asm.fa`, just fill in these fields and re-run puzzler.


> The script gave a warning: HapHiC for Fungus with 14 chrs failed, trying: 10 

If HapHic does not succeed with the specified chromosome numbers, it will re-run it with +/- 4 chromsomes. If you specified n = 14 chromosomes, `puzzler` will re-run HapHic with 10 - 18 chromoomes. If all of those fail, it will run scaffolding with YAHS instead. 

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

Please make a git issue with any problems, no matter how small. Otherwise, email Justin at heritabilities [@] gmail.com 

<!-- TOC --><a name="citation"></a>
## Citation

Please ensure you cite the developers of software within `puzzler`:

* [apptainer/singularity](https://github.com/apptainer/apptainer): https://doi.org/10.1371/journal.pone.0177459
* [hifiasm v0.25.0-r726](https://github.com/chhylp123/hifiasm): https://doi.org/10.1038/s41592-020-01056-5
* [purge_dups v1.2.6](https://github.com/dfguan/purge_dups): https://doi.org/10.1093/bioinformatics/btaa025
* [minimap2 v2.29-r1283](https://github.com/lh3/minimap2): https://doi.org/10.1093/bioinformatics/bty191
* [bwa-mem2 v2.2.1](https://github.com/bwa-mem2/bwa-mem2): https://doi.org/10.1109/IPDPS.2019.00041
* [samblaster v0.1.26](https://github.com/GregoryFaust/samblaster): https://doi.org/10.1093/bioinformatics/btu314
* [samtools v1.21](https://github.com/samtools/samtools): https://doi.org/10.1093/gigascience/giab008
* [HapHiC v1.0.6](https://github.com/zengxiaofei/HapHiC): https://doi.org/10.1038/s41477-024-01755-3 and https://doi.org/10.1038/s41477-019-0487-8
* [juicer v1.2](https://github.com/aidenlab/juicer): https://doi.org/10.1016/j.cels.2016.07.002

If using chromosome-renaming: 
* [merothon v0.4.2](https://github.com/merondun/merothon)
* [seqkit v2.10.0](https://bioinf.shenwei.me/seqkit/): https://doi.org/10.1002/imt2.191

If using optional software:
* [genomescope2 v2.0](https://github.com/tbenavi1/genomescope2.0): https://doi.org/10.1038/s41467-020-14998-3

If using, please ensure you cite the developers of software within `puzzle_quality`:
* [busco v5.8.3](https://busco.ezlab.org/): https://doi.org/10.1093/bioinformatics/btv351
* [yak v0.1-r69-dirty](https://github.com/lh3/yak): https://doi.org/10.1038/s41592-020-01056-5 (same as hifiasm)
* [blobtools v1.1.1](https://blobtools.readme.io/docs/what-is-blobtools): https://f1000research.com/articles/6-1287/v1
* [assembly_stats v0.1.2](https://github.com/MikeTrizna/assembly_stats): https://zenodo.org/record/3968774

<!-- TOC --><a name="changelog"></a>
## Changelog

`puzzler`

**v1.7**: merge puzzler & puzzle_quality. 
**v1.6**: add accomodation for specifying no reference ("NA" in `samples.tsv`). 
