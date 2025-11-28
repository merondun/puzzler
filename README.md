![Puzzler](/examples/figs/logo.png)

Check-point aware containerized shell pipeline for high-throughput genome assembly.

> What's it for?!

Scalable genome assembly. As simple and portable as possible: only requires conda install or `.sif` and a table with file paths. 

Puzzler v1.9 [![DOI](https://zenodo.org/badge/891638219.svg)](https://doi.org/10.5281/zenodo.15733730) [![Anaconda-Server Badge](https://anaconda.org/heritabilities/puzzler/badges/version.svg)](https://anaconda.org/heritabilities/puzzler)


<!-- TOC start (generated with https://github.com/derlin/bitdowntoc) -->

- [Installation](#installation)
- [Workflow](#workflow)
- [Quick Start](#quick-start)
- [Details](#details)
- [Visual Workflow](#visual)
- [Detailed Steps](#steps)
- [Outputs](#outputs)
- [FAQ](#faq)
- [Contact](#contact)
- [Citation](#citation)
- [Changelog](#changelog)

<!-- TOC end --> 

<!-- TOC --><a name="installation"></a>
## Installation

Install with either [conda](https://anaconda.org/heritabilities/puzzler) or [apptainer](https://apptainer.org/docs/admin/main/installation.html). 

The workhorse script `puzzler` is a check-point aware bash script: `puzzler --sample Fungus --map samples.tsv`. 

installation:  

**With Conda** :snake: 

This includes all dependencies and the `puzzler` binary. 

```bash
mamba create -n puzzler -c bioconda -c heritabilities -c conda-forge -c hcc puzzler -y
# add --channel-priority flexible if needed
mamba activate puzzler
```

:wrench: Note: this conda build depends on internal packaging of [HapHiC](https://github.com/zengxiaofei/HapHiC) v1.0.6 and bundled juicer/juicertools files. Please report issues with these steps here first so I can determine if they stem from my conda packaging or the external software.

**With apptainer** (safest). 

All dependencies are within the `.sif`. 

1) Download the container `.sif`:

```
apptainer pull --arch amd64 library://merondun/default/puzzler:v1.9

# Depending on your architecture, you might need to update your apptainer libraries:
apptainer remote add --no-login SylabsCloud cloud.sycloud.io
apptainer remote use SylabsCloud
apptainer pull --arch amd64 library://merondun/default/puzzler:v1.9
```
Note that Puzzler v1.8 and v1.9 have the same required dependencies. 

2) Download the `puzzler` shell script and add to path. `puzzler` automatically binds the directories listed in `samples.tsv`` and runs Apptainer. No manual commands required!

```
wget https://raw.githubusercontent.com/merondun/puzzler/main/bin/puzzler
chmod +x puzzler
INSTALL_PATH=$(dirname "$(realpath puzzler)")
export PATH="$PATH:$INSTALL_PATH"
```

Afterwards: 

Test: 

```bash
puzzler -h
Usage: puzzler -s sample -m samples.tsv [OPTIONS]

Options:
  -s, --sample SAMPLE   Sample name (required)
  -m, --map FILE        Path to .tsv/.csv map file (required)
  --threads t           Number of threads (optional; default 16)
  --mem MEM             Memory allocation (optional; default 128)
  --no_purge            Skip purge duplicates step entirely (optional)
  --no_juice            Skip juicer file creation entirely (optional; not recommended!)
  -v, --version         Show version and exit
  -h, --help            Show help and exit

  Required --map Structure:
  The provided map file (e.g., samples.txt) must contain the following columns in this order:
  RUNTIME CONTAINER WD HIFI HIC_R1 HIC_R2 NUM_CHRS REFERENCE HOM_COV BLOB_DB BUSCO_LINEAGE BUSCO_DB
  For optional columns (REFERENCE - BUSCO_DB), write NA if undesired.
```

<!-- TOC --><a name="workflow"></a>
## Workflow

The pipeline creates a collapsed completely *de novo* primary genome assembly using **both HiFi and HiC data**. The workflow is: 

:pushpin: **`puzzler -s $sample -m samples.tsv`**

1) [Hifiasm](https://github.com/chhylp123/hifiasm) assembly using HiFi + HiC reads.
2) Purging of haplotigs using [purge_dups](https://github.com/dfguan/purge_dups). Typically with sufficient HiFi data for diploid organisms, `hifiasm` adequately purges most coverage-based diplotigs. **This `puzzler` implementation is not coverage based** - only sequence similarity based, which works well in many species (even polyploids) and will not require re-mapping HiFi reads. This step can be omitted by passing the argument `--no_purge`. 
3) Scaffolding with [HapHiC](https://github.com/zengxiaofei/HapHiC), or [YAHS](https://github.com/c-zhou/yahs) if HapHic fails. 
4) Creating [Juicebox](https://github.com/aidenlab/Juicebox) manual curation inputs (`.hic`, `.assembly`), or skipping if `--no_juice` is given. 

:mag: **Manual Curation in Juicebox** [(see brief curation for 11 genomes here)](https://www.youtube.com/watch?v=rMUiNqZwEpA)

:pushpin: rerun: **`puzzler -s $sample -m samples.tsv`**

5) Finalizing assembly: extract juicer assembly, assign chromosome names based on a related reference (only scaffolds with >Chr or >chr!), reverse complement according to reverse, using [merothon](https://github.com/merondun/merothon).
6) QC: Re-map HiC for final assembly check with [HapHiC](https://github.com/zengxiaofei/HapHiC). 
7) QC: [yak](https://github.com/lh3/yak) base quality check. 
8) QC: [busco](https://busco.ezlab.org/) check. (optional; if database path specified)
9) QC: [blobtools](https://blobtools.readme.io/docs/what-is-blobtools) contaminant check against uniprot. (optional; if database path specified)
10) QC: [assembly_stats](https://github.com/MikeTrizna/assembly_stats) assembly statistic metrics. 


<!-- TOC --><a name="quick-start"></a>
## Quick Start

`puzzler` only requires two inputs: `--sample` and `--map`. 

Please see a small [24 Mb fungus example](/docs/Fungus_Example.md) and an example tutorial with [N = 12 species](/docs/All_N12_Genomes.md), with a small haploid fungal trial dataset (runs in 40 minutes with 8 cores and 64 Gb RAM) available on [Zenodo](https://doi.org/10.5281/zenodo.17756307). 

**1. Make pipeline map file** 

:boom: *Most important!*

Prepare a `samples.tsv` which outlines all necessary pipeline components. An example template is found in [/examples/samples.tsv](/examples/samples.tsv)

Below is an excerpt showing both conda and apptainer runtimes. Columns `Reference` and later are optional. 

| sample | runtime   | container                                        | wd                                                    | hifi                                                         | hic_r1                                                       | hic_r2                                                       | num_chrs | reference                                                    | hom_cov | blob_database                                             | busco_lineage  | busco_database                                             |
| ------ | --------- | ------------------------------------------------ | ----------------------------------------------------- | ------------------------------------------------------------ | ------------------------------------------------------------ | ------------------------------------------------------------ | -------- | ------------------------------------------------------------ | ------- | --------------------------------------------------------- | -------------- | ---------------------------------------------------------- |
| Fly    | conda     | NA                                               | /90daydata/coffea_pangenome/puzzler_trials/assemblies | /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Fly.HiFi.fastq.gz | /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Fly.HiC.R1.fastq.gz | /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Fly.HiC.R2.fastq.gz | 7        | /90daydata/coffea_pangenome/puzzler_trials/raw_data/references/GCF_028408465.1_idAnaLude1.1_genomic.fna | 44      | /90daydata/coffea_pangenome/puzzler_trials/blob_downloads | diptera_odb10  | /90daydata/coffea_pangenome/puzzler_trials/busco_downloads |
| Beaver | apptainer | /home/justin.merondun/apptainer/puzzler_v1.8.sif | /90daydata/coffea_pangenome/puzzler_trials/assemblies | /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Beaver.HiFi.fastq.gz | /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Beaver.HiC.R1.fastq.gz | /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Beaver.HiC.R2.fastq.gz | 20       | /90daydata/coffea_pangenome/puzzler_trials/raw_data/references/GCF_047511655.1_mCasCan1.hap1v2_genomic.fna | NA      | NA                                                        | mammalia_odb10 | /90daydata/coffea_pangenome/puzzler_trials/busco_downloads |
| Fungus | apptainer | /home/justin.merondun/apptainer/puzzler_v1.8.sif | /90daydata/coffea_pangenome/puzzler_trials/assemblies | /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Fungus.HiFi.fastq.gz | /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Fungus.HiC.R1.fastq.gz | /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Fungus.HiC.R2.fastq.gz | 14       | NA                                                           | NA      | NA                                                        | fungi_odb10    | /90daydata/coffea_pangenome/puzzler_trials/busco_downloads |

:page_with_curl: Map file descriptions (:exclamation: **Use full paths!!!**)

* **sample:** Sample ID, all assembly work will be saved in `$WD/$SAMPLE`.
* **runtime:** Either "apptainer" or "conda". Puzzler will automatically `--bind` necessary paths for apptainer.
* **container:** If conda, "NA", otherwise path to the apptainer `.sif`. Presumably works with Singularity (untested). 
* **wd:** Path to working directory to store all files.
* **hifi:** Path to HiFi reads.
* **hic_r1:** Path to HiC R1.
* **hic_r2:** Path to HiC R2.
* **chromosomes:** Chromosome number or best guess. Pipeline will attempt ±4 your estimate.

***OPTIONAL columns*** 
*Specify "NA" and the script will skip respective components.*

* **reference:** Path to related species genome for chromosome naming. 

:bomb: Scaffolds will be renamed to the closest syntenic chromosome **using the reference scaffold naming convention**. 

*Please ensure that this reference has appropriately named transferable chromosome labels, as only scaffold names containing 'Chr' or 'chr' will be renamed!*

For NCBI RefSeq genomes, this typically does the trick: `sed -i -e 's/>.*chromosome />chr/g; s/,.*//g' $reference.fna`

* **hom_cov:** Homozygous peak coverage, used for `hifiasm` if provided, otherwise write "NA".

This value is the homozygous peak coverage identified from k-mer coverage in the HiFi library. I prefer to quickly run genomescope2, where the `*_linear_plot.png` indicates the left peak with `kcov:`, which you can multiple by ploidy to get `--hom_cov`. 

If you do not run `genomescope2`, I **highly** recommend you inspect the hifiasm histogram within the log file to ensure it matches the provided homozygous peak (within `$WD/$SAMPLE/$SAMPLE.hifiasm.log`, line: `[M::purge_dups] homozygous read coverage threshold: X`).

For more details about this homozygous coverage and for a full example, see the [Fungus example](/docs/Fungus_Example.md). 

* **blob_database:** Directory to save all blobtools databases. *Note that this will add SIGNIFICANT time to your run, the database to download is > 200 Gb, and large genomes will run for > 24 hours even with 64 cores*. 
* **busco_lineage:** Busco odb10 version lineage.
* **busco_database:** Directory to save busco dbs.


<!-- TOC --><a name="details"></a>
## Details

:fire: `puzzler` is checkpoint-based and will skip already-completed tasks based on file existence and successful step completion (files with a size > 0 and `$WD/$SAMPLE/hifiasm.complete`). You can therefore trick the script by creating or copying manual files into the appropriate directories. 

For instance, if you submit the command `puzzler --sample Fungus --map samples.tsv` and you had already completed hifiasm, purge duplicates, and haphic steps, your `$WD/Fungus` directory will look like this:

```
Jun 16 14:02 01_hifiasm
Jun 16 17:34 02_purge_dups
Jun 16 14:02 03_haphic
Jun 16 13:44 hifiasm.complete
Jun 16 13:44 purge_dups.complete
Jun 16 13:44 align_hic.complete
Jun 16 13:44 scaffolding.complete
```

And since you miss the juicer outputs, the script will create the juicebox `.hic` files:

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
RUNTIME: apptainer
CONTAINER: /home/justin.merondun/apptainer/puzzler_v1.9.sif 
WD: /90daydata/coffea_pangenome/puzzler_trials/assemblies 
HIFI: /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Fungus.HiFi.fastq.gz
HIC_R1: /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Fungus.HiC.R1.fastq.gz
HIC_R2: /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Fungus.HiC.R2.fastq.gz
NUMBER CHRS: 14
REFERENCE: /90daydata/coffea_pangenome/puzzler_trials/raw_data/references/GCA_964106605.1_gdRhiKalk1.hap1.1_genomic.fna
HOM_COV: NA
BLOB_DB: /90daydata/coffea_pangenome/puzzler_trials/blob_downloads
BUSCO_LINEAGE: fungi_odb10
BUSCO_DB: /90daydata/coffea_pangenome/puzzler_trials/busco_downloads
Cores Requested: 24
Cores Available: 24
RAM Requested: 384
Memory Available: 2143.8 GB
PUZZLER command: apptainer exec --bind /90daydata/coffea_pangenome/puzzler_trials:/90daydata/coffea_pangenome/puzzler_trials --bind /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads:/90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads --bind /90daydata/coffea_pangenome/puzzler_trials/raw_data/references:/90daydata/coffea_pangenome/puzzler_trials/raw_data/references --bind /home/justin.merondun/apptainer:/home/justin.merondun/apptainer /home/justin.merondun/apptainer/puzzler_v1.8.sif
=======================================================================

~~~~ [+0m] Skipping assembly for Fungus: /90daydata/coffea_pangenome/puzzler_trials/assemblies/primary_asm/Fungus.fa exists ~~~~
~~~~ [+0m] Skipping purge for Fungus: /90daydata/coffea_pangenome/puzzler_trials/assemblies/Fungus/02_purge_dups/p_ctg.purged.fa exists ~~~~
~~~~ [+0m] Skipping HiC alignment for Fungus: /90daydata/coffea_pangenome/puzzler_trials/assemblies/Fungus/03_haphic/filtered.bam exists ~~~~
~~~~ [+0m] Skipping HapHiC for Fungus: /90daydata/coffea_pangenome/puzzler_trials/assemblies/Fungus/03_haphic/haphic/04.build/scaffolds.fa exists ~~~~
~~~~ [+1m] Creating .hic file for juicebox for Fungus, with reference alignment  ~~~~
```

:x: If `puzzler` fails, it will print an error like below. This will tell you which directory the command failed in, and the command which caused the failure. 

```
❌ Command failed in in /project/90daydata/coffea_pangenome/puzzler_trials/blob_downloads/nt: "${PUZZLER} update_blastdb.pl --force_ftp --num_threads ${t} --decompress nt > nt_check_database.log 2>&1" (line 677)
```

If you reach an error like this, please navigate to the specified directory (e.g. `/project/90daydata/coffea_pangenome/puzzler_trials/blob_downloads/nt`) and inspect the most recent log file for errors! 

:one: **Step 1:** `[sbatch] puzzler --sample Fungus --map samples.tsv`

This script will check for these files, and create them if they do not exist in these directories:

| File                                               | Completion?                      | Description                  | If missing, run      |
| -------------------------------------------------- | -------------------------------- | ---------------------------- | -------------------- |
| $WD/$SAMPLE/01_hifiasm/asm.hic.p_ctg.gfa           | $WD/$SAMPLE/hifiasm.complete     | Primary contig assembly      | hifiasm              |
| $WD/$SAMPLE/02_purge_dups/p_ctg.purged.fa          | $WD/$SAMPLE/purge_dups.complete  | Purged contig assembly       | purge_dups           |
| $WD/$SAMPLE/03_haphic/filtered.bam                 | $WD/$SAMPLE/align_hic.complete   | Filtered re-mapped HiC reads | bwa mem2             |
| $WD/$SAMPLE/03_haphic/haphic/04.build/scaffolds.fa | $WD/$SAMPLE/scaffolding.complete | HapHiC scaffolded assembly   | haphic pipeline/YAHS |
| $WD/primary_asm/juicer_files/$SAMPLE_JBAT.hic      | $WD/$SAMPLE/juicer.complete      | Juicebox input file          | juicer pre           |

<br/>

Note that the script will also output a `.complete` file for each step, in case the step didn't finish completely (e.g. ran out of time). **If you are manually copying in files to skip steps:** (e.g. you do a separate `purge_dups` process and copy in your own purged fasta, then make sure you create this file e.g. `touch $WD/$SAMPLE/purge_dups.complete` in addition to your `p_ctg.purged.fa`)

```
$WD/$SAMPLE/hifiasm.complete
$WD/$SAMPLE/purge_dups.complete
$WD/$SAMPLE/align_hic.complete
$WD/$SAMPLE/scaffolding.complete
$WD/$SAMPLE/juicer.complete
```

:mag: **Manual curation:** Juicebox! 


Open Juicebox, and drag the `.hic` file into the window. Import the `.assembly` file using `Assembly > Import Map Assembly`. Make any adjustments if necessary (only about 50% of my genomes need it), and then export the file with `Assembly > Export Assembly`. 

See example videos of manual curation on the 11 assembled videos on Youtube [here](https://www.youtube.com/watch?v=rMUiNqZwEpA). 

This will create e.g. `$SAMPLE_JBAT.review.assembly`. Copy that file as-is to `$WD/juicer_files/$SAMPLE_JBAT.review.assembly`. 

<br/>

:two: **Step 2:** resubmit `[sbatch] puzzler --sample Fungus --map samples.tsv`

This script will check for these files, and create them if they do not exist in these directories. The script will not start if you haven't added the Juicebox `.review` file or if you do not have a related reference genome. 

| File                                                      | Completion?                        | Description                     | If missing, run |
| --------------------------------------------------------- | ---------------------------------- | ------------------------------- | --------------- |
| $WD/primary_asm/juicer_files/$SAMPLE_JBAT.review.assembly |                                    | Juicebox output file            | NOTHING!        |
| $WD/$SAMPLE/05_postjuicebox/final_asm.fa                  |                                    | Renamed scaffold file           | seqkit renaming |
| $WD/$SAMPLE/06_realign_hic_hifi/final_asm.filtered.bam    | $WD/$SAMPLE/qc_align_hic.complete  | Final HiC re-mapped file        | bwa mem2        |
| $WD/primary_asm/stats/$SAMPLE.pdf                         |                                    | Final contact map               | haphic plot     |
| $WD/$SAMPLE/07_busco_yak_blob/sr.qv.txt                   | $WD/$SAMPLE/yak.complete           | YAK base quality file           | yak             |
| $WD/primary_asm/stats/$SAMPLE.busco.txt                   | $WD/$SAMPLE/busco.complete         | BUSCO stats                     | busco           |
| $WD/$SAMPLE/06_realign_hic_hifi/asm.hifi.bam              | $WD/$SAMPLE/qc_align_hifi.complete | Final HiFi re-mapped file       | minimap2        |
| $WD/$SAMPLE/07_busco_yak_blob/blast.out                   | $WD/$SAMPLE/blastn.complete        | Run blastn                      | blastn          |
| $WD/primary_asm/stats/$SAMPLE.blob.stats.txt              |                                    | Blobtools contaminants          | blobtools       |
| $WD/primary_asm/stats/$SAMPLE.summary.txt                 |                                    | Summarize assembly statistics   | assembly_stats  |

<br/>

:one: **Step 1(B):** `[sbatch] puzzler --sample Fungus --map samples.tsv --no_juice`

If you are feeling reckless, you can run the entire process in one command without juicebox manual curation. **I do not recommend this, as juicebox is essential, at least for fixing obvious misassemblies or merged chromosomes**. The command above would run the entire process from `hifiasm` to `assembly_stats`, using only the scaffolded assembly instead of the juicebox curated assembly. 

:one: **Step 1(C):** `[sbatch] puzzler --sample Fungus --map samples.tsv --no_purge`

If you would like to skip the `purge_dups` steps entirely, add this flag. It can of course be combined with `--no_juice`. 


<br/>

<!-- TOC --><a name="visual"></a>
## Visual Workflow

Below is an exhaustive workflow outline documenting the inputs and outputs for each `puzzler` step. Only three steps are required from raw reads to assembly summary statistics. 

![workflow](/examples/figs/Workflow.png)

<!-- TOC --><a name="steps"></a>
## Step Breakdown 

**1.**: `$sample/01_hifiasm/` 

Assembly. Default settings. 

```bash
hifiasm --primary -t ${t} -o asm --h1 ${HIC_R1} --h2 ${HIC_R2} ${HIFI}
```
If provided in `samples.tsv`, `--hom-cov` is given.

**2.**: `$sample/02_purge_dups/` 

Purge duplicates, skip if `--no_purge` is given. Default settings, no coverage map is given (only purge on sequence similarity).

```bash
split_fa pri.init.fa > pri.split.fa
minimap2 -t ${t} -xasm5 -DP pri.split.fa pri.split.fa 2> ${SAMPLE}.minimap.purge.log | gzip -c > pri.split.self.paf.gz
purge_dups pri.split.self.paf.gz > pri.dups.bed 2> ${SAMPLE}.purge.log
get_seqs pri.dups.bed pri.init.fa 2> "${SAMPLE}.getseqs.log"
```

**3.**: `$sample/03_haphic/` 

Align HiC reads. Default filtering recommended by HapHiC.

```bash
bwa-mem2 index ${WD}/${SAMPLE}/02_purge_dups/p_ctg.purged.fa > ${SAMPLE}.alignment.indexing.hic.log 2>&1
bwa-mem2 mem -5SP -t ${t} ${WD}/${SAMPLE}/02_purge_dups/p_ctg.purged.fa ${HIC_R1} ${HIC_R2} | \
  samblaster | \
  samtools view - -@ ${t} -S -h -b -F 3340 | \
  ${FILTER_BAM} - 1 --nm 3 --threads ${t} | \
  samtools view - -b -@ ${t} -o filtered.bam; } > ${SAMPLE}.alignment.hic.log 2>&1
```

**3A.**: `$sample/03_haphic/` 

Attempt [HapHiC](https://github.com/zengxiaofei/HapHiC) with specified `$NUM_CHRS` from `samples.tsv`. If fails, try +/- 4 from that number. 

Default recommended 2 rounds of correction. Test `--max_inflation` up to 10 which accomodates a huge variation of genomes. 

```bash
haphic pipeline ${WD}/${SAMPLE}/02_purge_dups/p_ctg.purged.fa ../filtered.bam ${NUM_CHRS} --correct_nrounds 2 --max_inflation 10.0 --threads ${t} --processes ${t} 2> ../${SAMPLE}.haphic.log
```

**3B.**: `$sample/03_haphic/` 

If `3A` scaffolding fails, scaffold with YAHS. Default settings. 

```bash
samtools faidx ${WD}/${SAMPLE}/02_purge_dups/p_ctg.purged.fa
yahs ${WD}/${SAMPLE}/02_purge_dups/p_ctg.purged.fa filtered.bam > yahs.log 2>&1
```

**4.**: `$sample/04_juicer/`

Prepare juicer input files (skipped if `--no_juice`). Recommended settings from developers. Mashmap uses a segment length of 10 Kb and Percent Identity of 85% to assign chromosome names. 

```bash
# This step is skipped if no reference is given 
mashmap -r ${REFERENCE} -q ${WD}/${SAMPLE}/02_purge_dups/p_ctg.purged.fa -t ${t} -s 10000 --perc_identity 85 -o asm_to_ref.paf 2> mashmap.ref.log
haphic refsort ${WD}/${SAMPLE}/03_haphic/haphic/04.build/scaffolds.raw.agp asm_to_ref.paf > alignment.agp 2> refsort.log
# If no reference, start here 
juicer pre \
  -a -q 1 \
  -o haphic_JBAT \
  ${WD}/${SAMPLE}/03_haphic/filtered.bam \
  alignment.agp \
  ${WD}/${SAMPLE}/02_purge_dups/p_ctg.purged.fa.fai > haphic_JBAT.log 2>&1
grep PRE_C_SIZE haphic_JBAT.log | \
  awk '{print $2" "$3}' > chrom.sizes
${RUN_JUICERTOOLS} pre \
  -r 5000000,4000000,3000000,2000000,1500000,1000000,750000,500000,250000,100000,50000 \
  haphic_JBAT.txt \
  haphic_JBAT.hic \
  chrom.sizes > juicer_pre.log 2>&1
```

**5A.**: `$sample/05_postjuicebox/`

Extract post-juicer assembly. Default settings.

```bash
juicer post \
  -o haphic-post_JBAT \
  ${WD}/juicer_files/${SAMPLE}_JBAT.review.assembly \
  ${WD}/${SAMPLE}/04_juicer/haphic_JBAT.liftover.agp \
  ${WD}/${SAMPLE}/02_purge_dups/p_ctg.purged.fa 2> ${SAMPLE}.juicer.post.log
```

**5B.**: `$sample/05_postjuicebox/`

Rename scaffolds according to syntenic chromosomes in reference, extract final assembly.

```bash
mashmap -r ${REFERENCE} -q post_juicer_asm.fa -t ${t} -s 10000 --perc_identity 85 -o asmpost_to_paf.paf 2> mashmap.postjuicer.log
samtools faidx post_juicer_asm.fa
samtools faidx ${REFERENCE}
map_chromosomes --paf asmpost_to_paf.paf --fai post_juicer_asm.fa.fai --out map.txt --min_size 0.1 &> mapping_renaming.log
... 
[lots of bash logic ensuring the scaffolds are renamed, rename unlocalized scaffolds (e.g. secondary hits to chr1 as chr1_unloc1)]
...
seqkit replace --line-width 0 -p "(.*)" -r "{kv}" -k chromosome_naming_map.txt orient.fa > haphic_renamed_unord.fa 2> seqkit_renaming.log
samtools faidx haphic_renamed_unord.fa
xargs samtools faidx haphic_renamed_unord.fa < sorted_chr.txt > final_asm.fa
```

**6.**: `$sample/06_realign_hic_hifi/`

Re-align HiC to final assembly and create contact plot. Same settings used for filtering as above. 

```
bwa-mem2 index final_asm.fa > alignment.indexing.final_hic.log 2>&1
samtools faidx final_asm.fa

bwa-mem2 mem -5SP -t ${t} final_asm.fa ${HIC_R1} ${HIC_R2} | \
samblaster | \
  samtools view - -@ ${t} -S -h -b -F 3340 | \
  ${FILTER_BAM} - 1 --nm 3 --threads ${t} --remove-dup | \
  samtools view - -b -@ ${t} -o final_asm.filtered.bam; > alignment.final_hic.log 2>&1
mock_agp_file.py final_asm.fa > final_asm.agp
haphic plot --threads ${t} final_asm.agp final_asm.filtered.bam --bin_size 1000 --min_len 1 2> haphic_plot.log
```

**7.**: `$sample/07_busco_yak_blob/`

Final QC. Run YAK, and then optionally run BUSCO and BlobTools if db paths are given in `samples.tsv`. 

Developer tutorial settings for [YAK](https://github.com/lh3/yak) and [BlobTools](https://blobtoolkit.genomehubs.org/blobtools2/blobtools2-tutorials/getting-started-with-blobtools2/), BUSCO default settings.

```bash
# YAK
yak count -b37 -t ${t} -o ccs.yak ${HIFI} > yak.count.log 2>&1
yak qv -t ${t} -p -K3.2g -l100k ccs.yak ${WD}/primary_asm/${SAMPLE}.fa > sr.qv.txt 2> yak.qv.log

# BUSCO
busco --download ${BUSCO_LINEAGE} --download_path ${BUSCO_DB} > ${BUSCO_DB}/${BUSCO_LINEAGE}_busco_download.log 2>&1
busco -i ${WD}/primary_asm/${SAMPLE}.fa \
  -l ${BUSCO_DB}/lineages/${BUSCO_LINEAGE} \
  -m genome \
  -c ${t} \
  -o ${SAMPLE} \
  -f \
  --offline > busco.log 2>&1

# Align HiFi reads for Blobtools
minimap2 -t ${t} -ax map-hifi \
  ${WD}/primary_asm/${SAMPLE}.fa ${HIFI} 2> ${SAMPLE}.hifi.minimap.log | \
  samtools sort --threads ${t} -o asm.hifi.bam; > hifi.alignment.log 2>&1
samtools index -c --threads ${t} asm.hifi.bam

# Prepare BlobTools DB
cd ${BLOB_DB}/nt
NT_LOCK_FILE="${BLOB_DB}/download_nt.lock"
update_blastdb.pl --force_ftp --num_threads ${t} --decompress nt > nt_check_database.log 2>&1
rm "$NT_LOCK_FILE"

# Run BlobTools
blastn -task megablast \
  -db ${BLOB_DB}/nt/nt \
  -query ${WD}/primary_asm/${SAMPLE}.fa \
  -outfmt "6 qseqid staxids bitscore std" \
  -max_target_seqs 10 \
  -max_hsps 1 \
  -evalue 1e-25 \
  -num_threads ${t} \
  -out blast.out > blast.log 2>&1
blobtools create -i ${WD}/primary_asm/${SAMPLE}.fa -b ${WD}/${SAMPLE}/06_realign_hic_hifi/asm.hifi.bam -t blast.out \
  --db ${BLOB_DB}/data/nodesDB.txt -o ${SAMPLE}
blobtools view -i ${SAMPLE}.blobDB.json
blobtools plot -i ${SAMPLE}.blobDB.json --format pdf 
```

**8.**: `primary_asm/stats`

Summarize final assembly. 

```bash
assembly_stats ${WD}/primary_asm/${SAMPLE}.fa
...
[bash logic to compile all optional QC information, if run]
...
echo -e "\e[46m~~~ $(elapsed) Your final assembly for ${SAMPLE} is: ${WD}/primary_asm/${SAMPLE}.fa ~~~\e[0m"
echo -e "\e[46m~~~ $(elapsed) Your final assembly stats for ${SAMPLE} are in: ${WD}/primary_asm/stats/${SAMPLE}.summary.txt ~~~\e[0m"
```


<!-- TOC --><a name="faq"></a>
## FAQ

> Do I need to specify the Hi-C enzyme motifs for HapHiC?

[HapHiC](https://github.com/zengxiaofei/HapHiC) is very good at recognizing fragments without specification, I haven't needed to do this - although I have mainly tested with Arima kits. 

___

> Do I need to trim my reads? 

Current HiFi libraries are predominantly highly accurate and adapter free. Integrate adapter trimming yourself if you would like, perhaps whilst you concatenate your libraries into a single file for Puzzler, using [HiFiAdapterFilt](https://github.com/sheinasim-USDA/HiFiAdapterFilt) or something similar, although I have found this will remove no reads or bases with current libraries. I have also not found this to be necessary for Hi-C reads with internal trials using HapHiC, but of course, feel free to do so. 

___

> I can't launch the container, "FATAL:   container creation failed: failed to resolve session directory /var/apptainer/mnt/session: lstat /var/apptainer: no such file or directory"

I see this warning if you try to run apptainer on a login node, which is typically a setting set by IT. You can try modifying some paths IF NECESSARY, but I really recommend testing on a compute node or submitting to a job scheduler. 

Potentially try:

```bash
export APPTAINER_TMPDIR=/your/tmp/dir
export APPTAINER_CACHEDIR=$HOME/.apptainer_cache
mkdir -p "$APPTAINER_TMPDIR"
apptainer exe $SIF hifiasm
```

___

> I receive an error for downloading busco / taxdump / nt databases!

Puzzler will create a `.lock` file when downloading the respective databases for busco and blobtools so that it will not overwrite or interupt on-going downloads. If you submitted Puzzler in parallel, other jobs will exit if they encounter this lock file. This is so that the resources will be freed up, otherwise they might hang for many hours while the `nt` database downloads. Once the databases are finished downloading, simply resubmit those jobs and Puzzler will take back up where needed. 

```bash
~~~ [+0m] Another instance already dl-ing taxdump db. Resubmit puzzler when finished, or rm /project/90daydata/coffea_pangenome/puzzler_trials/blob_downloads/download_taxdump.lock! Exiting ~~~
```

Depending on your HPC environment, you might run into issues about the compute nodes timing out, or the busco server's being unreachable. This will require manual inspection on your part, either for busco:

```bash
PUZZLER="apptainer exec /home/justin.merondun/apptainer/puzzler_v1.9.sif"
# for conda, busco is already on $PATH! 
BUSCO_DB=/90daydata/coffea_pangenome/puzzler_trials/busco_downloads
BUSCO_LINEAGE=fungi_odb10
cd ${BUSCO_DB}/../
${PUZZLER} busco --download ${BUSCO_LINEAGE} --download_path ${BUSCO_DB} 
```

Refseq taxdump (blobtools requirement):

```bash
BLOB_DB=/90daydata/coffea_pangenome/puzzler_trials/blob_downloads
cd ${BLOB_DB}
mkdir -p data
wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz -P data/
tar zxf data/taxdump.tar.gz -C data/ nodes.dmp names.dmp
```

Or Refseq nt (blobtools requirement): 

```bash
PUZZLER="apptainer exec /home/justin.merondun/apptainer/puzzler_v1.8.sif"
BLOB_DB=/90daydata/coffea_pangenome/puzzler_trials/blob_downloads
cd ${BLOB_DB}/nt
${PUZZLER} update_blastdb.pl --force_ftp --decompress nt
```

___

> I encountered the warning: Multiple scaffolds corresponding to a single Chr for ${SAMPLE} INSPECT!

This means that there are multiple scaffolds/chromosomes in your draft which correspond to a single chromosome in the reference genome. The script will automatically rename the duplicates into e.g. `chr1` `chr1_unloc1` `chr1_unloc2` according to length. 

You should probably inspect `$WD/$SAMPLE/05_postjuicebox/chromosome_naming_map.txt` and to ensure that you are happy with this naming scheme.

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

The script will automatically rename these chr1 and chr1_unloc1, but you could also modify this manually and resubmit:

```
head ${WD}/${SAMPLE}/05_postjuicebox/chromosome_naming_map.txt
scaffold_1      chr1    +       42.46%  35386063
scaffold_2      chr1_unloc1  +       41.33%  9847223
scaffold_3      chr11   +       47.76%  28452307
scaffold_4      chr2    +       39.75%  31132341
```

If you did any manual edits and want to re-generate `final_asm.fa`, just fill in these fields and re-run puzzler.

___

> The script gave a warning: HapHiC for Fungus with 14 chrs failed, trying: 10 

If HapHic does not succeed with the specified chromosome numbers, it will re-run it with +/- 4 chromsomes. If you specified n = 14 chromosomes, `puzzler` will re-run HapHic with 10 - 18 chromoomes. If all of those fail, it will run scaffolding with YAHS instead. **Please ensure that the script is entirely finished running and isn't continuing attempts with HapHiC +/- 4 chrs before resubmitting it! Strange errors may occur if you submit the sample $sample simultaneously**. 

<!-- TOC --><a name="outputs"></a>
## Outputs 

Relative to the `$WD` path, the outputs will be: 

* The final assembly will be found in: `$WD/primary_asm/$SAMPLE.fa`.

* The final assembly stats will be in: `$WD/primary_asm/stats/$SAMPLE.summary.txt`, e.g.: 

* The final contact map for assessment will be in: `$WD/primary_asm/stats/$SAMPLE.pdf` and should look something like this: 

![Final_Map](/examples/figs/contact_example.png)


Each directory within `$WD/$SAMPLE` corresponds to an assembly step, where all logs for all steps are saved in `.log` files. Please inspect these before opening an issue to see what the problem might be. 

A successful run should look like this within `$WD/$SAMPLE/` (excluding BlobTools, which dramatically increases runtime!):

```
01_hifiasm/
02_purge_dups/
03_haphic/
04_juicer/
05_postjuicebox/
06_realign_hic_hifi/
07_busco_yak_blob/
hifiasm.complete
purge_dups.complete
align_hic.complete
scaffolding.complete
juicer.complete
qc_align_hic.complete
yak.complete
busco.complete
```

<!-- TOC --><a name="contact"></a>
## Contact

Please make a git issue with any problems, no matter how small! I will respond as quickly as possible and it could help others. Otherwise, email Justin at heritabilities [@] gmail.com 

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
* [YAHS v1.2.2 (if HapHiC failed)](https://github.com/c-zhou/yahs): https://doi.org/10.1093/bioinformatics/btac808

If using chromosome-renaming: 
* [merothon v0.4.2](https://github.com/merondun/merothon)
* [seqkit v2.10.0](https://bioinf.shenwei.me/seqkit/): https://doi.org/10.1002/imt2.191
* [mashmap v3.1.3](https://github.com/marbl/MashMap): https://doi.org/10.1093/bioinformatics/btad512

If assessing assembly quality / contaminants:
* [busco v5.8.3](https://busco.ezlab.org/): https://doi.org/10.1093/bioinformatics/btv351
* [yak v0.1-r69-dirty](https://github.com/lh3/yak): https://doi.org/10.1038/s41592-020-01056-5 (same as hifiasm)
* [blobtools v1.1.1](https://blobtools.readme.io/docs/what-is-blobtools): https://doi.org/10.12688/f1000research.12232.1
* [blastn v2.16.0+](https://sequenceserver.com/blog/citing-ncbi-blast/): https://doi.org/10.1186/1471-2105-10-421
* [NCBI RefSeq nt database](https://www.ncbi.nlm.nih.gov/refseq/publications/): https://doi.org/10.1093/nar/gkae1038
* [assembly_stats v0.1.2](https://github.com/MikeTrizna/assembly_stats): https://zenodo.org/record/3968774

If using optional software:
* [genomescope2 v2.0](https://github.com/tbenavi1/genomescope2.0): https://doi.org/10.1038/s41467-020-14998-3



<!-- TOC --><a name="changelog"></a>
## Changelog

**v1.9**: 
- added optional `--no_purge` to skip purge dups. 
- Created a complete conda install for puzzler including a HapHic v1.0.6 and juicer jar files bundle. 
- Removed inherent SLURM blocks. If users want to submit directly `sbatch puzzler -s $sample -m $sample_file.tsv`, they can add the slurm block to the script, and find it with `which puzzler` if installed with conda. 
- Changed smaller scaffold naming scheme from e.g. Chr1 (largest), Chr1A, Chr1B.. etc, to Chr1, Chr1_unloc1, Chr1_unloc2, in case more than 26 unlocalized scaffolds exist.  

**v1.8**: switch minimap2 to mashmap for draft-reference mapping, magnitudes faster and works with >500mb chromosomes. Add time counter.

**v1.7**: merge puzzler & puzzle_quality. 

**v1.6**: add accomodation for specifying no reference ("NA" in `samples.tsv`). 
