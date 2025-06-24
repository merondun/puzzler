![Puzzler](/examples/figs/logo.png)

Simple check-point aware shell script wrapper for high-throughput genome assembly.

> What's it for?!

Scalable genome assembly. As simple and portable as possible: only requires the `.sif` and a table with file paths. 

Optimized for container-capable SLURM resources. 

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

Only a containerized release is supported because we shouldn't waste our time with software install, and [apptainer](https://apptainer.org/docs/admin/main/installation.html) is is straightforward to install, even without root. 

The singularity / apptainer `.sif` can be yanked with: 
`apptainer pull --arch amd64 library://merondun/default/puzzler:v1.7` 

Depending on your architecture, you might need to update your apptainer libraries:

```
apptainer remote add --no-login SylabsCloud cloud.sycloud.io
apptainer remote use SylabsCloud
apptainer pull --arch amd64 library://merondun/default/puzzler:v1.7
```

The workhorse script `puzzler` is simply a check-point aware bash script. It can be submitted directly with slurm, e.g. `sbatch puzzler --sample Fungus --map samples.tsv`. 

For installation:  

1) Download the container `.sif`
2) Download the `puzzler` shell script, ideally add to path:

```
git clone git@github.com:merondun/puzzler.git
cd puzzler
./setup.sh
"Installation complete! You can now run:
  puzzler --sample $SAMPLE --map_file examples/samples.tsv

# Refresh environment if needed
source ~/.bashrc
```

Or, just grab the `.sh` file and add to your path: 

```
wget https://raw.githubusercontent.com/merondun/puzzler/main/bin/puzzler
chmod +x puzzler
INSTALL_PATH=$(dirname "$(realpath puzzler)")
grep -qxF "export PATH=\$PATH:$INSTALL_PATH" ~/.bashrc || echo "export PATH=\$PATH:$INSTALL_PATH" >> ~/.bashrc
export PATH="$PATH:$INSTALL_PATH"
```

3) Modify your SLURM/HPC settings within the `puzzler` script to accomodate your scheduler. 

```
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

**Conda install**

1b) As a last resort, create conda environment from `/apptainer/environment.yml`, and install [HapHiC](https://github.com/zengxiaofei/HapHiC) and [assembly_stats](https://github.com/MikeTrizna/assembly_stats). If you go this route, please troubleshoot software and inspect the `.logs` before posting issues. 

:exclamation: If you go the conda route, you must make sure that `juicer pre` is available on path (included in HapHic installation), and you will need to modify path to `java -Xmx${MEM}G -jar /opt/HapHiC/utils/juicer_tools.1.9.9_jcuda.0.8.jar` to the correct path within `puzzler`. Simply search for that line in the script and replace the `/opt/` path with the path to your `.jarfile`. 


<!-- TOC --><a name="workflow"></a>
## Workflow

The pipeline creates a collapsed completely *de novo* primary genome assembly using both HiFi and HiC data. The workflow is: 

:pushpin: **`sbatch puzzler`**

1) [Hifiasm](https://github.com/chhylp123/hifiasm) assembly using HiFi + HiC reads.
2) Purging of haplotigs using [purge_dups](https://github.com/dfguan/purge_dups). Typically with sufficient HiFi data for diploid organisms, `hifiasm` adequately purges most coverage-based diplotigs. This `puzzler` implementation is not coverage based - only sequence similarity based, which works well in many species with variable ploidy. 
3) Scaffolding with [HapHiC](https://github.com/zengxiaofei/HapHiC), or [YAHS](https://github.com/c-zhou/yahs) if HapHic fails. 
4) Creating [Juicebox](https://github.com/aidenlab/Juicebox) manual curation inputs (`.hic`, `.assembly`). 

:mag: **Manual Curation in Juicebox** 

:pushpin: **`sbatch puzzler`**

5) Finalizing assembly: extract juicer assembly, assign chromosome names based on a related reference, reverse complement according to reverse, using [merothon](https://github.com/merondun/merothon).
6) QC: Re-map HiC for final assembly check with [HapHiC](https://github.com/zengxiaofei/HapHiC). 
7) QC: [yak](https://github.com/lh3/yak) base quality check. 
8) QC: [busco](https://busco.ezlab.org/) check.
9) QC: [blobtools](https://blobtools.readme.io/docs/what-is-blobtools) contaminant check against uniprot. 
10) QC: [assembly_stats](https://github.com/MikeTrizna/assembly_stats) assembly statistic metrics. 


<!-- TOC --><a name="quick-start"></a>
## Quick Start

`puzzler` only requires two inputs: `--sample` and `--map`. Please see a [fungus example](/docs/Fungus_Example.md) and an example with [N = 11 species](/docs/All_N11_Genomes.md). 

**1. Make pipeline map file** 

:boom: *Most important!*

Prepare a `samples.tsv` which outlines all necessary pipeline components. An example template is found in [/examples/samples.tsv](/examples/samples.tsv)

| sample | runtime   | container                                        | wd                                                    | hifi                                                         | hic_r1                                                       | hic_r2                                                       | num_chrs | reference                                                    | hom_cov | blob_database                                             | busco_lineage  | busco_database                                             |
| ------ | --------- | ------------------------------------------------ | ----------------------------------------------------- | ------------------------------------------------------------ | ------------------------------------------------------------ | ------------------------------------------------------------ | -------- | ------------------------------------------------------------ | ------- | --------------------------------------------------------- | -------------- | ---------------------------------------------------------- |
| Beaver | apptainer | /home/justin.merondun/apptainer/puzzler_v1.7.sif | /90daydata/coffea_pangenome/puzzler_trials/assemblies | /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Beaver.HiFi.fastq.gz | /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Beaver.HiC.R1.fastq.gz | /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Beaver.HiC.R2.fastq.gz | 20       | /90daydata/coffea_pangenome/puzzler_trials/raw_data/references/GCF_047511655.1_mCasCan1.hap1v2_genomic.fna | NA      | /90daydata/coffea_pangenome/puzzler_trials/blob_downloads | mammalia_odb10 | /90daydata/coffea_pangenome/puzzler_trials/busco_downloads |
| Crane  | apptainer | /home/justin.merondun/apptainer/puzzler_v1.7.sif | /90daydata/coffea_pangenome/puzzler_trials/assemblies | /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Crane.HiFi.fastq.gz | /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Crane.HiC.R1.fastq.gz | /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Crane.HiC.R2.fastq.gz | 40       | /90daydata/coffea_pangenome/puzzler_trials/raw_data/references/GCA_964106855.1_bGruGru1.hap1.1_genomic.fna | 48      | /90daydata/coffea_pangenome/puzzler_trials/blob_downloads | aves_odb10     | /90daydata/coffea_pangenome/puzzler_trials/busco_downloads |
| Fish   | apptainer | /home/justin.merondun/apptainer/puzzler_v1.7.sif | /90daydata/coffea_pangenome/puzzler_trials/assemblies | /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Fish.HiFi.fastq.gz | /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Fish.HiC.R1.fastq.gz | /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Fish.HiC.R2.fastq.gz | 24       | NA                                                           | NA      | NA                                                        | NA             | NA                                                         |

:page_with_curl: Map file descriptions (:exclamation: Use full paths!!!)

* **sample:** Sample ID, all assembly work will be saved in `$WD/$SAMPLE`.
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
* **blob_database:** Directory to save all blobtools databases. *Note that this will add SIGNIFICANT time to your run, the database to download is > 200 Gb, and large genomes will run for > 24 hours with 64 cores*. 
* **busco_lineage:** Busco odb10 version lineage.
* **busco_database:** Directory to save busco dbs.

***Homozygous peak coverage*** is the homozygous peak covewrage identified from k-mer coverage in the HiFi library. I prefer to quickly run genomescope2, where the `*_linear_plot.png` indicates the left peak with `kcov:`, which you can multiple by ploidy to get `--hom_cov`. 

If you do not run `genomescope2`, I **highly** recommend you inspect the hifiasm log file (in `$WD/$SAMPLE/$SAMPLE.hifiasm.log`).

For more details about this homozygous coverage and for a full example, see the [Fungus example](/docs/Fungus_Example.md). 


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
CONTAINER: /home/justin.merondun/apptainer/puzzler_v1.7.sif 
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
=======================================================================

~~~~ Skipping assembly for Fungus: /90daydata/coffea_pangenome/puzzler_trials/assemblies/primary_asm/Fungus.fa exists ~~~~
~~~~ Skipping purge for Fungus: /90daydata/coffea_pangenome/puzzler_trials/assemblies/Fungus/02_purge_dups/p_ctg.purged.fa exists ~~~~
~~~~ Skipping HiC alignment for Fungus: /90daydata/coffea_pangenome/puzzler_trials/assemblies/Fungus/03_haphic/filtered.bam exists ~~~~
~~~~ Skipping HapHiC for Fungus: /90daydata/coffea_pangenome/puzzler_trials/assemblies/Fungus/03_haphic/haphic/04.build/scaffolds.fa exists ~~~~
~~~~ Creating .hic file for juicebox for Fungus, with reference alignment  ~~~~
```

:x: If `puzzler` fails, it will print an error like below. This will tell you which directory the command failed in, and the command which caused the failure. 

```
❌ Command failed in in /project/90daydata/coffea_pangenome/puzzler_trials/blob_downloads/nt: "${PUZZLER} update_blastdb.pl --force_ftp --num_threads ${t} --decompress nt > nt_check_database.log 2>&1" (line 677)
```

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

This will create e.g. `$SAMPLE_JBAT.review.assembly`. Maintain that file name, and copy it to `$WD/juicer_files/$SAMPLE_JBAT.review.assembly`. 

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
| $WD/$SAMPLE/07_busco_yak_blob/diamond_1mb.out             | $WD/$SAMPLE/blastx.complete        | Run diamond blastx              | diamond         |
| $BLOB_DB/uniprot/reference_proteomes.dmnd                 | $BLOB_DB/uniprot_db.complete       | Create diamond uniprot database | diamond         |
| $WD/primary_asm/stats/$SAMPLE.blob.stats.txt              |                                    | Blobtools contaminants          | blobtools       |
| $WD/primary_asm/stats/$SAMPLE.summary.txt                 |                                    | Summarize assembly statistics   | assembly_stats  |

<br/>

:one: **Step 1(B):** `[sbatch] puzzler --sample Fungus --map samples.tsv --no_juice`

If you are feeling reckless, you can run the entire process in one command without juicebox manual curation. **I do not recommend this, as juicebox is essential, at least for fixing obvious misassemblies or merged chromosomes**. The command above would run the entire process from `hifiasm` to `assembly_stats`, using only the scaffolded assembly instead of the juicebox curated assembly. 


<br/>

<!-- TOC --><a name="visual"></a>
## Visual Workflow

Below is an exhaustive workflow outline documenting the inputs and outputs for each `puzzler` step. Only three steps are required from raw reads to assembly summary statistics. 

![workflow](/examples/figs/Workflow.png)

<!-- TOC --><a name="faq"></a>
## FAQ

> Do I need to specify the Hi-C enzyme motifs for HapHiC?

[HapHiC](https://github.com/zengxiaofei/HapHiC) is very good at recognizing fragments without specification, I haven't needed to do this - although I have mainly tested with Arima kits. 

___

> I can't launch the container, "FATAL:   container creation failed: failed to resolve session directory /var/apptainer/mnt/session: lstat /var/apptainer: no such file or directory"

I see this warning if you try to run apptainer on a login node, which is typically a setting set by IT. You can try modifying some paths IF NECESSARY, but I really recommend testing on a compute node or submitting to a job scheduler. 

Potentially try:

```
export APPTAINER_TMPDIR=/your/tmp/dir
export APPTAINER_CACHEDIR=$HOME/.apptainer_cache
mkdir -p "$APPTAINER_TMPDIR"
apptainer exe $SIF hifiasm
```

___

> I receive an error for downloading busco / taxdump / nt databases!

Puzzler will create a `.lock` file when downloading the respective databases for busco and blobtools so that it will not overwrite or interupt on-going downloads. If you submitted Puzzler in parallel, other jobs will exit if they encounter this lock file. This is so that the resources will be freed up, otherwise they might hang for many hours while the `nt` database downloads. Once the databases are finished downloading, simply resubmit those jobs and Puzzler will take back up where needed. 

```
~~~ [+0m] Another instance already dl-ing taxdump db. Resubmit puzzler when finished, or rm /project/90daydata/coffea_pangenome/puzzler_trials/blob_downloads/download_taxdump.lock! Exiting ~~~
```

Depending on your HPC environment, you might run into issues about the compute nodes timing out, or the busco server's being unreachable. This will require manual inspection on your part, either for busco:

```
PUZZLER="apptainer exec /home/justin.merondun/apptainer/puzzler_v1.8.sif"
BUSCO_DB=/90daydata/coffea_pangenome/puzzler_trials/busco_downloads
BUSCO_LINEAGE=fungi_odb10
cd ${BUSCO_DB}/../
${PUZZLER} busco --download ${BUSCO_LINEAGE} --download_path ${BUSCO_DB} 
```

Refseq taxdump (blobtools requirement):

```
BLOB_DB=/90daydata/coffea_pangenome/puzzler_trials/blob_downloads
cd ${BLOB_DB}
mkdir -p data
wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz -P data/
tar zxf data/taxdump.tar.gz -C data/ nodes.dmp names.dmp
```

Or Refseq nt (blobtools requirement): 

```
PUZZLER="apptainer exec /home/justin.merondun/apptainer/puzzler_v1.8.sif"
BLOB_DB=/90daydata/coffea_pangenome/puzzler_trials/blob_downloads
cd ${BLOB_DB}/nt
${PUZZLER} update_blastdb.pl --force_ftp --decompress nt
```

___

> I encountered the warning: Multiple scaffolds corresponding to a single Chr for ${SAMPLE} INSPECT!

This means that there are multiple scaffolds/chromosomes in your draft which correspond to a single chromosome in the reference genome. The script will automatically rename the duplicates into e.g. `chr1` `chr1A` `chr1B` according to length. 

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

The script will automatically rename these chr1 and chr1A, but you could also modify this manually and resubmit:

```
head ${WD}/${SAMPLE}/05_postjuicebox/chromosome_naming_map.txt
scaffold_1      chr1    +       42.46%  35386063
scaffold_2      chr1A  +       41.33%  9847223
scaffold_3      chr11   +       47.76%  28452307
scaffold_4      chr2    +       39.75%  31132341
```

If you did any manual edits and want to re-generate `final_asm.fa`, just fill in these fields and re-run puzzler.

___

> The script gave a warning: HapHiC for Fungus with 14 chrs failed, trying: 10 

If HapHic does not succeed with the specified chromosome numbers, it will re-run it with +/- 4 chromsomes. If you specified n = 14 chromosomes, `puzzler` will re-run HapHic with 10 - 18 chromoomes. If all of those fail, it will run scaffolding with YAHS instead. 

<!-- TOC --><a name="outputs"></a>
## Outputs 

Relative to the `$WD` path, the outputs will be: 

* The final assembly will be found in: `$WD/primary_asm/$SAMPLE.fa`.

* The final assembly stats will be in: `$WD/primary_asm/stats/$SAMPLE.summary.txt`, e.g.: 

* The final contact map for assessment will be in: `$WD/primary_asm/stats/$SAMPLE.pdf` and should look something like this: 

![Final_Map](/examples/figs/contact_example.png)


Each directory within `$WD/$SAMPLE` corresponds to an assembly step, where all logs for all steps are saved in `.log` files. Please inspect these before opening an issue to see what the problem might be. 

A successful run should look like this within `$WD/$SAMPLE/`:

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

`puzzler`

**v1.8**: switch minimap2 to mashmap for draft-reference mapping, magnitudes faster and works with >500mb chromosomes. Add time counter.

**v1.7**: merge puzzler & puzzle_quality. 

**v1.6**: add accomodation for specifying no reference ("NA" in `samples.tsv`). 
