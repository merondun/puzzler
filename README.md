![Puzzler](/examples/figs/logo.png)

Simple pipeline for assembling genomes from Hifi and HiC data. 

> What's it for?!

Genome assembly, primarily designed for SLURM sources on SciNet architecture. 

<!-- TOC start (generated with https://github.com/derlin/bitdowntoc) -->

- [Installation](#installation)
- [Workflow](#workflow)
- [Quick Start](#quick-start)
- [Details](#details)
- [Outputs ](#outputs)
- [Contact](#contact)

<!-- TOC end --> 

<!-- TOC --><a name="installation"></a>
## Installation

Currently, only a containerized release is supported. This includes all software necessary. 

This singularity / apptainer `.sif` container can be pulled from `singularity pull --arch amd64 library://merondun/default/puzzler:latest` 

The (large, ~4Gb) container contains many common tools including [hifiasm v0.20.0-r639](https://github.com/chhylp123/hifiasm), [purge_dups v1.2.5](https://github.com/dfguan/purge_dups), [gfatools v0.4-r214-dirty](https://github.com/lh3/gfatools), [minimap2 v2.28-r1209](https://github.com/lh3/minimap2), [samblaster v0.1.26](https://github.com/GregoryFaust/samblaster), [samtools v1.18](https://github.com/samtools/samtools), [bwa v0.7.18-r1243-dirty](https://bio-bwa.sourceforge.net/), and [HapHiC v1.0.6](https://github.com/zengxiaofei/HapHiC). It contains a lot of useful assembly and phasing tools, as I also use this for tools in development and as a 'stable' back-up for software compatability. 

As a last resort, you can create a conda environment from the `environment.yml` file in `/apptainer/`, and install [HapHiC](https://github.com/zengxiaofei/HapHiC). 

<!-- TOC --><a name="workflow"></a>
## Workflow

The pipeline creates a primary genome assembly using both HiFi and HiC data. The workflow is: 

:pushpin: **SCRIPT 1: `puzzler_asm`**

1) [Hifiasm](https://github.com/chhylp123/hifiasm) assembly using HiFi + HiC reads.
2) Purging of haplotigs using [purge_dups](https://github.com/dfguan/purge_dups), may require some manual adjustment based on contiguity, although if you provided accurate homozygous coverage peak information to hifiasm (`$HOM_COV`), I have found I never need to tweak this, especially for diploids and autopolyploids. Note that this implementation **does not** purge haplotigs according to coverage. I have designed this with polyploid assembly in mind, so this only removes based on sequence content! Based on numerous sensitivities, this seems to be preferable for the species assayed so far. 
3) Scaffolding with [HapHiC](https://github.com/zengxiaofei/HapHiC). 
4) Creating [Juicebox](https://github.com/aidenlab/Juicebox) manual curation inputs (`.hic`, `.assembly`). 

:mag: **Manual Curation in Juicebox** 

:pushpin: **SCRIPT 2: `puzzler_post`**

5) Finalizing assembly, assign chromosome names based on a related reference. 
6) Re-map HiC for final assembly check. 

<!-- TOC --><a name="quick-start"></a>
## Quick Start

**1. Clone** 

Clone this repo and pull the `.sif` container somewhere (see above). The two workhorse scripts are simply shell scripts, so you can add them to your bashrc with `./setup.sh`. 

You should see:

```
./setup.sh 
Adding Puzzler Pipeline to PATH...
Installation complete! You can now run:
  puzzler_asm --sample --map_file
  puzzler_post --sample --map_file
puzzler_asm
Error: --sample argument is required.
```

Otherwise, simply make the `/bin/puzzler_asm` and `/bin/puzzler_post` accessible and executable on the path.

**2. Make pipeline map file** 

:boom: *Most important!*

Prepare a `samples.csv` which outlines all necessary pipeline components. An example template is found in [/examples/samples.csv](/examples/samples.csv)

```
sample,Container,WD,ploidy,chromosomes,hom_cov,hifi,hic_r1,hic_r2,reference
HART001,/project/coffea_pangenome/Software/Merondun/apptainers/puzzler_v1.1.sif,/project/coffea_pangenome/Artocarpus/Assemblies/20250101_JustinAssemblies,2,28,68,/project/coffea_pangenome/Artocarpus/Concatenated_Reads/HART001.HiFi.fastq.gz,/project/coffea_pangenome/Artocarpus/Concatenated_Reads/HART001.HiC.R1.fastq.gz,/project/coffea_pangenome/Artocarpus/Concatenated_Reads/HART001.HiC.R2.fastq.gz,/project/coffea_pangenome/Artocarpus/Concatenated_Reads/ASM2540343.fa
```

:page_with_curl: Map file descriptions (:exclamation: Use full paths!!!)

* **sample:** Sample ID, all assembly work will be saved in $WD/$sample
* **container:** Path to the apptainer `.sif`
* **wd:** Path to working directory to store all files 
* **ploidy:** Ploidy of organism
* **chromosomes:** Number of chromosomes
* **hom_cov:** Homozygous peak coverage (see below) 
* **hifi:** Path to HiFi reads 
* **hic_r1:** Path to HiC R1 
* **hic_r2:** Path to HiC R2
* **reference:** Path to related species genome for chromosome naming. Scaffolds will be renamed to the closest syntenic chromosome using their naming convention. 

***Homozygous peak coverage*** is the left-most peak identified from k-mer coverage in the HiFi library. I prefer to quickly run genomescope2, where the `*_linear_plot.png` indicates the peak with `kcov:`, which you should manually confirm. You can also run hifiasm initially and check the hifiasm log (in $WD/logs/$SAMPLE.hifiasm.log). For a triploid with 43Gb of HiFi data and an estimated genome size of 800Mb, this left-most peak was around 16x.

<!-- TOC --><a name="details"></a>
## Details

:fire: `puzzler_asm` and `puzzler_post` are checkpoint-based, so they will skip already-completed tasks based on file existence (files with a size > 0). You can therefore trick the script by creating or copying manual files into the appropriate directories. 


:one: **Script 1:** `puzzler_asm`

This script will check for these files, and create them if they do not exist in these directories:

| Directory                                       | File                    | Description                  | If missing, run |
| ----------------------------------------------- | ----------------------- | ---------------------------- | --------------- |
| ${WD}/${SAMPLE}                                 | ${SAMPLE}.hic.p_ctg.gfa | Primary contig assembly      | hifiasm         |
| ${WD}/${SAMPLE}/01_scaffolding                  | all.purged.fa           | Purged contig assembly       | purge_dups      |
| ${WD}/${SAMPLE}/01_scaffolding                  | filtered.MQ1.bam        | Filtered re-mapped HiC reads | bwa mem         |
| ${WD}/${SAMPLE}/01_scaffolding/haphic/04.build/ | scaffolds.fa            | HapHiC scaffolded assembly   | haphic pipeline |
| ${WD}/logs/juicer/                              | ${SAMPLE}_JBAT.hic      | Juicebox input file          | juicer pre      |
<br/>

:mag: **Manual curation:** Juicebox! 

<br/>

Open Juicebox, and drag the `.hic` file into the window. Import the `.assembly` file using `Assembly > Import Map Assembly`. Make any adjustments if necessary (only about 50% of my genomes need it), and then export the file with `Assembly > Export Assembly`. 

This will create e.g. `${SAMPLE}.pri-MQ1_JBAT.review.assembly`. Maintain that file name, and copy it to ${WD}/${SAMPLE}/01_scaffolding. 

<br/><br/>

:one: **Script 2:** `puzzler_post`

This script will check for these files, and create them if they do not exist in these directories. The script will not start if you haven't added the Juicebox `.review` file. 

| Directory                      | File                           |                            | If missing, run |
| ------------------------------ | ------------------------------ | -------------------------- | --------------- |
| ${WD}/${SAMPLE}/01_scaffolding | ${SAMPLE}_JBAT.review.assembly | Juicebox output file       | NOTHING!        |
| ${WD}/${SAMPLE}/01_scaffolding | map.txt                        | Scaffold renaming file map | minimap2        |
| ${WD}/${SAMPLE}/01_scaffolding | haphic_renamed.fa              | Renamed scaffold file      | seqkit renaming |
| ${WD}/${SAMPLE}/01_scaffolding | haphic_renamed.filtered.bam    | Final HiC re-mapped file   | bwa mem         |
| ${WD}/logs/contact_maps/       | ${SAMPLE}.pdf                  | Final contact map          | haphic plot     |

<br/>
:bomb: :warning: **If you then encounter this warning during** `puzzler_post`: 

`~~~~ Multiple scaffolds corresponding to a single Chr for ${SAMPLE} INSPECT!  ~~~~` 

You must stop and inspect `${WD}/${SAMPLE}/01_scaffolding/02_orienting/hap_newID_map.txt` and deal with the duplicate 'Chr', because you have multiple scaffolds corresponding to a single reference chromosome. Rename them something unique (e.g. Chr1 and Chr1A) in column 2.

For example, in this file:

```
head ${WD}/${SAMPLE}/01_scaffolding/02_orienting/hap_newID_map.txt
scaffold_1      chr1    +       42.46%  35386063
scaffold_2      chr1  +       41.33%  9847223
scaffold_3      chr11   +       47.76%  28452307
scaffold_4      chr2    +       39.75%  31132341
```

We have two scaffolds with a high similarity (~42%) and a large amount of aligned bases (35 Mb and 9.8 Mb) to reference chr1. 

Simply rename the second one to e.g. 

```
head ${WD}/${SAMPLE}/01_scaffolding/02_orienting/hap_newID_map.txt
scaffold_1      chr1    +       42.46%  35386063
scaffold_2      chr1A  +       41.33%  9847223
scaffold_3      chr11   +       47.76%  28452307
scaffold_4      chr2    +       39.75%  31132341
```

And then run: 

```
cd $WD/$SAMPLE/02_$ITHapHiC/
seqkit replace --line-width 0 -p "(.*)" -r "kv" -k ${WD}/${SAMPLE}/01_scaffolding/02_orienting/hap_newID_map.txt ${WD}/${SAMPLE}/01_scaffolding/02_orienting/orient.fa | \
    seqkit sort --line-width 0 -n > ${WD}/${SAMPLE}/01_scaffolding/haphic_renamed.fa
```

Then, re-submit the script: `puzzler_post --sample HART001 --map-file samples.csv` and the script will identify the `haphic_renamed.fa`. 

<!-- TOC --><a name="outputs"></a>
## Outputs 

Relative to the `wd` path, the outputs will be: 

* The final assembly will be found in: `${WD}/primary_asm/${SAMPLE}.pri.fa`.

* The final contact map for assessment will be in: `${WD}/logs/contact_maps/${SAMPLE}.pdf` and should look something like this: 

![Final_Map](/examples/figs/contact_example.png)


Within `$WD/logs` you can find: <br/>

`*.hifiasm.log`: hifiasm logs <br/>
`/contact_maps/`: final pdf contact maps <br/>
`/haphic/`: haphic logs <br/>
`/juicer/`: juicer files <br/>


<!-- TOC --><a name="contact"></a>
## Contact

Either make an issue or send an email to Justin at heritabilities [@] gmail.com 
