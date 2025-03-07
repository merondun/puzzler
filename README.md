![Puzzler](/examples/figs/logo.png)

Simple pipeline for assembling genomes from Hifi and HiC data, stable release using bash scripts. Currently recommended to `ln -s` scripts as they are still in progress. Final versions will be wrapped into commands and ported with e.g. conda. 

> What's it for?!

Genome assembly, primarily designed for SLURM sources on SciNet architecture. 

<!-- TOC start (generated with https://github.com/derlin/bitdowntoc) -->

- [Installation](#installation)
- [Workflow](#workflow)
- [Quick Start](#quick-start)
- [Outputs ](#outputs)
- [Contact](#contact)

<!-- TOC end --> 

<!-- TOC --><a name="installation"></a>
## Installation

A singularity / apptainer `.sif` container can be pulled from `singularity pull --arch amd64 library://merondun/default/puzzler:latest` 

**You will need access to BUSCO5 for the assembly statistics phase. This install is difficult to port into a container, so ensure it is available on your command line!** 

There is also a definition file for the build in `/apptainer/`. It requires many common tools including [hifiasm v0.20.0-r639](https://github.com/chhylp123/hifiasm), [purge_dups v1.2.5](https://github.com/dfguan/purge_dups), [gfatools v0.4-r214-dirty](https://github.com/lh3/gfatools), [minimap2 v2.28-r1209](https://github.com/lh3/minimap2), [samblaster v0.1.26](https://github.com/GregoryFaust/samblaster), [samtools v1.18](https://github.com/samtools/samtools), [bwa v0.7.18-r1243-dirty](https://bio-bwa.sourceforge.net/), and [HapHiC v1.0.6](https://github.com/zengxiaofei/HapHiC).

As a last resort, you can create a conda environment from the `environment.yml` file in `/apptainer/`, and install [HapHiC](https://github.com/zengxiaofei/HapHiC) and BUSCO5, and you are good to go. 

<!-- TOC --><a name="workflow"></a>
## Workflow

The pipeline creates both a final collapsed primary assembly (pri) and haplotype-phased assembly (hap) following these steps:

1) [Hifiasm](https://github.com/chhylp123/hifiasm) assembly using HiFi + HiC reads, creating a number of haplotypes corresponding to inferred ploidy.
2) Purging of haplotigs using [purge_dups](https://github.com/dfguan/purge_dups), may require some manual adjustment based on contiguity. Note that this implementation **does not** purge haplotigs according to coverage. This is designed for polyploid genomes, so this only removes based on sequence content! Based on numerous sensitivities, this seems to be preferable for the species assayed at least so far. 
3) Re-scaffolding with [HapHiC](https://github.com/zengxiaofei/HapHiC). 
4) Manual curation with juicebox. 
5) Assembly statistics (BUSCO, [contiguity](https://github.com/MikeTrizna/assembly_stats)). **Note that BUSCO is set to use `embryophyta_odb10` as the database, so please modify that within the `Snakefile` is necessary. 

***Potential future additions***

6) Repeat annotation with [EarlGrey](https://github.com/TobyBaril/EarlGrey).
7) Gene annotation using any available intraspecific RNAseq data using [BRAKER](https://github.com/Gaius-Augustus/BRAKER) and the soft-masked genome.

<!-- TOC --><a name="quick-start"></a>
## Quick Start

**1. Clone** 

Clone this repo and save the `.sif` container somewhere (see above). 

**2. Make sample map** 

Prepare a `samples.csv` indicating sample name, ploidy, number of chromosomes, homozygous peak coverage, reads, and a related species to use for assigning chromosome names:

```
sample,ploidy,chromosomes,hom_cov,hifi,hic_r1,hic_r2,reference
HART001,2,28,68,/project/coffea_pangenome/Artocarpus/Concatenated_Reads/HART001.HiFi.fastq.gz,/project/coffea_pangenome/Artocarpus/Concatenated_Reads/HART001.HiC.R1.fastq.gz,/project/coffea_pangenome/Artocarpus/Concatenated_Reads/HART001.HiC.R2.fastq.gz,/project/coffea_pangenome/Artocarpus/Concatenated_Reads/ASM2540343.fa
```

**3. Copy basefiles** 

Copy the `00_Assembly.sh` script into your wd, and modify these lines to indicate the high level assembly directory, your samples.csv file, and the puzzler `.sif` file. 

```
WD=/project/coffea_pangenome/Artocarpus/Assemblies/20250101_JustinAssemblies
SAMPLE_FILE="/project/coffea_pangenome/Artocarpus/Assemblies/20250101_JustinAssemblies/samples.csv"
PUZZLER="apptainer exec /project/coffea_pangenome/Software/Merondun/apptainers/puzzler_v1.1.sif"
```

:exclamation: `00_Assembly.sh` only requires a positional argument for sample ID which matches `samples.csv`, e.g. `sbatch 00_Assembly.sh HART001`. 

:exclamation: These scripts will skip already-completed tasks based on file existence (files with a size > 0). You can therefore trick the script by creating or copying manual files into the appropriate directories. 


**Under the hood steps for `00_Assembly.sh`:**

Within: `${WD}/${SAMPLE}`

***If doesn't exist :arrow_right: then:***

`${SAMPLE}.hic.hap1.p_ctg.gfa` :arrow_right: hifiasm assembly

\
\

The script then checks for these steps for both primary assembly `$IT="pri"` and haplotype assembly `$IT="hap"`

Within: `${WD}/${SAMPLE}/02_${IT}HapHiC`

***If doesn't exist :arrow_right: then:***

`all.purged.fa` :arrow_right: purge duplicates \
`chr.divergence.txt` :arrow_right: estimate haplotype divergence (only for $IT="hap") \
`filtered.MQ1.bam` :arrow_right: align HiC reads \
`01_haphicMQ1/04.build/scaffolds.fa` :arrow_right: run HapHiC \
`${WD}/logs/juicer/${SAMPLE}.${IT}-MQ1_JBAT.hic` :arrow_right: create juicer files \

\
:exclamation: The pipeline STOPs and will create Juicebox manual curation files `logs/juicer/*_MQ1_JBAT.hic` and `logs/juicer/*_MQ1_JBAT.assembly`. Load those into Juicebox, perform any edits if necessary, export with `Assembly > Export Assembly` and add the `MQ1_JBAT.review.assembly` to the sample's assembly directory: `${WD}/${SAMPLE}/02_${IT}HapHiC`

\
\

**4. Post-manual finalization**

The script `01_PostJuicebox.sh` then takes up provided that e.g. `HART001.pri_JBAT.review.assembly` is within `HART001/02_priHapHiC`.

This script will identify which scaffolds correspond to which chromosomes from a related fasta, ensure there are no duplicate chromosomes identified, and rename and ensure strand-alignment. It will then re-map the HiC reads and create the final hic contact map pdf. 


The script checks for these steps for both primary assembly `$IT="pri"` and haplotype assembly `$IT="hap"`

Within: `${WD}/${SAMPLE}/02_${IT}HapHiC`

***If doesn't exist :arrow_right: then:***

`map.txt` :arrow_right: align haphic scaffold fasta with other species to get scaffold ~ chr map \
`haphic_renamed.fa` :arrow_right: rename scaffolds and ensure strands are in alignment with reference \
`pg_renamed.filtered.bam` :arrow_right: re-align HiC to final assembly \
`${WD}/logs/contact_maps/${SAMPLE}.${IT}.pdf` :arrow_right: create final contact map pdf \

\


:exclamation: If you then encounter this warning: 

`~~~~ Multiple scaffolds corresponding to a single Chr for HART038 pri, INSPECT!  ~~~~` 

You must stop and inspect `${WD}/${SAMPLE}/02_${IT}HapHiC/02_orienting/hap_newID_map.txt` and deal with the duplicate 'Chr', because you have multiple scaffolds corresponding to a single reference chromosome. Rename them something unique (e.g. Chr1 and Chr1A) in column 2, and then re-run:

```
cd ${WD}/${SAMPLE}/02_${IT}HapHiC/
seqkit replace --line-width 0 -p "(.*)" -r "{kv}" -k 02_orienting/hap_newID_map.txt 02_orienting/orient.fa | \
    seqkit sort --line-width 0 -n > haphic_renamed.fa
ln -s haphic_renamed.fa pg_renamed.fa
```

Then, re-submit the script: `sbatch 01_PostJuicebox.sh HART001` 

<!-- TOC --><a name="outputs"></a>
## Outputs 

Relative to the `wd` path, the assemblies will be:
```
grep '>' ${WD}/joint_scaffold/HART001.pri.fa
>Chr01
>Chr02
```

and

```
grep '>' ${WD}/joint_scaffold/HART001.hap.fa
>HART001#1#Chr01
>HART001#1#Chr02
...
>HART001#2#Chr01
>HART001#2#Chr02
```

Within `${WD}/logs` you can find:

`*.hifiasm.log`: hifiasm logs 

`/contact_maps/`: final pdf contact maps

`/haphic/`: haphic logs 

`/juicer/`: juicer files


<!-- TOC --><a name="contact"></a>
## Contact

Either make an issue or send an email to Justin at heritabilities [@] gmail.com 