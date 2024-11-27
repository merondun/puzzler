![Puzzler](/examples/figs/logo.png)

Simple pipeline for assembling genomes from Hifi and HiC data, implemented initially with nextflow but ported to Snakemake. 

> What's it for?!

Genome assembly, primarily designed for SLURM sources on SciNet architecture. 

## Installation

A singularity / apptainer `.sif` container can be pulled from `singularity pull --arch amd64 library://merondun/default/puzzler:latest` 

There is also a definition file for the build in `/apptainer/`. It requires many common tools including [hifiasm v0.20.0-r639](https://github.com/chhylp123/hifiasm), [purge_dups v1.2.5](https://github.com/dfguan/purge_dups), [gfatools v0.4-r214-dirty](https://github.com/lh3/gfatools), [minimap2 v2.28-r1209](https://github.com/lh3/minimap2), [samblaster v0.1.26](https://github.com/GregoryFaust/samblaster), [samtools v1.18](https://github.com/samtools/samtools), [bwa v0.7.18-r1243-dirty](https://bio-bwa.sourceforge.net/), [HapHiC v1.0.6](https://github.com/zengxiaofei/HapHiC), and [ragtag v.2.1.0](https://github.com/malonge/RagTag).

## Workflow

The pipeline follows these steps:

1) [Hifiasm](https://github.com/chhylp123/hifiasm) assembly using HiFi + HiC reads, creating a number of haplotypes corresponding to inferred ploidy from flow cytometry / kmer estimation.
   - If no HiC data is provided, the script will skip steps 3 - 4! Check the slurm / screen output to ensure the script recognized HiC reads. 
2) Purging of haplotigs using [purge_dups](https://github.com/dfguan/purge_dups), may require some manual adjustment based on contiguity. Note that this implementation **does not** purge haplotigs according to coverage. This is designed for polyploid genomes, so this only removes based on sequence content! Based on numerous sensitivities, this seems to be preferable for the species assayed at least so far. 
3) Re-scaffolding with [HapHiC](https://github.com/zengxiaofei/HapHiC). 
4) Final scaffolding to a known chromosome reference to ensure proper orientation and naming of chromosomes with [RagTag](https://github.com/malonge/RagTag). Note that this could dramatically improve your scaffold N50, particularly if you have a very low N50 from poor HiC data to begin with. Note that RagTag will not break any assemblies, so it won't impact inversions, etc. 

***Potential future additions***

5) Assembly statistics (BUSCO, [contiguity](https://pypi.org/project/assembly-stats/)).
6) Repeat annotation with [EarlGrey](https://github.com/TobyBaril/EarlGrey).
7) Gene annotation using any available intraspecific RNAseq data using [BRAKER](https://github.com/Gaius-Augustus/BRAKER) and the soft-masked genome.

## Quick Start

**1.** Clone this repo and save the `.sif` container somewhere. 

**2.** Prepare a `samples.csv` indicating sample name, ploidy, and reads:

```
sample,ploidy,hifi,hic_r1,hic_r2
HART001,2,/project/coffea_pangenome/Artocarpus/Assemblies/20241115_JustinAssemblies/HART001/HART001.HiFi.fastq.gz,/project/coffea_pangenome/Artocarpus/Assemblies/20241115_JustinAssemblies/HART001/HART001.HiC.R1.fastq.gz,/project/coffea_pangenome/Artocarpus/Assemblies/20241115_JustinAssemblies/HART001/HART001.HiC.R2.fastq.gz
```

**3** If using snakemake (recommended):

:snake: Ensure snakemake is installed. 

:exclamation: Modify these lines in the `snakemake.sh` file to match your system **IMPORTANT**.

```
# Where the assemblies are stored 
WORKDIR="/project/coffea_pangenome/Artocarpus/Assemblies/20241115_JustinAssemblies"
# The puzzler.sif
CONTAINER="/project/coffea_pangenome/Software/Merondun/apptainers/puzzler.sif"

...
# A closely-related species so you can orient your chromosomes similarly 
reference: "/project/coffea_pangenome/Artocarpus/WholeGenomeAlignments/fastas/ASM2540343.fa"
```

Now you can simply submit the script using a positional argument, matching the first line of your `samples.csv`:  `sbatch -J assembly_${SAMPLE} snakemake.sh ${SAMPLE}`.

The script will create a directory within the pipeline directory for `sample_runs/${SAMPLE}` where any logs, configs, and the actual Snakefile are stored. 



**3.b** If using the nextflow pipeline:

Modify the `nextflow.config` file. **IMPORTANT**. 

:exclamation: Nextflow will submit jobs via slurm for you, so it must be aware of your HPC syntax. Modify partitions accordingly. 

:exclamation: You must modify these 4 lines in this file to full paths on your system. 

```
    wd = "/project/coffea_pangenome/Artocarpus/Assemblies/20241115_JustinAssemblies" # Where your genomes will be saved
    workDir = "/project/coffea_pangenome/Artocarpus/tmp_nf" # tmp directory for all the nextflow intermediate steps
    reference = "/project/coffea_pangenome/Artocarpus/WholeGenomeAlignments/fastas/ASM2540343.fa" # the reference genome you want to orient your chromosomes against
    samples = "${params.wd}/samples.csv" # the samples file, as above
```
**4.** Submit. You can run this in a screen job, but I generally just prefer a single core compute job so I can track easily with a slurm output log:

```
#!/bin/bash
#SBATCH --job-name=puzzler_NF
#SBATCH --time=100:00:00
#SBATCH --cpus-per-task=1
#SBATCH --partition=long

# Navigate to the nextflow directory
cd /home/justin.merondun/puzzler/pipeline

# Run, resume if previously started
nextflow run main.nf -profile slurm -resume -with-apptainer /home/justin.merondun/apptainer/puzzler.sif
```

## Contact

Either make an issue or send an email to Justin at heritabilities [@] gmail.com 
