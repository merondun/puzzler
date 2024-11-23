![Puzzler](/examples/figs/logo.png)

Simple nextflow pipeline for assembly genomes from Hifi and HiC data. 

> What's it for?!

Genome assembly, primarily designed for SLURM sources on SciNet architecture. 

## Installation

A singularity / apptainer `.sif` container can be pulled from `singularity pull --arch amd64 library://merondun/default/puzzler:latest` 

There is also a definition file for the build in `/apptainer/`. It requires many common tools including [hifiasm v0.20.0-r639](https://github.com/chhylp123/hifiasm), [purge_dups v1.2.5](https://github.com/dfguan/purge_dups), [gfatools v0.4-r214-dirty](https://github.com/lh3/gfatools), [minimap2 v2.28-r1209](https://github.com/lh3/minimap2), [samblaster v0.1.26](https://github.com/GregoryFaust/samblaster), [samtools v1.18](https://github.com/samtools/samtools), [bwa v0.7.18-r1243-dirty](https://bio-bwa.sourceforge.net/), and [HapHiC v1.0.6](https://github.com/zengxiaofei/HapHiC), [ragtag v.2.1.0](https://github.com/malonge/RagTag).

## Workflow

The Nextflow DSL2 workflow follows these steps:

1) [Hifiasm](https://github.com/chhylp123/hifiasm) assembly using HiFi + HiC reads, creating a number of haplotypes corresponding to inferred ploidy from flow cytometry / kmer estimation.
2) Purging of haplotigs using [purge_dups](https://github.com/dfguan/purge_dups), may require some manual adjustment based on contiguity. Also performs first round of BUSCO checks.  
3) Re-scaffolding with [HapHiC](https://github.com/zengxiaofei/HapHiC). 
4) Final scaffolding to a known chromosome reference to ensure proper orientation and naming of chromosomes with [RagTag](https://github.com/malonge/RagTag).

***Potential future additions***

5) Assembly statistics (BUSCO, [contiguity](https://pypi.org/project/assembly-stats/)).
6) Repeat annotation with [EarlGrey](https://github.com/TobyBaril/EarlGrey).
7) Gene annotation using any available intraspecific RNAseq data using [BRAKER](https://github.com/Gaius-Augustus/BRAKER) and the soft-masked genome.

## tl;dr

**1.** Clone this repo and save the `.sif` container somewhere. 

**2.** Prepare a `samples.csv` indicating sample name, ploidy, and reads:

```
sample,ploidy,hifi,hic_r1,hic_r2
HART001,2,/project/coffea_pangenome/Artocarpus/Assemblies/20241115_JustinAssemblies/HART001/HART001.HiFi.fastq.gz,/project/coffea_pangenome/Artocarpus/Assemblies/20241115_JustinAssemblies/HART001/HART001.HiC.R1.fastq.gz,/project/coffea_pangenome/Artocarpus/Assemblies/20241115_JustinAssemblies/HART001/HART001.HiC.R2.fastq.gz
```

**3.** Modify the `nextflow.config` file. **IMPORTANT**. 

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
