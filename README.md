![Puzzler](/examples/figs/logo.png)

Simple nextflow pipeline for assembly genomes from Hifi and HiC data. 

> What's it for?!

Genome assembly, primarily designed for SLURM sources on SciNet architecture. 

## Installation

Puzzler currently requires access to BUSCO5, but all other dependencies will be installable via `conda` soon.

## Workflow

The Nextflow DSL2 workflow follows these steps:

1) [Hifiasm](https://github.com/chhylp123/hifiasm) assembly using HiFi + HiC reads, creating a number of haplotypes corresponding to inferred ploidy from flow cytometry / kmer estimation.
2) Purging of haplotigs using [purge_dups](https://github.com/dfguan/purge_dups), may require some manual adjustment based on contiguity. Also performs first round of BUSCO checks.  
3) Re-scaffolding with [HapHiC](https://github.com/zengxiaofei/HapHiC). 
4) Final scaffolding to a known chromosome reference to ensure proper orientation and naming of chromosomes with [RagTag](https://github.com/malonge/RagTag).
5) Assembly statistics (BUSCO, [contiguity](https://pypi.org/project/assembly-stats/)).
6) Repeat annotation with [EarlGrey](https://github.com/TobyBaril/EarlGrey).
7) Gene annotation using any available intraspecific RNAseq data using [BRAKER](https://github.com/Gaius-Augustus/BRAKER) and the soft-masked genome.
