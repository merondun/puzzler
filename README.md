![Puzzler](/examples/figs/logo.png)

Simple pipeline for assembling genomes from Hifi and HiC data, stable release in :snake:make  

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

There is also a definition file for the build in `/apptainer/`. It requires many common tools including [hifiasm v0.20.0-r639](https://github.com/chhylp123/hifiasm), [purge_dups v1.2.5](https://github.com/dfguan/purge_dups), [gfatools v0.4-r214-dirty](https://github.com/lh3/gfatools), [minimap2 v2.28-r1209](https://github.com/lh3/minimap2), [samblaster v0.1.26](https://github.com/GregoryFaust/samblaster), [samtools v1.18](https://github.com/samtools/samtools), [bwa v0.7.18-r1243-dirty](https://bio-bwa.sourceforge.net/), [HapHiC v1.0.6](https://github.com/zengxiaofei/HapHiC), and [ragtag v.2.1.0](https://github.com/malonge/RagTag).

<!-- TOC --><a name="workflow"></a>
## Workflow

The pipeline follows these steps:

1) [Hifiasm](https://github.com/chhylp123/hifiasm) assembly using HiFi + HiC reads, creating a number of haplotypes corresponding to inferred ploidy from flow cytometry / kmer estimation.
   - If no HiC data is provided, the script will skip steps 3 - 4! Check the slurm / screen output to ensure the script recognized HiC reads. 
2) Purging of haplotigs using [purge_dups](https://github.com/dfguan/purge_dups), may require some manual adjustment based on contiguity. Note that this implementation **does not** purge haplotigs according to coverage. This is designed for polyploid genomes, so this only removes based on sequence content! Based on numerous sensitivities, this seems to be preferable for the species assayed at least so far. 
3) Re-scaffolding with [HapHiC](https://github.com/zengxiaofei/HapHiC). 
4) Final scaffolding to a known chromosome reference to ensure proper orientation and naming of chromosomes with [RagTag](https://github.com/malonge/RagTag). Note that this could dramatically improve your scaffold N50, particularly if you have a very low N50 from poor HiC data to begin with. Note that RagTag will not break any assemblies, so it won't impact inversions, etc. 
5) Assembly statistics (BUSCO, [contiguity](https://github.com/MikeTrizna/assembly_stats)). **Note that BUSCO is set to use `embryophyta_odb10` as the database, so please modify that within the `Snakefile` is necessary. 

***Potential future additions***

6) Repeat annotation with [EarlGrey](https://github.com/TobyBaril/EarlGrey).
7) Gene annotation using any available intraspecific RNAseq data using [BRAKER](https://github.com/Gaius-Augustus/BRAKER) and the soft-masked genome.

<!-- TOC --><a name="quick-start"></a>
## Quick Start

**1.** Clone this repo and save the `.sif` container somewhere (see above). 

**2.** Prepare a `samples.csv` indicating sample name, ploidy, and reads:

```
sample,ploidy,hifi,hic_r1,hic_r2
HART001,2,/project/coffea_pangenome/Artocarpus/Assemblies/20241115_JustinAssemblies/HART001/HART001.HiFi.fastq.gz,/project/coffea_pangenome/Artocarpus/Assemblies/20241115_JustinAssemblies/HART001/HART001.HiC.R1.fastq.gz,/project/coffea_pangenome/Artocarpus/Assemblies/20241115_JustinAssemblies/HART001/HART001.HiC.R2.fastq.gz
```

**3** If using snakemake (recommended):

:snake: Ensure snakemake is installed. 

:exclamation: Puzzler requires a `samples.csv` in the submission directory, as well as 4 self explanatory arugments:

```
./snakemake.sh 
=======================================================================
__________ ____ _______________________.____     _____________________ 
\______   \    |   \____    /\____    /|    |    \_   _____/\______   \
 |     ___/    |   / /     /   /     / |    |     |    __)_  |       _/
 |    |   |    |  / /     /_  /     /_ |    |___  |        \ |    |   \
 |____|   |______/ /_______ \/_______ \|_______ \/_______  / |____|_  /
                           \/        \/        \/        \/         \/ 
=======================================================================
Error: Missing required arguments
Usage: sbatch ./snakemake.sh --sample <sample_name> --reference <reference_path> --wd <working_directory> --container <container_path>
Example: sbatch ./snakemake.sh --sample HART001 --reference /path/to/close_species_reference.fasta --wd /path/to/workdir --container /path/to/container.sif
```
**--sample** corresponds to the first column of `samples.csv`, and will run the pipeline on that sample 
**--reference** a closely-related species `.fa`, the assembly will have similarly oriented chromosomes
**--wd** directory where files will be saved 
**--container** path to the puzzler `.sif`

Now you can simply submit the script using a positional argument, matching the first line of your `samples.csv`:  

```
sbatch -J stats_HART001 snakemake.sh --sample HART001 --reference --wd /project/coffea_pangenome/Artocarpus/Assemblies/20241115_JustinAssemblies --container /project/coffea_pangenome/Software/Merondun/apptainers/puzzler_v1.1.sif
````

The script will create a directory within the pipeline directory for `sample_runs/${SAMPLE}` where any logs, configs, and the actual Snakefile are stored. 

:bomb: `snakemake.sh` understands `--dry-run`, so you can run `snakemake.sh --dry-run` on a login node and it will tell you which processes are completed, and which are yet to run:

```
./snakemake.sh --sample HART001 --reference /project/coffea_pangenome/Artocarpus/WholeGenomeAlignments/fastas/ASM2540343.fa --wd /project/coffea_pangenome/Artocarpus/Assemblies/20241115_JustinAssemblies --container /project/coffea_pangenome/Software/Merondun/apptainers/puzzler_v1.1.sif
=======================================================================
__________ ____ _______________________.____     _____________________ 
\______   \    |   \____    /\____    /|    |    \_   _____/\______   \
 |     ___/    |   / /     /   /     / |    |     |    __)_  |       _/
 |    |   |    |  / /     /_  /     /_ |    |___  |        \ |    |   \
 |____|   |______/ /_______ \/_______ \|_______ \/_______  / |____|_  /
                           \/        \/        \/        \/         \/ 
=======================================================================
Starting assembly pipeline for HART001
Ploidy: 2
HiFi reads: /project/coffea_pangenome/Artocarpus/Assemblies/20241115_JustinAssemblies/HART001/HART001.HiFi.fastq.gz
HiC reads: /project/coffea_pangenome/Artocarpus/Assemblies/20241115_JustinAssemblies/HART001/HART001.HiC.R1.fastq.gz, /project/coffea_pangenome/Artocarpus/Assemblies/20241115_JustinAssemblies/HART001/HART001.HiC.R2.fastq.gz
Assuming unrestricted shared filesystem usage.
host: ceres24-compute-20.ceres.scinet.usda.gov
Building DAG of jobs...
Using shell: /usr/bin/bash
Provided cores: 1 (use --cores to define parallelism)
Rules claiming more threads will be scaled down.
Conda environments: ignored
Job stats:
job                   count
------------------  -------
all                       1
summarize_assembly        1
total                     2

Select jobs to execute...
Execute 1 jobs...

[Wed Nov 27 14:18:51 2024]
localrule summarize_assembly:
    input: /project/coffea_pangenome/Artocarpus/Assemblies/20241115_JustinAssemblies/Primary_Assemblies/HART001.HifiasmHifiHiC-PurgeDups-HapHiC-RagTag.fa, /project/coffea_pangenome/Artocarpus/Assemblies/20241115_JustinAssemblies/HART001/HiCHiFi/03_HapHiC/04.build/scaffolds.fa, /project/coffea_pangenome/Artocarpus/Assemblies/20241115_JustinAssemblies/Primary_Assemblies/busco/HART001/short_summary.specific.embryophyta_odb10.HART001.txt
    output: /project/coffea_pangenome/Artocarpus/Assemblies/20241115_JustinAssemblies/Primary_Assemblies/assembly_summary_HART001.txt
    jobid: 8
    reason: Missing output files: /project/coffea_pangenome/Artocarpus/Assemblies/20241115_JustinAssemblies/Primary_Assemblies/assembly_summary_HART001.txt
    wildcards: sample=HART001
    resources: tmpdir=/local/bgfs/justin.merondun/13666519

Activating singularity image /project/coffea_pangenome/Software/Merondun/apptainers/puzzler_v1.1.sif
[Wed Nov 27 14:19:02 2024]
Finished job 8.
1 of 2 steps (50%) done
Select jobs to execute...
Execute 1 jobs...

[Wed Nov 27 14:19:02 2024]
localrule all:
    input: /project/coffea_pangenome/Artocarpus/Assemblies/20241115_JustinAssemblies/Primary_Assemblies/HART001.HifiasmHifiHiC-PurgeDups-HapHiC-RagTag.fa, /project/coffea_pangenome/Artocarpus/Assemblies/20241115_JustinAssemblies/Primary_Assemblies/busco/HART001/short_summary.specific.embryophyta_odb10.HART001.txt, /project/coffea_pangenome/Artocarpus/Assemblies/20241115_JustinAssemblies/Primary_Assemblies/haplotypes/HART001.hap1.purged.fa, /project/coffea_pangenome/Artocarpus/Assemblies/20241115_JustinAssemblies/Primary_Assemblies/haplotypes/HART001.hap2.purged.fa, /project/coffea_pangenome/Artocarpus/Assemblies/20241115_JustinAssemblies/Primary_Assemblies/assembly_summary_HART001.txt
    jobid: 0
    reason: Input files updated by another job: /project/coffea_pangenome/Artocarpus/Assemblies/20241115_JustinAssemblies/Primary_Assemblies/assembly_summary_HART001.txt
    resources: tmpdir=/local/bgfs/justin.merondun/13666519

[Wed Nov 27 14:19:02 2024]
Finished job 0.
2 of 2 steps (100%) done
```

snakemake has identified that the outputs from the previous steps exist, so it will only run assembly statistics. 

<!-- TOC --><a name="outputs"></a>
## Outputs 

Relative to the `wd` path, you can find these output files. Note that you can also tell snakemake to pick-up from certain steps if you copy these files into the appropriate directories. This can be useful if you want to iron out differences with e.g. purge dups, or had run hifiasm separately, you can just copy the output files to the folders that snakemake expects them: 

```
# Provided that you ran `snakemake.sh` HART001 as your sample
# Step 1 (hifiasm): 
HART001/HiCHiFi/HART001.hic.p_ctg.gfa
# Step 2 (purge_dups): 
HART001/HiCHiFi/02_PurgeDups/purged.fa
# Step 3 (HapHiC):
HART001/HiCHiFi/03_HapHiC/04.build/scaffolds.fa
# Step 4 (RagTag - final assembly!):
Primary_Assemblies/HART001.HifiasmHifiHiC-PurgeDups-HapHiC-RagTag.fa
# Step 5 (purge_dups - each haplotype, diploid here):
Primary_Assemblies/haplotypes/HART001.hap1.purged.fa
Primary_Assemblies/haplotypes/HART001.hap2.purged.fa
# Step 6 (BUSCO): 
Primary_Assemblies/busco/HART001/short_summary.specific.embryophyta_odb10.HART001.txt
# Step 7 (Statistics):
Primary_Assemblies/assembly_summary_HART001.txt
```


**3b.** (Deprecated) Nextflow pipeline:

Modify the `nextflow.config` file. **IMPORTANT**. 

:exclamation: Nextflow will submit jobs via slurm for you, so it must be aware of your HPC syntax. Modify partitions accordingly. 

:exclamation: You must modify these 4 lines in this file to full paths on your system. 

```
    wd = "/project/coffea_pangenome/Artocarpus/Assemblies/20241115_JustinAssemblies" # Where your genomes will be saved
    workDir = "/project/coffea_pangenome/Artocarpus/tmp_nf" # tmp directory for all the nextflow intermediate steps
    reference = "/project/coffea_pangenome/Artocarpus/WholeGenomeAlignments/fastas/ASM2540343.fa" # the reference genome you want to orient your chromosomes against
    samples = "${params.wd}/samples.csv" # the samples file, as above
```
Submit. You can run this in a screen job, but I generally just prefer a single core compute job so I can track easily with a slurm output log:

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

<!-- TOC --><a name="contact"></a>
## Contact

Either make an issue or send an email to Justin at heritabilities [@] gmail.com 