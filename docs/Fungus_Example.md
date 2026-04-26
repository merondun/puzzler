# Needlecast fungus Genome Assembly with **`puzzler`**

[Needlecast fungus](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_964106605.1/) has a small genome (25 Mb) and 14 chromosomes, so it's good for a quick tutorial. 

First install puzzler (this tutorial currently followed v1.9.2 from conda). Install is large and with mamba this took around 10 minutes. The only difference from the current version is Merqury has replaced YAK, and FCS has replaced blobtools. 

```
mamba create -n puzzler -c conda-forge -c bioconda -c hcc -c heritabilities puzzler -y
```

<!-- TOC start (generated with https://github.com/derlin/bitdowntoc) -->

- [File Preparation](#file-preparation)
- [Homozygous Coverage Estimates](#optional-determine-homozygous-coverage)
- [Tutorial Overview](#tutorial-overview)
- [Recommended Workflow: Juicer & Reference](#recommended-workflow-juicer-and-reference)
- [Juicer & No Reference](#juicer-with-no-reference)
- [No Juicer & Reference](#no-juicer-but-with-reference)
- [No Juicer & No Reference](#no-juicer-and-no-reference)
- [Outputs](#outputs)
- [Alignments](#alignments)

<!-- TOC end --> 

## File Preparation

### Reads

First, grab the [HiFi](https://www.ncbi.nlm.nih.gov/sra/ERX12134039[accn]) and [HiC](https://www.ncbi.nlm.nih.gov/sra/ERX12138301[accn]) data, which are subset from those SRA links and archived on the [puzzler Zenodo repository](https://zenodo.org/records/17756307) to 100K HiFi reads and 10M HiC paired reads. 

I am going to run this on a compute node with 8 cores and 48 GB RAM (BUSCO tends to run out of memory with ~32GB). 

First, set your WD and download the subset data.

```bash

WD=/project/coffea_pangenome/Artocarpus/puzzler_trials
cd ${WD}
wget https://zenodo.org/records/17756307/files/Fungus.HiFi.fastq.gz
wget https://zenodo.org/records/17756307/files/Fungus.HiC.R1.fastq.gz
wget https://zenodo.org/records/17756307/files/Fungus.HiC.R2.fastq.gz
```

### Simplify Reference Chr Names

I will name the chromosomes according to the published reference assembly. I will clean the chromosome names so that they are e.g. ">chr1" intead of ">OZ066564.1":

```bash
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/964/106/605/GCA_964106605.1_gdRhiKalk1.hap1.1/GCA_964106605.1_gdRhiKalk1.hap1.1_genomic.fna.gz
gunzip GCA_964106605.1_gdRhiKalk1.hap1.1_genomic.fna.gz

# Check headers
grep '>' GCA_964106605.1_gdRhiKalk1.hap1.1_genomic.fna | head -n 2
>OZ066564.1 Rhizosphaera kalkhoffii genome assembly, chromosome: 1
>OZ066565.1 Rhizosphaera kalkhoffii genome assembly, chromosome: 2

# Clean headers 
sed -i 's/>.*chromosome: />chr/g' GCA_964106605.1_gdRhiKalk1.hap1.1_genomic.fna
grep '>' GCA_964106605.1_gdRhiKalk1.hap1.1_genomic.fna | head -n 2
>chr1
>chr2
```

### Sample input file:

Puzzler relies on a `.tsv` tab-separated metadata file which organizes all the sample read-specific information. You can either create this in excel, or for single-genome studies, simply fill out the required info below and the `printf` command will create the necessary file: 

```

# --- REQUIRED inputs (must be provided) ---
SAMPLE=Fungus
RUNTIME=conda
CONTAINER=NA
WD=/project/coffea_pangenome/Artocarpus/puzzler_trials
HIFI=$(realpath Fungus.HiFi.fastq.gz)
HIC_R1=$(realpath Fungus.HiC.R1.fastq.gz)
HIC_R2=$(realpath Fungus.HiC.R2.fastq.gz)
NUM_CHRS=14
REF=$(realpath GCA_964106605.1_gdRhiKalk1.hap1.1_genomic.fna)

# --- Optional inputs (defaults set to your example, can be overridden) ---
HOM_COV=NA
FCS_DB=NA
FCS_TAXID=NA
BUSCO_LINEAGE=fungi_odb10
BUSCO_DB=${WD}/busco_downloads

# --- Write header + one row to TSV ---
printf "sample\truntime\tcontainer\twd\thifi\thic_r1\thic_r2\tnum_chrs\treference\thom_cov\tfcs_database\tfcs_taxid\tbusco_lineage\tbusco_database\n" > ${WD}/samples.tsv
printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
  "$SAMPLE" "$RUNTIME" "$CONTAINER" "$WD" "$HIFI" "$HIC_R1" "$HIC_R2" "$NUM_CHRS" "$REF" "$HOM_COV" "$FCS_DB" "$FCS_TAXID" "$BUSCO_LINEAGE" "$BUSCO_DB" >> ${WD}/samples.tsv

```

Details:

:heavy_exclamation_mark:***Use full paths***!! 

* **sample:** Sample ID, all assembly work will be saved in `$WD/Fungus`.
* **runtime:** Either "apptainer", "singularity", or "conda". If other runtime, ensure `$runtime exec puzzler.sif` works.
* **container:** Path to the apptainer `.sif`. If all software available on path, simply write 'conda'. 
* **wd:** Path to working directory to store all files.
* **hifi:** Path to HiFi reads.
* **hic_r1:** Path to HiC R1.
* **hic_r2:** Path to HiC R2.
* **chromosomes:** Number of chromosomes, or best guess. The pipeline will attempt +/- 4 your estimate if unknown.

***OPTIONAL columns***: *Specify "NA" if not needed! The script will skip respective components if "NA".*

* **reference:** Path to related species genome for chromosome naming. Scaffolds will be renamed to the closest syntenic chromosome **using their scaffold naming convention**.
* **hom_cov:** Homozygous peak coverage.
* **fcs_database:** Directory to save the ~500gb fcs database.
* **fcs_taxonid:** NCBI taxon ID from [here](https://www.ncbi.nlm.nih.gov/taxonomy/). 
* **busco_lineage:** Busco odb10/odb12 version lineage.
* **busco_database:** Directory to save busco dbs.

## OPTIONAL Determine homozygous coverage

We can determine homozygous coverage peaks using our HiFi reads which will also gives a sense of HiFi quality. We want enough HiFi data that we have a clean peak that is easily separated from the errors. Necessary software is included with `puzzler` from conda v1.9.1 and up.

```bash
SAMPLE=Fungus
WD=/project/coffea_pangenome/Artocarpus/puzzler_trials
HIFI=$(realpath Fungus.HiFi.fastq.gz)

# Kmer size
K=31

mkdir -p ${WD}/${SAMPLE}/00_Genomescope
cd ${WD}/${SAMPLE}/00_Genomescope
FastK -v -t1 -k${K} -M70 -T8 -NFastK_Table_${SAMPLE} ${HIFI}
Histex -G FastK_Table_${SAMPLE} > ${SAMPLE}.hist
genomescope2 --input ${SAMPLE}.hist --output . --kmer_length ${K} --name_prefix ${SAMPLE}
rm FastK* # to save space because this file can be quite large 
```

Will output this plot:

![genomescope2](/examples/figs/Fungus_linear_plot.png)

The left-peak corresponds to around 18, so we could set `--hom_cov` to 32 if it was diploid, although honestly hifiasm is quite good if your coverage is sufficient. 

**Instead of specifying `--hom_cov` directly, I will just let use `hifiasm` determine appropriate levels**, and check the `.log`.

## Tutorial Overview

I will pursue 4 analytical forks for the Fungus sample to show the full pipeline breadth:

1) With `juicer` manual curation.
2) Same with (1), except without a `$REFERENCE` relative species.
3) Without `juicer` manual curation. 
4) Same with (3), except without a `$REFERENCE` relative species.

:stop_sign: Analyses (1) and (2) will stop at the manual curation step, giving us the `.hic` and `.assembly` file within `$WD/juicer_files` to load into juicebox. 

:rabbit: Analyses (3) and (4) go all the way to the end, one-command genome assembly and QC. **This is not recommended, as manual curation is essential for checking misassemblies!** 

_____

### Inputs

**Puzzler script:**

Navigate to `$WD` and ensure `puzzler` is on path:

```bash
WD=/project/coffea_pangenome/Artocarpus/puzzler_trials
cd ${WD}
puzzler -h 
Usage (v2.0.0): puzzler -s sample -m samples.tsv [OPTIONS]

Options:
  -s, --sample SAMPLE   Sample name, corresponding to the first column in the map file (required)
  -m, --map FILE        Path to .tsv/.csv map file (required)
  --threads t           Number of threads (optional; default 48)
  --mem MEM             Memory allocation (optional; default 512)
  --no_purge            Skip purge duplicates step entirely (optional)
  --no_juice            Skip juicer file creation entirely (optional; not recommended!)
  --yahs                Run scaffolding with YAHS directly instead of HapHiC (optional)
  -v, --version         Show version and exit
  -h, --help            Show help and exit

  Required --map Structure:
  The provided map file (e.g., samples.tsv) must contain the following columns in this order:
  SAMPLE RUNTIME CONTAINER WD HIFI HIC_R1 HIC_R2 NUM_CHRS REFERENCE HOM_COV FCS_DB FCS_TAXID BUSCO_LINEAGE BUSCO_DB
  For optional columns (REFERENCE - BUSCO_DB), write NA if undesired.
```

Good to go. 

:bulb: The most important part is to prepare the map file `samples.tsv` with these columns, **in this specific order.**

For this example with 4 forks, I will add 3 more rows where the only difference is `$SAMPLE` (where files are stored) and `$REFERENCE`. 

| sample              | runtime | container | wd                                                  | hifi                                                         | hic_r1                                                       | hic_r2                                                       | num_chrs | reference                                                    | hom_cov | fcs_database | fcs_taxid | busco_lineage | busco_database                                               |
| ------------------- | ------- | --------- | --------------------------------------------------- | ------------------------------------------------------------ | ------------------------------------------------------------ | ------------------------------------------------------------ | -------- | ------------------------------------------------------------ | ------- | ------------ | --------- | ------------- | ------------------------------------------------------------ |
| Fungus              | conda   | NA        | /project/coffea_pangenome/Artocarpus/puzzler_trials | /project/coffea_pangenome/Artocarpus/puzzler_trials/Fungus.HiFi.fastq.gz | /project/coffea_pangenome/Artocarpus/puzzler_trials/Fungus.HiC.R1.fastq.gz | /project/coffea_pangenome/Artocarpus/puzzler_trials/Fungus.HiC.R2.fastq.gz | 14       | /project/coffea_pangenome/Artocarpus/puzzler_trials/GCA_964106605.1_gdRhiKalk1.hap1.1_genomic.fna | NA      | NA           | NA        | fungi_odb10   | /project/coffea_pangenome/Artocarpus/puzzler_trials/busco_downloads |
| Fungus_NoRef        | conda   | NA        | /project/coffea_pangenome/Artocarpus/puzzler_trials | /project/coffea_pangenome/Artocarpus/puzzler_trials/Fungus.HiFi.fastq.gz | /project/coffea_pangenome/Artocarpus/puzzler_trials/Fungus.HiC.R1.fastq.gz | /project/coffea_pangenome/Artocarpus/puzzler_trials/Fungus.HiC.R2.fastq.gz | 14       | NA                                                           | NA      | NA           | NA        | fungi_odb10   | /project/coffea_pangenome/Artocarpus/puzzler_trials/busco_downloads |
| Fungus_NoJuice      | conda   | NA        | /project/coffea_pangenome/Artocarpus/puzzler_trials | /project/coffea_pangenome/Artocarpus/puzzler_trials/Fungus.HiFi.fastq.gz | /project/coffea_pangenome/Artocarpus/puzzler_trials/Fungus.HiC.R1.fastq.gz | /project/coffea_pangenome/Artocarpus/puzzler_trials/Fungus.HiC.R2.fastq.gz | 14       | /project/coffea_pangenome/Artocarpus/puzzler_trials/GCA_964106605.1_gdRhiKalk1.hap1.1_genomic.fna | NA      | NA           | NA        | fungi_odb10   | /project/coffea_pangenome/Artocarpus/puzzler_trials/busco_downloads |
| Fungus_NoRefNoJuice | conda   | NA        | /project/coffea_pangenome/Artocarpus/puzzler_trials | /project/coffea_pangenome/Artocarpus/puzzler_trials/Fungus.HiFi.fastq.gz | /project/coffea_pangenome/Artocarpus/puzzler_trials/Fungus.HiC.R1.fastq.gz | /project/coffea_pangenome/Artocarpus/puzzler_trials/Fungus.HiC.R2.fastq.gz | 14       | NA                                                           | NA      | NA           | NA        | fungi_odb10   | /project/coffea_pangenome/Artocarpus/puzzler_trials/busco_downloads |

## Recommended Workflow Juicer and Reference

### Step 1: Draft Assembly

Begin: 

```bash
puzzler -s Fungus -m samples.tsv --threads 8 --mem 48

=======================================================================
__________ ____ _______________________.____     _____________________ 
\______   \    |   \____    /\____    /|    |    \_   _____/\______   \
 |     ___/    |   / /     /   /     / |    |     |    __)_  |       _/
 |    |   |    |  / /     /_  /     /_ |    |___  |        \ |    |   \
 |____|   |______/ /_______ \/_______ \|_______ \/_______  / |____|_  /
                           \/        \/        \/        \/         \/ 
=======================================================================
                         ~~~ v1.9.2 ~~~

=======================================================================
Parameters for sample: Fungus 
RUNTIME: conda
CONTAINER: /project/coffea_pangenome/Artocarpus/puzzler_trials/NA (ignore if using conda)
WD: /project/coffea_pangenome/Artocarpus/puzzler_trials 
HIFI: /project/coffea_pangenome/Artocarpus/puzzler_trials/Fungus.HiFi.fastq.gz
HIC_R1: /project/coffea_pangenome/Artocarpus/puzzler_trials/Fungus.HiC.R1.fastq.gz
HIC_R2: /project/coffea_pangenome/Artocarpus/puzzler_trials/Fungus.HiC.R2.fastq.gz
NUMBER CHRS: 14
REFERENCE: /project/coffea_pangenome/Artocarpus/puzzler_trials/GCA_964106605.1_gdRhiKalk1.hap1.1_genomic.fna
HOM_COV: NA
BLOB_DB: NA
BUSCO_LINEAGE: fungi_odb10
BUSCO_DB: /project/coffea_pangenome/Artocarpus/puzzler_trials/busco_downloads
Cores Requested: 8
Cores Available: 8
RAM Requested: 48
Memory Available: 651.3 GB
PUZZLER command:  (blank if using conda)
=======================================================================

~~~ [+0m] Checking software availability, this will take about 1 minute ~~~
~~~ [+0m] Software check complete ~~~
~~~ [+0m] Assembling genome for Fungus ~~~
~~~ [+0m] Running juicer, script will stop after .hic files created ~~~
~~~ [+0m] Running purge_dups step ~~~
~~~ [+0m] Scaffolding with HapHiC, using YAHS only as a fallback ~~~
~~~ [+0m] Starting hifiasm assembly for Fungus ~~~
~~~ [+11m] Starting Purge_Dups for Fungus ~~~
~~~ [+11m] Mapping HiC reads to Fungus draft ~~~
~~~ [+14m] Running HapHiC for Fungus  ~~~

❌ Command failed in in /project/coffea_pangenome/Artocarpus/puzzler_trials/Fungus/03_haphic/haphic: "${PUZZLER} haphic pipeline ${WD}/${SAMPLE}/02_purge_dups/p_ctg.purged.fa ../filtered.bam ${NUM_CHRS} --correct_nrounds 2 --max_inflation 10.0 --threads ${t} --processes ${t} 2> ../${SAMPLE}.haphic.log" (line 350)

~~~ [+17m] HapHiC for Fungus with 14 chrs failed, trying: 10 ~~~

❌ Command failed in in /project/coffea_pangenome/Artocarpus/puzzler_trials/Fungus/03_haphic/haphic: "${PUZZLER} haphic pipeline ${WD}/${SAMPLE}/02_purge_dups/p_ctg.purged.fa ../filtered.bam ${CHR_ATTEMPT} --correct_nrounds 2 --max_inflation 10.0 --threads ${t} --processes ${t} 2> ../${SAMPLE}.haphic_N_chrs${CHR_ATTEMPT}.log" (line 373)

~~~ [+19m] HapHiC for Fungus with 14 chrs failed, trying: 11 ~~~

❌ Command failed in in /project/coffea_pangenome/Artocarpus/puzzler_trials/Fungus/03_haphic/haphic: "${PUZZLER} haphic pipeline ${WD}/${SAMPLE}/02_purge_dups/p_ctg.purged.fa ../filtered.bam ${CHR_ATTEMPT} --correct_nrounds 2 --max_inflation 10.0 --threads ${t} --processes ${t} 2> ../${SAMPLE}.haphic_N_chrs${CHR_ATTEMPT}.log" (line 373)

~~~ [+22m] HapHiC for Fungus with 14 chrs failed, trying: 12 ~~~

❌ Command failed in in /project/coffea_pangenome/Artocarpus/puzzler_trials/Fungus/03_haphic/haphic: "${PUZZLER} haphic pipeline ${WD}/${SAMPLE}/02_purge_dups/p_ctg.purged.fa ../filtered.bam ${CHR_ATTEMPT} --correct_nrounds 2 --max_inflation 10.0 --threads ${t} --processes ${t} 2> ../${SAMPLE}.haphic_N_chrs${CHR_ATTEMPT}.log" (line 373)

~~~ [+25m] HapHiC for Fungus with 14 chrs failed, trying: 13 ~~~

❌ Command failed in in /project/coffea_pangenome/Artocarpus/puzzler_trials/Fungus/03_haphic/haphic: "${PUZZLER} haphic pipeline ${WD}/${SAMPLE}/02_purge_dups/p_ctg.purged.fa ../filtered.bam ${CHR_ATTEMPT} --correct_nrounds 2 --max_inflation 10.0 --threads ${t} --processes ${t} 2> ../${SAMPLE}.haphic_N_chrs${CHR_ATTEMPT}.log" (line 373)

~~~ [+27m] HapHiC for Fungus with 14 chrs failed, trying: 14 ~~~

❌ Command failed in in /project/coffea_pangenome/Artocarpus/puzzler_trials/Fungus/03_haphic/haphic: "${PUZZLER} haphic pipeline ${WD}/${SAMPLE}/02_purge_dups/p_ctg.purged.fa ../filtered.bam ${CHR_ATTEMPT} --correct_nrounds 2 --max_inflation 10.0 --threads ${t} --processes ${t} 2> ../${SAMPLE}.haphic_N_chrs${CHR_ATTEMPT}.log" (line 373)

~~~ [+30m] HapHiC for Fungus with 14 chrs failed, trying: 15 ~~~

❌ Command failed in in /project/coffea_pangenome/Artocarpus/puzzler_trials/Fungus/03_haphic/haphic: "${PUZZLER} haphic pipeline ${WD}/${SAMPLE}/02_purge_dups/p_ctg.purged.fa ../filtered.bam ${CHR_ATTEMPT} --correct_nrounds 2 --max_inflation 10.0 --threads ${t} --processes ${t} 2> ../${SAMPLE}.haphic_N_chrs${CHR_ATTEMPT}.log" (line 373)

~~~ [+33m] HapHiC for Fungus with 14 chrs failed, trying: 16 ~~~

❌ Command failed in in /project/coffea_pangenome/Artocarpus/puzzler_trials/Fungus/03_haphic/haphic: "${PUZZLER} haphic pipeline ${WD}/${SAMPLE}/02_purge_dups/p_ctg.purged.fa ../filtered.bam ${CHR_ATTEMPT} --correct_nrounds 2 --max_inflation 10.0 --threads ${t} --processes ${t} 2> ../${SAMPLE}.haphic_N_chrs${CHR_ATTEMPT}.log" (line 373)

~~~ [+35m] HapHiC for Fungus with 14 chrs failed, trying: 17 ~~~

❌ Command failed in in /project/coffea_pangenome/Artocarpus/puzzler_trials/Fungus/03_haphic/haphic: "${PUZZLER} haphic pipeline ${WD}/${SAMPLE}/02_purge_dups/p_ctg.purged.fa ../filtered.bam ${CHR_ATTEMPT} --correct_nrounds 2 --max_inflation 10.0 --threads ${t} --processes ${t} 2> ../${SAMPLE}.haphic_N_chrs${CHR_ATTEMPT}.log" (line 373)

~~~ [+38m] HapHiC for Fungus with 14 chrs failed, trying: 18 ~~~

❌ Command failed in in /project/coffea_pangenome/Artocarpus/puzzler_trials/Fungus/03_haphic/haphic: "${PUZZLER} haphic pipeline ${WD}/${SAMPLE}/02_purge_dups/p_ctg.purged.fa ../filtered.bam ${CHR_ATTEMPT} --correct_nrounds 2 --max_inflation 10.0 --threads ${t} --processes ${t} 2> ../${SAMPLE}.haphic_N_chrs${CHR_ATTEMPT}.log" (line 373)

~~~ [+40m] HapHiC for Fungus failed, scaffolding with YAHS instead ~~~
~~~ [+40m] Creating .hic file for juicebox for Fungus ~~~
~~~ [+41m] Post curation assembly file missing for Fungus: Run Juicebox & place in /project/coffea_pangenome/Artocarpus/puzzler_trials/juicer_files/Fungus_JBAT.review.assembly ~~~
```

:alarm_clock: This is a good time to check that all of your paths and parameters look appropriate. 

 We should check `Fungus/01_hifiasm/Fungus.hifiasm.log` to ensure an appropriate `--hom_cov` was default selected by hifiasm:

```bash
 1: ****************************************************************************************************> 905745
 2: ****************************************************************************************************> 63784
 3: ********************************** 12624
 4: ************** 5273
 5: ******** 2871
 6: **** 1643
 7: *** 1091
 8: ** 842
 9: ** 677
10: * 383
11: * 306
12: * 298
13: * 413
14: * 443
15: ** 577
16: ** 904
17: *** 1281
18: ***** 1752
19: ****** 2310
20: ********* 3375
21: ************ 4347
22: **************** 5891
23: ********************** 8012
24: **************************** 10458
25: *********************************** 12942
26: ***************************************** 14919
27: ************************************************** 18171
28: ********************************************************* 20927
29: *************************************************************** 22993
30: ********************************************************************** 25613
31: ******************************************************************************* 28887
32: ************************************************************************************* 31331
33: ******************************************************************************************* 33257
34: ************************************************************************************************ 35302
35: **************************************************************************************************** 36699
36: ************************************************************************************************** 36086
37: *********************************************************************************************** 34944
38: ********************************************************************************************** 34616
39: *************************************************************************************** 31801
40: ************************************************************************************ 30938
41: ******************************************************************************* 28950
42: ************************************************************************** 27087
43: ****************************************************************** 24120
44: ********************************************************** 21313
45: *********************************************** 17428
46: **************************************** 14650
47: ********************************** 12535
48: ***************************** 10638
49: ************************* 9016
50: ******************** 7335
51: *************** 5624
52: ************ 4336
53: ********* 3417
54: ******* 2709
55: ***** 1883
56: **** 1415
57: *** 1273
58: *** 981
59: ** 821
60: ** 719
...
94: * 215
95: * 208
[M::ha_hist_line]  rest: ************************ 8842
[M::ha_analyze_count] left: none
[M::ha_analyze_count] right: none
[M::ha_pt_gen] peak_hom: 35; peak_het: -1
```

`hifiasm` automatically selected 35, while `genomescope2` estimated left peak at 18.1 (if assuming diploid, `--hom_cov` = 36), so they both agree and we can leave it as-is. 

:exclamation: **❌ Command failed:**

The output indicates that puzzler could not run [HapHiC](https://github.com/zengxiaofei/HapHiC) with the specified number of chromosomes (`$NUM_CHR` from `samples.tsv`). Puzzler will automatically atttempt HapHiC with ±4 chromosomes to see if HapHiC will instead run with that number (it tried 10-18 here). If that fails, puzzler will instead run [YAHS](https://github.com/c-zhou/yahs), which in my experience has always run successfully.

As long as you see:

`~~~ [+40m] Creating .hic file for juicebox for Fungus ~~~`

Afterwards, then the run was successful. 

Of the many genomes I have run so far, Fungus is the only genome not successfully scaffolded with HapHiC, maybe having to do with the small chromosome sizes. 

**Runtime:**

It completed the assembly in 41 minutes.

The pipeline will also run `haphic refsort`, provided a related species is provided, which will align your draft scaffolds in order of the reference genome for juicebox curation. **This is not reference-guided assembly** and will not otherwise modify the scaffolds. Please see the [HapHiC](https://github.com/zengxiaofei/HapHiC) `refsort` description for more details. 

If you did not provide a `$REFERENCE` file, it will still output the `.hic` and `.assembly` files, just without running `haphic refsort`. 

### Step 2: Juicebox Curation

Now we load the `$WD/juicer_files/Fungus_JBAT.hic` and `$WD/juicer_files/Fungus_JBAT.assembly` file into juicebox. 

Perform our edits. Pull in the `.hic` file, import the `.assembly` file with `Assembly > Import Map Assembly`. There's very minimal edits on this genome, I just merge two blocks which originate from one chromosome: 

To watch the video, please download the [Fungus_Ref.mp4](https://zenodo.org/records/15693026/files/Fungus_Ref.mp4?download=1) file from the manuscript's [Zenodo link](https://doi.org/10.5281/zenodo.15693025). 

Afterwards, save the file with `Assembly > Export Assembly`, and save the  `Fungus_JBAT.review.assembly` here: `$WD/juicer_files/Fungus_JBAT.review.assembly`, **retaining the default file name!**

### Step 3: Chromosome Renaming, Remapping HiC

Afterwards, resubmit the same `puzzler` command and the script will assign chromosome names, output a final HiC contact map on the final assembly, and output assembly statistics. **I skip blobtools, because that will take several hours to blast against the nt database...**. 

Resubmit the same command, and the script will take up where it left off: 

```bash
puzzler -s Fungus -m samples.tsv --threads 8 --mem 48

=======================================================================
__________ ____ _______________________.____     _____________________ 
\______   \    |   \____    /\____    /|    |    \_   _____/\______   \
 |     ___/    |   / /     /   /     / |    |     |    __)_  |       _/
 |    |   |    |  / /     /_  /     /_ |    |___  |        \ |    |   \
 |____|   |______/ /_______ \/_______ \|_______ \/_______  / |____|_  /
                           \/        \/        \/        \/         \/ 
=======================================================================
                         ~~~ v1.9.2 ~~~

=======================================================================
Parameters for sample: Fungus 
RUNTIME: conda
CONTAINER: /project/coffea_pangenome/Artocarpus/puzzler_trials/NA (ignore if using conda)
WD: /project/coffea_pangenome/Artocarpus/puzzler_trials 
HIFI: /project/coffea_pangenome/Artocarpus/puzzler_trials/Fungus.HiFi.fastq.gz
HIC_R1: /project/coffea_pangenome/Artocarpus/puzzler_trials/Fungus.HiC.R1.fastq.gz
HIC_R2: /project/coffea_pangenome/Artocarpus/puzzler_trials/Fungus.HiC.R2.fastq.gz
NUMBER CHRS: 14
REFERENCE: /project/coffea_pangenome/Artocarpus/puzzler_trials/GCA_964106605.1_gdRhiKalk1.hap1.1_genomic.fna
HOM_COV: NA
BLOB_DB: NA
BUSCO_LINEAGE: fungi_odb10
BUSCO_DB: /project/coffea_pangenome/Artocarpus/puzzler_trials/busco_downloads
Cores Requested: 8
Cores Available: 8
RAM Requested: 64
Memory Available: 739.2 GB
PUZZLER command:  (blank if using conda)
=======================================================================

~~~ [+0m] Checking software availability, this will take about 1 minute ~~~
~~~ [+0m] Software check complete ~~~
~~~ [+0m] Assembling genome for Fungus ~~~
~~~ [+0m] Running juicer, script will stop after .hic files created ~~~
~~~ [+0m] Running purge_dups step ~~~
~~~ [+0m] Skipping HapHic and scaffolding with YAHS directly ~~~
~~~ [+0m] Skipping hifiasm for Fungus: /project/coffea_pangenome/Artocarpus/puzzler_trials/Fungus/01_hifiasm/asm.hic.p_ctg.gfa exists ~~~
~~~ [+0m] Skipping purge for Fungus: /project/coffea_pangenome/Artocarpus/puzzler_trials/Fungus/02_purge_dups/p_ctg.purged.fa exists ~~~
~~~ [+0m] Skipping HiC alignment for Fungus: /project/coffea_pangenome/Artocarpus/puzzler_trials/Fungus/03_haphic/filtered.bam exists ~~~
~~~ [+0m] Skipping HapHiC for Fungus: /project/coffea_pangenome/Artocarpus/puzzler_trials/Fungus/03_haphic/haphic/04.build/scaffolds.fa exists ~~~
~~~ [+0m] Skipping juicer HiC file creation for Fungus: /project/coffea_pangenome/Artocarpus/puzzler_trials/juicer_files/Fungus_JBAT.hic exists ~~~
~~~ [+0m] Extracting post-curation assembly and mapping to reference for Fungus ~~~
~~~ [+0m] Renaming chromosomes for Fungus ~~~
~~~ [+0m] Scaffold sanity check passed for renaming, proceeding! ~~~
~~~ [+0m] Single scaffolds corresponding to a single Chr for Fungus ~~~
~~~ [+0m] Assessing genome quality for Fungus ~~~
~~~ [+0m] Creating final HiC bam for Fungus ~~~
~~~ [+3m] Creating final HiC contact map for Fungus ~~~
~~~ [+3m] Creating final contact map for Fungus, plotting named chromosomes ~~~
~~~ [+3m] Running YAK on Fungus ~~~
~~~ [+4m] BUSCO lineage dataset already exists, skipping ~~~
~~~ [+4m] Running BUSCO for Fungus using lineage: fungi_odb10 ~~~
~~~ [+25m] Skipping blobtools for Fungus, not desired ~~~
~~~ [+25m] Summarizing Assembly for Fungus ~~~
~~~ [+25m] Your final assembly for Fungus is: /project/coffea_pangenome/Artocarpus/puzzler_trials/primary_asm/Fungus.fa ~~~
~~~ [+25m] Your final assembly stats for Fungus are in: /project/coffea_pangenome/Artocarpus/puzzler_trials/primary_asm/stats/Fungus.summary.txt ~~~
```

Which only ran for about 25 minutes. 

The script skips all the previously run steps, because we have the `.complete` files in `$WD/Fungus`, as well as any already-completed step-specific files (e.g. `$WD/juicer_files/$SAMPLE_JBAT.review.assembly` for post-curation assembly and renaming). 

```
drwxr-s---. 2  4096 Jun 26 15:21 01_hifiasm
drwxr-s---. 2  4096 Jun 26 15:47 02_purge_dups
drwxr-s---. 2  4096 Jun 26 15:47 03_haphic
drwxr-s---. 2  4096 Jun 26 15:48 04_juicer
drwxr-s---. 2  4096 Jun 26 16:04 05_postjuicebox
drwxr-s---. 2  4096 Jun 26 16:07 06_realign_hic_hifi
drwxr-s---. 2  4096 Jun 26 16:19 07_busco_yak_blob
-rw-r-----. 1     0 Jun 26 15:23 align_hic.complete
-rw-r-----. 1     0 Jun 26 16:20 busco.complete
-rw-r-----. 1     0 Jun 26 15:21 hifiasm.complete
-rw-r-----. 1     0 Jun 26 15:48 juicer.complete
-rw-r-----. 1     0 Jun 26 15:22 purge_dups.complete
-rw-r-----. 1     0 Jun 26 16:06 qc_align_hic.complete
-rw-r-----. 1     0 Jun 26 15:47 scaffolding.complete
-rw-r-----. 1     0 Jun 26 16:08 yak.complete
```

Which outputs a final HiC contact map at `$WD/primary_asm/stats/$SAMPLE.pdf`: 

![final hic contacts](/examples/figs/Fungus_Final.png)

And the final assembly stats here:

```
cat primary_asm/stats/Fungus.summary.txt
Sample  SizeBP  Sequences       Contigs Gaps    ContigN50       ScafN50 YAK_CV  YAK_QV  WithinChrsBP    PropWithinChrs  GapsInChrs      BUSCO_Complete  BUSCO_singlecopy
Fungus   24787743        22      23     1        1776619         1806824        0.999   72.007   23726189       0.9572  1       99.1    99.1
```
___

## Juicer with No Reference

This time, since I know that HapHiC won't succeed, I will skip straight ahead to using YAHS with `--yahs`.

Submit: `puzzler -s Fungus_NoRef -m samples.tsv --threads 8 --mem 32 --yahs` 

```
=======================================================================
__________ ____ _______________________.____     _____________________ 
\______   \    |   \____    /\____    /|    |    \_   _____/\______   \
 |     ___/    |   / /     /   /     / |    |     |    __)_  |       _/
 |    |   |    |  / /     /_  /     /_ |    |___  |        \ |    |   \
 |____|   |______/ /_______ \/_______ \|_______ \/_______  / |____|_  /
                           \/        \/        \/        \/         \/ 
=======================================================================
                         ~~~ v1.9.2 ~~~

=======================================================================
Parameters for sample: Fungus_NoRef 
RUNTIME: conda
CONTAINER: /project/coffea_pangenome/Artocarpus/puzzler_trials/NA (ignore if using conda)
WD: /project/coffea_pangenome/Artocarpus/puzzler_trials 
HIFI: /project/coffea_pangenome/Artocarpus/puzzler_trials/Fungus.HiFi.fastq.gz
HIC_R1: /project/coffea_pangenome/Artocarpus/puzzler_trials/Fungus.HiC.R1.fastq.gz
HIC_R2: /project/coffea_pangenome/Artocarpus/puzzler_trials/Fungus.HiC.R2.fastq.gz
NUMBER CHRS: 14
REFERENCE: NA
HOM_COV: NA
BLOB_DB: NA
BUSCO_LINEAGE: fungi_odb10
BUSCO_DB: /project/coffea_pangenome/Artocarpus/puzzler_trials/busco_downloads
Cores Requested: 8
Cores Available: 8
RAM Requested: 32
Memory Available: 2025.1 GB
PUZZLER command:  (blank if using conda)
=======================================================================

~~~ [+0m] Checking software availability, this will take about 1 minute ~~~
~~~ [+0m] Software check complete ~~~
~~~ [+0m] Assembling genome for Fungus_NoRef ~~~
~~~ [+0m] Running juicer, script will stop after .hic files created ~~~
~~~ [+0m] Will run purge_dups step ~~~
~~~ [+0m] Skipping HapHic and scaffolding with YAHS directly ~~~
~~~ [+0m] Starting hifiasm assembly for Fungus_NoRef ~~~
~~~ [+9m] Starting Purge_Dups for Fungus_NoRef ~~~
~~~ [+9m] Mapping HiC reads to Fungus_NoRef draft ~~~
~~~ [+13m] Running YAHS directly on Fungus_NoRef: skipping HapHiC ~~~
~~~ [+13m] Creating .hic file for juicebox for Fungus_NoRef ~~~
~~~ [+14m] Post curation assembly file missing for Fungus_NoRef: Run Juicebox & place in /project/coffea_pangenome/Artocarpus/puzzler_trials/juicer_files/Fungus_NoRef_JBAT.review.assembly ~~~
```

Ran in 14 minutes. 

Perform juicebox edits....

To watch the video, please download the [Fungus_NoRef.mp4](https://zenodo.org/records/15693026/files/Fungus_NoRef.mp4?download=1) file from the manuscript's [Zenodo link](https://doi.org/10.5281/zenodo.15693025). 

**Save this file in this location:** `$WD/juicer_files/Fungus_NoRef_JBAT.review.assembly`, simply re-run:

```
=======================================================================
__________ ____ _______________________.____     _____________________ 
\______   \    |   \____    /\____    /|    |    \_   _____/\______   \
 |     ___/    |   / /     /   /     / |    |     |    __)_  |       _/
 |    |   |    |  / /     /_  /     /_ |    |___  |        \ |    |   \
 |____|   |______/ /_______ \/_______ \|_______ \/_______  / |____|_  /
                           \/        \/        \/        \/         \/ 
=======================================================================
                         ~~~ v1.9.2 ~~~

=======================================================================
Parameters for sample: Fungus_NoRef 
RUNTIME: conda
CONTAINER: /project/coffea_pangenome/Artocarpus/puzzler_trials/NA (ignore if using conda)
WD: /project/coffea_pangenome/Artocarpus/puzzler_trials 
HIFI: /project/coffea_pangenome/Artocarpus/puzzler_trials/Fungus.HiFi.fastq.gz
HIC_R1: /project/coffea_pangenome/Artocarpus/puzzler_trials/Fungus.HiC.R1.fastq.gz
HIC_R2: /project/coffea_pangenome/Artocarpus/puzzler_trials/Fungus.HiC.R2.fastq.gz
NUMBER CHRS: 14
REFERENCE: NA
HOM_COV: NA
BLOB_DB: NA
BUSCO_LINEAGE: fungi_odb10
BUSCO_DB: /project/coffea_pangenome/Artocarpus/puzzler_trials/busco_downloads
Cores Requested: 8
Cores Available: 8
RAM Requested: 48
Memory Available: 601.0 GB
PUZZLER command:  (blank if using conda)
=======================================================================

~~~ [+0m] Checking software availability, this will take about 1 minute ~~~
~~~ [+1m] Software check complete ~~~
~~~ [+1m] Assembling genome for Fungus_NoRef ~~~
~~~ [+1m] Running juicer, script will stop after .hic files created ~~~
~~~ [+1m] Running purge_dups step ~~~
~~~ [+1m] Skipping HapHic and scaffolding with YAHS directly ~~~
~~~ [+1m] Skipping hifiasm for Fungus_NoRef: /project/coffea_pangenome/Artocarpus/puzzler_trials/Fungus_NoRef/01_hifiasm/asm.hic.p_ctg.gfa exists ~~~
~~~ [+1m] Skipping purge for Fungus_NoRef: /project/coffea_pangenome/Artocarpus/puzzler_trials/Fungus_NoRef/02_purge_dups/p_ctg.purged.fa exists ~~~
~~~ [+1m] Skipping HiC alignment for Fungus_NoRef: /project/coffea_pangenome/Artocarpus/puzzler_trials/Fungus_NoRef/03_haphic/filtered.bam exists ~~~
~~~ [+1m] Skipping HapHiC for Fungus_NoRef: /project/coffea_pangenome/Artocarpus/puzzler_trials/Fungus_NoRef/03_haphic/haphic/04.build/scaffolds.fa exists ~~~
~~~ [+1m] Creating .hic file for juicebox for Fungus_NoRef ~~~
~~~ [+1m] Skipping initial alignment between reference and draft, /project/coffea_pangenome/Artocarpus/puzzler_trials/Fungus_NoRef/04_juicer/alignment.agp exists ~~~]
~~~ [+2m] No reference provided for Fungus_NoRef: simply extracting assembly to: /project/coffea_pangenome/Artocarpus/puzzler_trials/primary_asm/Fungus_NoRef.fa ~~~
~~~ [+2m] No reference provided for Fungus_NoRef: no chromosome re-naming, so final assembly already: /project/coffea_pangenome/Artocarpus/puzzler_trials/primary_asm/Fungus_NoRef.fa ~~~
~~~ [+2m] Assessing genome quality for Fungus_NoRef ~~~
~~~ [+2m] Creating final HiC bam for Fungus_NoRef ~~~
~~~ [+5m] Creating final HiC contact map for Fungus_NoRef ~~~
~~~ [+5m] Creating final contact map for Fungus_NoRef, no ref specified, plotting scaffolds > 1mb ~~~
~~~ [+6m] Running YAK on Fungus_NoRef ~~~
~~~ [+6m] BUSCO lineage dataset already exists, skipping ~~~
~~~ [+6m] Running BUSCO for Fungus_NoRef using lineage: fungi_odb10 ~~~
~~~ [+32m] Skipping blobtools for Fungus_NoRef, not desired ~~~
~~~ [+32m] Summarizing Assembly for Fungus_NoRef ~~~
~~~ [+32m] Your final assembly for Fungus_NoRef is: /project/coffea_pangenome/Artocarpus/puzzler_trials/primary_asm/Fungus_NoRef.fa ~~~
~~~ [+32m] Your final assembly stats for Fungus_NoRef are in: /project/coffea_pangenome/Artocarpus/puzzler_trials/primary_asm/stats/Fungus_NoRef.summary.txt ~~~
```

With only 8 cores, the second step only took 32 minutes - most of that was BUSCO. 

___

## No Juicer but with Reference

Just add the `--no_juice` flag to skip juicer. Not recommended! 

I ran this with 64 GB RAM because BUSCO failed once. 

```
=======================================================================
__________ ____ _______________________.____     _____________________ 
\______   \    |   \____    /\____    /|    |    \_   _____/\______   \
 |     ___/    |   / /     /   /     / |    |     |    __)_  |       _/
 |    |   |    |  / /     /_  /     /_ |    |___  |        \ |    |   \
 |____|   |______/ /_______ \/_______ \|_______ \/_______  / |____|_  /
                           \/        \/        \/        \/         \/ 
=======================================================================
                         ~~~ v1.9.2 ~~~

=======================================================================
Parameters for sample: Fungus_NoJuice 
RUNTIME: conda
CONTAINER: /project/coffea_pangenome/Artocarpus/puzzler_trials/NA (ignore if using conda)
WD: /project/coffea_pangenome/Artocarpus/puzzler_trials 
HIFI: /project/coffea_pangenome/Artocarpus/puzzler_trials/Fungus.HiFi.fastq.gz
HIC_R1: /project/coffea_pangenome/Artocarpus/puzzler_trials/Fungus.HiC.R1.fastq.gz
HIC_R2: /project/coffea_pangenome/Artocarpus/puzzler_trials/Fungus.HiC.R2.fastq.gz
NUMBER CHRS: 14
REFERENCE: /project/coffea_pangenome/Artocarpus/puzzler_trials/GCA_964106605.1_gdRhiKalk1.hap1.1_genomic.fna
HOM_COV: NA
BLOB_DB: NA
BUSCO_LINEAGE: fungi_odb10
BUSCO_DB: /project/coffea_pangenome/Artocarpus/puzzler_trials/busco_downloads
Cores Requested: 8
Cores Available: 8
RAM Requested: 64
Memory Available: 1888.0 GB
PUZZLER command:  (blank if using conda)
=======================================================================

~~~ [+0m] Checking software availability, this will take about 1 minute ~~~
~~~ [+1m] Software check complete ~~~
~~~ [+1m] Assembling genome for Fungus_NoJuice ~~~
~~~ [+1m] Skipping juicer, no manual curation (not recommended!) ~~~
~~~ [+1m] Running purge_dups step ~~~
~~~ [+1m] Skipping HapHic and scaffolding with YAHS directly ~~~
~~~ [+1m] Starting hifiasm assembly for Fungus_NoJuice ~~~
~~~ [+10m] Starting Purge_Dups for Fungus_NoJuice ~~~
~~~ [+10m] Mapping HiC reads to Fungus_NoJuice draft ~~~
~~~ [+13m] Running YAHS directly on Fungus_NoJuice: skipping HapHiC ~~~
~~~ [+14m] Skipping juicer HiC file creation for Fungus_NoJuice: --no_juice given, not recommended! ~~~
~~~ [+14m] Skipping juicer extraction for Fungus_NoJuice, --no_juice requested! (not recommended!) ~~~
~~~ [+14m] Renaming chromosomes for Fungus_NoJuice ~~~
~~~ [+14m] Scaffold sanity check passed for renaming, proceeding! ~~~
~~~ [+14m] Multiple scaffolds corresponding to single Chr for Fungus_NoJuice, Renaming them e.g. Chr1_unloc1, Chr1_unloc2.. ~~~
~~~ [+14m] Assessing genome quality for Fungus_NoJuice ~~~
~~~ [+14m] Creating final HiC bam for Fungus_NoJuice ~~~
~~~ [+18m] Creating final HiC contact map for Fungus_NoJuice ~~~
~~~ [+18m] Creating final contact map for Fungus_NoJuice, plotting named chromosomes ~~~
~~~ [+21m] Running YAK on Fungus_NoJuice ~~~
~~~ [+22m] BUSCO lineage dataset already exists, skipping ~~~
~~~ [+22m] Running BUSCO for Fungus_NoJuice using lineage: fungi_odb10 ~~~
~~~ [+42m] Skipping blobtools for Fungus_NoJuice, not desired ~~~
~~~ [+42m] Summarizing Assembly for Fungus_NoJuice ~~~
~~~ [+43m] Your final assembly for Fungus_NoJuice is: /project/coffea_pangenome/Artocarpus/puzzler_trials/primary_asm/Fungus_NoJuice.fa ~~~
~~~ [+43m] Your final assembly stats for Fungus_NoJuice are in: /project/coffea_pangenome/Artocarpus/puzzler_trials/primary_asm/stats/Fungus_NoJuice.summary.txt ~~~
```

Completed in 43 minutes. 

___

## No Juicer and No Reference

Just add the `--no_juice` flag to skip juicer and call the `--sample/-s` with "NA" in the reference column. Not recommended! 

```
puzzler -s Fungus_NoRefNoJuice -m samples.tsv --threads 8 --mem 48 --yahs --no_juice

=======================================================================
__________ ____ _______________________.____     _____________________ 
\______   \    |   \____    /\____    /|    |    \_   _____/\______   \
 |     ___/    |   / /     /   /     / |    |     |    __)_  |       _/
 |    |   |    |  / /     /_  /     /_ |    |___  |        \ |    |   \
 |____|   |______/ /_______ \/_______ \|_______ \/_______  / |____|_  /
                           \/        \/        \/        \/         \/ 
=======================================================================
                         ~~~ v1.9.2 ~~~

=======================================================================
Parameters for sample: Fungus_NoRefNoJuice 
RUNTIME: conda
CONTAINER: /project/coffea_pangenome/Artocarpus/puzzler_trials/NA (ignore if using conda)
WD: /project/coffea_pangenome/Artocarpus/puzzler_trials 
HIFI: /project/coffea_pangenome/Artocarpus/puzzler_trials/Fungus.HiFi.fastq.gz
HIC_R1: /project/coffea_pangenome/Artocarpus/puzzler_trials/Fungus.HiC.R1.fastq.gz
HIC_R2: /project/coffea_pangenome/Artocarpus/puzzler_trials/Fungus.HiC.R2.fastq.gz
NUMBER CHRS: 14
REFERENCE: NA
HOM_COV: NA
BLOB_DB: NA
BUSCO_LINEAGE: fungi_odb10
BUSCO_DB: /project/coffea_pangenome/Artocarpus/puzzler_trials/busco_downloads
Cores Requested: 8
Cores Available: 8
RAM Requested: 48
Memory Available: 739.2 GB
PUZZLER command:  (blank if using conda)
=======================================================================

~~~ [+0m] Checking software availability, this will take about 1 minute ~~~
~~~ [+0m] Software check complete ~~~
~~~ [+0m] Assembling genome for Fungus_NoRefNoJuice ~~~
~~~ [+0m] Skipping juicer, no manual curation (not recommended!) ~~~
~~~ [+0m] Running purge_dups step ~~~
~~~ [+0m] Skipping HapHic and scaffolding with YAHS directly ~~~
~~~ [+0m] Starting hifiasm assembly for Fungus_NoRefNoJuice ~~~
~~~ [+9m] Starting Purge_Dups for Fungus_NoRefNoJuice ~~~
~~~ [+9m] Mapping HiC reads to Fungus_NoRefNoJuice draft ~~~
~~~ [+11m] Running YAHS directly on Fungus_NoRefNoJuice: skipping HapHiC ~~~
~~~ [+11m] Skipping juicer HiC file creation for Fungus_NoRefNoJuice: --no_juice given, not recommended! ~~~
~~~ [+11m] No reference provided for Fungus_NoRefNoJuice: simply extracting assembly to: /project/coffea_pangenome/Artocarpus/puzzler_trials/primary_asm/Fungus_NoRefNoJuice.fa ~~~
~~~ [+11m] Skipping juicer extraction for Fungus_NoRefNoJuice, --no_juice requested! (not recommended!) ~~~
~~~ [+11m] No reference provided for Fungus_NoRefNoJuice: no chromosome re-naming, so final assembly already: /project/coffea_pangenome/Artocarpus/puzzler_trials/primary_asm/Fungus_NoRefNoJuice.fa ~~~
~~~ [+11m] Assessing genome quality for Fungus_NoRefNoJuice ~~~
~~~ [+11m] Creating final HiC bam for Fungus_NoRefNoJuice ~~~
~~~ [+13m] Creating final HiC contact map for Fungus_NoRefNoJuice ~~~
~~~ [+13m] Creating final contact map for Fungus_NoRefNoJuice, no ref specified, plotting scaffolds > 1mb ~~~
~~~ [+14m] Running YAK on Fungus_NoRefNoJuice ~~~
~~~ [+15m] BUSCO lineage dataset already exists, skipping ~~~
~~~ [+15m] Running BUSCO for Fungus_NoRefNoJuice using lineage: fungi_odb10 ~~~
~~~ [+35m] Skipping blobtools for Fungus_NoRefNoJuice, not desired ~~~
~~~ [+35m] Summarizing Assembly for Fungus_NoRefNoJuice ~~~
~~~ [+35m] Your final assembly for Fungus_NoRefNoJuice is: /project/coffea_pangenome/Artocarpus/puzzler_trials/primary_asm/Fungus_NoRefNoJuice.fa ~~~
~~~ [+35m] Your final assembly stats for Fungus_NoRefNoJuice are in: /project/coffea_pangenome/Artocarpus/puzzler_trials/primary_asm/stats/Fungus_NoRefNoJuice.summary.txt ~~~
```

This fork ran in 35 minutes including BUSCO.

___

## Outputs

Within `$WD/primary_asm/stats/`, we have our summary statistics.

The samples which have a `$REFERENCE` have 2 additional columns: `WithinChrsBp` and `PropWithinChrs`, showing how much of the assembly is within named chromosomes:

```bash
awk 'NR == 1 || FNR > 1' Fungus.summary.txt Fungus_NoJuice.summary.txt 
Sample  SizeBP  Sequences       Contigs Gaps    ContigN50       ScafN50 YAK_CV  YAK_QV  WithinChrsBP    PropWithinChrs  GapsInChrs      BUSCO_Complete  BUSCO_singlecopy
Fungus   24787743        22      23     1        1776619         1806824        0.999   72.007   23726189       0.9572  1       99.1    99.1
Fungus_NoJuice   24787643        23      23     0        1776619         1776619        0.999   72.007   23548816       0.9500  0       99.1    99.1
```

| Sample         | SizeBP   | Sequences | Contigs | Gaps | ContigN50 | ScafN50 | YAK_CV | YAK_QV | WithinChrsBP | PropWithinChrs | GapsInChrs | BUSCO_Complete | BUSCO_singlecopy |
| -------------- | -------- | --------- | ------- | ---- | --------- | ------- | ------ | ------ | ------------ | -------------- | ---------- | -------------- | ---------------- |
| Fungus         | 24787743 | 22        | 23      | 1    | 1776619   | 1806824 | 0.999  | 72.007 | 23726189     | 0.9572         | 1          | 99.1           | 99.1             |
| Fungus_NoJuice | 24787643 | 23        | 23      | 0    | 1776619   | 1776619 | 0.999  | 72.007 | 23548816     | 0.95           | 0          | 99.1           | 99.1             |

While the assemblies with no `$REFERENCE` do not have those columns:

```bash
awk 'NR == 1 || FNR > 1' Fungus_NoRef.summary.txt Fungus_NoRefNoJuice.summary.txt
Sample  SizeBP  Sequences       Contigs Gaps    ContigN50       ScafN50 YAK_CV  YAK_QV  BUSCO_Complete  BUSCO_singlecopy
Fungus_NoRef     24787743        22      23     1        1776619         1806824        0.999   72.007  99.5    99.5
Fungus_NoRefNoJuice      24787643        23      23     0        1776619         1776619        0.999   72.007  99.5    99.5
```

| Sample              | SizeBP   | Sequences | Contigs | Gaps | ContigN50 | ScafN50 | YAK_CV | YAK_QV | BUSCO_Complete | BUSCO_singlecopy |
| ------------------- | -------- | --------- | ------- | ---- | --------- | ------- | ------ | ------ | -------------- | ---------------- |
| Fungus_NoRef        | 24787743 | 22        | 23      | 1    | 1776619   | 1806824 | 0.999  | 72.007 | 99.5           | 99.5             |
| Fungus_NoRefNoJuice | 24787643 | 23        | 23      | 0    | 1776619   | 1776619 | 0.999  | 72.007 | 99.5           | 99.5             |

Overall, all assemblies look good. 

![fungus_assemblies](/examples/figs/Fungus_Assemblies_HiC.png)

## Alignments

I will perform a quick WGA and create a dotplot between the assemblies. Since the no `$REFERENCE` assemblies do not have named chromosomes, I will just align any sequences > 200 Kb:

```bash
for i in Fungus Fungus_NoRef Fungus_RefNoJuice Fungus_NoRefNoJuice RhiKalk1.hap1; do 

    echo "Extracting chromosomes for ${i}"
    samtools faidx ${i}.fa
    #egrep 'chr|Chr' ${i}.fa.fai | egrep -v 'SUPER|unloc' | awk '{print $1}' > $i.tmp 
    awk '$2 > 2e5' $i.fa.fai | awk '{print $1}' > $i.tmp
    seqtk subseq ${i}.fa $i.tmp > ../../alignments/fastas/${i}.fa

done
```

I will compare each assembly to the published reference:

```bash
for i in Fungus Fungus_NoRef Fungus_RefNoJuice Fungus_NoRefNoJuice; do 

    echo "Alignment for ${i}"
    minimap2 -t ${t} -x asm20 --secondary=no -c -D --max-chain-skip 100 --max-chain-iter 1000 --frag yes ${i}.fa RhiKalk1.hap1.fa > fungus_pafs/${i}.paf
    paf2dotplot.R fungus_pafs/${i}.paf -r 1e4 -m 1e4 

done 
```

Alignments:

![fungus_alignments](/examples/figs/Fungus_Assemblies_Alignments.png)