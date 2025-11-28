# Genome Assembly for N=12 species with **`puzzler`**

I will create Puzzler assemblies for 12 species which have HiC and HiFi data from the NCBI:

| Sample      | Species                     | Species Name                 | Bioproject                                                   | Link                                                         |
| ----------- | --------------------------- | ---------------------------- | ------------------------------------------------------------ | ------------------------------------------------------------ |
| Rhodie      | Rhododendron                | Rhododendron vialii          | PRJNA971245                                                  | https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_030253575.1/ |
| Strawberry  | Strawberry                  | Fragaria x ananassa          | PRJNA1216395                                                 | https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_048987885.1/ |
| Fungus      | Rhizosphaera needle cast    | Rhizosphaera kalkhoffii      | PRJEB76005                                                   | https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_964106605.1/ |
| Fly         | Mexican Fruit Fly           | Anastrepha ludens            | PRJNA803324                                                  | https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_028408465.1/ |
| Squirt      | Solitary Sea Squirt         | Ascidia mentula              | PRJEB58134                                                   | https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_947561715.1/ |
| Fish        | Antarctic lanternfish       | Electrona antarctica         | PRJEB61834                                                   | https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_951216825.1/ |
| Stickleback | Ninespine Stickleback       | Pungitius pungitius          | PRJEB60134                                                   | https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_949316345.1/ |
| Toad        | Eastern Narrow-mouthed Toad | Gastrophryne carolinensis    | PRJNA923362                                                  | https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_027917425.1/ |
| Frog        | Guiana Rocket Frog          | Anomaloglossus baeobatrachus | PRJNA1226914                                                 | https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_048569485.1/ |
| Crane       | Eurasian Crane              | Grus grus                    | [PRJEB75999](https://www.ncbi.nlm.nih.gov/bioproject/PRJEB75999/) | https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_964106855.1/ |
| Beaver      | North American Beaver       | Castor canadensis            | PRJNA1197684                                                 | https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_047511655.1/ |
| Whale       | True's beaked whale         | Mesoplodon mirus             | PRJEB79280                                                   | https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_964341445.1/ |

# Acquire Data:

## Reference Fastas

Download the completed platinum assembly, and for simple chromosome renaming, strip all information before and after the chrID:

```bash
# Beaver
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/047/511/655/GCF_047511655.1_mCasCan1.hap1v2/GCF_047511655.1_mCasCan1.hap1v2_genomic.fna.gz
gunzip GCF_047511655.1_mCasCan1.hap1v2_genomic.fna.gz
sed -i 's/>.*chromosome />chr/g' GCF_047511655.1_mCasCan1.hap1v2_genomic.fna
sed -i 's/,.*//g' GCF_047511655.1_mCasCan1.hap1v2_genomic.fna

# Crane 
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/964/106/855/GCA_964106855.1_bGruGru1.hap1.1/GCA_964106855.1_bGruGru1.hap1.1_genomic.fna.gz
gunzip GCA_964106855.1_bGruGru1.hap1.1_genomic.fna.gz
sed -i 's/>.*chromosome: />chr/g' GCA_964106855.1_bGruGru1.hap1.1_genomic.fna

# Fish (Lanternfish)
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/951/216/825/GCA_951216825.1_fEleAnt2.1/GCA_951216825.1_fEleAnt2.1_genomic.fna.gz
gunzip GCA_951216825.1_fEleAnt2.1_genomic.fna.gz
sed -i 's/>.*chromosome: />chr/g' GCA_951216825.1_fEleAnt2.1_genomic.fna

# Fly
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/028/408/465/GCF_028408465.1_idAnaLude1.1/GCF_028408465.1_idAnaLude1.1_genomic.fna.gz
gunzip GCF_028408465.1_idAnaLude1.1_genomic.fna.gz
sed -i 's/>.*chromosome />chr/g' GCF_028408465.1_idAnaLude1.1_genomic.fna
sed -i 's/,.*//g' GCF_028408465.1_idAnaLude1.1_genomic.fna

# Frog
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/048/569/485/GCA_048569485.1_aAnoBae1.hap1/GCA_048569485.1_aAnoBae1.hap1_genomic.fna.gz
gunzip GCA_048569485.1_aAnoBae1.hap1_genomic.fna.gz
sed -i 's/>.*chromosome />chr/g' GCA_048569485.1_aAnoBae1.hap1_genomic.fna
sed -i 's/,.*//g' GCA_048569485.1_aAnoBae1.hap1_genomic.fna
sed -i 's/ /_/g' GCA_048569485.1_aAnoBae1.hap1_genomic.fna

# Rhodie (Rhododendron)
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/030/253/575/GCF_030253575.1_ASM3025357v1/GCF_030253575.1_ASM3025357v1_genomic.fna.gz
gunzip GCF_030253575.1_ASM3025357v1_genomic.fna.gz
sed -i 's/>.*chromosome />chr/g' GCF_030253575.1_ASM3025357v1_genomic.fna
sed -i 's/,.*//g' GCF_030253575.1_ASM3025357v1_genomic.fna

# Squirt
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/947/561/715/GCA_947561715.1_kaAscMent1.1/GCA_947561715.1_kaAscMent1.1_genomic.fna.gz
gunzip GCA_947561715.1_kaAscMent1.1_genomic.fna.gz
sed -i 's/>.*chromosome: />chr/g' GCA_947561715.1_kaAscMent1.1_genomic.fna

# Stickleback
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/949/316/345/GCF_949316345.1_fPunPun2.1/GCF_949316345.1_fPunPun2.1_genomic.fna.gz
gunzip GCF_949316345.1_fPunPun2.1_genomic.fna.gz
sed -i 's/>.*chromosome />chr/g' GCF_949316345.1_fPunPun2.1_genomic.fna
sed -i 's/,.*//g' GCF_949316345.1_fPunPun2.1_genomic.fna
sed -i 's/ /_/g' GCF_949316345.1_fPunPun2.1_genomic.fna

# Toad
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/027/917/425/GCA_027917425.1_aGasCar1.pri/GCA_027917425.1_aGasCar1.pri_genomic.fna.gz
gunzip GCA_027917425.1_aGasCar1.pri_genomic.fna.gz
sed -i 's/>.*chromosome />chr/g' GCA_027917425.1_aGasCar1.pri_genomic.fna
sed -i 's/,.*//g' GCA_027917425.1_aGasCar1.pri_genomic.fna

# Whale
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/964/341/445/GCA_964341445.1_mMesMir1.hap1.1/GCA_964341445.1_mMesMir1.hap1.1_genomic.fna.gz
gunzip GCA_964341445.1_mMesMir1.hap1.1_genomic.fna.gz
sed -i 's/>.*chromosome: />chr/g' GCA_964341445.1_mMesMir1.hap1.1_genomic.fna

# Fungus (Rhizosphaera)
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/964/106/605/GCA_964106605.1_gdRhiKalk1.hap1.1/GCA_964106605.1_gdRhiKalk1.hap1.1_genomic.fna.gz
gunzip GCA_964106605.1_gdRhiKalk1.hap1.1_genomic.fna.gz
sed -i 's/>.*chromosome: />chr/g' GCA_964106605.1_gdRhiKalk1.hap1.1_genomic.fna

# Strawberry
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/048/987/885/GCA_048987885.1_NARO_FxaKMR_1.0/GCA_048987885.1_NARO_FxaKMR_1.0_genomic.fna.gz
gunzip GCA_048987885.1_NARO_FxaKMR_1.0_genomic.fna.gz
sed -i 's/>.*chromosome />chr/g' GCA_048987885.1_NARO_FxaKMR_1.0_genomic.fna
sed -i 's/,.*//g' GCA_048987885.1_NARO_FxaKMR_1.0_genomic.fna
```

**Note that this is preference only**. Otherwise, your chromosomes would be named e.g.  `OZ066564.1 Rhizosphaera kalkhoffii genome assembly, chromosome: 1` 

```
grep '>' GCA_964106605.1_gdRhiKalk1.hap1.1_genomic.fna
>OZ066564.1 Rhizosphaera kalkhoffii genome assembly, chromosome: 1
>OZ066565.1 Rhizosphaera kalkhoffii genome assembly, chromosome: 2
>OZ066566.1 Rhizosphaera kalkhoffii genome assembly, chromosome: 3
>OZ066567.1 Rhizosphaera kalkhoffii genome assembly, chromosome: 4
>OZ066568.1 Rhizosphaera kalkhoffii genome assembly, chromosome: 5
>OZ066569.1 Rhizosphaera kalkhoffii genome assembly, chromosome: 6
>OZ066570.1 Rhizosphaera kalkhoffii genome assembly, chromosome: 7
>OZ066571.1 Rhizosphaera kalkhoffii genome assembly, chromosome: 8
>OZ066572.1 Rhizosphaera kalkhoffii genome assembly, chromosome: 9
>OZ066573.1 Rhizosphaera kalkhoffii genome assembly, chromosome: 10
>OZ066574.1 Rhizosphaera kalkhoffii genome assembly, chromosome: 11
>OZ066575.1 Rhizosphaera kalkhoffii genome assembly, chromosome: 12
>OZ066576.1 Rhizosphaera kalkhoffii genome assembly, chromosome: 13
>OZ066577.1 Rhizosphaera kalkhoffii genome assembly, chromosome: 14
>OZ066578.1 Rhizosphaera kalkhoffii genome assembly, organelle: mitochondrion
```

After stripping:

```
sed -i 's/>.*chromosome: />chr/g' GCA_964106605.1_gdRhiKalk1.hap1.1_genomic.fna

grep '>' GCA_964106605.1_gdRhiKalk1.hap1.1_genomic.fna
>chr1
>chr2
>chr3
>chr4
>chr5
>chr6
>chr7
>chr8
>chr9
>chr10
>chr11
>chr12
>chr13
>chr14
>OZ066578.1 Rhizosphaera kalkhoffii genome assembly, organelle: mitochondrion
```



## Reads

Beaver

```bash
#!/bin/bash

#SBATCH --time=4-00:00:00   
#SBATCH --nodes=1  
#SBATCH --ntasks-per-node=1
#SBATCH --partition=ceres

WD=/90daydata/coffea_pangenome/puzzler_trials

mkdir Castor_canadensis
cd Castor_canadensis
wget https://genomeark.s3.amazonaws.com/species/Castor_canadensis/mCasCan1/genomic_data/pacbio_hifi/m84091_240920_000814_s4.hifi_reads.bc1018.fastq.gz
wget https://genomeark.s3.amazonaws.com/species/Castor_canadensis/mCasCan1/genomic_data/pacbio_hifi/m84091_240927_164339_s1.hifi_reads.bc1018.fastq.gz
wget https://genomeark.s3.amazonaws.com/species/Castor_canadensis/mCasCan1/genomic_data/arima/mCasCan1_L1_S5_R1_001.fastq.gz
wget https://genomeark.s3.amazonaws.com/species/Castor_canadensis/mCasCan1/genomic_data/arima/mCasCan1_L1_S5_R2_001.fastq.gz
```


Crane

```bash
#!/bin/bash

#SBATCH --time=4-00:00:00   
#SBATCH --nodes=1  
#SBATCH --ntasks-per-node=1
#SBATCH --partition=ceres

WD=/90daydata/coffea_pangenome/puzzler_trials

mkdir Grus_grus
cd Grus_grus
wget https://genomeark.s3.amazonaws.com/species/Grus_grus/bGruGru1/genomic_data/pacbio_hifi/m84093_240219_121437_s3.ccs.bc2028.bam
wget "https://genomeark.s3.amazonaws.com/species/Grus_grus/bGruGru1/genomic_data/arima/bGruGru1_48542_3-4%237.cram"

pbindex m84093_240219_121437_s3.ccs.bc2028.bam
bam2fastq m84093_240219_121437_s3.ccs.bc2028.bam -o crane
samtools fastq -1 crane_R1.fastq.gz -2 crane_R2.fastq.gz -s crane_unpaired.fastq.gz -0 crane_other.fastq.gz "bGruGru1_48542_3-4%237.cram"
```

Fish (Lanternfish)

```bash
#!/bin/bash

#SBATCH --time=4-00:00:00   
#SBATCH --nodes=1  
#SBATCH --ntasks-per-node=1
#SBATCH --partition=ceres

WD=/90daydata/coffea_pangenome/puzzler_trials

mkdir Electrona_antarctica
cd Electrona_antarctica
wget https://genomeark.s3.amazonaws.com/species/Electrona_antarctica/fEleAnt2/genomic_data/pacbio_hifi/m64016e_220218_220054.ccs.bc1018_BAK8B_OA--bc1018_BAK8B_OA.bam
wget https://genomeark.s3.amazonaws.com/species/Electrona_antarctica/fEleAnt2/genomic_data/pacbio_hifi/m64094e_230228_152759.ccs.bc1016_BAK8B_OA--bc1016_BAK8B_OA.bam
wget https://genomeark.s3.amazonaws.com/species/Electrona_antarctica/fEleAnt2/genomic_data/pacbio_hifi/m64097e_210809_145236.ccs.bc1012_BAK8A_OA--bc1012_BAK8A_OA.bam
wget https://genomeark.s3.amazonaws.com/species/Electrona_antarctica/fEleAnt2/genomic_data/pacbio_hifi/m64230e_210724_235140.ccs.bc1012_BAK8A_OA--bc1012_BAK8A_OA.bam
wget "https://genomeark.s3.amazonaws.com/species/Electrona_antarctica/fEleAnt3/genomic_data/arima/fEleAnt3_40063_1%233_R1.fastq.gz"
wget "https://genomeark.s3.amazonaws.com/species/Electrona_antarctica/fEleAnt3/genomic_data/arima/fEleAnt3_40063_1%233_R2.fastq.gz"

pbindex m64016e_220218_220054.ccs.bc1018_BAK8B_OA--bc1018_BAK8B_OA.bam
bam2fastq m64016e_220218_220054.ccs.bc1018_BAK8B_OA--bc1018_BAK8B_OA.bam -o fish1
pbindex m64094e_230228_152759.ccs.bc1016_BAK8B_OA--bc1016_BAK8B_OA.bam
bam2fastq m64094e_230228_152759.ccs.bc1016_BAK8B_OA--bc1016_BAK8B_OA.bam -o fish2
pbindex m64097e_210809_145236.ccs.bc1012_BAK8A_OA--bc1012_BAK8A_OA.bam
bam2fastq m64097e_210809_145236.ccs.bc1012_BAK8A_OA--bc1012_BAK8A_OA.bam -o fish3
pbindex m64230e_210724_235140.ccs.bc1012_BAK8A_OA--bc1012_BAK8A_OA.bam
bam2fastq m64230e_210724_235140.ccs.bc1012_BAK8A_OA--bc1012_BAK8A_OA.bam -o fish4
```

Fly

```bash
#!/bin/bash

#SBATCH --time=4-00:00:00   
#SBATCH --nodes=1  
#SBATCH --ntasks-per-node=1
#SBATCH --partition=ceres

WD=/90daydata/coffea_pangenome/puzzler_trials

mkdir Anastrepha_ludens
cd Anastrepha_ludens
module load sratoolkit/3.2.0
# FLY
#HIFI SRR30045406   HIC SRR18053126
prefetch SRR18053126 --max-size 0
fasterq-dump --split-files SRR18053126

prefetch SRR30045406 --max-size 0
fasterq-dump --split-files SRR30045406
```

Frog

```
# Frog
wget https://genomeark.s3.amazonaws.com/species/Anomaloglossus_baeobatrachus/aAnoBae1/genomic_data/pacbio_hifi/m54306Ue_220202_203517.hifi_reads.fastq.gz
wget https://genomeark.s3.amazonaws.com/species/Anomaloglossus_baeobatrachus/aAnoBae1/genomic_data/pacbio_hifi/m64055e_220214_185145.hifi_reads.fastq.gz
wget https://genomeark.s3.amazonaws.com/species/Anomaloglossus_baeobatrachus/aAnoBae1/genomic_data/pacbio_hifi/m64334e_220209_203922.hifi_reads.fastq.gz
wget https://genomeark.s3.amazonaws.com/species/Anomaloglossus_baeobatrachus/aAnoBae1/genomic_data/pacbio_hifi/m64334e_220211_073330.hifi_reads.fastq.gz
wget https://genomeark.s3.amazonaws.com/species/Anomaloglossus_baeobatrachus/aAnoBae1/genomic_data/pacbio_hifi/m64334e_220212_183416.hifi_reads.fastq.gz
wget https://genomeark.s3.amazonaws.com/species/Anomaloglossus_baeobatrachus/aAnoBae1/genomic_data/pacbio_hifi/m64334e_220214_041926.hifi_reads.fastq.gz
wget https://genomeark.s3.amazonaws.com/species/Anomaloglossus_baeobatrachus/aAnoBae1/genomic_data/arima/aAnoBae1_1.fastq.gz
wget https://genomeark.s3.amazonaws.com/species/Anomaloglossus_baeobatrachus/aAnoBae1/genomic_data/arima/aAnoBae1_2.fastq.gz

```


Rhodie (Rhododendron)

```bash
#!/bin/bash

#SBATCH --time=4-00:00:00   
#SBATCH --nodes=1  
#SBATCH --ntasks-per-node=1
#SBATCH --partition=ceres

WD=/90daydata/coffea_pangenome/puzzler_trials
module load sratoolkit/3.2.0

mkdir Rhododendron_vialii
cd Rhododendron_vialii

# RHODIE
#HIFI SRR24501949   HIC  SRR24501947
prefetch SRR24501949 --max-size 0
fasterq-dump --split-files SRR24501949

prefetch SRR24501947 --max-size 0
fasterq-dump --split-files SRR24501947
```

Sea squirt

```bash
#!/bin/bash

#SBATCH --time=4-00:00:00   
#SBATCH --nodes=1  
#SBATCH --ntasks-per-node=1
#SBATCH --partition=ceres

WD=/90daydata/coffea_pangenome/puzzler_trials

mkdir Ascidia_mentula
cd Ascidia_mentula
wget https://genomeark.s3.amazonaws.com/species/Ascidia_mentula/kaAscMent1/genomic_data/pacbio_hifi/m64221e_210705_191039.ccs.bc1015_BAK8B_OA--bc1015_BAK8B_OA.bam
wget https://genomeark.s3.amazonaws.com/species/Ascidia_mentula/kaAscMent1/genomic_data/arima/kaAscMent1_45817_1%232.cram
samtools fastq -1 squirt_R1.fastq.gz -2 squirt_R2.fastq.gz -s squirt_unpaired.fastq.gz -0 squit_other.fastq.gz kaAscMent1_45817_1%232.cram

pbindex m64221e_210705_191039.ccs.bc1015_BAK8B_OA--bc1015_BAK8B_OA.bam
bam2fastq m64221e_210705_191039.ccs.bc1015_BAK8B_OA--bc1015_BAK8B_OA.bam -o squirt
```

Stickleback

```bash
#!/bin/bash

#SBATCH --time=4-00:00:00   
#SBATCH --nodes=1  
#SBATCH --ntasks-per-node=1
#SBATCH --partition=ceres

WD=/90daydata/coffea_pangenome/puzzler_trials

module load miniconda
source activate assembly 

mkdir Pungitius_pungitius
cd Pungitius_pungitius
wget https://genomeark.s3.amazonaws.com/species/Pungitius_pungitius/fPunPun2/genomic_data/pacbio_hifi/m64178e_230122_111455.ccs.bc1008_BAK8A_OA--bc1008_BAK8A_OA.bam
wget https://genomeark.s3.amazonaws.com/species/Pungitius_pungitius/fPunPun2/genomic_data/pacbio_hifi/m64178e_230122_111455.ccs.bc1008_BAK8A_OA--bc1008_BAK8A_OA.bam.pbi

wget https://genomeark.s3.amazonaws.com/species/Pungitius_pungitius/fPunPun2/genomic_data/arima/fPunPun2_46625_3%233.cram

samtools fastq -1 stickleback_R1.fastq.gz -2 stickleback_R2.fastq.gz -s stickleback_unpaired.fastq.gz -0 stickleback_other.fastq.gz fPunPun2*cram
bam2fastq m64178e_230122_111455.ccs.bc1008_BAK8A_OA--bc1008_BAK8A_OA.bam -o stickleback
```

Toad

```bash
#!/bin/bash

#SBATCH --time=4-00:00:00   
#SBATCH --nodes=1  
#SBATCH --ntasks-per-node=1
#SBATCH --partition=ceres

WD=/90daydata/coffea_pangenome/puzzler_trials

mkdir Gastrophryne_carolinensis
cd Gastrophryne_carolinensis
wget https://genomeark.s3.amazonaws.com/species/Gastrophryne_carolinensis/aGasCar1/genomic_data/pacbio_hifi/m64055e_211104_144808.hifi_reads.fastq.gz
wget https://genomeark.s3.amazonaws.com/species/Gastrophryne_carolinensis/aGasCar1/genomic_data/pacbio_hifi/m64330e_211108_203923.hifi_reads.fastq.gz
wget https://genomeark.s3.amazonaws.com/species/Gastrophryne_carolinensis/aGasCar1/genomic_data/pacbio_hifi/m64330e_211110_060034.hifi_reads.fastq.gz
wget https://genomeark.s3.amazonaws.com/species/Gastrophryne_carolinensis/aGasCar1/genomic_data/pacbio_hifi/m64334e_211119_061646.hifi_reads.fastq.gz
wget https://genomeark.s3.amazonaws.com/species/Gastrophryne_carolinensis/aGasCar1/genomic_data/arima/aGasCar1_S6_L002_R1_001.fastq.gz
wget https://genomeark.s3.amazonaws.com/species/Gastrophryne_carolinensis/aGasCar1/genomic_data/arima/aGasCar1_S6_L002_R2_001.fastq.gz

cat m64055e_211104_144808.hifi_reads.fastq.gz m64330e_211108_203923.hifi_reads.fastq.gz m64330e_211110_060034.hifi_reads.fastq.gz m64334e_211119_061646.hifi_reads.fastq.gz > toad.HiFi.fastq.gz
```

Whale

```bash
#!/bin/bash

#SBATCH --time=4-00:00:00   
#SBATCH --nodes=1  
#SBATCH --ntasks-per-node=1
#SBATCH --partition=ceres

WD=/90daydata/coffea_pangenome/puzzler_trials

mkdir Mesoplodon_mirus
cd Mesoplodon_mirus
wget https://genomeark.s3.amazonaws.com/species/Mesoplodon_mirus/mMesMir1/genomic_data/pacbio_hifi/m84047_240717_113325_s2.ccs.bc2062.bam
wget https://genomeark.s3.amazonaws.com/species/Mesoplodon_mirus/mMesMir1/genomic_data/pacbio_hifi/m84093_240514_112417_s1.ccs.bc2001.bam
wget "https://genomeark.s3.amazonaws.com/species/Mesoplodon_mirus/mMesMir1/genomic_data/arima/mMesMir1_49110_3-4%232.cram"
samtools fastq -1 whale_R1.fastq.gz -2 whale_R2.fastq.gz -s whale_unpaired.fastq.gz -0 whale_other.fastq.gz "mMesMir1_49110_3-4%232.cram"

pbindex m84047_240717_113325_s2.ccs.bc2062.bam
bam2fastq m84047_240717_113325_s2.ccs.bc2062.bam -o whale1
pbindex m84093_240514_112417_s1.ccs.bc2001.bam
bam2fastq m84093_240514_112417_s1.ccs.bc2001.bam -o whale2
```

Fungus

```bash
#!/bin/bash

#SBATCH --time=4-00:00:00   
#SBATCH --nodes=1  
#SBATCH --ntasks-per-node=1
#SBATCH --partition=ceres

WD=/90daydata/coffea_pangenome/puzzler_trials

mkdir Rhizosphaera_kalkhoffii
cd Rhizosphaera_kalkhoffii
module load sratoolkit/3.2.0

# FUNGUS
#HIFI ERR12760830          HIC ERR12765203
prefetch ERR12760830
fasterq-dump --split-files ERR12760830
head -n 100000 ERR12760830.fastq | gzip -c > ../concat_reads/Fungus.HiFi.fastq.gz

prefetch --max-size 0 ERR12765203
fasterq-dump --split-files ERR12765203
head -n 1000000 ERR12765203_1.fastq | gzip -c > ../concat_reads/Fungus.HiC.R1.fastq.gz
head -n 1000000 ERR12765203_2.fastq | gzip -c > ../concat_reads/Fungus.HiC.R2.fastq.gz
```

[Strawberry](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_048987885.1/) 

```
#!/bin/bash

#SBATCH --time=4-00:00:00   
#SBATCH --nodes=1  
#SBATCH --ntasks-per-node=1
#SBATCH --partition=ceres

WD=/90daydata/coffea_pangenome/puzzler_trials

mkdir Strawberry
cd Strawberry
module load sratoolkit/3.2.0
# STRAWBERRY
#HIFI SRX27493331   HIC SRX27493332
prefetch SRX27493331 --max-size 0
fasterq-dump --split-files SRX27493331

prefetch SRX27493332 --max-size 0
fasterq-dump --split-files SRX27493332
```



# Create Drafts -> Juicer

I will skip blobtools because that will double to triple the assembly time. I recommend running this (just set `blob_db` to a directory, but expect long runtimes! 

Here is my `samples.tsv`:

| sample      | runtime   | container                                        | wd                                                      | hifi                                                         | hic_r1                                                       | hic_r2                                                       | num_chrs | reference                                                    | hom_cov | blob_database | busco_lineage        | busco_database                                             |
| ----------- | --------- | ------------------------------------------------ | ------------------------------------------------------- | ------------------------------------------------------------ | ------------------------------------------------------------ | ------------------------------------------------------------ | -------- | ------------------------------------------------------------ | ------- | ------------- | -------------------- | ---------------------------------------------------------- |
| Beaver      | apptainer | /home/justin.merondun/apptainer/puzzler_v1.8.sif | /90daydata/coffea_pangenome/puzzler_trials/modest_128gb | /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Beaver.HiFi.fastq.gz | /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Beaver.HiC.R1.fastq.gz | /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Beaver.HiC.R2.fastq.gz | 20       | /90daydata/coffea_pangenome/puzzler_trials/raw_data/references/GCF_047511655.1_mCasCan1.hap1v2_genomic.fna | NA      | NA            | mammalia_odb10       | /90daydata/coffea_pangenome/puzzler_trials/busco_downloads |
| Crane       | apptainer | /home/justin.merondun/apptainer/puzzler_v1.8.sif | /90daydata/coffea_pangenome/puzzler_trials/modest_128gb | /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Crane.HiFi.fastq.gz | /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Crane.HiC.R1.fastq.gz | /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Crane.HiC.R2.fastq.gz | 40       | /90daydata/coffea_pangenome/puzzler_trials/raw_data/references/GCA_964106855.1_bGruGru1.hap1.1_genomic.fna | NA      | NA            | aves_odb10           | /90daydata/coffea_pangenome/puzzler_trials/busco_downloads |
| Fish        | apptainer | /home/justin.merondun/apptainer/puzzler_v1.8.sif | /90daydata/coffea_pangenome/puzzler_trials/modest_128gb | /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Fish.HiFi.fastq.gz | /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Fish.HiC.R1.fastq.gz | /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Fish.HiC.R2.fastq.gz | 24       | /90daydata/coffea_pangenome/puzzler_trials/raw_data/references/GCA_951216825.1_fEleAnt2.1_genomic.fna | NA      | NA            | actinopterygii_odb10 | /90daydata/coffea_pangenome/puzzler_trials/busco_downloads |
| Fly         | apptainer | /home/justin.merondun/apptainer/puzzler_v1.8.sif | /90daydata/coffea_pangenome/puzzler_trials/modest_128gb | /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Fly.HiFi.fastq.gz | /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Fly.HiC.R1.fastq.gz | /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Fly.HiC.R2.fastq.gz | 7        | /90daydata/coffea_pangenome/puzzler_trials/raw_data/references/GCF_028408465.1_idAnaLude1.1_genomic.fna | NA      | NA            | diptera_odb10        | /90daydata/coffea_pangenome/puzzler_trials/busco_downloads |
| Frog        | apptainer | /home/justin.merondun/apptainer/puzzler_v1.8.sif | /90daydata/coffea_pangenome/puzzler_trials/modest_128gb | /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Frog.HiFi.fastq.gz | /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Frog.HiC.R1.fastq.gz | /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Frog.HiC.R2.fastq.gz | 12       | /90daydata/coffea_pangenome/puzzler_trials/raw_data/references/GCA_048569485.1_aAnoBae1.hap1_genomic.fna | NA      | NA            | tetrapoda_odb10      | /90daydata/coffea_pangenome/puzzler_trials/busco_downloads |
| Fungus      | apptainer | /home/justin.merondun/apptainer/puzzler_v1.8.sif | /90daydata/coffea_pangenome/puzzler_trials/modest_128gb | /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Fungus.HiFi.fastq.gz | /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Fungus.HiC.R1.fastq.gz | /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Fungus.HiC.R2.fastq.gz | 14       | /90daydata/coffea_pangenome/puzzler_trials/raw_data/references/GCA_964106605.1_gdRhiKalk1.hap1.1_genomic.fna | NA      | NA            | fungi_odb10          | /90daydata/coffea_pangenome/puzzler_trials/busco_downloads |
| Rhodie      | apptainer | /home/justin.merondun/apptainer/puzzler_v1.8.sif | /90daydata/coffea_pangenome/puzzler_trials/modest_128gb | /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Rhodie.HiFi.fastq.gz | /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Rhodie.HiC.R1.fastq.gz | /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Rhodie.HiC.R2.fastq.gz | 13       | /90daydata/coffea_pangenome/puzzler_trials/raw_data/references/GCF_030253575.1_ASM3025357v1_genomic.fna | NA      | NA            | embryophyta_odb10    | /90daydata/coffea_pangenome/puzzler_trials/busco_downloads |
| Squirt      | apptainer | /home/justin.merondun/apptainer/puzzler_v1.8.sif | /90daydata/coffea_pangenome/puzzler_trials/modest_128gb | /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Squirt.HiFi.fastq.gz | /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Squirt.HiC.R1.fastq.gz | /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Squirt.HiC.R2.fastq.gz | 9        | /90daydata/coffea_pangenome/puzzler_trials/raw_data/references/GCA_947561715.1_kaAscMent1.1_genomic.fna | NA      | NA            | metazoa_odb10        | /90daydata/coffea_pangenome/puzzler_trials/busco_downloads |
| Stickleback | apptainer | /home/justin.merondun/apptainer/puzzler_v1.8.sif | /90daydata/coffea_pangenome/puzzler_trials/modest_128gb | /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Stickleback.HiFi.fastq.gz | /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Stickleback.HiC.R1.fastq.gz | /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Stickleback.HiC.R2.fastq.gz | 21       | /90daydata/coffea_pangenome/puzzler_trials/raw_data/references/GCF_949316345.1_fPunPun2.1_genomic.fna | NA      | NA            | actinopterygii_odb10 | /90daydata/coffea_pangenome/puzzler_trials/busco_downloads |
| Toad        | apptainer | /home/justin.merondun/apptainer/puzzler_v1.8.sif | /90daydata/coffea_pangenome/puzzler_trials/modest_128gb | /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Toad.HiFi.fastq.gz | /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Toad.HiC.R1.fastq.gz | /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Toad.HiC.R2.fastq.gz | 11       | /90daydata/coffea_pangenome/puzzler_trials/raw_data/references/GCA_027917425.1_aGasCar1.pri_genomic.fna | NA      | NA            | tetrapoda_odb10      | /90daydata/coffea_pangenome/puzzler_trials/busco_downloads |
| Whale       | apptainer | /home/justin.merondun/apptainer/puzzler_v1.8.sif | /90daydata/coffea_pangenome/puzzler_trials/modest_128gb | /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Whale.HiFi.fastq.gz | /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Whale.HiC.R1.fastq.gz | /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Whale.HiC.R2.fastq.gz | 21       | /90daydata/coffea_pangenome/puzzler_trials/raw_data/references/GCA_964341445.1_mMesMir1.hap1.1_genomic.fna | NA      | NA            | mammalia_odb10       | /90daydata/coffea_pangenome/puzzler_trials/busco_downloads |
| Strawberry  | apptainer | /home/justin.merondun/apptainer/puzzler_v1.8.sif | /90daydata/coffea_pangenome/puzzler_trials/modest_128gb | /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Strawberry.HiFi.fastq.gz | /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Strawberry.HiC.R1.fastq.gz | /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Strawberry.HiC.R2.fastq.gz | 28       | /90daydata/coffea_pangenome/puzzler_trials/raw_data/references/GCA_048987885.1_NARO_FxaKMR_1.0_genomic.fna | NA      | NA            | embryophyta_odb10    | /90daydata/coffea_pangenome/puzzler_trials/busco_downloads |

You can then submit an assembly for each. Since I am on a SLURM HPC, I create a sbatch script:

```
#!/bin/bash

#SBATCH --time=5-00:00:00   
#SBATCH --nodes=1  
#SBATCH --cpus-per-task=48
#SBATCH --mem=512Gb

SAMPLE=$1

puzzler -s ${SAMPLE} -m samples.tsv --threads 48 --mem 48
```

Extract the sample name from the samples file, and submit each sample using xargs: 

```bash
cut -f1 samples.tsv | sed '1d' | xargs -I {} echo sbatch -J asm_{} Submit_Puzzler.sh {}
```

And then just remove the `echo` to submit in parallel. 

All samples will run until juicer file creation:

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
RUNTIME: conda
CONTAINER: /90daydata/coffea_pangenome/puzzler_trials/rerun_juicer/NA 
WD: /90daydata/coffea_pangenome/puzzler_trials/rerun_juicer 
HIFI: /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Fungus.HiFi.fastq.gz
HIC_R1: /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Fungus.HiC.R1.fastq.gz
HIC_R2: /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Fungus.HiC.R2.fastq.gz
NUMBER CHRS: 14
REFERENCE: /90daydata/coffea_pangenome/puzzler_trials/raw_data/references/GCA_964106605.1_gdRhiKalk1.hap1.1_genomic.fna
HOM_COV: NA
BLOB_DB: NA
BUSCO_LINEAGE: fungi_odb10
BUSCO_DB: /90daydata/coffea_pangenome/puzzler_trials/busco_downloads
Cores Requested: 8
Cores Available: 8
RAM Requested: 128
Memory Available: 2238.7 GB
PUZZLER command: 
=======================================================================

~~~ [+0m] Checking software availability, this will take about 1 minute ~~~
~~~ [+0m] Software check complete ~~~
~~~ [+0m] Assembling genome for Fungus ~~~
~~~ [+0m] Running juicer, script will stop after .hic files created ~~~
~~~ [+0m] Starting hifiasm assembly for Fungus ~~~
~~~ [+6m] Starting Purge_Dups for Fungus ~~~
~~~ [+6m] Mapping HiC reads to Fungus draft ~~~
~~~ [+8m] Running HapHiC for Fungus  ~~~

❌ Command failed in in /90daydata/coffea_pangenome/puzzler_trials/rerun_juicer/Fungus/03_haphic/haphic: "${PUZZLER} haphic pipeline ${WD}/${SAMPLE}/02_purge_dups/p_ctg.purged.fa ../filtered.bam ${NUM_CHRS} --correct_nrounds 2 --max_inflation 10.0 --threads ${t} --processes ${t} 2> ../${SAMPLE}.haphic.log" (line 298)

~~~ [+9m] HapHiC for Fungus with 14 chrs failed, trying: 10 ~~~

❌ Command failed in in /90daydata/coffea_pangenome/puzzler_trials/rerun_juicer/Fungus/03_haphic/haphic: "${PUZZLER} haphic pipeline ${WD}/${SAMPLE}/02_purge_dups/p_ctg.purged.fa ../filtered.bam ${CHR_ATTEMPT} --correct_nrounds 2 --max_inflation 10.0 --threads ${t} --processes ${t} 2> ../${SAMPLE}.haphic_N_chrs${CHR_ATTEMPT}.log" (line 321)

~~~ [+11m] HapHiC for Fungus with 14 chrs failed, trying: 11 ~~~

❌ Command failed in in /90daydata/coffea_pangenome/puzzler_trials/rerun_juicer/Fungus/03_haphic/haphic: "${PUZZLER} haphic pipeline ${WD}/${SAMPLE}/02_purge_dups/p_ctg.purged.fa ../filtered.bam ${CHR_ATTEMPT} --correct_nrounds 2 --max_inflation 10.0 --threads ${t} --processes ${t} 2> ../${SAMPLE}.haphic_N_chrs${CHR_ATTEMPT}.log" (line 321)

~~~ [+13m] HapHiC for Fungus with 14 chrs failed, trying: 12 ~~~

❌ Command failed in in /90daydata/coffea_pangenome/puzzler_trials/rerun_juicer/Fungus/03_haphic/haphic: "${PUZZLER} haphic pipeline ${WD}/${SAMPLE}/02_purge_dups/p_ctg.purged.fa ../filtered.bam ${CHR_ATTEMPT} --correct_nrounds 2 --max_inflation 10.0 --threads ${t} --processes ${t} 2> ../${SAMPLE}.haphic_N_chrs${CHR_ATTEMPT}.log" (line 321)

~~~ [+14m] HapHiC for Fungus with 14 chrs failed, trying: 13 ~~~

❌ Command failed in in /90daydata/coffea_pangenome/puzzler_trials/rerun_juicer/Fungus/03_haphic/haphic: "${PUZZLER} haphic pipeline ${WD}/${SAMPLE}/02_purge_dups/p_ctg.purged.fa ../filtered.bam ${CHR_ATTEMPT} --correct_nrounds 2 --max_inflation 10.0 --threads ${t} --processes ${t} 2> ../${SAMPLE}.haphic_N_chrs${CHR_ATTEMPT}.log" (line 321)

~~~ [+16m] HapHiC for Fungus with 14 chrs failed, trying: 14 ~~~

❌ Command failed in in /90daydata/coffea_pangenome/puzzler_trials/rerun_juicer/Fungus/03_haphic/haphic: "${PUZZLER} haphic pipeline ${WD}/${SAMPLE}/02_purge_dups/p_ctg.purged.fa ../filtered.bam ${CHR_ATTEMPT} --correct_nrounds 2 --max_inflation 10.0 --threads ${t} --processes ${t} 2> ../${SAMPLE}.haphic_N_chrs${CHR_ATTEMPT}.log" (line 321)

~~~ [+17m] HapHiC for Fungus with 14 chrs failed, trying: 15 ~~~

❌ Command failed in in /90daydata/coffea_pangenome/puzzler_trials/rerun_juicer/Fungus/03_haphic/haphic: "${PUZZLER} haphic pipeline ${WD}/${SAMPLE}/02_purge_dups/p_ctg.purged.fa ../filtered.bam ${CHR_ATTEMPT} --correct_nrounds 2 --max_inflation 10.0 --threads ${t} --processes ${t} 2> ../${SAMPLE}.haphic_N_chrs${CHR_ATTEMPT}.log" (line 321)

~~~ [+19m] HapHiC for Fungus with 14 chrs failed, trying: 16 ~~~

❌ Command failed in in /90daydata/coffea_pangenome/puzzler_trials/rerun_juicer/Fungus/03_haphic/haphic: "${PUZZLER} haphic pipeline ${WD}/${SAMPLE}/02_purge_dups/p_ctg.purged.fa ../filtered.bam ${CHR_ATTEMPT} --correct_nrounds 2 --max_inflation 10.0 --threads ${t} --processes ${t} 2> ../${SAMPLE}.haphic_N_chrs${CHR_ATTEMPT}.log" (line 321)

~~~ [+20m] HapHiC for Fungus with 14 chrs failed, trying: 17 ~~~

❌ Command failed in in /90daydata/coffea_pangenome/puzzler_trials/rerun_juicer/Fungus/03_haphic/haphic: "${PUZZLER} haphic pipeline ${WD}/${SAMPLE}/02_purge_dups/p_ctg.purged.fa ../filtered.bam ${CHR_ATTEMPT} --correct_nrounds 2 --max_inflation 10.0 --threads ${t} --processes ${t} 2> ../${SAMPLE}.haphic_N_chrs${CHR_ATTEMPT}.log" (line 321)

~~~ [+22m] HapHiC for Fungus with 14 chrs failed, trying: 18 ~~~

❌ Command failed in in /90daydata/coffea_pangenome/puzzler_trials/rerun_juicer/Fungus/03_haphic/haphic: "${PUZZLER} haphic pipeline ${WD}/${SAMPLE}/02_purge_dups/p_ctg.purged.fa ../filtered.bam ${CHR_ATTEMPT} --correct_nrounds 2 --max_inflation 10.0 --threads ${t} --processes ${t} 2> ../${SAMPLE}.haphic_N_chrs${CHR_ATTEMPT}.log" (line 321)

~~~ [+23m] HapHiC for Fungus failed, scaffolding with YAHS instead ~~~
~~~ [+23m] Creating .hic file for juicebox for Fungus ~~~
~~~ [+24m] Post curation assembly file missing for Fungus: Run Juicebox & place in /90daydata/coffea_pangenome/puzzler_trials/rerun_juicer/juicer_files/Fungus_JBAT.review.assembly ~~~
```

We see that HapHic was unable to obtain an assembly with the specified number of chromosomes, so it fell back on YAHS for scaffolding. In reality - this seems to be an issue with scaffolding small chromosomes in this specific fungal genome - as juicer clearly reveals 14 chromosome structures. 

So at this time, you will need to do juicebox manual curation. 

# Juicer Curation

Manual curation files and tutorials (very basic and quick) for n=11 assemblies can be found either on Zenodo or on [Youtube](https://www.youtube.com/watch?v=rMUiNqZwEpA). 

Save the Juicebox curation files as default, with `Assembly > Export Assembly`, which will generate e.g. `Stickleback_JBAT.review.assembly`. Pull that back into `$WD/juicer_files/`. 

# Post Juicebox Extraction & QC

Simply resubmit `puzzler`: 

```bash
cut -f1 samples.tsv | sed '1d' | xargs -I {} echo sbatch -J asm_{} Submit_Puzzler.sh {}
```

You should see an output like:

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
Parameters for sample: Frog 
RUNTIME: apptainer
CONTAINER: /home/justin.merondun/apptainer/puzzler_v1.8.sif 
WD: /90daydata/coffea_pangenome/puzzler_trials/assemblies 
HIFI: /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Frog.HiFi.fastq.gz
HIC_R1: /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Frog.HiC.R1.fastq.gz
HIC_R2: /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Frog.HiC.R2.fastq.gz
NUMBER CHRS: 12
REFERENCE: /90daydata/coffea_pangenome/puzzler_trials/raw_data/references/GCA_048569485.1_aAnoBae1.hap1_genomic.fna
HOM_COV: NA
BLOB_DB: NA
BUSCO_LINEAGE: tetrapoda_odb10
BUSCO_DB: /90daydata/coffea_pangenome/puzzler_trials/busco_downloads
Cores Requested: 64
Cores Available: 64
RAM Requested: 512
Memory Available: 2238.4 GB
PUZZLER command: apptainer exec --bind /90daydata/coffea_pangenome/puzzler_trials:/90daydata/coffea_pangenome/puzzler_trials --bind /90daydata/coffea_pangenome/puzzler_trials/assemblies:/90daydata/coffea_pangenome/puzzler_trials/assemblies --bind /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads:/90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads --bind /90daydata/coffea_pangenome/puzzler_trials/raw_data/references:/90daydata/coffea_pangenome/puzzler_trials/raw_data/references --bind /home/justin.merondun/apptainer:/home/justin.merondun/apptainer /home/justin.merondun/apptainer/puzzler_v1.8.sif
=======================================================================

~~~ [+0m] Assembling genome for Frog ~~~
~~~ [+0m] Running juicer, script will stop after .hic files created ~~~
~~~ [+0m] Skipping hifiasm for Frog: /90daydata/coffea_pangenome/puzzler_trials/assemblies/Frog/01_hifiasm/asm.hic.p_ctg.gfa exists ~~~
~~~ [+0m] Skipping purge for Frog: /90daydata/coffea_pangenome/puzzler_trials/assemblies/Frog/02_purge_dups/p_ctg.purged.fa exists ~~~
~~~ [+0m] Skipping HiC alignment for Frog: /90daydata/coffea_pangenome/puzzler_trials/assemblies/Frog/03_haphic/filtered.bam exists ~~~
~~~ [+0m] Skipping HapHiC for Frog: /90daydata/coffea_pangenome/puzzler_trials/assemblies/Frog/03_haphic/haphic/04.build/scaffolds.fa exists ~~~
~~~ [+0m] Skipping juicer HiC file creation for Frog: /90daydata/coffea_pangenome/puzzler_trials/assemblies/juicer_files/Frog_JBAT.hic exists ~~~
~~~ [+0m] Extracting post-curation assembly and mapping to reference for Frog ~~~
~~~ [+62m] Renaming chromosomes for Frog ~~~
~~~ [+62m] Scaffold sanity check passed for renaming, proceeding! ~~~
~~~ [+63m] Multiple scaffolds corresponding to single Chr for Frog, Renaming them e.g. Chr1A, Chr1B.. ~~~
~~~ [+66m] Assessing genome quality for Frog ~~~
~~~ [+66m] Creating final HiC bam for Frog ~~~
~~~ [+672m] Creating final contact map for Frog, plotting named chromosomes ~~~
~~~ [+705m] Running YAK on Frog ~~~
```

Finished! 

## Statistics on Reference Assemblies

Of course, we can also extract the same statistics from the reference fastas (e.g. BUSCO, N50). I will just copy the reference fastas in `$WD/primary_asm/$SAMPLE_Ref.fa`, and then submit `puzzler` with duplicated rows, except appending "_Ref":

| sample          | runtime   | container                                        | wd                                                    | hifi                                                         | hic_r1                                                       | hic_r2                                                       | num_chrs | reference                                                    | hom_cov | blob_database | busco_lineage        | busco_database                                             |
| --------------- | --------- | ------------------------------------------------ | ----------------------------------------------------- | ------------------------------------------------------------ | ------------------------------------------------------------ | ------------------------------------------------------------ | -------- | ------------------------------------------------------------ | ------- | ------------- | -------------------- | ---------------------------------------------------------- |
| Frog_Ref        | apptainer | /home/justin.merondun/apptainer/puzzler_v1.8.sif | /90daydata/coffea_pangenome/puzzler_trials/assemblies | /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Frog.HiFi.fastq.gz | /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Frog.HiC.R1.fastq.gz | /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Frog.HiC.R2.fastq.gz | 12       | /90daydata/coffea_pangenome/puzzler_trials/raw_data/references/GCA_048569485.1_aAnoBae1.hap1_genomic.fna | NA      | NA            | tetrapoda_odb10      | /90daydata/coffea_pangenome/puzzler_trials/busco_downloads |
| Beaver_Ref      | apptainer | /home/justin.merondun/apptainer/puzzler_v1.8.sif | /90daydata/coffea_pangenome/puzzler_trials/assemblies | /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Beaver.HiFi.fastq.gz | /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Beaver.HiC.R1.fastq.gz | /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Beaver.HiC.R2.fastq.gz | 20       | /90daydata/coffea_pangenome/puzzler_trials/raw_data/references/GCF_047511655.1_mCasCan1.hap1v2_genomic.fna | NA      | NA            | mammalia_odb10       | /90daydata/coffea_pangenome/puzzler_trials/busco_downloads |
| Stickleback_Ref | apptainer | /home/justin.merondun/apptainer/puzzler_v1.8.sif | /90daydata/coffea_pangenome/puzzler_trials/assemblies | /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Stickleback.HiFi.fastq.gz | /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Stickleback.HiC.R1.fastq.gz | /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Stickleback.HiC.R2.fastq.gz | 21       | /90daydata/coffea_pangenome/puzzler_trials/raw_data/references/GCF_949316345.1_fPunPun2.1_genomic.fna | NA      | NA            | actinopterygii_odb10 | /90daydata/coffea_pangenome/puzzler_trials/busco_downloads |

Submit with:

```bash
cut -f1 samples_Ref.tsv | sed '1d' | xargs -I {} echo sbatch -J asm_{} Submit_Puzzler.sh {}
```

Example output from Frog:

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
Parameters for sample: Frog_Ref 
RUNTIME: apptainer
CONTAINER: /home/justin.merondun/apptainer/puzzler_v1.8.sif 
WD: /90daydata/coffea_pangenome/puzzler_trials/assemblies 
HIFI: /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Frog.HiFi.fastq.gz
HIC_R1: /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Frog.HiC.R1.fastq.gz
HIC_R2: /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Frog.HiC.R2.fastq.gz
NUMBER CHRS: 12
REFERENCE: /90daydata/coffea_pangenome/puzzler_trials/raw_data/references/GCA_048569485.1_aAnoBae1.hap1_genomic.fna
HOM_COV: NA
BLOB_DB: NA
BUSCO_LINEAGE: tetrapoda_odb10
BUSCO_DB: /90daydata/coffea_pangenome/puzzler_trials/busco_downloads
Cores Requested: 64
Cores Available: 64
RAM Requested: 512
Memory Available: 2228.6 GB
PUZZLER command: apptainer exec --bind /90daydata/coffea_pangenome/puzzler_trials:/90daydata/coffea_pangenome/puzzler_trials --bind /90daydata/coffea_pangenome/puzzler_trials/assemblies:/90daydata/coffea_pangenome/puzzler_trials/assemblies --bind /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads:/90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads --bind /90daydata/coffea_pangenome/puzzler_trials/raw_data/references:/90daydata/coffea_pangenome/puzzler_trials/raw_data/references --bind /home/justin.merondun/apptainer:/home/justin.merondun/apptainer /home/justin.merondun/apptainer/puzzler_v1.8.sif
=======================================================================

~~~ [+0m] Skipping assembly for Frog_Ref: /90daydata/coffea_pangenome/puzzler_trials/assemblies/primary_asm/Frog_Ref.fa exists ~~~
~~~ [+0m] Assessing genome quality for Frog_Ref ~~~
~~~ [+0m] Creating final HiC bam for Frog_Ref ~~~
```

## Assembly Summaries

Merge the output `*summary.txt` files: 

```bash
awk 'NR == 1 || FNR > 1' *summary.txt > Summaries.tsv
```

And visualize in R:

```
#### Summarize assemblies
setwd('/90daydata/coffea_pangenome/puzzler_trials/assemblies/primary_asm/stats/')
library(GGally)
library(tidyverse)
library(RColorBrewer)
library(ggpubr)
library(meRo) #devtools::install_github('merondun/meRo')
library(gt)
library(Biostrings)

# Import metadata / order
md <- read_tsv('/90daydata/coffea_pangenome/puzzler_trials/assemblies/metadata.tsv')
label <- md %>% mutate(lab = paste0(Species,': ',Time)) %>% select(Sample,lab,Phylo_Order) %>% arrange(Phylo_Order)

dats <- NULL
# Count chromosomes with 'N' (gaps) for both draft and reference 
for (sp in unique(md$Sample)) { 
  
  cat('Working on sample: ',sp,'\n')
  #Draft
  fa <- readDNAStringSet(paste0('/90daydata/coffea_pangenome/puzzler_trials/assemblies/primary_asm/',sp,'.fa'))
  gapless_chrs <- tibble(
    chr = names(fa),
    has_N = vcountPattern("N", fa) > 0,
    length = width(fa)) 
  gc <- gapless_chrs %>% mutate(Sample = sp)
  dats <- rbind(dats,gc)
  rm(fa,gc,gapless_chrs)
  
  # Reference 
  fa <- readDNAStringSet(paste0('/90daydata/coffea_pangenome/puzzler_trials/assemblies/primary_asm/',sp,'_Ref.fa'))
  gapless_chrs <- tibble(
    chr = names(fa),
    has_N = vcountPattern("N", fa) > 0,
    length = width(fa)) %>% filter(grepl('chr|Chr',chr))
  gc <- gapless_chrs %>% mutate(Sample = paste0(sp,'_Ref'))
  dats <- rbind(dats,gc)
}

# Match the reference and puzzler chromosomes
compiled <- dats %>% filter(!grepl('scaff|organ|mito|_unloc',chr)) %>%
  mutate(
    Species = gsub("_Ref", "", Sample),
    ref_type = ifelse(grepl("_Ref", Sample), "Ref", "")
  ) %>%
  select(chr, Species, ref_type, has_N, length) %>%
  pivot_wider(
    names_from = ref_type,
    values_from = c(has_N, length),
    names_glue = "{.value}{ifelse(ref_type == 'Ref', '_Ref', '')}") 
write.table(compiled,'20251124_Gapless_Summaries.txt',quote=F,sep='\t',row.names=F)
compiled <- read_tsv('20251124_Gapless_Summaries.txt')

# First, set Puzzler chrs to gapped if they are less than 3 Mb from Reference size, then count number of gapless and gapped for puzzler and ref 
compstats <- compiled %>% 
  mutate(dlength = (length-length_Ref)/1e6,
         has_N = ifelse(abs(dlength) > 1,TRUE,has_N)) %>% 
  pivot_longer(c(has_N,has_N_Ref)) %>% 
  mutate(value = gsub('TRUE','Gapped',gsub('FALSE','Gapless',value))) %>% 
  replace_na(list(value='Missing')) %>% 
  group_by(Species,name) %>% 
  count(value) %>% ungroup %>% 
  complete(Species,name,value, fill = list(n = 0)) %>% 
  mutate(Assembly = ifelse(grepl('Ref',name),'Reference','Puzzler'))

# Plot 
gap_in <- left_join(compstats %>% dplyr::rename(Sample=Species),md %>% select(Sample,Name,Species,Phylo_Order)) %>% 
  mutate(splab = paste0(Name,'\n',Species)) %>% arrange(Phylo_Order)
gap_in$splab <- factor(gap_in$splab,levels=rev(unique(gap_in$splab)))
gap_plot <- gap_in %>% 
  ggplot(aes(y=splab,x=n,fill=value))+
  geom_col(width=1,col='black',position=position_dodge(width=0.7))+
  geom_text(aes(label=n),position=position_dodge(width=1),hjust=-0.5)+
  scale_fill_manual('Chromosome',values=c('black','white','grey'))+
  facet_grid(.~Assembly,scales='free')+
  coord_cartesian(xlim=c(0,43))+
  xlab('Number of Gapless Chromosomes\nRequiring That Puzzler Chrs Are ± 1 Mb From Reference Size')+
  theme_bw()
gap_plot
ggsave('20251126_GapCounts.pdf',gap_plot,height=7,width=7)

# Add metadata information
s <- read_tsv('Summaries.tsv')
sl <- s %>% 
  pivot_longer(!c(Sample,Species,Name,Assembly,Phylo_Order))
levels <- unique(sl$name) 
sl$name <- factor(sl$name,levels=levels)
slf <- left_join(sl,label %>% arrange(Phylo_Order) %>% select(-Phylo_Order))

# Estimate fold change from Reference
fc <- slf %>%
  group_by(Sample, Species, name) %>%
  mutate(
    ref_value = value[Assembly == "Ref"],
    ref_diff = value - ref_value, 
    max_x = max(value),
    # compute % relative to reference
    value_pct = 100 * value / ref_value,
    delta_pct = value_pct - 100
  ) %>%
  ungroup()

fc %>% group_by(Assembly,name) %>% sum_stats(delta_pct) %>%   
  mutate(lab = paste0(round(mean,1),' (',round(min,1),'–',round(max,1),')')) %>% 
  select(Assembly,name,lab) %>% 
  pivot_wider(names_from = name,values_from = lab)

fc_plot %>% 
  filter(Assembly != 'Ref' & grepl('Size|Seq|N50|BUSCO_Complete',name)) %>% 
  ggplot(aes(y=name,x=mean,xmin=pmax(conf_low,0),xmax=conf_high,fill=Assembly))+
  geom_col(position=position_dodge(width=0.8))+
  geom_errorbar(position=position_dodge(width=0.8))+
  scale_fill_manual(values = c("cadetblue", "goldenrod1")) +
  theme_bw()

# More summaries
s %>% pivot_longer(c(Gaps,ContigN50,ScafN50,YAK_CV,YAK_QV,BUSCO_Complete,BUSCO_singlecopy)) %>% 
  group_by(Assembly,name) %>% 
  sum_stats(value) %>% 
  data.frame

# labs for diff between puzzler and ref 
labs <- fc %>% filter(grepl('ContigN50|ScafN50|SizeBP|BUSCO_Complete',name))
#labs$Species = factor(labs$Species,levels=unique(slf$Species))
slf$Species = factor(slf$Species,levels=rev(unique(slf$Species)))

# Plot
all <- slf %>%
  filter(grepl('ContigN50|ScafN50|SizeBP|BUSCO_Complete',name)) %>% 
  ggplot(aes(y = Species, fill = Assembly, x = value, shape = Assembly)) +
  #geom_bar(stat = "identity", position = position_dodge(width = 0.9), col = "grey50") +
  geom_point(position = position_dodge(width = 0.8),stroke=0.25)+
  geom_text(data = labs, aes(x = max_x *1.01, y = Species, label = signif(ref_diff,2)), 
            position = position_dodge(width = 0.8), hjust = 0.5,size=1) +  # Adjust label position dynamically
  scale_fill_manual(values = brewer.pal(3,'Paired')) +
  scale_shape_manual(values=c(21,22,25))+
  facet_wrap(name ~ ., scales = "free_x",nrow=1) +  # Allows individual x-axis scaling per facet
  theme_bw() +
  theme(
    strip.text = element_text(size = 12),  # Increase facet label size
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels for readability
    panel.spacing = unit(2, "lines")  # Increase spacing between facets
  )
all
ggsave('20251128_Assembly_Summaries.pdf',all,height=3.5,width=10)

# Summarize assembly time vs inputs (64 cores)
asm_stats <- s %>% filter(Assembly == 'Ref') %>% select(Sample,Species,SizeBP,ContigN50,ScafN50,Gaps,PropWithinChrs)
tl <- left_join(md %>% select(-SizeBP),asm_stats) %>% 
  dplyr::rename(hifi = `HiFi Bases`, hic = `HiC Bases`, length = `HiFi Length`, cpus = `Peak CPUs`, mem = `Peak RAM`) %>% 
  mutate(time_gb = SizeBP / Minutes,
         hifi_bp = hifi / SizeBP,
         hic_bp = hic / SizeBP,
         ContigN50 = ContigN50,
         length = length / 1e3) %>% 
  select(Sample, time_gb, Heterozygosity, hifi_bp, hic_bp, SizeBP, hifi, hic,ContigN50,ScafN50,Gaps,PropWithinChrs,length)
tl
tl$Sample = factor(tl$Sample,levels=unique(slf$Sample))

# COVERAGE: Define scaling factor
scale_factor <- max(tl$ContigN50) / max(tl$hifi_bp, tl$hic_bp)
tl <- tl %>%
  mutate(hifi_scaled = hifi_bp * scale_factor,
         hic_scaled = hic_bp * scale_factor)
input_in <- left_join(tl,md %>% select(Sample,Name,Species,Phylo_Order)) %>% 
  mutate(splab = paste0(Name,'\n',Species)) %>% arrange(Phylo_Order)
input_in$splab <- factor(input_in$splab,levels=unique(input_in$splab))
inputs <- ggplot(input_in, aes(x = splab)) +
  geom_col(aes(y = hifi_scaled, group = 1), fill = "dodgerblue",width=0.35,position=position_nudge(x=-0.1)) +
  geom_col(aes(y = hic_scaled, group = 1), fill = "orange",width=0.35,position=position_nudge(x=0.1)) +
  geom_point(aes(y = ContigN50,size=Heterozygosity), color = "black") +
  #geom_line(aes(y = hifi_scaled, group = 1), color = "dodgerblue", linetype = "dashed") +
  #geom_line(aes(y = hic_scaled, group = 1), color = "orange", linetype = "dotted") +
  scale_y_continuous(
    name = "Contig N50 (Mb)",
    sec.axis = sec_axis(~./scale_factor, name = "Coverage (Input Bases / Assembly Size)", breaks = scales::pretty_breaks(5))) +
  xlab("") +
  theme_bw(base_size = 9) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.y.right = element_text(color = "grey20"),
    axis.text.y.right = element_text(color = "grey20"))
inputs

ggsave('20251128_InputsPlot.pdf',inputs,height=4,width=5.5)

input_in %>% pivot_longer(c(hifi_bp,hic_bp, length)) %>% group_by(name) %>% summarize(min = min(value),mean=mean(value),max=max(value))

tl %>% sum_stats(time_gb)
# # A tibble: 1 × 9
# mean   min   max    sd    se median   iqr conf_low conf_high
# <dbl> <dbl> <dbl> <dbl> <dbl>  <dbl> <dbl>    <dbl>     <dbl>
#   1  1.75 0.535  3.77 0.921 0.266   1.67  1.03     1.16      2.33

###### Plot correlation pairs: assembly quality 
for_pairs <- tl %>% select(Sample,ScafN50,ContigN50,Gaps,PropWithinChrs,Heterozygosity,SizeBP,hifi_bp,hic_bp,length)
for_pairs$Sample <- factor(for_pairs$Sample)
pdf('20251128_GGpairs_Correlations_Reference.pdf',height=7,width=9)
for_pairs %>% ggpairs + theme_bw()
dev.off()

###### Assembly times by step
steps <- read_tsv('/90daydata/coffea_pangenome/puzzler_trials/rerun_assemblies/Assembly_Times.tsv')
steps$Species <- factor(steps$Species,levels=unique(steps$Species))

# calculate time per step not total
steps <- steps %>%
  group_by(Species) %>%
  mutate(
    Process = row_number(),
    Step_Time = if_else(Process == 1, Time, Time - lag(Time))
  )

steps %>% ungroup %>% group_by(Step) %>% summarize(min = min(Step_Time),mean = mean(Step_Time),max=max(Step_Time))
steps %>% group_by(Step) %>% sum_stats(Step_Time)
steps %>% 
  select(-Sample,-Phylo_Order) %>% 
  pivot_wider()

```

# Alignments

Whole genome alignment between chromosomes. I will repeat this for both manual curation genomes and `--no_juice` genomes. 

Just grab the named chromosomes for alignment:

```bash
#!/bin/bash

#SBATCH --time=1-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=10
#SBATCH --mem=64Gb
#SBATCH --partition=atlas
#SBATCH --account=coffea_pangenome

module load apptainer
MAP_FILE=/90daydata/coffea_pangenome/puzzler_trials/assemblies/samples.tsv

mkdir -p pafs_juicer

for SAMPLE in $(cat Samples.list); do

# Read either tsv or csv and assign variables from map file
IFS=$'\t,' read -r _ RUNTIME CONTAINER WD HIFI HIC_R1 HIC_R2 NUM_CHRS REFERENCE HOM_COV BLOB_DB BUSCO_LINEAGE BUSCO_DB < <(
    awk -F'[\t,]' -v sample="$SAMPLE" '$1 == sample {print $0}' "${MAP_FILE}"
)

PUZZLER="apptainer exec ${CONTAINER}"

# Only grab scaffolds > 1mb and not 'unlocalized' superscaffolds for reference
$PUZZLER samtools faidx ${REFERENCE}
egrep 'chr|Chr' ${REFERENCE}.fai | awk '$2 > 1e6' | egrep -v 'SUPER|super|unloc' | cut -f1 > ${SAMPLE}_Ref.tmp
$PUZZLER seqtk subseq ${REFERENCE} ${SAMPLE}_Ref.tmp > fastas/${SAMPLE}_Ref.fa

# Only grab scaffolds > 1mb and not 'unlocalized' superscaffolds for draft 
$PUZZLER samtools faidx ${WD}/primary_asm/${SAMPLE}.fa
egrep 'chr|Chr' ${WD}/primary_asm/${SAMPLE}.fa.fai | awk '$2 > 1e6' | egrep -v 'SUPER|super|unloc' | cut -f1 > ${SAMPLE}.tmp
$PUZZLER seqtk subseq ${WD}/primary_asm/${SAMPLE}.fa ${SAMPLE}.tmp > fastas/${SAMPLE}.fa

echo "WGA for ${SAMPLE}"
${PUZZLER} mashmap -r fastas/${SAMPLE}_Ref.fa -q fastas/${SAMPLE}.fa -t 10 -s 10000 --perc_identity 99 -o pafs_juicer/${SAMPLE}.paf 2> pafs_juicer/${SAMPLE}.mashmap.log
${PUZZLER} paf2dotplot.R pafs_juicer/${SAMPLE}.paf -r 1e6 -m 1e4 -p 8 -c 1

done

rm *tmp
```

Plot in [R](https://dwinter.github.io/pafr/reference/dotplot.html) .

```R
# install.packages("devtools")
# library(devtools)
# devtools::install_github("dwinter/pafr")
library(pafr)
library(viridis)
library(tidyverse)
library(stringr)
library(ggpubr)

md <- read_tsv('/90daydata/coffea_pangenome/puzzler_trials/assemblies/metadata.tsv')
label <- md %>% mutate(lab = paste0(Name,'\n',Species,' (',Time,')')) %>% select(Sample,lab,Phylo_Order) %>% arrange(Phylo_Order)

setwd('/90daydata/coffea_pangenome/puzzler_trials/alignments/pafs_juicer/')
pafs = list.files(path = '.', pattern = '*paf$')
counter = 0
for (species in unique(label$Sample)) { counter = counter + 1

cat('Working on : ',species,'\n')
paf_file <- read_paf(paste0(species,'.paf'))
qchr <- data.frame(chr=unique(paf_file$qname)) %>% 
  mutate(chrs=ifelse(grepl('chrX|chrY',chr),99,chr),
         chrs=as.numeric(gsub("\\D", "", chrs))) %>% arrange(chrs,chr) %>% pull(chr)
tchr <- data.frame(chr=unique(paf_file$tname)) %>% 
  mutate(chrs=ifelse(grepl('chrX|chrY',chr),99,chr),
         chrs=as.numeric(gsub("\\D", "", chrs))) %>% arrange(chrs,chr) %>% pull(chr)
d <- dotplot(paf_file,
             order_by = 'provided', ordering=list(qchr,tchr) ,
             label_seqs = FALSE,alignment_colour = 'blue',xlab = '',ylab='',dashes = FALSE,line_size=0.5) +
  ggtitle(label %>% filter(Sample == species) %>% pull(lab))+
  theme_classic(base_size=12) 
d
assign(paste0('p',counter),d)
}

#pdf('../20250623_Dotplots-NoLines.pdf',height=15,width=25)
#png('../20250623_Dotplots.png',height=7,width=10,res=300,units='in')
png('../20251119_Dotplots-NoHashes-Juicer.png',height=12,width=17,res=300,units='in')
#pdf('../20251119_Dotplots-NoHashes-Juicer.pdf',height=12,width=17)
ggarrange(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12)
dev.off()

```

# Purge Duplicates Sensitivity

```
#!/bin/bash

########### EDIT THIS BLOCK WITH SLURM & APPTAINER/SINGULARITY SETTINGS ############
#SBATCH --time=2-00:00:00   
#SBATCH --nodes=1  
#SBATCH --cpus-per-task=32
#SBATCH --mem=128Gb
#SBATCH --partition=ceres
#SBATCH --account=coffea_pangenome

module load apptainer &> /dev/null || true
SINGULARITY_TMPDIR=$APPTAINER_TMPDIR

########### EDIT THIS BLOCK WITH SLURM & APPTAINER/SINGULARITY SETTINGS ############

# Default values
t=32
MEM=128
SAMPLE=""
MAP_FILE=""
VERSION="v1.8"
NO_JUICE="FALSE"

# Display help message
function show_help {
    echo "Usage: puzzler -s sample -m samples.tsv [OPTIONS]"
    echo ""
    echo "Options:"
    echo "  -s, --sample SAMPLE   Sample name (required)"
    echo "  -m, --map FILE        Path to .tsv/.csv map file (required)"
    echo "  --threads t           Number of threads (optional; default 64)"
    echo "  --mem MEM             Memory allocation (optional; default 512)"
    echo "  --no_juice            Skip juicer file creation entirely (optional; not recommended!)"
    echo "  -v, --version         Show version and exit"
    echo "  -h, --help            Show help and exit"
    echo ""
    echo "  Required --map Structure:"
    echo "  The provided map file (e.g., samples.txt) must contain the following columns in this order:"
    echo "  RUNTIME CONTAINER WD HIFI HIC_R1 HIC_R2 NUM_CHRS REFERENCE HOM_COV BLOB_DB BUSCO_LINEAGE BUSCO_DB"
    echo "  For optional columns (REFERENCE - BUSCO_DB), write "NA" if undesired."
    echo ""

    exit 0
}

########## start checks #############

# Parse command-line arguments
while [[ "$#" -gt 0 ]]; do
    case "$1" in
        -h|--help) show_help ;;  
        -v|--version) echo "puzzler version: $VERSION"; exit 0 ;;
        -s|--sample) SAMPLE="$2"; shift ;;
        -m|--map) MAP_FILE="$2"; shift ;;
        --no_juice) NO_JUICE="TRUE" ;;
        --threads) t="$2"; shift ;;
        --mem) MEM="$2"; shift ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
    shift
done

# Validate required arguments
if [ -z "$SAMPLE" ]; then
    echo "Error: --sample argument is required."
    exit 1
fi

if [ -z "$MAP_FILE" ]; then
    echo "Error: --map argument is required."
    exit 1
fi

# Read either tsv or csv and assign variables from map file
IFS=$'\t,' read -r _ RUNTIME CONTAINER WD HIFI HIC_R1 HIC_R2 NUM_CHRS REFERENCE HOM_COV BLOB_DB BUSCO_LINEAGE BUSCO_DB < <(
    awk -F'[\t,]' -v sample="$SAMPLE" '$1 == sample {print $0}' "${MAP_FILE}"
)

# Check if mandatory fields (RUNTIME - NUM_CHRS) are empty
if [[ -z "$RUNTIME" || -z "$CONTAINER" || -z "$WD" || -z "$HIFI" || -z "$HIC_R1" || -z "$HIC_R2" || -z "$NUM_CHRS" ]]; then
    echo -e "\e[41m~~~ ERROR: Missing required fields (RUNTIME - NUM_CHRS). Exiting. ~~~\e[0m"
    exit 1
fi

# Set REFERENCE - BUSCO_DB to "NA" if empty, otherwise extract full path
HOM_COV="${HOM_COV:-NA}"
BUSCO_LINEAGE="${BUSCO_LINEAGE:-NA}"
for option in REFERENCE BLOB_DB BUSCO_DB; do
    [[ "${!option}" != "NA" && -n "${!option}" ]] && declare "$option=$(realpath "${!option}")" || declare "$option=NA"
done

# In case relative paths were provided... extract full paths
CONTAINER=$(realpath ${CONTAINER})
WD=$(realpath ${WD})
HIFI=$(realpath ${HIFI})
HIC_R1=$(realpath ${HIC_R1})
HIC_R2=$(realpath ${HIC_R2})

if [ "$RUNTIME" = "conda" ]; then
    PUZZLER=""  # No runtime needed for Conda
    HAPHIC_PATH=$(which haphic)
    JARDIR=$(dirname "$HAPHIC_PATH")/utils
    JARFILE="${JARDIR}/juicer_tools.1.9.9_jcuda.0.8.jar"
    RUN_JUICERTOOLS="${PUZZLER} java -Xmx16G -jar ${JARFILE}"
else
    # On some cluster architecture, you need to specify which directories to bind.... we will simply bind all paths. 
    binds=()

    for varname in CONTAINER WD HIFI HIC_R1 HIC_R2 REFERENCE BLOB_DB BUSCO_DB; do
        val="${!varname}"
        if [[ -n "$val" && ( -e "$val" || -d "$val" ) ]]; then
            dir=$(realpath "$val" | xargs dirname)
            binds+=("$dir")
        fi
    done

    # deduplicate
    IFS=$'\n' read -r -d '' -a unique_binds < <(printf "%s\n" "${binds[@]}" | sort -u && printf '\0')

    bind_flags=()
    for b in "${unique_binds[@]}"; do
        bind_flags+=(--bind "$b:$b")
    done

    PUZZLER="${RUNTIME} exec ${bind_flags[@]} ${CONTAINER}"
    JARFILE="/opt/HapHiC/utils/juicer_tools.1.9.9_jcuda.0.8.jar"
    RUN_JUICERTOOLS="${PUZZLER} java -Xmx16G -jar ${JARFILE}"
fi

if [ "$NO_JUICE" = "TRUE" ]; then
    JUICE_PRINT="Skipping juicer, no manual curation (not recommended!)"
else
    JUICE_PRINT="Running juicer, script will stop after .hic files created"
fi

cat << "EOF"

=======================================================================
__________ ____ _______________________.____     _____________________ 
\______   \    |   \____    /\____    /|    |    \_   _____/\______   \
 |     ___/    |   / /     /   /     / |    |     |    __)_  |       _/
 |    |   |    |  / /     /_  /     /_ |    |___  |        \ |    |   \
 |____|   |______/ /_______ \/_______ \|_______ \/_______  / |____|_  /
                           \/        \/        \/        \/         \/ 
=======================================================================

EOF

core_avail=$(nproc)
mem_avail=$(awk '/MemAvailable/ {printf "%.1f GB\n", $2/1024/1024}' /proc/meminfo)

echo -e "=======================================================================\nParameters for sample: ${SAMPLE} \nRUNTIME: ${RUNTIME}\nCONTAINER: ${CONTAINER} \nWD: ${WD} \nHIFI: ${HIFI}\nHIC_R1: ${HIC_R1}\nHIC_R2: ${HIC_R2}\nNUMBER CHRS: ${NUM_CHRS}\nREFERENCE: ${REFERENCE}\nHOM_COV: ${HOM_COV}\nBLOB_DB: ${BLOB_DB}\nBUSCO_LINEAGE: ${BUSCO_LINEAGE}\nBUSCO_DB: ${BUSCO_DB}\nCores Requested: ${t}\nCores Available: ${core_avail}\nRAM Requested: ${MEM}\nMemory Available: ${mem_avail}\nPUZZLER command: ${PUZZLER}\n=======================================================================\n"

set -euo pipefail

### MONITORING ###
exec 3>&2
trap 'echo -e "\n❌ \e[41mCommand failed in in $(pwd):\e[0m \"$BASH_COMMAND\" (line $LINENO)\n" >&2' ERR
SCRIPT_START_TIME=$SECONDS
elapsed() {
    local mins=$(( (SECONDS - SCRIPT_START_TIME) / 60 ))
    printf "[+%dm]" "$mins"
}
### END MONITORING ###

echo -e "\e[43m~~~ $(elapsed) Checking software availability, this will take about 1 minute ~~~\e[0m"

declare -A TOOL_CMDS=(
  [hifiasm]="hifiasm --version"
  [haphic]="haphic --version"
  [juicer]="juicer pre --version"
  [yak]="yak version"
  [busco]="busco --version"
  [asm_stats]="assembly_stats -h"
  [juicertools]="java -Xmx2G -jar ${JARFILE} pre --version"
)

for tool in "${!TOOL_CMDS[@]}"; do
  if ! ${PUZZLER} ${TOOL_CMDS[$tool]} > /dev/null 2>&1; then
    echo -e "\e[41m~~~ $tool not found, is runtime set and container exists? ~~~\e[0m"
    exit 1
  fi
done

echo -e "\e[43m~~~ $(elapsed) Software check complete ~~~\e[0m"

########## End checks #############

# Purge dups trial 
mkdir -p ${WD}/${SAMPLE}/02_purge_dups
echo -e "\e[43m~~~ $(elapsed) Starting Purge_Dups Trial for ${SAMPLE} ~~~\e[0m"
cd ${WD}/${SAMPLE}/02_purge_dups

# Purge duplicates
${PUZZLER} purge_dups pri.split.self.paf.gz > pri.dups.DEFAULT.bed 2> ${SAMPLE}.purge.DEFAULT.log

# Check if pri.dups.bed is empty
if [ ! -s pri.dups.DEFAULT.bed ]; then
    echo -e "\e[42m~~~ $(elapsed) Skipping purge for ${SAMPLE}: no duplicates found! ~~~\e[0m"
    cp pri.split.fa p_ctg.purged.DEFAULT.fa 
else
    # Continue with purge process if duplicates exist
    ${PUZZLER} get_seqs pri.dups.DEFAULT.bed pri.init.fa 2> "${SAMPLE}.getseqs.DEFAULT.log"
    mv purged.fa p_ctg.purged.DEFAULT.fa
fi

##### ALIGN HIC #####
if [ -s "${WD}/${SAMPLE}/03_haphic_nopurge/filtered.bam" ] && [ -f "${WD}/${SAMPLE}/align_hic.complete" ]; then
    echo -e "\e[42m~~~ $(elapsed) Skipping HiC alignment for ${SAMPLE}: ${WD}/${SAMPLE}/03_haphic_nopurge/filtered.bam exists ~~~\e[0m"
elif [ ! -s "${WD}/${SAMPLE}/02_purge_dups/p_ctg.purged.fa" ]; then
    echo -e "\e[41m~~~ $(elapsed) Skipping HiC alignment for ${SAMPLE}, missing input: ${WD}/${SAMPLE}/02_purge_dups/p_ctg.purged.fa ~~~\e[0m"
    exit 1 
else
    mkdir -p ${WD}/${SAMPLE}/03_haphic_nopurge
    echo -e "\e[43m~~~ $(elapsed) Mapping HiC reads to ${SAMPLE} draft ~~~\e[0m" 
    cd ${WD}/${SAMPLE}/03_haphic_nopurge
    # Index reference
    ${PUZZLER} bwa-mem2 index ${WD}/${SAMPLE}/02_purge_dups/pri.init.fa > ${SAMPLE}.alignment.indexing.hic.log 2>&1

    # Align Hi-C reads
    { ${PUZZLER} bwa-mem2 mem -5SP -t ${t} ${WD}/${SAMPLE}/02_purge_dups/pri.init.fa ${HIC_R1} ${HIC_R2} | \
        ${PUZZLER} samblaster | \
        ${PUZZLER} samtools view - -@ ${t} -S -h -b -F 3340 | \
        ${PUZZLER} filter_bam - 1 --nm 3 --threads ${t} | \
        ${PUZZLER} samtools view - -b -@ ${t} -o filtered.bam; } > ${SAMPLE}.alignment.hic.log 2>&1
    touch ${WD}/${SAMPLE}/align_hic.complete
fi

##### Also do the HapHiC step without purge dups... #####
if [ -s "${WD}/${SAMPLE}/03_haphic_nopurge/haphic/04.build/scaffolds.fa" ] && [ -f "${WD}/${SAMPLE}/scaffolding_nopurge.complete" ]; then
    echo -e "\e[42m~~~ $(elapsed) Skipping HapHiC for ${SAMPLE}: ${WD}/${SAMPLE}/03_haphic_nopurge/haphic/04.build/scaffolds.fa exists ~~~\e[0m"
elif [ ! -s "${WD}/${SAMPLE}/03_haphic_nopurge/filtered.bam" ]; then
    echo -e "\e[41m~~~ $(elapsed) Skipping HapHiC for ${SAMPLE}, missing input: ${WD}/${SAMPLE}/03_haphic/filtered.bam ~~~\e[0m"
    exit 1 
elif [[ -z "$NUM_CHRS" || "$NUM_CHRS" == "NA" ]]; then
    echo -e "\e[41m~~~ $(elapsed) Skipping HapHiC for ${SAMPLE}, must specify chromosome number estimate or guess ~~~\e[0m"
    exit 1 
else

    echo -e "\e[43m~~~ $(elapsed) Running HapHiC for ${SAMPLE}  ~~~\e[0m" 
    mkdir -p ${WD}/${SAMPLE}/03_haphic_nopurge/haphic
    set +euo pipefail
    # Clean up any previous runs, otherwise fails ... 
    rm -rf ${WD}/${SAMPLE}/03_haphic_nopurge/haphic/*

    cd ${WD}/${SAMPLE}/03_haphic_nopurge/haphic/
    ${PUZZLER} haphic pipeline ${WD}/${SAMPLE}/02_purge_dups/pri.init.fa ../filtered.bam ${NUM_CHRS} --correct_nrounds 2 --max_inflation 10.0 --threads ${t} --processes ${t} 2> ../${SAMPLE}.haphic.log
    cp ${WD}/${SAMPLE}/03_haphic_nopurge/haphic/01.cluster/HapHiC_cluster.log ${WD}/${SAMPLE}/03_haphic_nopurge/${SAMPLE}.haphic.cluster_WOULD_CONTAIN_ERRORS.log

    # IF HapHic failed with the specified number of chromosomes $NUM_CHRS, then loop from NUM_CHRS - 1, NUM_CHRS - 2, until successful 
    if [ ! -s "${WD}/${SAMPLE}/03_haphic_nopurge/haphic/04.build/scaffolds.fa" ]; then
    rm -rf ${WD}/${SAMPLE}/03_haphic_nopurge/haphic/*

    # Define range of chrs to test: Within 4 of NUM_CHRS
    START=$(( NUM_CHRS - 4 ))
    END=$(( NUM_CHRS + 4 ))

    # Prevent negative starting values (if NUM_CHRS < 4)
    if (( START < 1 )); then
        START=1
    fi

    # Loop through values within 4 of NUM_CHRS
    for (( CHR_ATTEMPT=START; CHR_ATTEMPT<=END; CHR_ATTEMPT++ )); do

        echo -e "\e[43m~~~ $(elapsed) HapHiC for ${SAMPLE} with ${NUM_CHRS} chrs failed, trying: ${CHR_ATTEMPT} ~~~\e[0m"
        mkdir -p ${WD}/${SAMPLE}/03_haphic_nopurge/haphic
        cd ${WD}/${SAMPLE}/03_haphic_nopurge/haphic
        
        ${PUZZLER} haphic pipeline ${WD}/${SAMPLE}/02_purge_dups/p_ctg.purged.fa ../filtered.bam ${CHR_ATTEMPT} --correct_nrounds 2 --max_inflation 10.0 --threads ${t} --processes ${t} 2> ../${SAMPLE}.haphic_N_chrs${CHR_ATTEMPT}.log

        # Check if scaffolds.fa exists
        if [ -s ${WD}/${SAMPLE}/03_haphic_nopurge/haphic/04.build/scaffolds.fa ]; then
            echo -e "\e[43m~~~ $(elapsed)  HapHiC completed successfully with NUM_CHRS=${CHR_ATTEMPT} ~~~\e[0m"
            break
        else
            rm -rf ${WD}/${SAMPLE}/03_haphic_nopurge/haphic/*
        fi
    done

    fi 

    # IF this still failed, then just scaffold with YAHS 
    if [ ! -s "${WD}/${SAMPLE}/03_haphic_nopurge/haphic/04.build/scaffolds.fa" ]; then

        echo -e "\e[43m~~~ $(elapsed) HapHiC for ${SAMPLE} failed, scaffolding with YAHS instead ~~~\e[0m"
        cd ${WD}/${SAMPLE}/03_haphic_nopurge
        rm -rf ${WD}/${SAMPLE}/03_haphic_nopurge/haphic/*
        mkdir -p ${WD}/${SAMPLE}/03_haphic_nopurge/haphic/04.build

        ${PUZZLER} samtools faidx ${WD}/${SAMPLE}/02_purge_dups/pri.init.fa
        ${PUZZLER} yahs ${WD}/${SAMPLE}/02_purge_dups/pri.init.fa filtered.bam > yahs.log 2>&1
        cp ${WD}/${SAMPLE}/03_haphic_nopurge/yahs.out_scaffolds_final.fa ${WD}/${SAMPLE}/03_haphic_nopurge/haphic/04.build/scaffolds.fa
        cp ${WD}/${SAMPLE}/03_haphic_nopurge/yahs.out_scaffolds_final.agp ${WD}/${SAMPLE}/03_haphic_nopurge/haphic/04.build/scaffolds.raw.agp
    fi 
    touch ${WD}/${SAMPLE}/scaffolding_nopurge.complete
    set -euo pipefail
fi

##### Summarize duplicate forks 
echo -e "\e[43m~~~ $(elapsed) Summarizing Assembly for ${SAMPLE} ~~~\e[0m"
cd ${WD}

OUTFILE="${SAMPLE}_assembly_fork_results.tsv"
echo -e "Sample\tFork\tSize\tSeqs\tCtgs\tScafN50\tContigN50" > "$OUTFILE"

# Define fork names and corresponding file paths
declare -A paths=(
  [Initial]="${WD}/${SAMPLE}/02_purge_dups/pri.init.fa"
  [Purge]="${WD}/${SAMPLE}/02_purge_dups/p_ctg.purged.fa"
  [Purge_Default]="${WD}/${SAMPLE}/02_purge_dups/p_ctg.purged.DEFAULT.fa"
  [NoPurge_HapHic]="${WD}/${SAMPLE}/03_haphic_nopurge/haphic/04.build/scaffolds.fa"
  [HapHic]="${WD}/${SAMPLE}/03_haphic/haphic/04.build/scaffolds.fa"
  [Final]="${WD}/primary_asm/${SAMPLE}.fa"
)

# Loop through each fork
for fork in Initial Purge Purge_Default NoPurge_HapHic HapHic Final; do
  FILE="${paths[$fork]}"
  STATS=$(${PUZZLER} assembly_stats "$FILE")

  SIZE=$(echo "$STATS" | tr ':' '\n' | grep -A 1 'total_bps' | tail -n1 | sed 's/,//g')
  SEQS=$(echo "$STATS" | tr ':' '\n' | grep -A 1 'sequence' | tail -n1 | sed 's/,//g')
  CTGS=$(echo "$STATS" | tr ':' '\n' | grep -A 1 'sequence' | head -n2 | tail -n1 | sed 's/,//g')
  SCAF_N50=$(echo "$STATS" | tr ':' '\n' | grep -A 1 'N50' | tail -n1 | sed 's/,//g')
  CONT_N50=$(echo "$STATS" | tr ':' '\n' | grep -A 1 'N50' | head -n2 | tail -n1 | sed 's/,//g')

  echo -e "${SAMPLE}\t${fork}\t${SIZE}\t${SEQS}\t${CTGS}\t${SCAF_N50}\t${CONT_N50}" >> "$OUTFILE"
done 
```

## Plot 

Plot:

```R
setwd('/90daydata/coffea_pangenome/puzzler_trials/rerun_assemblies/')
library(viridis)
library(tidyverse)
library(ggpubr)
library(stringr)
library(scales)

md <- read_tsv('/90daydata/coffea_pangenome/puzzler_trials/assemblies/metadata.tsv')
purge <- read_tsv('/90daydata/coffea_pangenome/puzzler_trials/rerun_assemblies/Purge_Duplicates_Trials_NoJuice_20251123.tsv')

pi <- left_join(purge,md %>% select(Sample,Species,Phylo_Order,Heterozygosity)) %>% 
  mutate(CN50 = round(ContigN50/1e6,1),
         SN50 = round(ScafN50/1e6,1),
         SizeGB = round(Size/1e9,2))
cols <- viridis(4)

pi_plot <- pi %>% select(Sample,Fork,Seqs,SizeGB,CN50,SN50)
pis <- left_join(pi_plot,md %>% select(Sample,Species))
pis$Species <- factor(pis$Species,levels=unique(md$Species))
pis$Fork <- factor(pis$Fork,levels=c('Initial','NoPurge_HapHic','Purge','Purge_Default','HapHic','Final','Reference'))

# Estimate fold change from Reference
pis_norm <- pis %>%
  group_by(Sample, Species) %>%
  # get reference values per Sample
  mutate(
    ref_Seqs  = Seqs[ Fork == "Reference" ],
    ref_Size  = SizeGB[ Fork == "Reference" ],
    ref_CN50  = CN50[ Fork == "Reference" ],
    ref_SN50  = SN50[ Fork == "Reference" ],
    # compute % relative to reference
    Seqs_pct  = 100 * Seqs / ref_Seqs,
    Size_pct  = 100 * SizeGB / ref_Size,
    CN50_pct  = 100 * CN50 / ref_CN50,
    SN50_pct  = 100 * SN50 / ref_SN50
  ) %>%
  ungroup()
pis_norm

pis_norm %>% 
  select(matches(c('Sample','Species','Fork','pct'))) %>% 
  pivot_longer(!c(Sample,Species,Fork)) %>% 
  group_by(Fork,name) %>% sum_stats(value) %>% 
  mutate(lab = paste0(round(mean,0),' (',round(min,0),'–',round(max,0),')')) %>% 
  select(Fork,name,lab) %>% 
  pivot_wider(names_from = name,values_from = lab)

# Contig/Scaffold N50
p1 <- pis %>%
  ggplot(aes(y=Fork,x=CN50))+
  geom_col(col='black',width=0.75,position=position_dodge(width=0.75),fill=cols[1])+
  facet_grid(.~Species,scales='free',
             labeller = labeller(Species = function(x) str_wrap(x, width = 10)))+
  scale_x_continuous(breaks = pretty_breaks(n = 2))+
  xlab('Contig N50 (Mb)')+ylab('')+
  ggtitle('Contig N50')+
  theme_bw(base_size=8)
p1
ggsave('20251124_Purge_Trials_CN50.pdf',p1,height=2.5,width=9)

# Scaffold N50
p2 <- pis %>%
  ggplot(aes(y=Fork,x=SN50))+
  geom_col(col='black',width=0.75,position=position_dodge(width=0.75),fill=cols[2])+
  facet_grid(.~Species,scales='free',
             labeller = labeller(Species = function(x) str_wrap(x, width = 10)))+
  scale_x_continuous(breaks = pretty_breaks(n = 2))+
  xlab('Scaffold N50 (Mb)')+ylab('')+
  ggtitle('Scaffold N50')+
  theme_bw(base_size=8)
p2
ggsave('20251124_Purge_Trials_SN50.pdf',p2,height=2.5,width=9)

# Number of sequences 
p3 <- pis %>%
  ggplot(aes(y=Fork,x=Seqs))+
  geom_col(col='black',width=0.75,position=position_dodge(width=0.75),fill=cols[3])+
  facet_grid(.~Species,scales='free',
             labeller = labeller(Species = function(x) str_wrap(x, width = 10)))+
  ggtitle('Sequences')+
  xlab('Sequences (N)')+ylab('')+
  scale_x_continuous(breaks = pretty_breaks(n = 2))+
  theme_bw(base_size=8)
p3
ggsave('20251124_Purge_Trials_Seqs.pdf',p3,height=2.5,width=9)

# Size assembly
p4 <- pis %>%
  ggplot(aes(y=Fork,x=SizeGB))+
  geom_col(col='black',width=0.75,position=position_dodge(width=0.75),fill=cols[4])+
  facet_grid(.~Species,scales='free',
             labeller = labeller(Species = function(x) str_wrap(x, width = 10)))+
  ggtitle('Assembly Size')+
  xlab('Size (Gb)')+ylab('')+
  scale_x_continuous(breaks = pretty_breaks(n = 2))+
  theme_bw(base_size=8)
p4
ggsave('20251124_Purge_Trials_Size.pdf',p4,height=2.5,width=9)

```





| sample | runtime   | container                                        | wd                                                    | hifi                                                         | hic_r1                                                       | hic_r2                                                       | num_chrs | reference                                                    | hom_cov | blob_database                                             | busco_lineage  | busco_database                                             |
| ------ | --------- | ------------------------------------------------ | ----------------------------------------------------- | ------------------------------------------------------------ | ------------------------------------------------------------ | ------------------------------------------------------------ | -------- | ------------------------------------------------------------ | ------- | --------------------------------------------------------- | -------------- | ---------------------------------------------------------- |
| Fly    | conda     | NA                                               | /90daydata/coffea_pangenome/puzzler_trials/assemblies | /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Fly.HiFi.fastq.gz | /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Fly.HiC.R1.fastq.gz | /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Fly.HiC.R2.fastq.gz | 7        | /90daydata/coffea_pangenome/puzzler_trials/raw_data/references/GCF_028408465.1_idAnaLude1.1_genomic.fna | 44      | /90daydata/coffea_pangenome/puzzler_trials/blob_downloads | diptera_odb10  | /90daydata/coffea_pangenome/puzzler_trials/busco_downloads |
| Beaver | apptainer | /home/justin.merondun/apptainer/puzzler_v1.8.sif | /90daydata/coffea_pangenome/puzzler_trials/assemblies | /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Beaver.HiFi.fastq.gz | /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Beaver.HiC.R1.fastq.gz | /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Beaver.HiC.R2.fastq.gz | 20       | /90daydata/coffea_pangenome/puzzler_trials/raw_data/references/GCF_047511655.1_mCasCan1.hap1v2_genomic.fna | NA      | NA                                                        | mammalia_odb10 | /90daydata/coffea_pangenome/puzzler_trials/busco_downloads |
| Fungus | apptainer | /home/justin.merondun/apptainer/puzzler_v1.8.sif | /90daydata/coffea_pangenome/puzzler_trials/assemblies | /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Fungus.HiFi.fastq.gz | /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Fungus.HiC.R1.fastq.gz | /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Fungus.HiC.R2.fastq.gz | 14       | NA                                                           | NA      | NA                                                        | fungi_odb10    | /90daydata/coffea_pangenome/puzzler_trials/busco_downloads |
