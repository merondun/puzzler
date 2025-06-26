# Genome Assembly for N=11 species with **`puzzler`**

I will create Puzzler assemblies for 11 species which have HiC and HiFi data from the NCBI:

| Sample      | Species                     | Species Name                 | Bioproject                                                   | Link                                                         |
| ----------- | --------------------------- | ---------------------------- | ------------------------------------------------------------ | ------------------------------------------------------------ |
| Rhodie      | Rhododendron                | Rhododendron vialii          | PRJNA971245                                                  | https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_030253575.1/ |
| Fungus      | Rhizosphaera needle cast    | Rhizosphaera kalkhoffii      | PRJEB76005                                                   | https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_964106605.1/ |
| Fly         | Mexican Fruit Fly           | Anastrepha ludens            | PRJNA803324                                                  | https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_028408465.1/ |
| Squirt      | Solitary Sea Squirt         | Ascidia mentula              | PRJEB58134                                                   | https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_947561715.1/ |
| Fish        | Antarctic lanternfish       | Electrona antarctica         | PRJEB61834                                                   | https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_951216825.1/ |
| Stickleback | Ninespine Stickleback       | Pungitius pungitius          | PRJEB60134                                                   | https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_949316345.1/ |
| Toad        | Eastern Narrow-mouthed Toad | Gastrophryne carolinensis    | PRJNA923362                                                  | https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_027917425.1/ |
| Frog        | Guiana Rocket Frog          | Anomaloglossus baeobatrachus | PRJNA1226914                                                 | https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_048569485.1/ |
| Crane       | Eurasian Crane              | Grus grus                    | PRJEB75999 | https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_964106855.1/ |
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

# Create Drafts -> Juicer

I will skip blobtools because that will double to triple the assembly time. I recommend running this (just set `blob_db` to a directory, but expect long runtimes! 

Here is my `samples.tsv`:

| sample      | runtime   | container                                        | wd                                                    | hifi                                                         | hic_r1                                                       | hic_r2                                                       | num_chrs | reference                                                    | hom_cov | blob_database | busco_lineage        | busco_database                                             |
| ----------- | --------- | ------------------------------------------------ | ----------------------------------------------------- | ------------------------------------------------------------ | ------------------------------------------------------------ | ------------------------------------------------------------ | -------- | ------------------------------------------------------------ | ------- | ------------- | -------------------- | ---------------------------------------------------------- |
| Beaver      | apptainer | /home/justin.merondun/apptainer/puzzler_v1.8.sif | /90daydata/coffea_pangenome/puzzler_trials/assemblies | /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Beaver.HiFi.fastq.gz | /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Beaver.HiC.R1.fastq.gz | /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Beaver.HiC.R2.fastq.gz | 20       | /90daydata/coffea_pangenome/puzzler_trials/raw_data/references/GCF_047511655.1_mCasCan1.hap1v2_genomic.fna | NA      | NA            | mammalia_odb10       | /90daydata/coffea_pangenome/puzzler_trials/busco_downloads |
| Crane       | apptainer | /home/justin.merondun/apptainer/puzzler_v1.8.sif | /90daydata/coffea_pangenome/puzzler_trials/assemblies | /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Crane.HiFi.fastq.gz | /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Crane.HiC.R1.fastq.gz | /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Crane.HiC.R2.fastq.gz | 40       | /90daydata/coffea_pangenome/puzzler_trials/raw_data/references/GCA_964106855.1_bGruGru1.hap1.1_genomic.fna | 48      | NA            | aves_odb10           | /90daydata/coffea_pangenome/puzzler_trials/busco_downloads |
| Fish        | apptainer | /home/justin.merondun/apptainer/puzzler_v1.8.sif | /90daydata/coffea_pangenome/puzzler_trials/assemblies | /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Fish.HiFi.fastq.gz | /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Fish.HiC.R1.fastq.gz | /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Fish.HiC.R2.fastq.gz | 24       | /90daydata/coffea_pangenome/puzzler_trials/raw_data/references/GCA_951216825.1_fEleAnt2.1_genomic.fna | 42      | NA            | actinopterygii_odb10 | /90daydata/coffea_pangenome/puzzler_trials/busco_downloads |
| Fly         | apptainer | /home/justin.merondun/apptainer/puzzler_v1.8.sif | /90daydata/coffea_pangenome/puzzler_trials/assemblies | /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Fly.HiFi.fastq.gz | /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Fly.HiC.R1.fastq.gz | /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Fly.HiC.R2.fastq.gz | 7        | /90daydata/coffea_pangenome/puzzler_trials/raw_data/references/GCF_028408465.1_idAnaLude1.1_genomic.fna | 44      | NA            | diptera_odb10        | /90daydata/coffea_pangenome/puzzler_trials/busco_downloads |
| Frog        | apptainer | /home/justin.merondun/apptainer/puzzler_v1.8.sif | /90daydata/coffea_pangenome/puzzler_trials/assemblies | /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Frog.HiFi.fastq.gz | /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Frog.HiC.R1.fastq.gz | /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Frog.HiC.R2.fastq.gz | 12       | /90daydata/coffea_pangenome/puzzler_trials/raw_data/references/GCA_048569485.1_aAnoBae1.hap1_genomic.fna | NA      | NA            | tetrapoda_odb10      | /90daydata/coffea_pangenome/puzzler_trials/busco_downloads |
| Fungus      | apptainer | /home/justin.merondun/apptainer/puzzler_v1.8.sif | /90daydata/coffea_pangenome/puzzler_trials/assemblies | /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Fungus.HiFi.fastq.gz | /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Fungus.HiC.R1.fastq.gz | /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Fungus.HiC.R2.fastq.gz | 14       | /90daydata/coffea_pangenome/puzzler_trials/raw_data/references/GCA_964106605.1_gdRhiKalk1.hap1.1_genomic.fna | NA      | NA            | fungi_odb10          | /90daydata/coffea_pangenome/puzzler_trials/busco_downloads |
| Rhodie      | apptainer | /home/justin.merondun/apptainer/puzzler_v1.8.sif | /90daydata/coffea_pangenome/puzzler_trials/assemblies | /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Rhodie.HiFi.fastq.gz | /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Rhodie.HiC.R1.fastq.gz | /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Rhodie.HiC.R2.fastq.gz | 13       | /90daydata/coffea_pangenome/puzzler_trials/raw_data/references/GCF_030253575.1_ASM3025357v1_genomic.fna | 60      | NA            | embryophyta_odb10    | /90daydata/coffea_pangenome/puzzler_trials/busco_downloads |
| Squirt      | apptainer | /home/justin.merondun/apptainer/puzzler_v1.8.sif | /90daydata/coffea_pangenome/puzzler_trials/assemblies | /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Squirt.HiFi.fastq.gz | /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Squirt.HiC.R1.fastq.gz | /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Squirt.HiC.R2.fastq.gz | 9        | /90daydata/coffea_pangenome/puzzler_trials/raw_data/references/GCA_947561715.1_kaAscMent1.1_genomic.fna | 74      | NA            | vertebrata_odb10     | /90daydata/coffea_pangenome/puzzler_trials/busco_downloads |
| Stickleback | apptainer | /home/justin.merondun/apptainer/puzzler_v1.8.sif | /90daydata/coffea_pangenome/puzzler_trials/assemblies | /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Stickleback.HiFi.fastq.gz | /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Stickleback.HiC.R1.fastq.gz | /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Stickleback.HiC.R2.fastq.gz | 21       | /90daydata/coffea_pangenome/puzzler_trials/raw_data/references/GCF_949316345.1_fPunPun2.1_genomic.fna | NA      | NA            | actinopterygii_odb10 | /90daydata/coffea_pangenome/puzzler_trials/busco_downloads |
| Toad        | apptainer | /home/justin.merondun/apptainer/puzzler_v1.8.sif | /90daydata/coffea_pangenome/puzzler_trials/assemblies | /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Toad.HiFi.fastq.gz | /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Toad.HiC.R1.fastq.gz | /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Toad.HiC.R2.fastq.gz | 11       | /90daydata/coffea_pangenome/puzzler_trials/raw_data/references/GCA_027917425.1_aGasCar1.pri_genomic.fna | 24      | NA            | tetrapoda_odb10      | /90daydata/coffea_pangenome/puzzler_trials/busco_downloads |
| Whale       | apptainer | /home/justin.merondun/apptainer/puzzler_v1.8.sif | /90daydata/coffea_pangenome/puzzler_trials/assemblies | /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Whale.HiFi.fastq.gz | /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Whale.HiC.R1.fastq.gz | /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Whale.HiC.R2.fastq.gz | 21       | /90daydata/coffea_pangenome/puzzler_trials/raw_data/references/GCA_964341445.1_mMesMir1.hap1.1_genomic.fna | 54      | NA            | mammalia_odb10       | /90daydata/coffea_pangenome/puzzler_trials/busco_downloads |

You can then submit an assembly for each. First I recommend ensuring the command looks right: 

```bash
cut -f1 samples.tsv | sed '1d' | xargs -I {} echo sbatch -J asm_{} puzzler -s {} -m samples.tsv 
```

And then just remove the `echo` to submit in parallel. 

All samples will run until juicer file creation:

```
cat slurm-17666544.out

=======================================================================
__________ ____ _______________________.____     _____________________
\______   \    |   \____    /\____    /|    |    \_   _____/\______   \
 |     ___/    |   / /     /   /     / |    |     |    __)_  |       _/
 |    |   |    |  / /     /_  /     /_ |    |___  |        \ |    |   \
 |____|   |______/ /_______ \/_______ \|_______ \/_______  / |____|_  /
                           \/        \/        \/        \/         \/
=======================================================================

=======================================================================
Parameters for sample: Stickleback
RUNTIME: apptainer
CONTAINER: /home/justin.merondun/apptainer/puzzler_v1.8.sif
WD: /project/90daydata/coffea_pangenome/puzzler_trials/assemblies
HIFI: /project/90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Stickleback.HiFi.fastq.gz
HIC_R1: /project/90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Stickleback.HiC.R1.fastq.gz
HIC_R2: /project/90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads/Stickleback.HiC.R2.fastq.gz
NUMBER CHRS: 21
REFERENCE: /project/90daydata/coffea_pangenome/puzzler_trials/raw_data/references/GCF_949316345.1_fPunPun2.1_genomic.fna
HOM_COV: NA
BLOB_DB: /project/90daydata/coffea_pangenome/puzzler_trials/blob_downloads
BUSCO_LINEAGE: actinopterygii_odb10
BUSCO_DB: /project/90daydata/coffea_pangenome/puzzler_trials/busco_downloads
Cores Requested: 32
Cores Available: 32
RAM Requested: 256
Memory Available: 362.8 GB
PUZZLER command: apptainer exec --bind /90daydata/coffea_pangenome/puzzler_trials:/90daydata/coffea_pangenome/puzzler_trials --bind /90daydata/coffea_pangenome/puzzler_trials/assemblies:/90daydata/coffea_pangenome/puzzler_trials/assemblies --bind /90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads:/90daydata/coffea_pangenome/puzzler_trials/raw_data/concat_reads --bind /90daydata/coffea_pangenome/puzzler_trials/raw_data/references:/90daydata/coffea_pangenome/puzzler_trials/raw_data/references --bind /home/justin.merondun/apptainer:/home/justin.merondun/apptainer /home/justin.merondun/apptainer/puzzler_v1.8.sif
=======================================================================

~~~~ Assembling genome for Stickleback ~~~~
~~~~ Running juicer, script will stop after .hic files created ~~~~
~~~~ Starting hifiasm assembly for Stickleback ~~~~
~~~~ Starting Purge_Dups for Stickleback ~~~~
~~~~ Mapping HiC reads to Stickleback draft ~~~~
~~~~ Running HapHiC for Stickleback  ~~~~
~~~~ Creating .hic file for juicebox for Stickleback ~~~~
~~~~ Genome size is 0GB, running default minimap2 command ~~~~
~~~~ Post curation assembly file missing for Stickleback: Run Juicebox & place in /project/90daydata/coffea_pangenome/puzzler_trials/assemblies/juicer_files/Stickleback_JBAT.review.assembly ~~~~
```

So at this time, you will need to do juicebox manual curation. 

# Juicer Curation

Manual curation files can be found either on Zenodo or on [Youtube](https://www.youtube.com/watch?v=rMUiNqZwEpA). 

Save the Juicebox curation files as default, with `Assembly > Export Assembly`, which will generate e.g. `Stickleback_JBAT.review.assembly`. Pull that back into `$WD/juicer_files/`. 

# Post Juicebox Extraction & QC

Simply resubmit `puzzler`: 

```bash
cut -f1 samples.tsv | sed '1d' | xargs -I {} echo sbatch -J asm_{} puzzler -s {} -m samples.tsv 
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
sed '1d' samples_ref.tsv | cut -f1 | xargs -I {} sbatch -J ref_quality_{} puzzler -s {} -m samples_ref.tsv
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

```R
#### Summarize assemblies
setwd('/90daydata/coffea_pangenome/puzzler_trials/assemblies/primary_asm/stats/')
library(tidyverse)
library(RColorBrewer)
library(ggpubr)
library(meRo) #devtools::install_github('merondun/meRo')
library(gt)
library(Biostrings)
library(meRo)

# Import metadata / order
md <- read_tsv('/90daydata/coffea_pangenome/puzzler_trials/assemblies/metadata.tsv')
md <- md %>% mutate(hr_min = str_glue("{puzzler_no_juice %/% 60}h:{str_pad(puzzler_no_juice %% 60, 2, pad = '0')}m"))
label <- md %>% mutate(lab = paste0(Species,': ',`Species Name`,' ',hr_min)) %>% select(Sample,lab,Phylo_Order) %>% arrange(Phylo_Order)

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
    length = width(fa)) 
  gc <- gapless_chrs %>% mutate(Sample = paste0(sp,'_Ref'))
  dats <- rbind(dats,gc)
}

# Match the reference and puzzler chromosomes
compiled <- dats %>% filter(!grepl('scaff|organ|mito',chr)) %>% filter(!grepl("[A-W]", chr)) %>% 
  mutate(
    Species = gsub("_Ref", "", Sample),
    ref_type = ifelse(grepl("_Ref", Sample), "Ref", "")
  ) %>%
  select(chr, Species, ref_type, has_N, length) %>%
  pivot_wider(
    names_from = ref_type,
    values_from = c(has_N, length),
    names_glue = "{.value}{ifelse(ref_type == 'Ref', '_Ref', '')}") 
write.table(compiled,'Gapless_Summaries.txt',quote=F,sep='\t',row.names=F)
compiled <- read_tsv('Gapless_Summaries.txt')

# First, set Puzzler chrs to gapped if they are less than 3 Mb from Reference size, then count number of gapless and gapped for puzzler and ref 
compstats <- compiled %>% 
  mutate(dlength = (length-length_Ref)/1e6,
         has_N = ifelse(abs(dlength) > 1,TRUE,has_N)) %>% 
  pivot_longer(c(has_N,has_N_Ref)) %>% 
  replace_na(list(value=TRUE)) %>% 
  group_by(Species,name) %>% 
  count(value) %>% ungroup %>% 
  complete(Species,name,value, fill = list(n = 0)) %>% 
  mutate(Assembly = ifelse(grepl('Ref',name),'Reference','Puzzler'),
         value = gsub('TRUE','Gapped',gsub('FALSE','Gapless',value)))

# Plot 
gap_plot <- compstats %>% 
  ggplot(aes(y=Species,x=n,fill=value))+
  geom_col(width=1,col='black',position=position_dodge(width=0.7))+
  geom_text(aes(label=n),position=position_dodge(width=1),hjust=-0.5)+
  scale_fill_manual('Chromosome',values=c('black','white'))+
  facet_grid(.~Assembly,scales='free')+
  coord_cartesian(xlim=c(0,43))+
  xlab('Number of Gapless Chromosomes\nRequiring That Puzzler Chrs Are ± 1 Mb From Reference Size')+
  theme_bw()
gap_plot
ggsave('GapCounts.pdf',gap_plot,height=6,width=7)

# Add metadata information
s <- read_tsv('Summaries.tsv')
sl <- s %>% 
  select(Sample,SizeBP,ContigN50,ScafN50,BUSCO_Complete) %>% 
  mutate(Species = gsub('_Ref','',Sample),
         Assembly = ifelse(grepl('_Ref',Sample),'Reference','Puzzler'),
         ScafN50 = ScafN50 / 1e6,
         ContigN50 = ContigN50 / 1e6,
         SizeBP = SizeBP / 1e6,) %>% 
  pivot_longer(!c(Sample,Species,Assembly))
levels <- unique(sl$name) 
sl$name <- factor(sl$name,levels=levels)
slf <- left_join(sl,label %>% dplyr::rename(Species=Sample)) %>% arrange(Phylo_Order)

labs <- slf %>% 
  group_by(Species,name) %>% 
  mutate(max_x = max(value),
         diff = signif(value[Assembly == 'Puzzler'] - value[Assembly == 'Reference'],2))
labs %>% ungroup %>% group_by(name) %>% sum_stats(diff)
# # A tibble: 4 × 8
# name            mean      sd     se median    iqr conf_low conf_high
# <fct>          <dbl>   <dbl>  <dbl>  <dbl>  <dbl>    <dbl>     <dbl>
#   1 SizeBP          8.79 113.    24.2     2.6  75.6   -41.4       59.0  
# 2 ContigN50       5.13  24.9    5.31    0.4   6.72   -5.91      16.2  
# 3 ScafN50         7.19  13.0    2.76    0.92  6.45    1.45      12.9  
# 4 BUSCO_Complete  0.2    0.531  0.113   0     0.175  -0.0354     0.435

#labs$Species = factor(labs$Species,levels=unique(slf$Species))
slf$Species = factor(slf$Species,levels=rev(unique(slf$Species)))

# Plot
all <- slf %>%
  ggplot(aes(y = Species, fill = Assembly, x = value)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), col = "grey50") +
  geom_text(data = labs, aes(x = max_x *1.2, y = Species, label = diff), hjust = 0.5,size=1) +  # Adjust label position dynamically
  scale_fill_manual(values = c("cadetblue", "goldenrod1")) +
  facet_wrap(name ~ ., scales = "free_x",nrow=1) +  # Allows individual x-axis scaling per facet
  theme_bw() +
  theme(
    strip.text = element_text(size = 12),  # Increase facet label size
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels for readability
    panel.spacing = unit(2, "lines")  # Increase spacing between facets
  )
all
ggsave('Assembly_Summaries.pdf',all,height=2.5,width=10)

# Summarize assembly time vs inputs (48 cores)
tl <- left_join(md,s) %>% 
  dplyr::rename(hifi = `HiFi Bases`, hic = `HiC Bases`, length = `HiFi Length`) %>% 
  mutate(SizeBP = SizeBP,
         time_gb = (SizeBP/1e6) / puzzler_no_juice,
         hifi_bp = hifi / SizeBP,
         hic_bp = hic / SizeBP,
         ContigN50 = ContigN50 / 1e6,
         length = length / 1e3) %>% 
  select(Sample, time_gb, hifi_bp, hic_bp, SizeBP, hifi, hic,ContigN50, length)
tl
tl$Sample = factor(tl$Sample,levels=unique(slf$Sample))

# COVERAGE: Define scaling factor
scale_factor <- max(tl$time_gb) / max(tl$hifi_bp, tl$hic_bp)
tl <- tl %>%
  mutate(hifi_scaled = hifi_bp * scale_factor,
         hic_scaled = hic_bp * scale_factor)
inputs <- ggplot(tl, aes(x = Sample)) +
  geom_point(aes(y = time_gb,size=ContigN50), color = "black") +
  geom_line(aes(y = hifi_scaled, group = 1), color = "dodgerblue", linetype = "dashed") +
  geom_line(aes(y = hic_scaled, group = 1), color = "orange", linetype = "dotted") +
  scale_y_continuous(
    name = "Assembly Speed (Minutes per Mb)",
    sec.axis = sec_axis(~./scale_factor, name = "Coverage (Input Bases / Assembly Size)", breaks = scales::pretty_breaks(5))) +
  xlab("") +
  theme_bw(base_size = 9) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.y.right = element_text(color = "grey20"),
    axis.text.y.right = element_text(color = "grey20"))

# LENGTH: Define scaling factor
scale_factor <- max(tl$time_gb) / max(tl$length)
tl <- tl %>%
  mutate(length_scaled = length * scale_factor)
lengths <- ggplot(tl, aes(x = Sample)) +
  geom_point(aes(y = time_gb,size=ContigN50), color = "black") +
  geom_line(aes(y = length_scaled, group = 1), color = "seagreen", linetype = "dashed") +
  scale_y_continuous(
    name = "Assembly Speed (Minutes per Mb)",
    sec.axis = sec_axis(~./scale_factor, name = "Coverage (Input Bases / Assembly Size)", breaks = scales::pretty_breaks(5))) +
  xlab("") +
  theme_bw(base_size = 9) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.y.right = element_text(color = "grey20"),
    axis.text.y.right = element_text(color = "grey20"))
lengths

# INPUTS x LENGTH: Define scaling factor
scale_factor <- max(tl$length) / max(tl$hifi_bp, tl$hic_bp)
tl <- tl %>%
  mutate(hifi_scaled = hifi_bp * scale_factor,
         hic_scaled = hic_bp * scale_factor)
inputs <- ggplot(tl, aes(x = Sample)) +
  geom_point(aes(y = length,size=ContigN50), color = "black") +
  geom_point(aes(y = hifi_scaled, group = 1), color = "dodgerblue", pch=4) +
  geom_point(aes(y = hic_scaled, group = 1), color = "orange", pch=8) +
  #geom_line(aes(y = hifi_scaled, group = 1), color = "dodgerblue", linetype = "dashed") +
  #geom_line(aes(y = hic_scaled, group = 1), color = "orange", linetype = "dotted") +
  scale_y_continuous(
    name = "HiFi Read Length (Kb)",
    sec.axis = sec_axis(~./scale_factor, name = "Coverage (Input Bases / Assembly Size)", breaks = scales::pretty_breaks(5))) +
  xlab("") +
  theme_bw(base_size = 9) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.y.right = element_text(color = "grey20"),
    axis.text.y.right = element_text(color = "grey20"))
inputs

# Only plot inputs
tl %>% 
  select(!c(time_gb,SizeBP,hifi,hic)) %>% 
  pivot_longer(c(hifi_bp,hic_bp,length)) %>% 
  ggplot(aes(y=ContigN50,x=value,col=Sample))+
  geom_point()+
  facet_grid(.~name,scales='free')+
  theme_bw()

ggsave('InputsPlot.pdf',inputs,height=2,width=4.5)

tl %>% sum_stats(time_gb)
# # A tibble: 1 × 7
# mean    sd    se median   iqr conf_low conf_high
# <dbl> <dbl> <dbl>  <dbl> <dbl>    <dbl>     <dbl>
# 1.86  1.11 0.335   1.63  1.59     1.12      2.61


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
md <- md %>% mutate(hr_min = str_glue("{puzzler_no_juice %/% 60}h:{str_pad(puzzler_no_juice %% 60, 2, pad = '0')}m"))
label <- md %>% mutate(lab = paste0(Species,': ',`Species Name`,' (',hr_min,')')) %>% select(Sample,lab,Phylo_Order) %>% arrange(Phylo_Order)

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
png('../20250623_Dotplots-NoHashes-NoJuicer.png',height=12,width=17,res=300,units='in')
ggarrange(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11)
dev.off()

```
