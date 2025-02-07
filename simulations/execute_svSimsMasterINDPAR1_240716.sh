#!/bin/bash
#SBATCH --nodes=1
#SBATCH --mem=10G
#SBATCH --time=02:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=nb.boev@mail.utoronto.ca
#SBATCH --job-name=MastersvSim
#SBATCH --output=MastersvSim%j.o
#SBATCH --error=MastersvSim%j.e
#SBATCH --array=0-25

# in this script, we need to push each chromosome through job arrays
module --force purge ; module load StdEnv/2020
source $HOME/virtenv/bin/activate
module load gcc/9.3.0 r/4.0.2
module load scipy-stack

# Use this version to allow for scratch space to be used instead
# dirs that need to be created before running:
# in scratch:/home/nboev/scratch/data/SimulatedSVs/CPTAC3/bf584e4c-f981-473e-b40a-c73d7f3699d2.wgs.sanger_raw_pindel.raw_somatic_mutation with filenamesIND  parameterIND  rawIND dirs
# in the sim dir: /home/nboev/projects/def-sushant/nboev/data/SimulatedSVs/CPTAC3/bf584e4c-f981-473e-b40a-c73d7f3699d2.wgs.sanger_raw_pindel.raw_somatic_mutation with logIND  processedIND dirs

#names=($(cat /home/nboev/projects/def-sushant/nboev/data/SimulatedSVs/20220422_3202_phased_SNV_INDEL_SV_bychrom/chromoX3.txt))
names=($(cat /home/nboev/projects/def-sushant/nboev/data/SimulatedSVs/HGSVC2_v2_integrated_callset/chromolong.txt))
#names=($(cat /home/nboev/projects/def-sushant/nboev/preprocess/HGSVC2_v2_integrated_callset/chromosXY123.txt))
#names=($(cat ./chromosontmissing.txt))
#names=($(cat ./chromos5.txt))
files=${names[${SLURM_ARRAY_TASK_ID}]}
echo "${files}"

mkdir /home/nboev/scratch/data/SimulatedSVs/$1
mkdir /home/nboev/scratch/data/SimulatedSVs/$1/$2
mkdir /home/nboev/scratch/data/SimulatedSVs/$1/$2/parameterIND
mkdir /home/nboev/scratch/data/SimulatedSVs/$1/$2/filenamesIND
mkdir /home/nboev/scratch/data/SimulatedSVs/$1/$2/rawIND
mkdir /home/nboev/projects/def-sushant/nboev/data/SimulatedSVs/$1
mkdir /home/nboev/projects/def-sushant/nboev/data/SimulatedSVs/$1/$2
mkdir /home/nboev/projects/def-sushant/nboev/data/SimulatedSVs/$1/$2/processedIND
mkdir /home/nboev/scratch/data/SimulatedSVs/$1/$2/parameterIND/${files}
mkdir /home/nboev/scratch/data/SimulatedSVs/$1/$2/rawIND/${files}
mkdir /home/nboev/projects/def-sushant/nboev/data/SimulatedSVs/$1/$2/processedIND/${files}

python generating_lensIND_240424.py \
/home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/SVlen/$3 \
${files} \
parameter_file \
$1 $2 \
>> /home/nboev/projects/def-sushant/nboev/data/SimulatedSVs/$1/$2/generating_lensIND_240424.py.txt

ls /home/nboev/scratch/data/SimulatedSVs/$1/$2/parameterIND/${files}/ \
> /home/nboev/scratch/data/SimulatedSVs/$1/$2/filenamesIND/filenames.${files}.txt

split -l 250 -d --additional-suffix=.txt \
/home/nboev/scratch/data/SimulatedSVs/$1/$2/filenamesIND/filenames.${files}.txt \
/home/nboev/scratch/data/SimulatedSVs/$1/$2/rawIND/${files}/${files}.split

splits=$(ls -1 /home/nboev/scratch/data/SimulatedSVs/$1/$2/rawIND/${files}/${files}.split* | wc -l)
echo "${splits}"

ls /home/nboev/scratch/data/SimulatedSVs/$1/$2/rawIND/${files}/${files}.split* \
> /home/nboev/scratch/data/SimulatedSVs/$1/$2/rawIND/${files}/splitnames.txt

cp execute_SurvivorTemplate_240716.sh \
/home/nboev/scratch/data/SimulatedSVs/$1/$2/rawIND/${files}/execute_Survivor_${files}.sh

sed -i "s/PROJECT/$1/g" /home/nboev/scratch/data/SimulatedSVs/$1/$2/rawIND/${files}/execute_Survivor_${files}.sh
sed -i "s/LOCATION/$2/g" /home/nboev/scratch/data/SimulatedSVs/$1/$2/rawIND/${files}/execute_Survivor_${files}.sh
sed -i "s/array=0-0/array=0-${splits}/g" /home/nboev/scratch/data/SimulatedSVs/$1/$2/rawIND/${files}/execute_Survivor_${files}.sh
sed -i "s/CHROM/${files}/g" /home/nboev/scratch/data/SimulatedSVs/$1/$2/rawIND/${files}/execute_Survivor_${files}.sh

sbatch /home/nboev/scratch/data/SimulatedSVs/$1/$2/rawIND/${files}/execute_Survivor_${files}.sh $1 $2


# How to run: change this location!
# Going to start with CPTAC files = 2024/07/16
# sbatch execute_svSimsMasterINDPAR1_240716.sh CPTAC3 bf584e4c-f981-473e-b40a-c73d7f3699d2.wgs.sanger_raw_pindel.raw_somatic_mutation bf584e4c-f981-473e-b40a-c73d7f3699d2.wgs.sanger_raw_pindel.raw_somatic_mutation.tsv
# sbatch execute_svSimsMasterINDPAR1_240716.sh CPTAC3 fcc8a59b-19b0-4b78-9ca1-f5b9e580f183.wgs.sanger_raw_pindel.raw_somatic_mutation fcc8a59b-19b0-4b78-9ca1-f5b9e580f183.wgs.sanger_raw_pindel.raw_somatic_mutation.tsv

# Running on the actual ONT dataset itself = 2024/07/22
# sbatch execute_svSimsMasterINDPAR1_240716.sh ONT_1000Gs_100FREEZE sv merged_GT.SVlength.csvfilt.csv
# Updated file = 2024/08/03
# sbatch execute_svSimsMasterINDPAR1_240716.sh ONT_1000Gs_100FREEZE sv 20240423_jasmine_intrasample_noBND_custom_suppvec_alphanumeric_header_JASMINE.SVlengthFILT.csv

# Running on the Belyeu2021 dataset = 2024/07/24
# sbatch execute_svSimsMasterINDPAR1_240716.sh Belyeu2021 sv mmc2.SVlength.csv

# Running on 20220422_3202_phased_SNV_INDEL_SV_bychrom dataset = 2024/07/25
# sbatch execute_svSimsMasterINDPAR1_240716.sh 20220422_3202_phased_SNV_INDEL_SV_bychrom SVTrue_typedeletion_resTrue SVTrue_typedeletion_resTrue.csv

# Getting DepMap sims = 2024/10/20
# sbatch execute_svSimsMasterINDPAR1_240716.sh DepMap CA922 CA922_allSV.csv
# sbatch execute_svSimsMasterINDPAR1_240716.sh DepMap DOTC24510 DOTC24510_allSV.csv
# sbatch execute_svSimsMasterINDPAR1_240716.sh DepMap FU97 FU97_allSV.csv
# sbatch execute_svSimsMasterINDPAR1_240716.sh DepMap HCC1428 HCC1428_allSV.csv
# sbatch execute_svSimsMasterINDPAR1_240716.sh DepMap HSKTC HSKTC_allSV.csv
# sbatch execute_svSimsMasterINDPAR1_240716.sh DepMap KMCH1 KMCH1_allSV.csv
# sbatch execute_svSimsMasterINDPAR1_240716.sh DepMap KP363T KP363T_allSV.csv
# sbatch execute_svSimsMasterINDPAR1_240716.sh DepMap KYSE30 KYSE30_allSV.csv
# sbatch execute_svSimsMasterINDPAR1_240716.sh DepMap KYSE410 KYSE410_allSV.csv
# sbatch execute_svSimsMasterINDPAR1_240716.sh DepMap KYSE450 KYSE450_allSV.csv
# sbatch execute_svSimsMasterINDPAR1_240716.sh DepMap NCIH748 NCIH748_allSV.csv
# sbatch execute_svSimsMasterINDPAR1_240716.sh DepMap SH10TC SH10TC_allSV.csv
# sbatch execute_svSimsMasterINDPAR1_240716.sh DepMap SNU119 SNU119_allSV.csv
# sbatch execute_svSimsMasterINDPAR1_240716.sh DepMap SNU251 SNU251_allSV.csv
# sbatch execute_svSimsMasterINDPAR1_240716.sh DepMap SUDHL1 SUDHL1_allSV.csv
# sbatch execute_svSimsMasterINDPAR1_240716.sh DepMap UPCISCC074 UPCISCC074_allSV.csv
# sbatch execute_svSimsMasterINDPAR1_240716.sh DepMap VMRCMELG VMRCMELG_allSV.csv


