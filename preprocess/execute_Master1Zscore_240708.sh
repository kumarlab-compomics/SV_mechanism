#!/bin/bash
#SBATCH --nodes=1
#SBATCH --mem=15G
#SBATCH --time=04:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=nb.boev@mail.utoronto.ca
#SBATCH --job-name=MastersZscore
#SBATCH --output=MasterZscore_%j.o
#SBATCH --error=MasterZscore_%j.e
#SBATCH --qos=privilege
#SBATCH --array=0-25

# remember to add back all the jobs later

# in this script, we need to push each chromosome through job arrays
module --force purge ; module load StdEnv/2020
source $HOME/virtenv/bin/activate
module load gcc/9.3.0 r/4.0.2
module load scipy-stack

# for now, we'll use the chromo list in the simulation directory
names=($(cat ./chromos.txt))

files=${names[${SLURM_ARRAY_TASK_ID}]}
echo "${files}"

mkdir /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/zscore

python adding_zscores_240610.py \
/home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/merge_FeatMatrix/${files}.full.csv.DEL.FeatMatrix.csv \
$1 $2 \
20240625.${files}.deletion_FeatMatrix \
${files} \
/home/nboev/projects/def-sushant/nboev/preprocess/$1/$3/merge_FeatMatrix/${files}.full.csv.DEL.FeatMatrix.csv \
>> /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/adding_zscores_240610.py.txt

python adding_zscores_240610.py \
/home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/merge_FeatMatrix/${files}.full.csv.INS.FeatMatrix.csv \
$1 $2 \
20240625.${files}.insertion_FeatMatrix \
${files} \
/home/nboev/projects/def-sushant/nboev/preprocess/$1/$3/merge_FeatMatrix/${files}.full.csv.INS.FeatMatrix.csv \
>> /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/adding_zscores_240610.py.txt


# Running on 1000 genomes project = 2024/09/11:
# sbatch execute_Master1Zscore_240708.sh 20220422_3202_phased_SNV_INDEL_SV_bychrom SVTrue_typedeletion_resTrue SVTrue_typedeletion_resTrueSIM


# Running for cptac validation files = 2024/07/19
# sbatch execute_Master1Zscore_240708.sh CPTAC3 bf584e4c-f981-473e-b40a-c73d7f3699d2.wgs.sanger_raw_pindel.raw_somatic_mutation bf584e4c-f981-473e-b40a-c73d7f3699d2.wgs.sanger_raw_pindel.SIM
# sbatch execute_Master1Zscore_240708.sh CPTAC3 fcc8a59b-19b0-4b78-9ca1-f5b9e580f183.wgs.sanger_raw_pindel.raw_somatic_mutation fcc8a59b-19b0-4b78-9ca1-f5b9e580f183.wgs.sanger_raw_pindel.SIM

# running on hgsvc2 files = 2024/07/23, noticed a problem with chr1
# sbatch execute_Master1Zscore_240708.sh HGSVC2_v2_integrated_callset sv svSim

# Running on Belyeu2021 data set = 2024/07/29
# sbatch execute_Master1Zscore_240708.sh Belyeu2021 sv svSim

# Running on the HGSVC3 dataset = 2024/08/20
# sbatch execute_Master1Zscore_240708.sh allHGSVC3_20240415_Freeze4_GRCh38 sv svSim

# Running on ONT dataset = 2024/09/11
# sbatch execute_Master1Zscore_240708.sh ONT_1000Gs_100FREEZE sv svSim

# Running on Depmap = 2024/10/21
# sbatch execute_Master1Zscore_240708.sh DepMap CA922 CA922svSim
# sbatch execute_Master1Zscore_240708.sh DepMap DOTC24510 DOTC24510svSim
# sbatch execute_Master1Zscore_240708.sh DepMap FU97 FU97svSim
# sbatch execute_Master1Zscore_240708.sh DepMap HCC1428 HCC1428svSim
# sbatch execute_Master1Zscore_240708.sh DepMap HSKTC HSKTCsvSim
# sbatch execute_Master1Zscore_240708.sh DepMap KMCH1 KMCH1svSim
# sbatch execute_Master1Zscore_240708.sh DepMap KP363T KP363TsvSim
# sbatch execute_Master1Zscore_240708.sh DepMap KYSE30 KYSE30svSim
# sbatch execute_Master1Zscore_240708.sh DepMap KYSE410 KYSE410svSim
# sbatch execute_Master1Zscore_240708.sh DepMap KYSE450 KYSE450svSim
# sbatch execute_Master1Zscore_240708.sh DepMap NCIH748 NCIH748svSim
# sbatch execute_Master1Zscore_240708.sh DepMap SH10TC SH10TCsvSim
# sbatch execute_Master1Zscore_240708.sh DepMap SNU119 SNU119svSim
# sbatch execute_Master1Zscore_240708.sh DepMap SNU251 SNU251svSim
# sbatch execute_Master1Zscore_240708.sh DepMap SUDHL1 SUDHL1svSim
# sbatch execute_Master1Zscore_240708.sh DepMap UPCISCC074 UPCISCC074svSim
# sbatch execute_Master1Zscore_240708.sh DepMap VMRCMELG VMRCMELGsvSim




