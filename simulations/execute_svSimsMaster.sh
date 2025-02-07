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

names=($(cat ../preprocess/chromos.txt))
files=${names[${SLURM_ARRAY_TASK_ID}]}

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

# Running on 20220422_3202_phased_SNV_INDEL_SV_bychrom dataset = 2024/07/25
# sbatch execute_svSimsMasterINDPAR1_240716.sh 20220422_3202_phased_SNV_INDEL_SV_bychrom SVTrue_typedeletion_resTrue SVTrue_typedeletion_resTrue.csv



