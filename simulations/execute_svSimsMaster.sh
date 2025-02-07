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

# Loading in virtual environment and modules
module --force purge ; module load StdEnv/2020
source $HOME/virtenv/bin/activate
module load gcc/9.3.0 r/4.0.2
module load scipy-stack

# To run 24 jobs in array, based on this list of chromosomes
names=($(cat ../preprocess/chromos.txt))
files=${names[${SLURM_ARRAY_TASK_ID}]}

# The input: 
# 1: Project name and location
# 2: The type of SVs to annotate (ie. sv or svSim). This is a subdirectory that holds the vcf/csv file
# 3: The name of the file to annotate 

# This script leverages "scratch" space on an HPC given many files will be generated. The following creates directories to hold the files produced in scratch + those in the project space (for long term storage)

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

# We push this python script to create SURVIVOR configuration file used for SV simulation generation
python generating_lensIND.py \
/home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/SVlen/$3 \
${files} \
parameter_file \
$1 $2 \
>> /home/nboev/projects/def-sushant/nboev/data/SimulatedSVs/$1/$2/generating_lensIND_240424.py.txt

# Getting a list of the number of configuration files for this chromosome. In this case, I have chosen to batch across 250 real SVs for simulations
ls /home/nboev/scratch/data/SimulatedSVs/$1/$2/parameterIND/${files}/ \
> /home/nboev/scratch/data/SimulatedSVs/$1/$2/filenamesIND/filenames.${files}.txt

split -l 250 -d --additional-suffix=.txt \
/home/nboev/scratch/data/SimulatedSVs/$1/$2/filenamesIND/filenames.${files}.txt \
/home/nboev/scratch/data/SimulatedSVs/$1/$2/rawIND/${files}/${files}.split

splits=$(ls -1 /home/nboev/scratch/data/SimulatedSVs/$1/$2/rawIND/${files}/${files}.split* | wc -l)
echo "${splits}"

ls /home/nboev/scratch/data/SimulatedSVs/$1/$2/rawIND/${files}/${files}.split* \
> /home/nboev/scratch/data/SimulatedSVs/$1/$2/rawIND/${files}/splitnames.txt

# We now have to alter some parameters within the execute_SurvivorTemplate_240716.sh script, adjusted for each chromosome + number of splts required
cp execute_SurvivorTemplate.sh \
/home/nboev/scratch/data/SimulatedSVs/$1/$2/rawIND/${files}/execute_Survivor_${files}.sh

sed -i "s/PROJECT/$1/g" /home/nboev/scratch/data/SimulatedSVs/$1/$2/rawIND/${files}/execute_Survivor_${files}.sh
sed -i "s/LOCATION/$2/g" /home/nboev/scratch/data/SimulatedSVs/$1/$2/rawIND/${files}/execute_Survivor_${files}.sh
sed -i "s/array=0-0/array=0-${splits}/g" /home/nboev/scratch/data/SimulatedSVs/$1/$2/rawIND/${files}/execute_Survivor_${files}.sh
sed -i "s/CHROM/${files}/g" /home/nboev/scratch/data/SimulatedSVs/$1/$2/rawIND/${files}/execute_Survivor_${files}.sh

# We push the sbatch command script for this particular chromosome
sbatch /home/nboev/scratch/data/SimulatedSVs/$1/$2/rawIND/${files}/execute_Survivor_${files}.sh $1 $2

# example: 
# sbatch execute_svSimsMaster.sh 20220422_3202_phased_SNV_INDEL_SV_bychrom SVTrue_typedeletion_resTrue SVTrue_typedeletion_resTrue.csv



