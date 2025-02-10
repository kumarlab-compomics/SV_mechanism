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
#SBATCH --array=0-24

# Loading in virtual environment and modules
module --force purge ; module load StdEnv/2020
source $HOME/virtenv/bin/activate
module load gcc/9.3.0 r/4.0.2
module load scipy-stack

# To run 24 jobs in array, based on this list of chromosomes
names=($(cat ./chromos.txt))
files=${names[${SLURM_ARRAY_TASK_ID}]}

# The input: 
# 1: Project name and location
# 2: The type of Real SVs to annotate. This is a subdirectory that holds the vcf/csv file
# 3: The type of Simulated SVs to annotate. This is a subdirectory that holds the vcf/csv file

# Creating a directory to hold the chromosome and SV type specific z-score file we create
mkdir /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/zscore

# Recall, we have separate merged files for insertions and deletions. 
# We push the python script, adding_zscores.py, to match the matching real-simulations to calculate the z-scores for each feature.
# Therefore, for example, we want to match the real insertions to the simulated insertions

python adding_zscores.py \
/home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/merge_FeatMatrix/${files}.full.csv.DEL.FeatMatrix.csv \
$1 $2 \
20240625.${files}.deletion_FeatMatrix \
${files} \
/home/nboev/projects/def-sushant/nboev/preprocess/$1/$3/merge_FeatMatrix/${files}.full.csv.DEL.FeatMatrix.csv \
>> /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/adding_zscores_240610.py.txt

python adding_zscores.py \
/home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/merge_FeatMatrix/${files}.full.csv.INS.FeatMatrix.csv \
$1 $2 \
20240625.${files}.insertion_FeatMatrix \
${files} \
/home/nboev/projects/def-sushant/nboev/preprocess/$1/$3/merge_FeatMatrix/${files}.full.csv.INS.FeatMatrix.csv \
>> /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/adding_zscores_240610.py.txt


# Example
# sbatch execute_Master1Zscore.sh HGSVC2_v2_integrated_callset sv svSim


