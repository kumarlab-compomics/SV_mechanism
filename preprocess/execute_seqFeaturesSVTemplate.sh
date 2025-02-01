#!/bin/bash
#SBATCH --nodes=1
#SBATCH --mem=300000
#SBATCH --time=00:30:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=nb.boev@mail.utoronto.ca
#SBATCH --job-name=seqFeaturesSVCHROM
#SBATCH --output=/home/nboev/projects/def-sushant/nboev/preprocess/PROJECT/LOCATION/seqFeaturesSV/log/seqFeaturesCHROM_%j.o
#SBATCH --error=/home/nboev/projects/def-sushant/nboev/preprocess/PROJECT/LOCATION/seqFeaturesSV/log/seqFeaturesCHROM_%j.e
#SBATCH --array=0-0

# Loading in virtual environment and modules
module --force purge ; module load StdEnv/2020
source $HOME/virtenv/bin/activate
module load scipy-stack

# For this chromosome, we have a job for each split
names=($(cat /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/SVlen/CHROM/splitnames.txt))
files=${names[${SLURM_ARRAY_TASK_ID}]}
filename=$(basename "$files")
filenamenoext=$(basename "${filename%.*}")

# The input (comes from execute_Master1.sh)
# 1: Project name and location
# 2: The type of SVs to annotate (ie. sv or svSim). This is a subdirectory that holds the vcf/csv file

# We now push this python script which will annotate the SV's sequence features
python adding_seqFeaturesSV_240118.py \
${files} \
$1 $2 \
${filenamenoext} \
CHROM \
>> /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/adding_seqFeaturesSV_240118.py.txt
