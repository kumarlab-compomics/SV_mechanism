#!/bin/bash
#SBATCH --nodes=1
#SBATCH --mem=300000
#SBATCH --time=00:45:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=nb.boev@mail.utoronto.ca
#SBATCH --job-name=RepliSeqCHROM
#SBATCH --output=/home/nboev/projects/def-sushant/nboev/preprocess/PROJECT/LOCATION/RepliSeq/log/RepliSeqCHROM_%j.o
#SBATCH --error=/home/nboev/projects/def-sushant/nboev/preprocess/PROJECT/LOCATION/RepliSeq/log/RepliSeqCHROM_%j.e
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

# We push this python script to annotate the replication timing of the SV's position
# Note: This requires a file which is generated using a RepliSeq experiment (H1 cell line, 16-stage downloaded from 4D) then processed with RepliSeq
  # See: https://github.com/CL-CHEN-Lab/RepliSeq 
# The processed file (tempS50.csv) is held in the /data/RepliSeq/H1/ directory
python adding_RepliSeq.py \
/home/nboev/projects/def-sushant/nboev/data/RepliSeq/H1/Multi-stage/tempS50.csv \
${files} \
$1 $2 \
${filenamenoext} \
CHROM \
>> /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/adding_RepliSeq_240517.py.txt
