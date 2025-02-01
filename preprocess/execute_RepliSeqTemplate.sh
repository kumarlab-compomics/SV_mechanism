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

module --force purge ; module load StdEnv/2020
source $HOME/virtenv/bin/activate
module load scipy-stack

# looping through files from the split
names=($(cat /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/SVlen/CHROM/splitnames.txt))
files=${names[${SLURM_ARRAY_TASK_ID}]}
echo "${files}"
filename=$(basename "$files")
filenamenoext=$(basename "${filename%.*}")

echo "${filename}"
echo "${filenamenoext}"


# I think we need to add something in the python file so they don't overwrite!
python adding_RepliSeq_240517.py \
/home/nboev/projects/def-sushant/nboev/data/RepliSeq/H1/Multi-stage/tempS50.csv \
${files} \
$1 $2 \
${filenamenoext} \
CHROM \
>> /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/adding_RepliSeq_240517.py.txt

