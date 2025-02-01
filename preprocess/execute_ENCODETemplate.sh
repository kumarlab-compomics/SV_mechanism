#!/bin/bash
#SBATCH --nodes=1
#SBATCH --mem=300000
#SBATCH --time=04:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=nb.boev@mail.utoronto.ca
#SBATCH --job-name=ENCODECHROM
#SBATCH --output=/home/nboev/projects/def-sushant/nboev/preprocess/PROJECT/LOCATION/ENCODE/log/ENCODECHROM_%j.o
#SBATCH --error=/home/nboev/projects/def-sushant/nboev/preprocess/PROJECT/LOCATION/ENCODE/log/ENCODECHROM_%j.e
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

python adding_epiFeaturesSV_240517.py \
${files} \
$1 $2 \
${filenamenoext} \
CHROM \
/home/nboev/projects/def-sushant/nboev/data/ENCODE/$3/ \
$3 \
>> /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/adding_epiFeaturesSV_240517.py.txt

python adding_epiFeaturesflanks_240517.py \
${files} \
$1 $2 \
${filenamenoext} \
CHROM \
/home/nboev/projects/def-sushant/nboev/data/ENCODE/$3/ \
2000 \
$3 \
>> /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/adding_epiFeaturesflanks_240517.py.txt
