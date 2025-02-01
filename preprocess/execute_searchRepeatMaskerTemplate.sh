#!/bin/bash
#SBATCH --nodes=1
#SBATCH --mem=300000
#SBATCH --time=04:30:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=nb.boev@mail.utoronto.ca
#SBATCH --job-name=searchRepeatMaskerCHROM
#SBATCH --output=/home/nboev/projects/def-sushant/nboev/preprocess/PROJECT/LOCATION/searchRepeatMasker/log/searchRepeatMaskerCHROM_%j.o
#SBATCH --error=/home/nboev/projects/def-sushant/nboev/preprocess/PROJECT/LOCATION/searchRepeatMasker/log/searchRepeatMaskerCHROM_%j.e
#SBATCH --array=0-0

# Loading in virtual environment and modules
# In this case, we need repeatmasker
module --force purge ; module load StdEnv/2020 gcc/9.3.0
source $HOME/virtenv/bin/activate
module load scipy-stack
module load repeatmasker

# For this chromosome, we have a job for each split
names=($(cat /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/SVlen/CHROM/splitnames.txt))
files=${names[${SLURM_ARRAY_TASK_ID}]}
filename=$(basename "$files")
filenamenoext=$(basename "${filename%.*}")

# The input (comes from execute_Master1.sh)
# 1: Project name and location
# 2: The type of SVs to annotate (ie. sv or svSim). This is a subdirectory that holds the vcf/csv file

# Make a directory to hold the RepeatMasker results, thereby to not overwrite the directories RepeatMasker will generate
mkdir /home/nboev/projects/def-sushant/nboev/preprocess/PROJECT/LOCATION/searchRepeatMasker/CHROM/${filenamenoext}

# We push a python script which sends the DEL or INS sequence through RepeatMasker
python adding_searchRepeatMasker.py \
${files} \
$1 $2 \
${filenamenoext} \
CHROM \
>> /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/adding_searchRepeatMasker_240521.py.txt

rm /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/searchRepeatMasker/CHROM/${filenamenoext}/${filenamenoext}.fa.tbl
rm /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/searchRepeatMasker/CHROM/${filenamenoext}/${filenamenoext}.fa.masked
rm /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/searchRepeatMasker/CHROM/${filenamenoext}/${filenamenoext}.fa.cat

# We then push a python script which will append the RepeatMasker output to the vcf/csv file
python adding_searchRepeatMaskermerges.py \
${files} \
$1 $2 \
${filenamenoext} \
CHROM \
/home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/searchRepeatMasker/CHROM/${filenamenoext}/${filenamenoext}.fa.out \
>> /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/adding_searchRepeatMaskermerges_240521.py.txt

rm /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/searchRepeatMasker/CHROM/${filenamenoext}/${filenamenoext}.fa
rm /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/searchRepeatMasker/CHROM/${filenamenoext}/${filenamenoext}.fa.out
