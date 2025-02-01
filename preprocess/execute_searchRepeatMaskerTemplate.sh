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

module --force purge ; module load StdEnv/2020 gcc/9.3.0
source $HOME/virtenv/bin/activate
module load scipy-stack
module load repeatmasker

# looping through files from the split
names=($(cat /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/SVlen/CHROM/splitnames.txt))
files=${names[${SLURM_ARRAY_TASK_ID}]}
echo "${files}"
filename=$(basename "$files")
filenamenoext=$(basename "${filename%.*}")

echo "${filename}"
echo "${filenamenoext}"

mkdir /home/nboev/projects/def-sushant/nboev/preprocess/PROJECT/LOCATION/searchRepeatMasker/CHROM/${filenamenoext}

python adding_searchRepeatMasker_240521.py \
${files} \
$1 $2 \
${filenamenoext} \
CHROM \
>> /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/adding_searchRepeatMasker_240521.py.txt

rm /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/searchRepeatMasker/CHROM/${filenamenoext}/${filenamenoext}.fa.tbl
rm /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/searchRepeatMasker/CHROM/${filenamenoext}/${filenamenoext}.fa.masked
rm /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/searchRepeatMasker/CHROM/${filenamenoext}/${filenamenoext}.fa.cat

python adding_searchRepeatMaskermerges_240521.py \
${files} \
$1 $2 \
${filenamenoext} \
CHROM \
/home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/searchRepeatMasker/CHROM/${filenamenoext}/${filenamenoext}.fa.out \
>> /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/adding_searchRepeatMaskermerges_240521.py.txt


rm /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/searchRepeatMasker/CHROM/${filenamenoext}/${filenamenoext}.fa
rm /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/searchRepeatMasker/CHROM/${filenamenoext}/${filenamenoext}.fa.out
