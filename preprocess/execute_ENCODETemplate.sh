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
# 3: Cell line name to use for epigenomic annotation (must lead to the location of the ENCODE files)

# We now push this python script which will annotate the SVs' epigenetic features
python adding_epiFeaturesSV.py \
${files} \
$1 $2 \
${filenamenoext} \
CHROM \
/home/nboev/projects/def-sushant/nboev/data/ENCODE/$3/ \
$3 \
>> /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/adding_epiFeaturesSV.py.txt

# We now push this python script which will annotate the SVs flanking (2000bp) epigenetic features
python adding_epiFeaturesflanks.py \
${files} \
$1 $2 \
${filenamenoext} \
CHROM \
/home/nboev/projects/def-sushant/nboev/data/ENCODE/$3/ \
2000 \
$3 \
>> /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/adding_epiFeaturesflanks.py.txt
