#!/bin/bash
#SBATCH --nodes=1
#SBATCH --mem=3G
#SBATCH --time=3-16:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=nb.boev@mail.utoronto.ca
#SBATCH --job-name=flankSeqCHROM
#SBATCH --output=/home/nboev/projects/def-sushant/nboev/preprocess/PROJECT/LOCATION/flankSeq/log/flankSeqCHROM_%j.o
#SBATCH --error=/home/nboev/projects/def-sushant/nboev/preprocess/PROJECT/LOCATION/flankSeq/log/flankSeqCHROM_%j.e
#SBATCH --array=0-0

module --force purge ; module load StdEnv/2020 gcc/9.3.0 blast+/2.13.0
source $HOME/virtenv/bin/activate
module load scipy-stack
module load samtools
module load blast/2.2.26
module load r/4.0.2

# looping through files from the split
names=($(cat /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/SVlen/CHROM/splitnames.txt))
files=${names[${SLURM_ARRAY_TASK_ID}]}
echo "${files}"
filename=$(basename "$files")
filenamenoext=$(basename "${filename%.*}")

echo "${filename}"
echo "${filenamenoext}"


# I think we need to add something in the python file so they don't overwrite!
python adding_flankSeq_240517.py \
${files} \
$1 $2 \
${filenamenoext} \
CHROM \
2000 \
>> /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/adding_flankSeq_240517.py.txt

# Adding seq Features
python adding_seqFeatures_240118.py \
/home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/flankSeq/CHROM/${filenamenoext}.with2000.flankseq_samtools.csv \
$1 $2 \
${filenamenoext} \
CHROM \
2000 \
>> /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/adding_seqFeatures_240118.py.txt


# Were doing Blast and DNA shape together

python adding_BlastDNAShape_240522.py \
/home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/flankSeq/CHROM/${filenamenoext}.with2000.flankseq_samtools.csv \
$1 $2 \
${filenamenoext} \
CHROM \
2000 \
>> /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/adding_BlastDNAShape_240522.py.txt

python adding_Blastmerges_240522.py \
/home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/flankSeq/CHROM/${filenamenoext}.with2000.flankseq_samtools.csv \
$1 $2 \
${filenamenoext} \
CHROM \
/home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/Blast/CHROM/premerge/${filenamenoext}pre_post_resultsungap.txt \
/home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/Blast/CHROM/premerge/${filenamenoext}pre_sv_resultsungap.txt \
/home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/Blast/CHROM/premerge/${filenamenoext}post_sv_resultsungap.txt \
>> /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/adding_Blastmerges_240522.py.txt


# cleaning up the merged files
rm /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/Blast/CHROM/premerge/${filenamenoext}pre_post_resultsungap.txt
rm /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/Blast/CHROM/premerge/${filenamenoext}pre_sv_resultsungap.txt
rm /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/Blast/CHROM/premerge/${filenamenoext}post_sv_resultsungap.txt
rm /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/Blast/CHROM/premerge/${filenamenoext}SV_testing.csv


