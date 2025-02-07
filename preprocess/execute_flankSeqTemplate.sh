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

# Loading in virtual environment and modules
module --force purge ; module load StdEnv/2020 gcc/9.3.0 blast+/2.13.0
source $HOME/virtenv/bin/activate
module load scipy-stack
module load samtools
module load blast/2.2.26
module load r/4.0.2

# For this chromosome, we have a job for each split
names=($(cat /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/SVlen/CHROM/splitnames.txt))
files=${names[${SLURM_ARRAY_TASK_ID}]}
filename=$(basename "$files")
filenamenoext=$(basename "${filename%.*}")

# The input (comes from execute_Master1.sh)
# 1: Project name and location
# 2: The type of SVs to annotate (ie. sv or svSim). This is a subdirectory that holds the vcf/csv file

# We push this python script which pulls the 2000bp flanking sequences (up and downstream) using faidx
# Note : This script requires an accessible copy of the reference genome
python adding_flankSeq.py \
${files} \
$1 $2 \
${filenamenoext} \
CHROM \
2000 \
>> /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/adding_flankSeq_240517.py.txt

# We push this python script which annotates the SV flanks with its GC content, complexity, flexibility and complexity
# Similar to: adding_seqFeaturesSV.py
python adding_seqFeatures.py \
/home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/flankSeq/CHROM/${filenamenoext}.with2000.flankseq_samtools.csv \
$1 $2 \
${filenamenoext} \
CHROM \
2000 \
>> /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/adding_seqFeatures_240118.py.txt

# We push this script which takes the 2000bp flanking sequences to :
# 1. Perform various blast alignments to describe local homology
# 2. Employ DNAShapeR to determine local and SV DNA shape features
python adding_BlastDNAShape.py \
/home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/flankSeq/CHROM/${filenamenoext}.with2000.flankseq_samtools.csv \
$1 $2 \
${filenamenoext} \
CHROM \
2000 \
>> /home/nboev/projects/def-sushant/nboev/preprocess/$1/$2/adding_BlastDNAShape_240522.py.txt

# We pusht his script to merge the homology results from Blast to the original vcf/csv file
python adding_Blastmerges.py \
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


